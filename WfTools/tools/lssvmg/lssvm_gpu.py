import numpy as np
import cupy as cp
import os
from os.path import expanduser, join
from sklearn.metrics.pairwise import rbf_kernel as rbf_kernel_cpu

def reshape (mat, shape=(-1,1)):
    if type(mat) is np.ndarray:
        return np.rollaxis(mat, 1).reshape(shape)
    else:
        return cp.rollaxis(mat, 1).reshape(shape)

def reshape_F (mat, shape=(-1,1)):
    if type(mat) is np.ndarray:
        return mat.reshape(shape, order="F")
    else:
        return mat.reshape(shape[::-1]).T

def to_gpu (arr):
    return cp.array(arr)

################################################################
### rbf gpu

def rbf_kernel (X, Y, gamma):
    K = euclidean_distances(X,Y)
    K *= -gamma
    cp.exp(K,out=K)
    return K


def row_norms (X):
    return cp.sum(X**2, axis = -1)

### this brings gpu memory issues
# def row_norms(X):
#     return cp.einsum('ij,ij->i', X, X)

def euclidean_distances (X, Y):
    XX = row_norms(X)[:, np.newaxis]    
    if X is Y:  # shortcut in the common case euclidean_distances(X, X)
        YY = XX.T
    else:
        YY = row_norms(Y)[np.newaxis, :]
        
    distances = cp.dot(X, Y.T)
    distances *= -2
    distances += XX
    distances += YY
    cp.maximum(distances, 0, out=distances)
    
    if X is Y:
        # Ensure that distances between vectors and themselves are set to 0.0.
        # This may not be the case due to floating point rounding errors.
        # distances.ravel()[::distances.shape[0] + 1] = 0.0
        cp.fill_diagonal(distances, 0)
    
    return distances

################################################################



class LSSVM_GPU ():
    def __init__ (self, gpu=0, gamma=0.5, lamb=2, rbf_gamma=0, metric='rmse'):
        if gpu == -1:
            gpu = 0
        
        os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"   # see issue #152
        os.environ["CUDA_VISIBLE_DEVICES"]=str(gpu)

#        cp.cuda.Device(gpu).use()
        self.gamma = gamma
        self.lamb = lamb
        self.rbf_gamma = rbf_gamma
        self.fitted = False
        self.metric = metric

    def get_params (self, deep=False):
        return {'gamma' : self.gamma, 'lamb' : self.lamb, 'rbf_gamma' : self.rbf_gamma}


    def set_params (self, **params):
        if params.get('gamma') is not None:
            self.gamma = params['gamma']
        if params.get('lamb') is not None:
            self.lamb = params['lamb']
        if params.get('rbf_gamma') is not None:
            self.rbf_gamma = params['rbf_gamma']
        return self

    def fit (self, X, Y):
        
        if self.rbf_gamma == 0: # rbf gamma guess
            self.rbf_gamma = 1./X.shape[1]
            
        if len(X) != len(Y):
            raise Exception("different num of samples in X and Y")
        self.X_train = X
        self.Y_train = Y
        X = to_gpu(X)
        Y = to_gpu(Y)
        n, outdim = Y.shape
        self.outdim = outdim

        K = rbf_kernel(X,X, self.rbf_gamma)
        H = cp.tile(K,(outdim,outdim))
        H.ravel()[::H.shape[1]+1] += 1./self.gamma
        P = cp.zeros((outdim*n,outdim))
        
        for i in range(outdim):
            id1 = n*i
            id2 = n*(i+1)
            H[id1:id2,id1:id2] += K*(float(outdim)/self.lamb)
            P[id1:id2,i] = cp.ones(n)

        del K
        eta = cp.linalg.solve(H,P)
        nu = cp.linalg.solve(H, reshape(Y))
        del H
        S = cp.dot(P.T, eta)
        b = cp.dot(cp.dot(cp.linalg.inv(S), eta.T), reshape(Y))
        alpha = nu - cp.dot(eta, b)
        alpha = reshape_F(alpha, (n,outdim))

        self.alpha = alpha.get()
        self.b = b.get()
        self.fitted = True
        return self
    
    def predict (self, X):
        if not self.fitted:
            raise Exception("Fit model before applying")
        
        n = len(self.X_train)
        n_test = len(X)

        X = to_gpu(X)
        X_train = to_gpu(self.X_train)
        alpha = to_gpu(self.alpha)
        b = to_gpu(self.b)
        
        K = rbf_kernel(X, X_train, self.rbf_gamma)
        Y_pred = cp.tile( cp.sum(cp.dot(K, alpha), axis=1).reshape((-1,1)), (1,self.outdim)) + cp.dot(K,alpha)*(float(self.outdim)/self.lamb) + cp.tile(b.T,(n_test,1))

        return Y_pred.get()

    def predict_cpu (self, X):
        if not self.fitted:
            raise Exception("Fit model before applying")
        
        n = len(self.X_train)
        n_test = len(X)

        X_train = self.X_train
        alpha = self.alpha
        b = self.b
        
        K = rbf_kernel_cpu(X, X_train, self.rbf_gamma)
        Y_pred = np.tile( np.sum(np.dot(K, alpha), axis=1).reshape((-1,1)), (1,self.outdim)) + np.dot(K,alpha)*(float(self.outdim)/self.lamb) + np.tile(b.T,(n_test,1))

        return Y_pred

    
    def score (self, X, Y):
        outdim = Y.shape[1]
        Y_pred = self.predict(X)
        if self.metric == 'rmse':
            return np.sqrt(np.mean((Y-Y_pred)**2))
        elif self.metric == 'mae':
            return np.mean(np.abs(Y-Y_pred))
        elif self.metric == 'r2':
            return -1*r2_score(Y, Y_pred, multioutput='variance_weighted')
        else:
            raise Exception('unknown metric')
