import numpy as np

def reshape (mat, shape=(-1,1)):
    return mat.reshape(shape, order="F")

from sklearn.metrics.pairwise import rbf_kernel
from sklearn.metrics import r2_score

class LSSVM ():
    def __init__ (self, gpu=-1, gamma=0.5, lamb=2, rbf_gamma=0, metric='rmse'):
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
        n, outdim = Y.shape
        self.outdim = outdim

        K = rbf_kernel(X,X, self.rbf_gamma)
        H = np.tile(K,(outdim,outdim)) + np.eye(outdim*n) / self.gamma
        P = np.zeros((outdim*n,outdim))

        for i in range(outdim):
            id1 = n*i
            id2 = n*(i+1)
            H[id1:id2,id1:id2] += K*(float(outdim)/self.lamb)
            P[id1:id2,i] = np.ones(n)

        eta = np.linalg.solve(H,P)
        nu = np.linalg.solve(H, reshape(Y))
        S = np.dot(P.T, eta)
        b = np.dot(np.dot(np.linalg.inv(S), eta.T), reshape(Y))
        alpha = nu - np.dot(eta, b)
        alpha = reshape(alpha, (n,outdim))
        
        self.alpha = alpha
        self.b = b
        self.fitted = True
        return self

    def predict (self, X):
        if not self.fitted:
            raise Exception("Fit model before applying")
        
        n = len(self.X_train)
        n_test = len(X)
        K = rbf_kernel(X, self.X_train, self.rbf_gamma)
        Y_pred = np.tile( reshape(np.sum(np.dot(K, self.alpha), axis=1)), (1,self.outdim)) + np.dot(K,self.alpha)*(float(self.outdim)/self.lamb) + np.tile(self.b.T,(n_test,1))

        return Y_pred

    def predict_cpu (self, X):
        return self.predict(X)
    
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

    
