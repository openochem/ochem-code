import numpy as np
import csv

from scipy.optimize import minimize
from sklearn.metrics import r2_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA

import bees

class base ():
    def __init__ (self, LSSVM, gpu=-1, center=True, pca=False, seed=0, metric='rmse'):
        self.do_center = center
        self.do_pca = pca
        self.seed = seed
        self.lssvm = LSSVM(gpu=gpu)
        self.rng = np.random.RandomState(self.seed)
        self.metric = metric
        self.device_str = "CPU" if gpu==-1 else "GPU"

    def load_train_data (self, path, ntargets):
        with open(path) as fh:
            reader = csv.reader(fh)
            data = []
            for row in reader:
                data.append(row)
        data = np.array(data, np.float)
        self.rng.shuffle(data)
            
        X = data[:,:-ntargets]
        Y = data[:,-ntargets:]

        if self.do_center:
            self.center_x = StandardScaler(with_std=True)
            self.center_y = StandardScaler(with_std=False)
            X = self.center_x.fit_transform(X)
            Y = self.center_y.fit_transform(Y)

        if self.do_pca:
            self.pca = PCA()
            X = self.pca.fit_transform(X)
            S = self.pca.singular_values_
            Shist = np.cumsum(S)
            Shist /= Shist[-1]
            self.use_comps = len(list(filter(lambda x: x, Shist/Shist[-1]<0.98)))
            X = X[:,:self.use_comps]

        self.X = X
        self.Y = Y
        self.train_mode = True
        self.iteration = 0
        del X,Y

    def load_test_data (self, path):
        with open(path) as fh:
            reader = csv.reader(fh)
            data = []
            for row in reader:
                data.append(row)
        X = np.array(data, np.float)

        if self.do_center:
            X = self.center_x.transform(X)

        if self.do_pca:
            X = self.pca.transform(X)
            X = X[:,:self.use_comps]

        self.X = X
        del X
        self.train_mode = False

    def cv_fit (self, nfold, shuffle=False):
        score = 0
        for train,test in KFold(nfold, shuffle=True, random_state=self.rng).split(self.X, self.Y):
            self.lssvm.fit(self.X[train], self.Y[train])
            score += self.lssvm.score(self.X[test], self.Y[test])**2*len(test)
        score = np.sqrt(score/len(test)/nfold)
        return score

    def train (self, nfold, glob_opt=False):
        self.lssvm.fit(self.X, self.Y)
        self.report()
        
        if glob_opt:
            self.optimize_bees(nfold)
            self.optimize_bfgs_polish(nfold)
        else:
            self.optimize_bfgs(nfold)

        self.lssvm.fit(self.X, self.Y)
        self.report()

    def apply (self, output_path):
        Y = self.lssvm.predict_cpu(self.X)
        if self.do_center:
            Y = self.center_y.inverse_transform(Y)
        with open(output_path, 'w') as fh:
            writer = csv.writer(fh, delimiter=',', quoting=csv.QUOTE_MINIMAL)
            for row in Y:
                writer.writerow(row)

    def report (self):
        self.iteration = self.iteration + 1
        if self.train_mode:
            print ("MESSAGE: epoch: %d train score: %.4g (%s) / %s" % (self.iteration, abs(self.lssvm.score(self.X, self.Y)), self.lssvm.metric, self.device_str))
            
            
    def optimize_bees (self, nfold):
        def fitness (vec):
            vec = 2**np.array(vec)
            gamma = vec[0]# if vec[0]>0.1 else 0.1
            lamb = vec[1]# if vec[1]>0.1 else 0.1
            rbf_gamma = vec[2]# if vec[0]>1e-6 else 1e-6

            self.lssvm.set_params(gamma=gamma, lamb=lamb, rbf_gamma=rbf_gamma)
            score = self.cv_fit(nfold, shuffle=True)
            return score
        H = bees.hive(fitness, np.log2([self.lssvm.gamma, self.lssvm.lamb, self.lssvm.rbf_gamma]),
                          g_range = [[-5, 10],
                                     [-5, 5],
                                     [-15, -2]],
                        range=[2,2,2], maxiter=30, g_density=10, density=4, ngood=5)
        res = H.run()
        self.lssvm.gamma = 2**res[0]
        self.lssvm.lamb = 2**res[1]
        self.lssvm.rbf_gamma = 2**res[2]
        return res

    def optimize_bfgs (self, nfold):
        def fitness (vec):
            self.lssvm.set_params(gamma=2**vec[0], lamb=2**vec[1], rbf_gamma=2**vec[2])
            score = self.cv_fit(nfold)
            return score
        
        solver_opts = {#'maxcor' : 3,
                       #'ftol' : 5e-4,
                       #'gtol' : 5e-4,
                       'disp': True}#,
                       #'maxiter' : 10}

        minimizer_opts = {'method' : 'L-BFGS-B',
                            'bounds' : [[-5, 10],
                                        [-5, 5],
                                        [-15, -2]]}
        res = minimize(fitness, np.log2([self.lssvm.gamma, self.lssvm.lamb, self.lssvm.rbf_gamma]),\
                            options=solver_opts, **minimizer_opts)
        self.lssvm.gamma = 2**res['x'][0]
        self.lssvm.lamb = 2**res['x'][1]
        self.lssvm.rbf_gamma = 2**res['x'][2]
        return res['x']

    def optimize_bfgs_polish (self, nfold):
        def fitness (vec):
            self.lssvm.set_params(gamma=vec[0], lamb=vec[1], rbf_gamma=vec[2])
            score = self.cv_fit(int(nfold*1.5))
            return score

        solver_opts = { 'disp': True }
            
        minimizer_opts = {'method' : 'L-BFGS-B',
                          'bounds' : [[0.25*self.lssvm.gamma,4.*self.lssvm.gamma],
                                      [0.25*self.lssvm.lamb,4.*self.lssvm.lamb],
                                      [0.5*self.lssvm.rbf_gamma,2.*self.lssvm.rbf_gamma]]}
        res= minimize(fitness, [self.lssvm.gamma, self.lssvm.lamb, self.lssvm.rbf_gamma], options=solver_opts, **minimizer_opts)

        self.lssvm.gamma = res['x'][0]
        self.lssvm.lamb = res['x'][1]
        self.lssvm.rbf_gamma = res['x'][2]
        return res['x']
