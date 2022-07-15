import numpy as np
import sys
from time import time
from pebble import ProcessPool
from concurrent.futures import TimeoutError
#from multiprocessing import Pool, TimeoutError

sys.path.append("models")

models_list = ['random_forest', 'ada_boost', 'extra_trees', 'linear', 'ada_boost_trees', 'gaussian_process', 'gradient_boosting', 'orthogonal_matching_pursuit_cv', 'bagging', 'huber', 'bagging_trees', 'kernel_ridge', 'bayesian_ridge', 'k_neighbors', 'ransac_trees', 'elastic_net_cv', 'lasso_lars_cv', 'ridge_cv', 'elastic_net', 'lasso_lars_ic', 'svm']

cv_models_list = ['orthogonal_matching_pursuit_cv', 'elastic_net_cv', 'lasso_lars_cv', 'ridge_cv']
from sklearn.model_selection import train_test_split

class full_sklearn ():
    def __init__ (self, rng, njobs, mlt, cv_only):
        self.cv_only = cv_only
        self.mlt = mlt
        self.rng = rng
        self.njobs = njobs
        if self.cv_only:
            self.models_list = cv_models_list
            self.timeout_coef = 10
        else:
            self.models_list = models_list
            self.timeout_coef = 3
        self.models = []
        for model_name in self.models_list:
            self.models.append(__import__(model_name).builder(self.rng, self.njobs, self.mlt))
        self.scores = [-1]*len(self.models)

    def worker (self, m, data):
        if self.cv_only:
            X,Y = data
            m1 = m.fit(X,Y)
            return m1, self.score(m1, X,Y)
        else:
            X_tr,X_v,Y_tr,Y_v = data
            m1 = m.fit(X_tr,Y_tr)
            return m1, self.score(m1, X_v,Y_v)
    
    def score (self, m, X, Y):
        Y_pred = m.predict(X)
        return np.sqrt(np.mean((Y-Y_pred)**2))

    def ensemble_score (self, data):
        if self.cv_only:
            X,Y = data
        else:
            X,Y = data[1],data[3]
        Y_pred = self.predict(X)
        return np.sqrt(np.mean((Y-Y_pred)**2))
    
    def fit (self, X, Y):
        start_moment = time()
        i = 0
        m = self.models[i]
        if self.cv_only:
            data = X,Y
        else:
            data = train_test_split(X,Y, test_size=0.2, random_state=self.rng)
            
        m,self.scores[i] = self.worker(m,data)
        rf_time = time()-start_moment
        print("MESSAGE: model {} trained with score: {}; time spent: {}s".format(self.models_list[0].upper(), round(self.scores[i], 4), round(rf_time, 1)))
        
        pool = ProcessPool(4)
        for j,m in enumerate(self.models[1:]):
            i = j+1
            try:
                start_moment = time()
                res = pool.schedule(self.worker, [m,data], timeout=rf_time*self.timeout_coef)
                self.models[i],self.scores[i] = res.result()
                print("MESSAGE: model {} trained with score: {}; time spent: {}s".format(self.models_list[i].upper(), round(self.scores[i], 4), round(time()-start_moment, 1)))
            except TimeoutError:
                print("MESSAGE: model {} stopped by timeout; time spent: {}s".format(self.models_list[i].upper(), round(time()-start_moment, 1)))
                if res.running():
                    res.cancel()
                self.scores[i] = -1
            except Exception as e:
                print("MESSAGE: model {} failed; reason: {}; time spent: {}s".format(self.models_list[i].upper(), e, round(time()-start_moment, 1)))
                self.scores[i] = -1
                if res.running():
                    res.cancel()
        print("MESSAGE: training finished; ensemble score: {}".format(self.ensemble_score(data)))
    
    def predict (self, X):
        result = None
        n = 0.
        for i,m in enumerate(self.models):
            sc = self.scores[i]
            if sc != -1:
                #if sc < 0.1:
                #    sc = 0.01
                sc = (1./sc)**5
                if result is None:
                    result = m.predict(X)*sc
                else:
                    result += m.predict(X)*sc
                n += sc
        return result/n

def builder (rng, njobs, mlt, cv_only):
    return full_sklearn(rng, njobs, mlt, cv_only)
