import csv
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA

from lssvm import lssvm

# params:
# init:      kernel (str), center (bool), whiten (bool), reduce_dim (bool), nfold (int)
# load_data: path (str), target_columns (int[]), test_prop (float)
# train:     optimize (bool), gamma (float), lamb (float), rbf_gamma (float), metric (str)
# apply:     output_path (str)



class model ():
    def __init__ (self, center=True, whiten=True, reduce_dim=True, nfold=3):
        self.center = center
        self.whiten = whiten
        self.reduce_dim = reduce_dim
        self.nfold = nfold
        
        self.y_avail = False
        self.test_avail = False
        

    def load_data (self, path, target_columns=[], test_prop=0):
        self.test_prop = test_prop
        self.test_avail = False
        self.y_avail = False

        with open(path) as fh:
            reader = csv.reader(fh)
            data = []
            for row in reader:
                data.append(row)
        data = np.array(data, np.float)
        ncols = data.shape[1]

        if len(target_columns)>0:
            rand_state = np.random.RandomState(0)
            rand_state.shuffle(data)
            self.y_avail = True
            
            self.target_columns = map(lambda x: x if x>=0 else x+ncols, target_columns)
            self.input_columns = filter(lambda x: x not in self.target_columns, np.arange(ncols))
            self.X = data[:,self.input_columns]
            self.Y = data[:,self.target_columns]
            
            if self.test_prop > 0:
                self.test_avail = True
                self.X, self.Xtest, self.Y, self.Ytest = train_test_split(self.X, self.Y, test_size=self.test_prop, shuffle=True)

        else:
            self.X = data

    def fit_preproc (self):

        if self.center:
            self.center_x = StandardScaler(with_std=True)
            self.X = self.center_x.fit_transform(self.X)
#            if self.y_avail:
            self.center_y = StandardScaler(with_std=False)
            self.Y = self.center_y.fit_transform(self.Y)

        if self.whiten:
            self.decomposer = PCA(whiten=False)#True)
            self.X = self.decomposer.fit_transform(self.X)
            if self.reduce_dim:
                S = self.decomposer.singular_values_
                Shist = np.cumsum(S)
                Shist /= Shist[-1]
                self.use_comps = len(filter(lambda x: x, Shist/Shist[-1]<0.98))
                self.X = self.X[:,:self.use_comps]

    def preproc (self, X, Y=None):
        if self.center:
            X = self.center_x.transform(X)
            if Y is not None:
                Y = self.center_y.transform(Y)
        if self.whiten:
            X = self.decomposer.transform(X)
            if self.reduce_dim:
                X = X[:,:self.use_comps]
        if Y is not None:
            return X,Y
        else:
            return X

    def postproc (self, Y):
        if self.center:
            Y = self.center_y.inverse_transform(Y)
        return Y

    def train (self, kernel='rbf', optimize=True, optimize_algo='bfgs_log', gamma=0.5, lamb=2.0, rbf_gamma=0, metric='rmse'):
        if self.y_avail:
            self.lssvm = lssvm(kernel=kernel, gamma=gamma, lamb=lamb, rbf_gamma=rbf_gamma, metric=metric)
            print "%d-fold cv score: %.4f (%s)" % (self.nfold, abs(self.lssvm.cv_fit(self.X, self.Y, self.nfold)), metric)
            self.lssvm.fit(self.X, self.Y)
            self.report()
            if optimize:
                self.lssvm.optimize(self.X, self.Y, self.nfold, optimize_algo[0])
                print "\nScores after hyper optimization"
                self.lssvm.cv_fit(self.X, self.Y, self.nfold)
                print "%d-fold cv score: %.4f (%s)" % (self.nfold, abs(self.lssvm.cv_fit(self.X, self.Y, self.nfold)), metric)
                self.lssvm.fit(self.X, self.Y)
                self.report()

    def report (self):
        print "train score: %.4g (%s)" % (abs(self.lssvm.score(self.X, self.Y)), self.lssvm.metric)
        if self.test_avail:
            X,Y = self.preproc(self.Xtest, self.Ytest)
            print "validation score: %.g (%s)" % (abs(self.lssvm.score(X,Y)), self.lssvm.metric)

    def apply (self, output_path):
        Xp = self.preproc(self.X)
        Yp = self.lssvm.predict(Xp)
        Yp = self.postproc(Yp)
        with open(output_path, 'w') as fh:
            writer = csv.writer(fh, delimiter=',', quoting=csv.QUOTE_MINIMAL)
            for row in Yp:
                writer.writerow(row)
