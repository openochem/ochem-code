import numpy as np
import csv
import sys
sys.path.append("models")

def subdict (d, keys):
    return dict((k, d[k]) for k in keys if k in d.keys())

def libsvm_parse (path, rng=None):
    with open(path) as fh:
        reader = csv.reader(fh, delimiter=" ")
        data = []
        for row in reader:
            data.append(row)

            
    def get_dims (row):
        xdim = np.max(list(map(lambda x: int(x.split(":")[0]), row[1:])))
        ydim = len(row[0].split(","))
        return xdim, ydim
    
    xdim,ydim = get_dims(data[0])
    X = np.zeros((len(data), xdim), np.float)
    Y = np.zeros((len(data), ydim), np.float)

    
    order = np.arange(len(data))
    if rng is not None:
        rng.shuffle(order)
        
    for i in order:
        row = data[i]
        y = row[0].split(",")
        for j in range(ydim):
            val = y[j]
            if val == "":
                val = 0 #None
            else:
                val = float(val)
            Y[i,j] = val
            
        for pos,val in map(lambda x: x.split(":"), row[1:]):
            X[i, int(pos)-1] = float(val)
    return X,Y

class base ():
    def __init__ (self, task):
        self.task = task
        self.seed = task['seed']
        self.rng = np.random.RandomState(self.seed)

    def build_model (self):
        njobs = -1 if self.task['parallelize'] else 1
        if self.task['base_model'] == 'all':
            self.model = __import__(self.task['base_model']).builder(self.rng, njobs, self.ntargets!=1, self.task['cv_only'])
            self.model.calc_score = self.score
        else:
            self.model = __import__(self.task['base_model']).builder(self.rng, njobs, self.ntargets!=1)
            
    def load_train_data (self):
        self.X,self.Y = libsvm_parse(self.task['train_data_file'], self.rng)
        self.train_mode = True

        self.ntargets = self.Y.shape[1]
        if self.ntargets == 1:
            self.Y = self.Y.ravel()

    def load_test_data (self):
        self.X,_ = libsvm_parse(self.task['apply_data_file'])
        self.train_mode = False

    def train (self):
        if self.train_mode:
            self.model.fit(self.X, self.Y)
            
            
    def apply (self):
        Y_pred = self.model.predict(self.X).reshape((-1, self.ntargets))
        with open(self.task['result_file'], 'w') as fh:
            writer = csv.writer(fh, delimiter=',', quoting=csv.QUOTE_MINIMAL)
            for row in Y_pred:
                writer.writerow(row)
        fh.close()

    def score (self, X, Y):
        Y_pred = self.model.predict(X)
        return np.sqrt(np.mean((Y-Y_pred)**2))
            
    def validate (self, valid_file):
        # for dev purpose
        X,Y = libsvm_parse(valid_file)
        return self.score(X, Y)
