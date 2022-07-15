import tarfile
import shutil
import pandas as pd
import numpy as np
import sys
import os
import joblib
import csv

#New version

from cnf import CNF
from load import read_task
task = read_task(sys.argv[1])
                
def pack_model (savefile, savedir):
    with tarfile.open(savefile, "w") as tar:
        tar.add(savedir)    
    shutil.rmtree(savedir, ignore_errors=True)

def unpack_model (savefile, savedir):
    with tarfile.open(savefile, "r") as tar:
        tar.extractall()


def write_csv (smiles, Y, file):
    with open(file, 'w') as fh:
        writer = csv.writer(fh, delimiter=',', quoting=csv.QUOTE_MINIMAL)
        for i in range(len(smiles)):
            writer.writerow([smiles[i]] + Y[i].tolist())

def get_X (data):
    S = np.array(data.get("smiles")).reshape((-1,1))
    desc_header = list(filter(lambda x: 'desc' in x, data.columns))
    D = np.array(data.get(desc_header), np.float32)
    X = np.hstack((S,D))
    return X
            
if __name__ == "__main__":
    if task['train_mode']:
        data = pd.read_csv(task['train_data_file'], sep=",")
        Y = np.array(data.get(task['task_names']), np.float32)
        X = get_X(data)
        m = CNF(name="CNF", task=task)

        W = np.array(task['weights']) if task['weights'] else None
        print("start fitting")
        m.fit(X, Y, W)        
        m.augment = 0 # augmenting? Still?
        Yp = m.predict(X)
        print(">>> Train RMSE: {}".format("%.4g" % np.sqrt(np.mean((Yp-Y)**2))))
        joblib.dump(m, os.path.join("model_dir", "model.pkl"))
        pack_model(task['model_file'], "model_dir")
    
    if not task['train_mode'] or os.path.exists(task['apply_data_file']):
        unpack_model(task['model_file'], "model_dir")
        m = joblib.load(os.path.join("model_dir", "model.pkl"))
        m.gpu = task['gpu']
        data = pd.read_csv(task['apply_data_file'], sep=",")
        X = get_X(data)
        m.augment = 0 # augmenting? Still?
        Yp = m.predict(X)
        write_csv(X[:,0], Yp, task['result_file'])

        try:
            Y = np.array(data.get(task['task_names']), np.float32)
            print(">>> Application RMSE: {}".format("%.4g" % np.sqrt(np.mean((Yp-Y)**2))))
        except:
            pass

'''
data = pd.read_csv(task['apply_data_file'], sep=",")
X = get_X(data)
Yp = m.predict(X)
Y = np.array(data.get(task['task_names']), np.float32)
np.sqrt(np.mean((Yp-Y)**2))

'''
