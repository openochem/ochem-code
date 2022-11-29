import tarfile
import shutil
import pandas as pd
import numpy as np
import random
import sys
import os
import joblib
import json
import csv
import tensorflow as tf

os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'

def pack_model (savefile, savedir):
    with tarfile.open(savefile, "w") as tar:
        tar.add(savedir)    
    shutil.rmtree(savedir, ignore_errors=True)

def unpack_model (savefile, savedir):
    with tarfile.open(savefile, "r") as tar:
        def is_within_directory(directory, target):
            
            abs_directory = os.path.abspath(directory)
            abs_target = os.path.abspath(target)
        
            prefix = os.path.commonprefix([abs_directory, abs_target])
            
            return prefix == abs_directory
        
        def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
        
            for member in tar.getmembers():
                member_path = os.path.join(path, member.name)
                if not is_within_directory(path, member_path):
                    raise Exception("Attempted Path Traversal in Tar File")
        
            tar.extractall(path, members, numeric_owner=numeric_owner) 
            
        
        safe_extract(tar)

from load import read_task
task = read_task(sys.argv[1])

random.seed(task['seed'])
tf.random.set_seed(task['seed'])
np.random.seed(task['seed'])
import deepchem as dc
from callback import ValidationCallback
from rdkit import Chem

_MODEL_DIR_ = "model_dir"

import glob
def remove_prev_checkpoints (model_dir):
    with open(model_dir+"/checkpoint", 'r') as fh:
        data = fh.read().split('\n')

    global tmp
    tmp = data

    ckpt_name_use = data[0].split('"')[1]
    for d in data[1:]:
        try:
            name = d.split('"')[1]
            if name != ckpt_name_use:
                for f in glob.glob(os.path.join(model_dir, name+'*')):
                    os.remove(f)
        except: pass
        

def write_csv (smiles, Y, file):
    with open(file, 'w') as fh:
        writer = csv.writer(fh, delimiter=',', quoting=csv.QUOTE_MINIMAL)
        for i in range(len(smiles)):
            y = Y[i].tolist()
            if np.isnan(y[0]):
                y = ["error"] *len(y)
            writer.writerow([smiles[i]] + y)

def get_X (data):
    if 'SMILES' in data.keys():
        name = 'SMILES'
    else:
        name = 'smiles'
    S = np.array(data.get(name)).tolist()
    desc_header = list(filter(lambda x: 'desc' in x, data.columns))
    D = np.array(data.get(desc_header), np.float32)
    #D = np.random.rand(len(S), 2) # for test
    if D.shape[1] == 0:
        D = None
    header = [name] + desc_header
    return S,D, header

def roc_auc_score (y_true, y_score, **kwargs):
    try:
        return dc.metrics.roc_auc_score(y_true, y_score, **kwargs)
    except:
        # fallback if got one label for a class
        # just return something feasible
        return dc.metrics.rms_score(y_true, y_score)
    
def fix_missed_values (Y):
    w = (~np.isnan(Y)).astype(np.float)
    Yfix = np.nan_to_num(Y.copy(), 0)
    return Yfix, w

def is_regression (Y):
    Y = np.nan_to_num(Y.copy(), 0)
    vals = np.sort(np.unique(Y))
    if len(vals) == 2 and vals[0] == 0 and vals[1] == 1:
        return False
    else:
        return True

class DeepchemModel ():
    def __init__ (self, task):
        self.task = task

    def set_mode (self, Y):
        if is_regression(Y):
            self.task['mode'] = 'regression'
            self.metric = dc.metrics.Metric(
                dc.metrics.rms_score, np.mean, mode=self.task['mode'])

        else:
            self.task['mode'] = 'classification'
            self.metric = dc.metrics.Metric(
                roc_auc_score, np.mean, mode=self.task['mode'])

        
    def fit_scale (self, Y):
        if self.task['mode'] == 'regression':
            self.mY = np.nanmean(Y, axis=0)
            Y = Y-self.mY
            self.stdY = np.nanstd(Y, axis=0)
            Y = Y/self.stdY
        
        return Y

    def scale (self, Y):
        if self.task['mode'] == 'regression':
            Y = (Y-self.mY)/self.stdY
        
        return Y

    def unscale (self, Y):
        if self.task['mode'] == 'regression':
            Y = Y*self.stdY+self.mY
            #Y = Y*torch.FloatTensor(self.stdY).to(Y.device)
            #Y = Y+torch.FloatTensor(self.mY).to(Y.device)
        
        return Y

    def featurize (self, S):
        featurizer = self.featurizer
        mols = [Chem.MolFromSmiles(smiles) for smiles in S]
        features = featurizer.featurize(mols)
        return features

    def transform (self, ds):
        try:
            transformer = self.transformer
        except:
            return ds

        return transformer.transform(ds)

    def set_max_size (self, S):
        if self.task['deepchem_model'] in ['DAG', 'ChemCeption']:
            n_atoms, n_chars = 0,0
            for s in S:
                nch = len(s)
                if nch > n_chars: n_chars = nch
                
                na = Chem.MolFromSmiles(s).GetNumAtoms()
                if na > n_atoms: n_atoms = na
                
            self.task['max_atoms'] = int(n_atoms*1.5)
            self.task['max_chars'] = int(n_chars*1.5)

        elif self.task['deepchem_model'] in ['TextCNN']:
            max_chars = np.max([len(s) for s in S])
            self.task['max_chars'] = int(max_chars*1.5)
    
    def fit (self, S, Y):
        self.set_mode(Y)
        self.set_max_size(S)
        Y = self.fit_scale(Y)
        
        self.dc_model_builder(**self.get_kwargs())

        Y,w = fix_missed_values(Y)
        dataset = dc.data.NumpyDataset(X=self.featurize(S), y=Y, ids=S, w=w)
        dataset = self.transform(dataset)

        train_split = int(self.task['train_proportion']*len(dataset))
        train_ids = np.arange(len(dataset))[:train_split]
        train_ds = dataset.select(train_ids)
        valid_ids  = np.arange(len(dataset))[train_split:]
        valid_ds = dataset.select(valid_ids)

        validation_interval = int(np.ceil(len(train_ds)/self.task['batch_size']))
        callback = ValidationCallback(
            valid_ds, validation_interval, [self.metric],
            save_dir=_MODEL_DIR_),

        self.model.fit(train_ds, nb_epoch=self.task['n_epochs'],
                       callbacks=callback)

    def predict (self, S):
        # deepchem requires targets even at prediction step, feeding zeros
        dataset = dc.data.NumpyDataset(X=self.featurize(S),
                                       y=np.zeros((len(S), len(self.task['task_names']))),
                                       ids=S)
        dataset = self.transform(dataset)

        return self.model.predict(dataset)

    def load (self):
        
        self.dc_model_builder(**self.get_kwargs(), model_dir=_MODEL_DIR_)
        self.model.restore()
        
################################################################
        
class GraphConvModel (DeepchemModel):
    def __init__ (self, task):
        super().__init__(task)
        self.featurizer = dc.feat.ConvMolFeaturizer()

    def dc_model_builder (self, **kwargs):
        self.model = dc.models.graph_models.GraphConvModel(**kwargs)
        self.model.optimizer.learning_rate = self.task['learning_rate']

    def get_kwargs (self):
        return { "n_tasks" : len(self.task['task_names']),
                 "batch_size" : self.task['batch_size'],
                 "graph_conv_layers" : self.task['graph_conv_layers'], # [64,64]
                 "dense_layer_size" : self.task['dense_layer_size'], # 128
                 "dropout" : self.task['dropout'],
                 "number_atom_features" : self.task['number_atom_features'], # 75
                 "mode" : self.task['mode']}

class DAGModel (DeepchemModel):
    def __init__ (self, task):
        super().__init__(task)
        self.featurizer = dc.feat.ConvMolFeaturizer()

    def dc_model_builder (self, **kwargs):
        self.transformer = dc.trans.DAGTransformer(max_atoms=self.task['max_atoms'])
        self.model = dc.models.graph_models.DAGModel(**kwargs)
        
    def get_kwargs (self):
        return { "n_tasks" : len(self.task['task_names']),
                 "max_atoms" : self.task['max_atoms'],
                 "batch_size" : self.task['batch_size'],
                 'n_atom_feat' : 75,
                 'n_graph_feat' : self.task['n_graph_feat'], #30
                 'n_outputs' :  self.task['n_outputs'], #30,
                 'layer_sizes' : self.task['layer_sizes'], #[100]
                 'layer_sizes_gather' : self.task['layer_sizes_gather'], #[100],
                 "dropout" : self.task['dropout'],
                 "learning_rate" : self.task['learning_rate'],
                 "mode" : self.task['mode'],
                 "use_queue" : False }

class WeaveModel (DeepchemModel):
    # pip install tensorflow-probability
    def __init__ (self, task):
        super().__init__(task)
        self.featurizer = dc.feat.WeaveFeaturizer(use_chirality=True)

    def dc_model_builder (self, **kwargs):
        self.model = dc.models.graph_models.WeaveModel(**kwargs)
        self.model.optimizer.learning_rate = self.task['learning_rate']

    def get_kwargs (self):
        return { "n_tasks" : len(self.task['task_names']),
                 'n_atom_feat' : 78, # 75 if use_chirality=False
                 'n_pair_feat' : 18, # 14 if use_chirality=False
                 'n_graph_feat' : self.task['n_graph_feat'], #128
                 'n_weave' : self.task['n_weave'], #2
                 'n_hidden' : self.task['n_hidden'], #50
                 "batch_size" : self.task['batch_size'],
                 "dropouts" : self.task['dropout'],
                 'use_queue' : False,
                 "mode" : self.task['mode']}


class MPNNModel (DeepchemModel):
    def __init__ (self, task):
        super().__init__(task)
        self.featurizer = dc.feat.WeaveFeaturizer(use_chirality=True)

    def get_kwargs (self):        
        return { "n_tasks" : len(self.task['task_names']),
                 'n_atom_feat' : 78, #75 if use_chirality=False
                 'n_pair_feat' : 18, #14 if use_chirality=False
                 'T' : 3,
                 'M' : 5,
                 "learning_rate" : self.task['learning_rate'],
                 "batch_size" : self.task['batch_size'],
                 "dropouts" : self.task['dropout'],
                 'use_queue' : False,
                 "mode" : self.task['mode'] }

    def dc_model_builder (self, **kwargs):
        self.model = dc.models.graph_models.MPNNModel(**kwargs)

class TextCNNModel (DeepchemModel):
    def __init__ (self, task):
        super().__init__(task)
        
        self.featurizer = dc.feat.RawFeaturizer(smiles=False)
        self.make_char_dict()
        
    def make_char_dict (self):
        ascii_lowercase = 'abcdefghijklmnopqrstuvwxyz'
        ascii_uppercase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        digits = '0123456789'
        punctuation = r"""!"#$%&'()*+,-./:;<=>?@[\]^_`{|}~"""
        tokens = list(ascii_lowercase) + \
                 list(ascii_uppercase) + list(digits) + list(punctuation)
        self.char_dict = {v:k for k,v in enumerate(tokens)}

    def get_kwargs (self):
        return { "n_tasks" : len(self.task['task_names']),
                 "char_dict" : self.char_dict,
                 "seq_length" : self.task['max_chars'],
                 "n_embeddings" : self.task['n_embeddings'], #75
                 "kernel_sizes" : self.task['kernel_sizes'], #[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20],
                 "num_filters" : self.task['num_filters'], #[100, 200, 200, 200, 200, 100, 100, 100, 100, 100, 160, 160]
                 "learning_rate" : self.task['learning_rate'],
                 "batch_size" : self.task['batch_size'],
                 "dropout" : self.task['dropout'],
                 'use_queue' : False,
                 "mode" : self.task['mode'] }
        
    def dc_model_builder (self, **kwargs):
        self.model = dc.models.TextCNNModel(**kwargs)

class ChemCeptionModel (DeepchemModel):
    # not working, including their own example
    def __init__ (self, task):
        super().__init__(task)
        pass

    def get_kwargs (self):
        return { "n_tasks" : len(self.task['task_names']),
                 "img_spec" : "engd", 
                 "mode" : self.task['mode'] }

    def dc_model_builder (self, **kwargs):
        self.featurizer = dc.feat.SmilesToImage(img_size=self.task['max_atoms'],
                                                img_spec='std')
        self.model = dc.models.ChemCeption(**kwargs)
        self.model.optimizer.learning_rate = self.task['learning_rate']

        
def model_builder (task):
    name = task['deepchem_model']
    if name == 'GraphConv':
        return GraphConvModel(task)
    elif name == "DAG":
        return DAGModel(task)
    elif name == "Weave":
        return WeaveModel(task)
    elif name == 'ChemCeption':
        return ChemCeptionModel(task)
    elif name == 'MPNN':
        return MPNNModel(task)
    elif name == 'TextCNN':
        return TextCNNModel(task)
        
if __name__ == "__main__":
    if task['train_mode']:
        shutil.rmtree(_MODEL_DIR_, ignore_errors=True)
        
        data = pd.read_csv(task['train_data_file'], sep=",")
        S,D,_ = get_X(data)
        Y = np.array(data.get(task['task_names']), np.float32)
        
        m = model_builder(task)
        m.fit(S, Y)
        joblib.dump(m.task, os.path.join(_MODEL_DIR_, "train_task.pkl"))
        remove_prev_checkpoints(_MODEL_DIR_)
        pack_model(task['model_file'], _MODEL_DIR_)

    if not task['train_mode'] or os.path.exists(task['apply_data_file']):
        unpack_model(task['model_file'], _MODEL_DIR_)
        train_task = joblib.load(os.path.join(_MODEL_DIR_, "train_task.pkl"))
        m = model_builder(train_task)
        m.load()
        
        data = pd.read_csv(task['apply_data_file'], sep=",")
        S,D,_ = get_X(data)
        Yp = m.predict(S)
        if train_task['mode'] == 'classification':
            Yp = Yp[:,:,1]
        else:
            if m.task['deepchem_model'] in ['TextCNN']:
                Yp = Yp[:,:,0]
            
        write_csv(S, Yp, task['result_file'])
