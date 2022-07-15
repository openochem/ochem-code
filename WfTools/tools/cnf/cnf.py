import uuid
import csv
import pandas as pd
import numpy as np
from rdkit import Chem
import time
#import tensorflow as tf
import tensorflow.compat.v1 as tf
import os
import shutil
import pickle
from time import sleep
from collections import Counter
import sys
import re
import math
from multiprocessing import Process, Lock, Queue
import random as python_random
from RAdam import RAdamOptimizer

# Bug fix version: correct augmentation, two Processes, significantly optimized speed

tf.disable_v2_behavior()

def uid ():
    return '{:x}'.format(int(uuid.uuid4()))

def scatter_mean (data, idx):
    result = []
    cur_x = None
    cur_count = 0
    cur_idx = idx[0]
    for i in range(len(data)):
        if cur_idx != idx[i]:
            result.append(cur_x/cur_count)
            cur_x = None
            cur_count = 0
        
        cur_idx = idx[i]
        
        if cur_x is None:
            cur_x = data[i].copy()
        else:
            cur_x += data[i]
        cur_count += 1
    if cur_count > 0:
        result.append(cur_x/cur_count)
    return np.array(result)


################################################################


################################################################

#global_tokens = set({u'\xae': 0, u'!': 1, u' ': 2, u'#': 3, u'%': 4, u'&': 5, u')': 6, u'(': 7, u'+': 8, u'*': 9, u'-': 10, u'/': 11, u'.': 12, u'1': 13, u'0': 14, u'3': 15, u'2': 16, u'5': 17, u'4': 18, u'7': 19, u'6': 20, u'9': 21, u'8': 22, u':': 23, u'=': 24, u'A': 25, u'@': 26, u'C': 27, u'B': 28, u'E': 29, u'D': 30, u'G': 31, u'F': 32, u'I': 33, u'H': 34, u'K': 35, u'J': 36, u'M': 37, u'L': 38, u'O': 39, u'N': 40, u'Q': 41, u'P': 42, u'S': 43, u'R': 44, u'U': 45, u'T': 46, u'W': 47, u'V': 48, u'Y': 49, u'X': 50, u'[': 51, u'Z': 52, u']': 53, u'\\': 54, u'a': 55, u'c': 56, u'e': 57, u'g': 58, u'i': 59, u'l': 60, u'o': 61, u'n': 62, u's': 63, u'r': 64, u'u': 65, u'\xa8': 66, u'|': 67})
global_tokens = set([' ', '!', '#', '%', '&', '(', ')', '*', '+', '-', '.', '/', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ':', '=', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '[', '\\', ']', 'a', 'c', 'e', 'g', 'i', 'l', 'n', 'o', 'r', 's', 'u', '|'])

def can_smile (sm):
    m = Chem.MolFromSmiles(sm)
    if m is None:
        return smil
    N =m.GetNumAtoms()
    if N==0:
        return sm
    try:
        t= Chem.MolToSmiles(m, canonical=True, isomericSmiles=True)
    except :
        return sm
    return t

def smiles_char_tokenizer (smi):
    smi = smi.replace("Cl", "L").replace("Br","R").replace("Pt","T")
    return [ch for ch in smi]
    
def smiles_atom_tokenizer (smi):
    pattern =  "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
    regex = re.compile(pattern)
    tokens = [token for token in regex.findall(smi)]
    return tokens

def smiles_tokenizer (smi, type="char"):
    return {"char":smiles_char_tokenizer, "atom":smiles_atom_tokenizer}[type](smi)


################################################################

def emptySmile (sm):
    augs = []
    #augs.append(sm) -- does not help: actually makes sense only for very small molecules
    return augs

def aug_smile (sm, nmax, exclude):
    m = Chem.MolFromSmiles(sm)

    augs = []

    if m is None:
        emptySmile(sm)
    
    N = m.GetNumAtoms()
    if N==0:
        emptySmile(sm)
    
    if nmax >= 2:
        max = nmax*2
    else:
        max = 1
    
    aids = list(range(N))
    np.random.shuffle(aids)
    m = Chem.RenumberAtoms(m,aids)

    for i in range(max):
        try:
            n = np.random.randint(0,N)
            aug = Chem.MolToSmiles(m, rootedAtAtom=n, doRandom=True, canonical=False, isomericSmiles=True)
            if aug not in augs and aug not in exclude:
                augs.append(aug)
        except:
            pass

    if len(augs) == 0:
        emptySmile(sm)
    
    return augs[:nmax]

def aug_smiles (smiles, nmax=10, exclude=[]):
    res = []
    idx = []

    for i,sm in enumerate(smiles):
        augs = aug_smile(sm, nmax=nmax, exclude=exclude)
        res += augs
        idx += [i]*len(augs)
  
    return np.array(res), idx

sess_config = tf.ConfigProto(allow_soft_placement=True, log_device_placement=False)
sess_config.gpu_options.allow_growth = True

def check_device (gpu):
    if gpu == -1:
        device_name = 'cpu:0'
        device_str = "CPU"
    else:
        device_name = "gpu:0"#+str(gpu)
        session =  tf.Session(config=sess_config)
        if True in (map(lambda x: device_name.upper() in x.name, session.list_devices())):
            #os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"
            #os.environ["CUDA_VISIBLE_DEVICES"]=str(gpu)
            try:
                from cupy.cuda import cudnn
                cudnn.getVersion()
                cudnn = True
            except:
                cudnn = False
        else:
            raise Exception ('Device {} not available, fallback to cpu'.format(device_name))
        session.close()
        if device_name != "cpu:0":
            if cudnn:
                device_str = "GPU+CuDNN"
            else:
                device_str = "GPU"
        else:
            device_str = "CPU"
    return device_name, device_str

################################################################

def tf_contrib_layer_norm(input_tensor, name=None):
  return tf.keras.layers.LayerNormalization(
      name=name, axis=-1, epsilon=1e-12, dtype=tf.float32
  )(input_tensor)


class CNF ():
    def __init__ (self, name=None, task={}):
        def get_value (name):
            val = task.get(name)
            if val is None: raise Exception(f"Parameter {name} is not defined")
            return val

        self.gpu = get_value('gpu')
        self.seed = get_value('seed')
        self.batch_size = get_value('batch_size')
        self.n_epochs = get_value('n_epochs')

        self.normalize = get_value('normalize')
        self.canonize = get_value('canonize')
        self.augment = get_value('augment')
        if self.augment != 0:
            self.canonize = False
        self.taxonomy_graph = get_value('graph')
        self.error_function = get_value('error_function')
        self.exp_coef = get_value('exp_coef')
        
        self.mlp_dims = get_value('mlp_dims')
        self.learning_rate = get_value('learning_rate')
        self.nfp_type = get_value('nfp_type')
        self.n_fp_layers = get_value('n_fp_layers')
        self.fp_layer_dim = get_value('fp_layer_dim')
        self.filter_size = get_value('filter_size')
        self.relu_type = get_value('relu_type')
        self.nfp_alpha = get_value('nfp_alpha')
        self.debug = get_value('debug')
        self.highway = get_value('highway')
        self.allow_overfit = get_value('allow_overfit')
        self.tokenizer_type = get_value('tokenizer_type')
        if self.tokenizer_type == 'old':
            self.tokenizer_type = 'char'
            self.preload_dictionary = True
        else:
            self.preload_dictionary = False
        
        self.dropout = get_value('dropout')
        if self.dropout == -1:
            self.dropout = 0.3
            self.adjust_dropout = True
        else:
            self.adjust_dropout = False

        self.early = get_value('early')

        self.name = uid() if name is None else name
        self.model_path = "model_dir/model.cpkt"
        if os.path.exists("model_dir"):
            shutil.rmtree("model_dir", ignore_errors=True)
        os.makedirs("model_dir")

        self.validation_interval = 5.
        if self.preload_dictionary:
            self.tokens = global_tokens
        else:
            self.tokens = set()
        
    def fill_dictionary (self, smiles):
        for sm in smiles:
            self.tokens |= set(smiles_tokenizer(sm, self.tokenizer_type))
        self.tokens = dict(zip(sorted(list(self.tokens)), list(range(len(self.tokens)))))
        print(">>> Using smiles chars: {}".format(' '.join(sorted(list(set(self.tokens))))))
        if not self.preload_dictionary:
            self.tokens["pad"] = len(self.tokens) # gpu/cpu embedding fix
        
    def embed_smiles (self, smiles):
        lmax = max([len(sm) for sm in smiles])
        def embed_smile (sm):
            pad_id = -1 if self.preload_dictionary else self.tokens["pad"] # gpu/cpu embedding fix
            tokens = [self.tokens.get(char, pad_id) for char in smiles_tokenizer(sm, self.tokenizer_type)]
            padsize = lmax-len(tokens)
            return np.array(tokens + [pad_id]*padsize, np.float)
        return np.array([embed_smile(sm) for sm in smiles], np.float)


    def fit (self, X, Y, W=None):
        np.random.seed(self.seed)
        python_random.seed(self.seed)
        tf.set_random_seed(self.seed)
        device_name, device_str = check_device(self.gpu)
        self.n_targets = Y.shape[1]

        m1 = np.ones_like(Y, np.float32)
        m2 = np.zeros_like(Y, np.float32)
        nan_places = np.vstack(np.where(np.isnan(Y))).T
        for j,k in nan_places:
            Y[j,k] = 0.0
            m1[j,k] = 0.0
            m2[j,k] = 1.0

        S = X[:,0]
        D = X[:,1:]
        self.desc_dim = D.shape[1]
        
        if W is None: # no weights
            w = np.ones_like(Y)
        elif len(W.shape) == 1: # class weights
            w = np.tile(W, (len(S), 1))
        else:
            w = W
   

        ################################################################

        def shuffle_set (dataset):
            idx = np.arange(len(dataset[0]))
            np.random.shuffle(idx)
            return [obj[idx] for obj in dataset]
        
        def getTestIds(len,early):
            test_idx = np.full(len, False, dtype=bool)
            test_idx[int((1-early)*len):] = True
            return test_idx

        def subset (dataset, idx):
            return [obj[idx] for obj in dataset]
        
        self.multiprocessing = False
                    
        if self.augment == 0:
            self.fill_dictionary(S)
            #S,D,Y,m1,m2,w = shuffle_set([S,D,Y,m1,m2,w])
            #test_idx = np.random.binomial(1,self.early,len(S)).astype("bool")
            test_idx  = getTestIds(len(S),self.early)
            X = self.embed_smiles(S)
            self.trainset = shuffle_set(subset( [S,X,D,Y,m1,m2,w], ~test_idx ))
            testset  = subset( [S,X,D,Y,m1,m2,w], test_idx )

        #elif self.augment < 0: # simulate pre-augmentation
        #    S,D,Y,m1,m2,w = shuffle_set([S,D,Y,m1,m2,w])
        #    test_idx = np.random.binomial(1,0.2,len(S)).astype("bool")
        #    SA_test, idxtestnew = aug_smiles(S[test_idx], nmax=-self.augment)
        #
        #    if self.allow_overfit:
        #        SA, idx = aug_smiles(S, nmax=-self.augment, exclude=SA_test)
        #        self.fill_dictionary(np.append(SA,SA_test))
        #        self.trainset = [SA, self.embed_smiles(SA)] + subset([D,Y,m1,m2,w], idx)
        #        self.trainset = shuffle_set(self.trainset)
        #    else:
        #        SA, idx = aug_smiles(S[~test_idx], nmax=-self.augment)
        #        self.fill_dictionary(np.append(SA,SA_test))
        #        self.trainset = [SA, self.embed_smiles(SA)] + subset(subset([D,Y,m1,m2,w], ~test_idx), idx)
        #        self.trainset = shuffle_set(self.trainset)
        #        
        #    testset  = [SA_test, self.embed_smiles(SA_test)] + subset( subset([D,Y,m1,m2,w], test_idx), idxtestnew)

        else: # >= 0
            #S,D,Y,m1,m2,w = shuffle_set([S,D,Y,m1,m2,w])
            #test_idx = np.random.binomial(1,self.early,len(S)).astype("bool")
            test_idx  = getTestIds(len(S),self.early)
            SA_test, idxtestnew = aug_smiles(S[test_idx], nmax=10) # increasing the size of augmentation set for validation to have a smoother estimation of errors
#            SA_test, idxtestnew = aug_smiles(S[test_idx], nmax=self.augment)

            self.multiprocessing = True

            ################################################################

            if self.multiprocessing:
                Q = Queue()
                
                def select_smiles (smiles, nmax, mols):
                    res = []
                    idx = []
                    for i,sm in enumerate(smiles):
                        augs = select_one_smiles(sm, nmax, mols)
                        res += augs
                        idx += [i]*len(augs)

                    return np.array(res), idx
                
                def select_one_smiles(sm, nmax, mols):
                    if len(mols[sm]) <= nmax: return mols[sm]
                    x = []
                    for i in range(2*nmax):
                        smile = mols[sm][np.random.randint(0, len(mols[sm]))]
                        if smile not in x: x.append(smile)
                        if len(x) == nmax: break
                    return x
                
                def get_ready_smiles():
                    while  Q.empty():
                        print ("MESSAGE: generating augmented SMILES")
                        sleep (1)
 
                    return pickle.loads(Q.get())
                
                def get_trainset_parallel ():
                    SA, idx = get_ready_smiles()
           
                    if self.allow_overfit:
                        self.trainset = [SA, self.embed_smiles(SA)] + subset([D,Y,m1,m2,w], idx)
                    else:
                        self.trainset = [SA, self.embed_smiles(SA)] + subset(subset([D,Y,m1,m2,w], ~test_idx), idx)
                   
                    self.trainset = shuffle_set(self.trainset)
                    return self.trainset
                
                def smile_continuous_generate(Q, ds, nmax, exclude=[]):
                    setsmi = set(ds)
                    exclude = set(exclude)
                    print("MESSAGE: start SMILES generation nmax =", nmax," records:",len(ds)," molecules:",len(setsmi)," excluded:",len(exclude))
                    mols = {}
                    for sm in setsmi:
                        mols.update({ sm: []} )
                    i = 0
                    j = 0
                    second = False
                    while True:
                        for sm in setsmi:
                            augs = aug_smile(sm, nmax*10, exclude=exclude)
                            if second and Q.empty():
                                Q.put(pickle.dumps(select_smiles(smiles = ds, nmax = nmax, mols=mols)))
                                print ("SMILES: ",i," were generated.")

                            for smile in augs:
                                i += 1
                                if math.log10(i) > j:
                                    j += 1
                                    print ("MESSAGE: SMILES: ",i," were generated.")
                                if smile not in mols[sm]:
                                    if len(mols[sm]) < nmax: mols[sm].append(smile)
                                    else:
                                        mols[sm][np.random.randint(0, len(mols[sm])) ] = smile

                        if not second:
                            second = True
                            Q.put(pickle.dumps(select_smiles(smiles = ds, nmax = nmax*10, mols=mols)))  # first time all SMILES are provided
                            print ("MESSAGE: SMILES: ",i," were generated.")
                         
                if self.allow_overfit:
                    P = Process(target = smile_continuous_generate, args=(Q, S, self.augment,SA_test, ), daemon=True)
                else:
                    P = Process(target = smile_continuous_generate, args=(Q, S[~test_idx], self.augment,SA_test, ), daemon=True) # some molecules can be duplicated -- we prevent for them to use same SMILES
    
                P.start()
                self.fill_dictionary(np.append(SA_test,get_ready_smiles()[0]))
            else:
                self.fill_dictionary(np.append(SA_test,aug_smiles(S,nmax=10)[0]))
             
            testset  = [SA_test, self.embed_smiles(SA_test)] + subset( subset([D,Y,m1,m2,w], test_idx), idxtestnew) 
            ################################################################

        def get_trainset ():
            if self.augment > 0:
                if self.allow_overfit:
                    SA, idx = aug_smiles(S, nmax=self.augment, exclude=SA_test)
                    self.trainset = [SA, self.embed_smiles(SA)] + subset([D,Y,m1,m2,w], idx)
                else:
                    SA, idx = aug_smiles(S[~test_idx], nmax=self.augment)
                    self.trainset = [SA, self.embed_smiles(SA)] + subset(subset([D,Y,m1,m2,w], ~test_idx), idx)
            
            self.trainset = shuffle_set(self.trainset) # always shuffling -- should not be worse!
            return self.trainset
        
        def get_testset ():
            return testset

        try:
            tf.set_random_seed(self.seed)
            with tf.variable_scope(self.name):
              with tf.device(device_name):
                self.model = self._build_model()
                self.sess = tf.Session(config=sess_config)
                self.sess.run(tf.global_variables_initializer())
                self._train(get_trainset_parallel if self.multiprocessing else get_trainset, get_testset, device_str)
                del self.trainset
                self.sess.close()
                self.sess = None
                self.model = None
            if self.multiprocessing:
                P.terminate()
            tf.reset_default_graph()
        except :
            if self.multiprocessing:
                P.terminate()
 
    def predict (self, X):
        S = X[:,0]
        D = X[:,1:]
        if self.augment != 0:
            S, idx = aug_smiles(S, nmax=10)
            D = D[idx]
        
        tf.set_random_seed(self.seed)
        np.random.seed(self.seed)
        
        result = []
        device_name, device_str = check_device(self.gpu)
        with tf.device(device_name):
            tf.reset_default_graph()
            
            sess = tf.Session(config=sess_config)
            saver = tf.train.import_meta_graph(self.model_path+".meta")
            saver.restore(sess, self.model_path)
            graph = tf.get_default_graph()
            
            train_mode = graph.get_tensor_by_name(self.name+"/train_mode:0")
            dropout_prob = graph.get_tensor_by_name(self.name+"/dropout_prob:0")
            input = graph.get_tensor_by_name(self.name+"/input:0")
            desc = graph.get_tensor_by_name(self.name+"/desc:0")
            prediction = graph.get_tensor_by_name(self.name+"/OUTPUT/BiasAdd:0")
            
            # apply
            X = self.embed_smiles(S)
            batch_size = self.batch_size
            i = 0
            
            while i < len(X):
                    if i+batch_size <= len(X):
                        Sb = np.array(S[i:i+batch_size])
                        Xb = X[i:i+batch_size]
                        Db = D[i:i+batch_size]
                    else:
                        Sb = np.array(S[i:])
                        Xb = X[i:]
                        Db = D[i:]
                    i += batch_size
                    result.append( sess.run([prediction], feed_dict={input: Xb, desc: Db, train_mode:False, dropout_prob:0})[0] )
        sess.close()
        tf.reset_default_graph()
        res = np.vstack(result)
        if self.augment != 0:
            res = scatter_mean(res, idx)
        return res
                

    def _build_model (self):
        train_mode = tf.placeholder(dtype=tf.bool, shape=[], name="train_mode")
        dropout_prob = tf.placeholder(tf.float32, name="dropout_prob")
    
        input = tf.placeholder(dtype=tf.int64, shape=[None,None], name="input")
        target = tf.placeholder(dtype=tf.float32,shape=[None, self.n_targets], name="target")

        nfp_function = {'0' : self.cnf0, '1' : self.cnf1, '2' : self.cnf2}[self.nfp_type]
        relu_function = {'crelu' : tf.nn.crelu, 'lrelu' : tf.nn.leaky_relu, 'relu' : tf.nn.relu, 'elu' : tf.nn.elu, 'swish' : tf.nn.swish}[self.relu_type]

        fp = nfp_function(train_mode, input, relu_function)
        desc = tf.placeholder(dtype=tf.float32, shape=[None, self.desc_dim], name="desc")
        fp = tf.concat([fp, desc], axis=1)
        
        for dim in self.mlp_dims:
            fp = tf.layers.dense(fp, dim, 
                                 activation=tf.nn.elu,
                                 kernel_initializer = tf.initializers.glorot_uniform()
                                 )
            fp = tf.layers.dropout(fp, rate=dropout_prob, training=train_mode, seed=self.seed)
        
        prediction = tf.layers.dense(fp, units=self.n_targets, activation=None, name="OUTPUT")
         
        mv1 = tf.placeholder(dtype=tf.float32, shape = [None, self.n_targets])
        mv2 = tf.placeholder(dtype=tf.float32, shape = [None, self.n_targets])
        target_mv = tf.add(tf.multiply(target, mv1), tf.multiply(prediction, mv2))    
        weights = tf.placeholder(dtype=tf.float32, shape=[None, self.n_targets])

        ################################################################

        if self.error_function == 'rmse':
            loss = tf.reduce_mean(tf.multiply(tf.squared_difference(target_mv, prediction), weights))
    
        elif self.error_function == 'taxonomy':
            exp_coef = self.exp_coef
            # regloss
            loss = tf.reduce_mean(tf.squared_difference(target_mv, prediction))
            # logloss
            loss += exp_coef*tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(labels=target_mv,logits=prediction))
            # graphloss
            reduce_coef = tf.constant( 1./(n_targets**2), np.float32)
            for edge in self.taxonomy_graph:
                i,j,w = edge
                loss += tf.multiply(tf.reduce_mean(tf.squared_difference(prediction[:,i], prediction[:,j])), w) * reduce_coef

        elif self.error_function == 'entropy':
            exp_coef = self.exp_coef
            # regloss
            loss = tf.reduce_mean(tf.squared_difference(target_mv, prediction))
            # logloss
            loss += exp_coef*tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(labels=target_mv,logits=prediction))

        elif self.error_function == 'multilabelnn':
            one = tf.constant(1, np.float32)
            about_one = tf.constant(1-0.0001, np.float32)
            about_zero = tf.constant(0.0001, np.float32) 
            
            loss = target_mv[:,0]-target_mv[:,0]
            for k in range(n_targets):
                for l in range(n_targets):
                    loss += ( target_mv[:,k]*(one-target_mv[:,l]) ) * tf.exp(prediction[:,l]-prediction[:,k])
                    loss /= tf.reduce_sum(about_zero+target_mv, axis=1)
            loss /= tf.reduce_sum(about_one-target_mv, axis=1)
            loss = tf.reduce_sum(loss)
        
        else:
            raise Exception("bad error function")
         
        ################################################################

        #optimizer = tf.train.AdamOptimizer(learning_rate=self.learning_rate)
        optimizer = RAdamOptimizer(learning_rate=self.learning_rate, beta1=0.9, beta2=0.999, weight_decay=0.0)
        update_ops = tf.get_collection(tf.GraphKeys.UPDATE_OPS)
        with tf.control_dependencies(update_ops):
            train_op_all = optimizer.minimize(
                loss=loss,
                global_step=tf.train.get_global_step())

        
        model = [prediction, input, desc, train_mode, dropout_prob, train_op_all, loss, target, mv1, mv2, weights]
        return model

     
    def _train (self, trainset_fn, testset_fn, device_str):
        [prediction, input, desc, train_mode, dropout_prob, train_op_all, loss, target, mv1, mv2, weights] = self.model
        start_moment = time.time()
        global epoch
        #dropout_inc = 0.1 if self.adjust_dropout and self.allow_overfit else 0
        dropout_inc = 0

        def run (dataset, train=True):
            S,X,D,Y,m1,m2,w = dataset
            
            if self.debug:
                dirname = os.path.split(os.path.split(self.model_path)[0])[0]+"/"
                fname = "debug_"+str(epoch)+"_"+ ("train" if train else "valid") + ".csv"
                np.savetxt(dirname+fname, np.hstack((S.reshape(-1,1),Y)), fmt='%s')
                cnt = Counter([can_smile(s) for s in S])
                np.savetxt(dirname+"counter_"+fname, np.array([['{0:<40}'.format(k+","), cnt[k]] for k in cnt.keys()]), fmt="%s\t%s")
            
            sum_loss = 0
            i = 0
            bs = self.batch_size
                                 
            while i < len(X):
                if i+bs <= len(X):
                    Xb,Yb,Db = X[i:i+bs], Y[i:i+bs], D[i:i+bs]
                    m1b, m2b = m1[i:i+bs], m2[i:i+bs]
                    wb = w[i:i+bs]
                else:
                    Xb,Yb,Db = X[i:], Y[i:], D[i:]
                    m1b, m2b = m1[i:], m2[i:]
                    wb = w[i:]
                i += bs
                if train:
                    feed_dict = { input : Xb,
                                  target : Yb,
                                  desc : Db,
                                  train_mode : True,
                                  dropout_prob : self.dropout+dropout_inc,
                                  mv1 : m1b, mv2 : m2b, weights: wb}
                    _,l = self.sess.run([train_op_all,loss], feed_dict=feed_dict)
                    sum_loss += l*len(Xb)
                else:
                    feed_dict = { input : Xb,
                                  target : Yb,
                                  desc: Db,
                                  train_mode : False,
                                  dropout_prob : 0.,
                                  mv1 : m1b, mv2 : m2b, weights: wb }
                    l = loss.eval(session=self.sess, feed_dict=feed_dict)
                    sum_loss += l*len(Xb)
            return np.sqrt(sum_loss/len(X))

        best_score = np.inf
        count = 0
        model_saved_flag = False
        early_tolerance = 0.1 # 10% to diagnose early stopping

        for epoch in range(self.n_epochs):
            if epoch > self.validation_interval and os.path.exists('stop'):
                os.remove('stop')
                break
            trainset = trainset_fn()
            run(trainset, train=True)
             
            if epoch % self.validation_interval == 0:
                val_score = run(testset_fn(), train=False)
                train_score = run(trainset, train=False)
                elapsed_time = time.time()-start_moment

                if self.adjust_dropout and epoch > 10:
                    self.dropout *= val_score/train_score
                    if self.dropout > 0.9:
                        self.dropout = 0.9

                basic_message = "MESSAGE: train score: {} / validation score: {} / at epoch: {} / {} / elapsed time: {}s / dropout: {} ".format(round(float(train_score), 3), round(float(val_score), 3), epoch, device_str, round(elapsed_time, 1), round(self.dropout, 3))
                    
                if  val_score < 0.99*best_score: # to prevent fluctuations
                    best_score = val_score
                    count = 0
                    print(basic_message+ " / left {}s / GOT IMPROVEMENT".format(round(elapsed_time/(epoch+1)*(self.n_epochs-epoch-1))))
                    saver = tf.train.Saver()
                    saver.save(self.sess, self.model_path)
                    model_saved_flag = True
                else:
                    epoch_left = round(self.n_epochs*early_tolerance) - count*self.validation_interval
                    count += 1
                    if epoch_left <= 0:
                        print(basic_message+" / EARLY STOPPED")
                        np.savetxt("earlystop", [epoch, self.n_epochs], fmt="%d")
                        return
                    else:
                        print(basic_message+" / left {}s - epochs {} to early stop ".format(round(elapsed_time*epoch_left/(epoch+1)),round(epoch_left)))
                        
        if not model_saved_flag:
            saver = tf.train.Saver()
            saver.save(self.sess, self.model_path)
    
    def highway_layer (self, x, size, activation, carry_bias=-1.0, N=10):
        W_T = tf.Variable(tf.truncated_normal([size, size], stddev=0.1), name="weight_transform")
        b_T = tf.Variable(tf.constant(carry_bias, shape=[size]), name="bias_transform")
 
        W = tf.Variable(tf.truncated_normal([size, size], stddev=0.1), name="weight")
        b = tf.Variable(tf.constant(0.1, shape=[size]), name="bias")

        for i in range(N):
            T = tf.sigmoid(tf.matmul(x, W_T) + b_T, name="transform_gate")
            H = activation(tf.matmul(x, W) + b, name="activation")
            C = tf.subtract(1.0, T, name="carry_gate")
            x = tf.add(tf.multiply(H, T), tf.multiply(x, C))
        return x

    def batch_vm (self, v,m):
        shape = tf.shape(v)
        A = tf.reshape(v, [-1, shape[-1]])
        C = tf.matmul(A, m)
        return tf.reshape(C, [shape[0], shape[1], tf.shape(m)[1]])

    def normalize_fn (self, output, train_mode):
        if self.normalize=='layer':
            output = tf_contrib_layer_norm(output)
        elif self.normalize=='batch':
            output = tf.layers.batch_normalization(output, training=train_mode)
        return output
    
    def cnf0 (self, train_mode, input, relu):
        filter_size = self.filter_size
        hiddim = self.fp_layer_dim
        if relu == tf.nn.crelu:
            hiddim2 = hiddim*2
        else:
            hiddim2 = hiddim

        if self.preload_dictionary:
            extdim = len(self.tokens)
            layers = tf.one_hot(input, extdim)
        else:
            extdim = 150 
            
            WE = tf.Variable( tf.random_uniform([len(self.tokens), extdim], -1.0, 1.0, seed=self.seed))
            mask_padding_zero_op = tf.scatter_update(WE, 
                                                     self.tokens.get("pad"), 
                                                     tf.zeros([extdim,]))
            with tf.control_dependencies([mask_padding_zero_op]):
                layers = tf.nn.embedding_lookup(WE, input)
        

        out = []
        
        if filter_size == -1:
            filter_size = 2
            
        for l in range(self.n_fp_layers):
            filters = tf.Variable(
                tf.truncated_normal(
                    dtype=tf.float32,
                    shape=[filter_size, extdim, extdim],
                    stddev=0.1,
                    seed=self.seed))
            H = tf.Variable(
                tf.truncated_normal(
                    shape=[extdim, hiddim],
                    seed=self.seed))
        
            layers = tf.nn.conv1d(
                layers, filters,
                stride=1,
                padding='SAME',
                data_format='NHWC')            
            activation = relu(self.batch_vm(layers, H))
            output = tf.reduce_sum(activation,axis=1)
            output.set_shape([None,hiddim2])
            output = self.normalize_fn(output, train_mode)
            out.append(output)

        fp = tf.concat(out, axis=1)
        fp.set_shape([None, hiddim2*self.n_fp_layers])
        fp = self.normalize_fn(fp, train_mode) #error?
                
        if self.highway>0:
            fp = self.highway_layer(fp, hiddim2*self.n_fp_layers, relu, N=self.highway)
        return fp

    def cnf1 (self, train_mode, input, relu):
        filter_size = self.filter_size
        hiddim = self.fp_layer_dim
        if relu == tf.nn.crelu:
            hiddim2 = hiddim*2
        else:
            hiddim2 = hiddim

        
        if self.preload_dictionary:
            extdim = len(self.tokens)
            input = tf.one_hot(input, extdim)
        else:
            extdim = 150
            WE = tf.Variable( tf.random_uniform([len(self.tokens), extdim], -1.0, 1.0, seed=self.seed))
            mask_padding_zero_op = tf.scatter_update(WE, 
                                                     self.tokens.get("pad"), 
                                                     tf.zeros([extdim,]))
            with tf.control_dependencies([mask_padding_zero_op]):
                input = tf.nn.embedding_lookup(WE, input)
            
        out = []
        for l in range(self.n_fp_layers):
            if l > 0:
                if self.filter_size == -1:
                    filter_size = l+1
                else:
                    filter_size = self.filter_size 
                filters = tf.Variable(
                    tf.truncated_normal(
                        dtype=tf.float32,
                        shape=[filter_size, extdim, extdim],
                        stddev=0.1,
                        seed=self.seed))
        
                layers = tf.nn.conv1d(
                    input, filters,
                    stride=1,
                    padding='SAME',
                    data_format='NHWC')
            else:
                layers = input
      
            H = tf.Variable(
                tf.truncated_normal(
                    shape=[extdim, hiddim],
                seed=self.seed))
            activation = relu(self.batch_vm(layers, H))
            output = tf.reduce_sum(activation,axis=1)
            output.set_shape([None,hiddim2])
            output = self.normalize_fn(output, train_mode)
            out.append(output)
            
        fp = tf.concat(out, axis=1)
        fp.set_shape([None, hiddim2*self.n_fp_layers])
        fp = self.normalize_fn(fp, train_mode)
        if self.highway>0:
            fp = self.highway_layer(fp, hiddim2*self.n_fp_layers, relu, N=self.highway)
        return fp

    def cnf2 (self, train_mode, input, relu):
        filter_size = self.filter_size
        hiddim = self.fp_layer_dim
        if relu == tf.nn.crelu:
            hiddim2 = hiddim*2
        else:
            hiddim2 = hiddim

        
        if self.preload_dictionary:
            extdim = len(self.tokens)
            input = tf.one_hot(input, extdim)
        else:
            extdim = 150
            WE = tf.Variable( tf.random_uniform([len(self.tokens), extdim], -1.0, 1.0, seed=self.seed))
            mask_padding_zero_op = tf.scatter_update(WE, 
                                                     self.tokens.get("pad"), 
                                                     tf.zeros([extdim,]))
            with tf.control_dependencies([mask_padding_zero_op]):
                input = tf.nn.embedding_lookup(WE, input)
        
        out = []
        
        for l in range(self.n_fp_layers):
            if l > 0:
                if self.filter_size == -1:
                    filter_size = l+1
                else:
                    filter_size = self.filter_size 
                filters = tf.Variable(
                    tf.truncated_normal(
                        dtype=tf.float32,
                        shape=[filter_size, extdim, extdim],
                        stddev=0.1,
                    seed=self.seed))
        
                layers = tf.nn.conv1d(
                    layers*self.nfp_alpha+input*(1-self.nfp_alpha), filters,
                    stride=1,
                    padding='SAME',
                    data_format='NHWC')
            else:
                layers = input
      
            H = tf.Variable(
                tf.truncated_normal(
                    shape=[extdim, hiddim],
                    seed=self.seed))
            activation = relu(self.batch_vm(layers, H))
            output = tf.reduce_sum(activation,axis=1)
            output.set_shape([None,hiddim2])
            output = self.normalize_fn(output, train_mode)
            out.append(output)
            
        fp = tf.concat(out, axis=1)
        fp.set_shape([None, hiddim2*self.n_fp_layers])
        fp = self.normalize_fn(fp, train_mode)
        if self.highway > 0:
            fp = self.highway_layer(fp, hiddim2*self.n_fp_layers, relu, N=self.highway)
        return fp
