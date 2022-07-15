# params: train_data_file, apply_data_file, target_name, result_file, gpu, train_mode


try:
    # python2
    import ConfigParser
except ImportError:
    # python3
    import configparser as ConfigParser

def check_option (config, section, name, default=None, type=str):
        if not config.has_section(section):
            raise Exception('No such config section: '+section) 
        
        if type==str:
            accessor = config.get
        elif type==int:
            accessor = config.getint
        elif type==float:
            accessor = config.getfloat
        elif type==bool:
            accessor = config.getboolean
        else:
            raise Exception('Internal exception #389423')
        
        if not config.has_option(section, name):
            if default is None:
                raise Exception('Mandatory config option '+name+' not found')
            else:
                print(f'Using default for parameter "{name}" ({section})\n')
                return default
        else:
            try:
                result = accessor(section, name)
                config.remove_option(section, name)
                return result
            except:
                raise Exception('Bad value for config option: '+name+', expected type: '+str(type))

import re
import sys
import os
import numpy as np
import csv

def read_weights (file):
    try:
        with open(file, 'r') as fh:
            w = fh.readline()
            w = re.sub(r"\s+", "", w)
            w = w.split(',')
            if w[-1] == "":
                w = w[:-1]
            w = list(map(float, w))
            w = (np.array(w)/np.sum(w)*len(w)).tolist()
            return w
    except:
        print("bad weights file")
        sys.exit(0)

###

def read_graph (file, targets):
    edges = []
    with open(file) as fh:
        reader = csv.reader(fh, delimiter=',')
        for row in reader:
            i,j = row[:2]
            if len(row) > 2:
                w_ij = row[2]
            else:
                w_ij = 1
            if i in targets and j in targets:
                edges.append([targets.index(i), targets.index(j), w_ij])
    if len(edges) > 0:
        return edges
    else:
        return False
    
def read_task (path):
    config = ConfigParser.ConfigParser()
    with open(path, 'r') as fh:
        config.readfp(fh)
    task = {}

    task['train_mode'] = check_option(config, 'Task', 'train_mode', type=bool, default=False)
    task['model_file'] = check_option(config, 'Task', 'model_file')
    task['batch_size'] = check_option(config, 'Details', 'batch_size', type=int, default=16)
    task['gpu'] = check_option(config, 'Details', 'gpu', type=int, default=-1)
 
    os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"
    if task['gpu'] != -1:
        os.environ["CUDA_VISIBLE_DEVICES"]=str(task['gpu'])
    else:
        os.environ["CUDA_VISIBLE_DEVICES"]=""
    
    if task['train_mode']:
        task['seed'] = check_option(config, 'Details', 'seed', type=int, default=0)
        task['learning_rate'] = check_option(config, 'Details', 'learning_rate', type=float, default=0.001)
        task['train_data_file'] = check_option(config, 'Task', 'train_data_file')
        task['apply_data_file'] = check_option(config, 'Task', 'apply_data_file', default="nonexistent")
        task['result_file'] = check_option(config, 'Task', 'result_file', default="nonexistent")

        task['weights'] = check_option(config, 'Task', 'weights', default=False)
        if task['weights']:
            task['weights'] = read_weights(task['weights'])

        task['error_function'] = check_option(config, 'Details', 'error_function', default='rmse')
        task['n_epochs'] = check_option(config, 'Details', 'n_epochs', type=int, default=100)
        task['fp_layer_dim'] = check_option(config, 'Details', 'fp_layer_dim', type=int, default=32)
        task['n_fp_layers'] = check_option(config, 'Details', 'n_fp_layers', type=int, default=3)
        task['normalize'] = check_option(config, 'Details', 'normalize', default="layer")
        task['nfp_type'] = check_option(config, 'Details', 'nfp_type', default='2')
        task['relu_type'] = check_option(config, 'Details', 'relu_type', default='relu')
        task['dropout'] = check_option(config, 'Details', 'dropout', type=float, default=0.3)
        task['tokenizer_type'] = check_option(config, 'Details', 'tokenizer_type', default='char')
        task['early'] = check_option(config, 'Details', 'early', type=float, default= 0.2)

        task['task_names'] = check_option(config, 'Task', 'task_names')
        task['task_names'] = list(map(lambda x: x.strip(), task['task_names'].split(',')))

        task['mlp_dims'] = check_option(config, 'Details', 'mlp_dims', default="128,64")
        task['mlp_dims'] = map(int, task['mlp_dims'].split(","))

        task['augment'] = check_option(config, 'Details', 'augment', type=int, default=0)
        task['augment'] = -task['augment'] # current values are just 0 and -1 -- on-line augmentation

        #almost not used, obsolete
        task['debug'] = check_option(config, 'Details', 'debug', type=bool, default=False)
        task['exp_coef'] = check_option(config, 'Details', 'exp_coef', type=float, default=0.1)
        task['graph'] = check_option(config, 'Task', 'graph_file', default=False)
        if task['graph']:
            task['graph'] = read_graph(task['graph'], task['task_names'])

        #not used options
        #task['highway'] = check_option(config, 'Details', 'highway', type=int, default=0)
        #task['allow_overfit'] = check_option(config, 'Details', 'allow_overfit', type=bool, default=False)
        #task['canonize'] = check_option(config, 'Details', 'canonize', type=bool, default=False)
        #task['filter_size'] = check_option(config, 'Details', 'filter_size', type=int, default=-1) # will use default settings, 
        #task['nfp_alpha'] = check_option(config, 'Details', 'nfp_alpha', type=float, default=0.5)
        
        task['nfp_alpha'] = 0.5
        task['filter_size'] = -1
        task['highway'] = 0
        task['allow_overfit'] = False
        task['canonize'] = False
    else:
        task['apply_data_file'] = check_option(config, 'Task', 'apply_data_file')
        task['result_file'] = check_option(config, 'Task', 'result_file')

    print ("\nNot used options in section Task:")
    print (config.options("Task"))
    print ("\nNot used options in section Details:")
    print (config.options("Details"),"\n")

    return task
