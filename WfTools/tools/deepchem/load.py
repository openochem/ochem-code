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
    
def read_task (path):
    config = ConfigParser.ConfigParser()
    with open(path, 'r') as fh:
        config.readfp(fh)
    task = {}

    task['train_mode'] = check_option(config, 'Task', 'train_mode', type=bool, default=False)
    task['model_file'] = check_option(config, 'Task', 'model_file')
    task['batch_size'] = check_option(config, 'Details', 'batch_size', type=int, default=32)
    task['gpu'] = check_option(config, 'Details', 'gpu', type=int, default=-1)
    os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"
    if task['gpu'] != -1:
        os.environ["CUDA_VISIBLE_DEVICES"]=str(task['gpu'])
    else:
        os.environ["CUDA_VISIBLE_DEVICES"]=""

    task['seed'] = check_option(config, 'Details', 'seed', type=int, default=42)
    
    if task['train_mode']:
        task['train_data_file'] = check_option(config, 'Task', 'train_data_file')
        task['apply_data_file'] = check_option(config, 'Task', 'apply_data_file', default="nonexistent")
        task['result_file'] = check_option(config, 'Task', 'result_file', default="nonexistent")

        task['task_names'] = check_option(config, 'Task', 'task_names')
        task['task_names'] = list(map(lambda x: x.strip(), task['task_names'].split(',')))
        
        task['class_weights'] = check_option(config, 'Task', 'class_weights', default=False)
        if task['class_weights']:
            task['class_weights'] = read_weights(task['class_weights'])

        task['regression'] = check_option(config, 'Task', 'regression', type=bool, default=True)

        task['n_epochs'] = check_option(config, 'Details', 'n_epochs', type=int, default=100)
        task['learning_rate'] = check_option(config, 'Details', 'learning_rate', type=float, default=0.001)
        task['dropout'] = check_option(config, 'Details', 'dropout', type=float, default=0.3)
        task['train_proportion'] = check_option(config, 'Details', 'train_proportion', type=float, default=0.9)
        task['n_best_nets'] = check_option(config, 'Details', 'n_best_nets', type=int, default=3)
        task['deepchem_model'] = check_option(config, 'Details', 'deepchem_model', default='GraphConv')

        if task['deepchem_model'] == 'GraphConv':
            task['graph_conv_layers'] = check_option(config, 'Details', 'graph_conv_layers', default="64,64")
            task['graph_conv_layers'] = list(map(int, task['graph_conv_layers'].split(",")))
            task['dense_layer_size'] = check_option(config, 'Details', 'dense_layer_size', type=int, default=128)
            task['number_atom_features'] = check_option(config, 'Details', 'number_atom_features', type=int, default=75)

        elif task['deepchem_model'] == 'DAG':
            task['n_graph_feat'] = check_option(config, 'Details', 'n_graph_feat_dag', type=int, default=30)
            task['n_outputs'] = check_option(config, 'Details', 'n_outputs', type=int, default=30)
            task['layer_sizes'] = check_option(config, 'Details', 'layer_sizes', default="100")
            task['layer_sizes'] = list(map(int, task['layer_sizes'].split(",")))
            task['layer_sizes_gather'] = check_option(config, 'Details', 'layer_sizes_gather', default="100")
            task['layer_sizes_gather'] = list(map(int, task['layer_sizes_gather'].split(",")))

        elif task['deepchem_model'] == 'Weave':
            task['n_graph_feat'] = check_option(config, 'Details', 'n_graph_feat_weave', type=int, default=128)
            task['n_weave'] = check_option(config, 'Details', 'n_weave', type=int, default=2)
            task['n_hidden'] = check_option(config, 'Details', 'n_hidden', type=int, default=50)

        elif task['deepchem_model'] == 'MPNN':
            pass

        elif task['deepchem_model'] == 'TextCNN':
            task['n_embeddings'] = check_option(config, 'Details', 'n_embeddings', type=int, default=75)
            
            task['num_filters'] = check_option(config, 'Details', 'num_filters', default="100,200,200,200,200,100,100,100,100,100,160,160")
            task['kernel_sizes'] = check_option(config, 'Details', 'kernel_sizes', default="1,2,3,4,5,6,7,8,9,10,15,20")
            task['num_filters'] = list(map(int, task['num_filters'].split(",")))
            task['kernel_sizes'] = list(map(int, task['kernel_sizes'].split(",")))
            
        elif task['deepchem_model'] == 'ChemCeption':
            raise Exception('ChemCeption is broken')
    else:
        task['task_names'] = check_option(config, 'Task', 'task_names')
        task['task_names'] = list(map(lambda x: x.strip(), task['task_names'].split(',')))
        task['server'] = check_option(config, 'Task', 'server', type=bool, default=False)

        task['apply_data_file'] = check_option(config, 'Task', 'apply_data_file')
        task['result_file'] = check_option(config, 'Task', 'result_file')

    print ("\nNot used options in section Task:")
    print (config.options("Task"))
    print ("\nNot used options in section Details:")
    print (config.options("Details"),"\n")
        
    return task

