import tensorflow as tf
import keras.backend as K
import numpy as np
import time
import os
import argparse
import json
from copy import deepcopy
from tensorflow_addons import optimizers
#from sklearn.model_selection import KFold
import kgcnn.training.schedule
import kgcnn.training.scheduler
import kgcnn.training.callbacks
from kgcnn.data.utils import save_json_file
from kgcnn.utils.models import get_model_class
from kgcnn.data.moleculenet import MoleculeNetDataset
from kgcnn.mol.encoder import OneHotEncoder
from kgcnn.hyper.hyper import HyperParameter

from tensorflow.keras.callbacks import EarlyStopping, ReduceLROnPlateau, CSVLogger
import configparser
import sys
import tarfile
from distutils.util import strtobool
# Standard imports
import json
import os
import pandas as pd
import numpy as np
import pickle
import tensorflow_addons as tfa

# Define the config
conf_name = sys.argv[1]
config = configparser.ConfigParser();
config.read(conf_name);

print("Load config file: ", conf_name);


def getConfig(section, attribute, default=""):
    try:
        return config[section][attribute];
    except:
        return default;


def SmilesOK(smi):
    try:
        mol = MolToSmiles(MolFromSmiles(smi))
        return True
    except:
        return False


# Categorical multiclass including Binary case "2"
def CCEmask(y_true, y_pred):
    y_true = tf.cast(y_true,
                     dtype=tf.float32)  # case of Classifiers # we may don't need it if we have the split between Reg/Class
    masked = tf.where(tf.math.is_nan(y_true), 0., 1.)
    return K.categorical_crossentropy(y_true * masked, y_pred * masked)


def BCEmask(y_true, y_pred):
    y_true = tf.cast(y_true,
                     dtype=tf.float32)  # case of Classifiers # we may don't need it if we have the split between Reg/Class
    masked = tf.where(tf.math.is_nan(y_true), 0., 1.)
    return K.binary_crossentropy(y_true * masked, y_pred * masked)


## multi sparse task
def RMSEmask(y_true, y_pred):
    # Compute the square error, and subsequently masked
    y_true = tf.cast(y_true, dtype=tf.float32)  # case of Classifiers
    masked = tf.where(tf.math.is_nan(y_true), 0., 1.)
    y_true_ = tf.where(tf.math.is_nan(y_true), 0., y_true)
    y_pred_ = tf.where(tf.math.is_nan(y_true), 0., y_pred)

    err = (y_true_ - y_pred_) * (y_true_ - y_pred_)

    # Select the qualifying values and mask
    sumsq = K.sum(masked * err, keepdims=True)
    num = tf.cast(K.sum(masked, keepdims=True), dtype=tf.float32)
    return K.sqrt(sumsq / (num + K.epsilon()))


def isClassifer(graph_labels_all):
    t = np.unique(graph_labels_all)
    if len(t) == 2 and np.all(t == [0, 1]):
        return True
    else:
        return False


target_list = ['target']
# Define standard parameters from config.cfg
TRAIN = getConfig("Task", "train_mode");
MODEL_FILE = getConfig("Task", "model_file");
TRAIN_FILE = getConfig("Task", "train_data_file");
APPLY_FILE = getConfig("Task", "apply_data_file", "train.csv");
RESULT_FILE = getConfig("Task", "result_file", "results.csv");
# Define details options from config.cfg
nbepochs = int(getConfig("Details", "nbepochs", 100));
seed = int(getConfig("Details", "seed", 101));
batch_size = int(getConfig("Details", "batch", 128));
earlystop_patience = int(getConfig("Details", "earlystop_patience", 30));
redlr_patience = int(getConfig("Details", "redlr_patience", 10));
redlr_factor = float(getConfig("Details", "redlr_factor", 0.8));
val_split = float(getConfig("Details", "reduce_lr_factor", 0.2));
gpu = int(getConfig("Details", "gpu", 0));

# Define network parameters
nn_embedding = int(getConfig("Details", "embedding_size", 32));
nn_lstmunits = int(getConfig("Details", "lstm_units", 16));
nn_denseunits = int(getConfig("Details", "dense_units", 16));
nn_return_proba = int(getConfig("Details", "return_proba", 0));
nn_loss = str(getConfig("Details", "nn_loss", "RMSE"));
nn_lr_start = float(getConfig("Details", "nn_lr_start", 1e-2));
output_dim = int(getConfig("Details", "output_dim", 1));

# model selection parameters for selecting which architecture to run
architecture_name = getConfig("Details", "architecture_name", 'AttFP');

print("Architecture selected:", architecture_name)

if architecture_name == "HamNet":
    hyper = { "model": {
            "class_name": "make_model",
            "module_name": "kgcnn.literature.HamNet",
            "config": {
                "name": "HamNet",
                "inputs": [
                    {"shape": [None, 41], "name": "node_attributes", "dtype": "float32", "ragged": True},
                    {"shape": [None, 11], "name": "edge_attributes", "dtype": "float32", "ragged": True},
                    {"shape": [None, 2], "name": "edge_indices", "dtype": "int64", "ragged": True},
                    {"shape": [None, 3], "name": "node_coordinates", "dtype": "float32", "ragged": True}
                ],
                "input_embedding": {"node": {"input_dim": 95, "output_dim": 64},
                                    "edge": {"input_dim": 5, "output_dim": 64}},
                "message_kwargs": {"units": 200, 
                                    "units_edge": 200,
                                    "rate": 0.2, "use_dropout": True},
                "fingerprint_kwargs": {"units": 200, "units_attend": 200, 
                                    "rate": 0.5, "use_dropout": True,
                                     "depth": 3},
                "gru_kwargs": {"units": 200},
                "verbose": 10, "depth": 3, 
                "union_type_node": "gru",
                "union_type_edge": "None",
                "given_coordinates": True,
                'output_embedding': 'graph',
                'output_mlp': {"use_bias": [True, False], "units": [200, 1],
                               "activation": ['relu', 'linear'],
                               "use_dropout": [True,  False],
                               "rate": [0.5, 0.0]}
            }
        },
        "training": {
            "fit": {
                "batch_size": 16, "epochs": 800, "validation_freq": 1, "verbose": 2,
                "callbacks": []
            },
            "compile": {
                "optimizer": {"class_name": "Addons>AdamW", "config": {"lr": 0.001,
                                                                       "weight_decay": 1e-05}},
                "loss": "mean_squared_error"
            },
            "cross_validation": {"class_name": "KFold",
                                 "config": {"n_splits": 5, "random_state": None, "shuffle": True}},
            "scaler": {"class_name": "StandardScaler",
                       "config": {"with_std": True, "with_mean": True, "copy": True}}
        },
        "data": {
            "dataset": {"class_name": "MoleculeNetDataset",
                        "config": {},
                        "methods": [{"set_attributes": {}}]
            },
            "data_unit": "mol/L"
        },
        "info": {
            "postfix": "",
            "kgcnn_version": "2.0.2"
        }
    }

if architecture_name == "PAiNN":
    hyper = {
        "model": {
            "class_name": "make_model",
            "module_name": "kgcnn.literature.PAiNN",
            "config": {
                "name": "PAiNN",
                "inputs": [
                    {"shape": [None], "name": "node_number", "dtype": "float32", "ragged": True},
                    {"shape": [None, 3], "name": "node_coordinates", "dtype": "float32", "ragged": True},
                    {"shape": [None, 2], "name": "range_indices", "dtype": "int64", "ragged": True}
                ],
                "input_embedding": {"node": {"input_dim": 95, "output_dim": 128}},
                "bessel_basis": {"num_radial": 20, "cutoff": 5.0, "envelope_exponent": 5},
                "pooling_args": {"pooling_method": "sum"}, "conv_args": {"units": 128, "cutoff": None},
                "update_args": {"units": 128}, "depth": 3, "verbose": 10,
                "output_embedding": "graph",
                "output_mlp": {"use_bias": [True, True], "units": [128, 1], "activation": ["swish", "linear"]}
            }
        },
        "training": {
            "fit": {
                "batch_size": 32, "epochs": 250, "validation_freq": 10, "verbose": 2,
                "callbacks": []
            },
            "compile": {
                "optimizer": {
                    "class_name": "Addons>MovingAverage", "config": {
                        "optimizer": {
                            "class_name": "Adam", "config": {
                                "learning_rate": {
                                    "class_name": "kgcnn>LinearWarmupExponentialDecay", "config": {
                                        "learning_rate": 0.001, "warmup_steps": 30.0, "decay_steps": 40000.0,
                                        "decay_rate": 0.01
                                    }
                                }, "amsgrad": True
                            }
                        },
                        "average_decay": 0.999
                    }
                },
                "loss": "mean_absolute_error",
            },
            "cross_validation": {"class_name": "KFold",
                                 "config": {"n_splits": 5, "random_state": None, "shuffle": True}},
            "scaler": {"class_name": "StandardScaler", "config": {"with_std": True, "with_mean": True, "copy": True}}
        },
        "data": {
            "dataset": {
                "class_name": "MoleculeNetDataset",
                "config": {},
                "methods": [
                    {"map_list": {"method": "set_range", "max_distance": 3, "max_neighbours": 10000}}
                ]
            },
            "data_unit": "mol/L"
        },
        "info": {
            "postfix": "",
            "kgcnn_version": "2.0.3"
        }
    }

if architecture_name == 'GCN':
    hyper = {
        "model": {
            "class_name": "make_model",
            "module_name": "kgcnn.literature.GCN",
            "config": {
                "name": "GCN",
                "inputs": [
                    {"shape": [None, 41], "name": "node_attributes", "dtype": "float32", "ragged": True},
                    {"shape": [None, 1], "name": "edge_weights", "dtype": "float32", "ragged": True},
                    {"shape": [None, 2], "name": "edge_indices", "dtype": "int64", "ragged": True}],
                "input_embedding": {"node": {"input_dim": 95, "output_dim": 100},
                                    "edge": {"input_dim": 10, "output_dim": 100}},
                "gcn_args": {"units": 200, "use_bias": True, "activation": "relu"},
                "depth": 5, "verbose": 10,
                "output_embedding": "graph",
                "output_mlp": {"use_bias": [True, True, False], "units": [200, 100, 1],
                               "activation": ["kgcnn>leaky_relu", "kgcnn>leaky_relu", "linear"]}
            }
        },
        "training": {
            "fit": {
                "batch_size": 32,
                "epochs": 800,
                "validation_freq": 10,
                "verbose": 2,
                "callbacks": [
                    {
                        "class_name": "kgcnn>LinearLearningRateScheduler", "config": {
                        "learning_rate_start": 1e-03, "learning_rate_stop": 5e-05, "epo_min": 250, "epo": 800,
                        "verbose": 0
                    }
                    }
                ]
            },
            "compile": {
                "optimizer": {"class_name": "Adam", "config": {"lr": 1e-03}}
            },
            "cross_validation": {"class_name": "KFold",
                                 "config": {"n_splits": 5, "random_state": None, "shuffle": True}}
        },
        "data": {
            "dataset": {
                "class_name": "MoleculeNetDataset",
                "config": {},
                "methods": [
                    {"map_list": {"method": "normalize_edge_weights_sym"}}
                ]
            },
            "data_unit": "mol/L"
        },
        "info": {
            "postfix": "",
            "kgcnn_version": "2.0.3"
        }
    }

if architecture_name == 'GAT':
    hyper = {
        "model": {
            "class_name": "make_model",
            "module_name": "kgcnn.literature.GAT",
            "config": {
                "name": "GAT",
                "inputs": [
                    {"shape": [None, 41], "name": "node_attributes", "dtype": "float32", "ragged": True},
                    {"shape": [None, 11], "name": "edge_attributes", "dtype": "float32", "ragged": True},
                    {"shape": [None, 2], "name": "edge_indices", "dtype": "int64", "ragged": True}
                ],
                "input_embedding": {
                    "node": {"input_dim": 95, "output_dim": 100},
                    "edge": {"input_dim": 8, "output_dim": 100}},
                "attention_args": {"units": 100, "use_bias": True, "use_edge_features": True,
                                   "use_final_activation": False, "has_self_loops": True},
                "pooling_nodes_args": {"pooling_method": "sum"},
                "depth": 4, "attention_heads_num": 10,
                "attention_heads_concat": False, "verbose": 10,
                "output_embedding": "graph",
                "output_mlp": {"use_bias": [True, True, False], "units": [200, 100, 1],
                               "activation": ["kgcnn>leaky_relu", "selu", "linear"]}
            }
        },
        "training": {
            "fit": {
                "batch_size": 32, "epochs": 500, "validation_freq": 2, "verbose": 2,
                "callbacks": [
                    {"class_name": "kgcnn>LinearLearningRateScheduler", "config": {
                        "learning_rate_start": 0.5e-03, "learning_rate_stop": 1e-05, "epo_min": 250, "epo": 500,
                        "verbose": 0}
                     }
                ]
            },
            "compile": {
                "optimizer": {"class_name": "Adam", "config": {"lr": 5e-03}},
                "loss": "mean_absolute_error"
            },
            "cross_validation": {"class_name": "KFold",
                                 "config": {"n_splits": 5, "random_state": None, "shuffle": True}}
        },
        "data": {
            "dataset": {
                "class_name": "MoleculeNetDataset",
                "config": {},
                "methods": [
                ]
            },
            "data_unit": "mol/L"
        },
        "info": {
            "postfix": "",
            "kgcnn_version": "2.0.3"
        }
    }

if architecture_name == 'GATv2':
    hyper = {
        "model": {
            "class_name": "make_model",
            "module_name": "kgcnn.literature.GATv2",
            "config": {
                "name": "GATv2",
                "inputs": [
                    {"shape": [None, 41], "name": "node_attributes", "dtype": "float32", "ragged": True},
                    {"shape": [None, 11], "name": "edge_attributes", "dtype": "float32", "ragged": True},
                    {"shape": [None, 2], "name": "edge_indices", "dtype": "int64", "ragged": True}
                ],
                "input_embedding": {
                    "node": {"input_dim": 95, "output_dim": 100},
                    "edge": {"input_dim": 8, "output_dim": 100}},
                "attention_args": {"units": 100, "use_bias": True, "use_edge_features": True,
                                   "use_final_activation": False, "has_self_loops": True},
                "pooling_nodes_args": {"pooling_method": "sum"},
                "depth": 4, "attention_heads_num": 10,
                "attention_heads_concat": False, "verbose": 10,
                "output_embedding": "graph",
                "output_mlp": {"use_bias": [True, True, False], "units": [200, 100, 1],
                               "activation": ["kgcnn>leaky_relu", "selu", "linear"]}
            }
        },
        "training": {
            "fit": {
                "batch_size": 32, "epochs": 500, "validation_freq": 2, "verbose": 2,
                "callbacks": [
                    {"class_name": "kgcnn>LinearLearningRateScheduler", "config": {
                        "learning_rate_start": 0.5e-03, "learning_rate_stop": 1e-05, "epo_min": 250, "epo": 500,
                        "verbose": 0}
                     }
                ]
            },
            "compile": {
                "optimizer": {"class_name": "Adam", "config": {"lr": 5e-03}},
                "loss": "mean_absolute_error"
            },
            "cross_validation": {"class_name": "KFold",
                                 "config": {"n_splits": 5, "random_state": None, "shuffle": True}}
        },
        "data": {
            "dataset": {
                "class_name": "MoleculeNetDataset",
                "config": {},
                "methods": []
            },
            "data_unit": "mol/L"
        },
        "info": {
            "postfix": "",
            "kgcnn_version": "2.0.3"
        }
    }

if architecture_name == 'Schnet':
    hyper = {
        "model": {
            "class_name": "make_model",
            "module_name": "kgcnn.literature.Schnet",
            "config": {
                "name": "Schnet",
                "inputs": [
                    {"shape": [None], "name": "node_number", "dtype": "float32", "ragged": True},
                    {"shape": [None, 3], "name": "node_coordinates", "dtype": "float32", "ragged": True},
                    {"shape": [None, 2], "name": "range_indices", "dtype": "int64", "ragged": True}
                ],
                "input_embedding": {"node": {"input_dim": 95, "output_dim": 64}},
                "output_embedding": "graph",
                'output_mlp': {"use_bias": [True, True], "units": [64, 1],
                               "activation": ['kgcnn>shifted_softplus', "linear"]},
                'last_mlp': {"use_bias": [True, True], "units": [128, 64],
                             "activation": ['kgcnn>shifted_softplus', 'kgcnn>shifted_softplus']},
                "interaction_args": {
                    "units": 128, "use_bias": True, "activation": "kgcnn>shifted_softplus", "cfconv_pool": "sum"
                },
                "node_pooling_args": {"pooling_method": "sum"},
                "depth": 4,
                "gauss_args": {"bins": 20, "distance": 4, "offset": 0.0, "sigma": 0.4}, "verbose": 10
            }
        },
        "training": {
            "cross_validation": {"class_name": "KFold",
                                 "config": {"n_splits": 5, "random_state": None, "shuffle": True}},
            "scaler": {"class_name": "StandardScaler", "config": {"with_std": True, "with_mean": True, "copy": True}},
            "fit": {
                "batch_size": 32, "epochs": 800, "validation_freq": 10, "verbose": 2,
                "callbacks": [
                    {"class_name": "kgcnn>LinearLearningRateScheduler", "config": {
                        "learning_rate_start": 0.0005, "learning_rate_stop": 1e-05, "epo_min": 100, "epo": 800,
                        "verbose": 0}
                     }
                ]
            },
            "compile": {
                "optimizer": {"class_name": "Adam", "config": {"lr": 0.0005}},
                "loss": "mean_absolute_error"
            }
        },
        "data": {
            "dataset": {
                "class_name": "MoleculeNetDataset",
                "config": {},
                "methods": [
                    {"map_list": {"method": "set_range", "max_distance": 4, "max_neighbours": 10000}}
                ]
            },
        },
        "info": {
            "postfix": "",
            "kgcnn_version": "2.0.3"
        }
    }

if architecture_name == 'GraphSAGE':
    hyper = {
        "model": {
            "class_name": "make_model",
            "module_name": "kgcnn.literature.GraphSAGE",
            "config": {
                "name": "GraphSAGE",
                "inputs": [
                    {"shape": [None, 41], "name": "node_attributes", "dtype": "float32", "ragged": True},
                    {"shape": [None, 11], "name": "edge_attributes", "dtype": "float32", "ragged": True},
                    {"shape": [None, 2], "name": "edge_indices", "dtype": "int64", "ragged": True}],
                "input_embedding": {
                    "node": {"input_dim": 95, "output_dim": 64},
                    "edge": {"input_dim": 32, "output_dim": 32}},
                "node_mlp_args": {"units": [64, 32], "use_bias": True, "activation": ["relu", "linear"]},
                "edge_mlp_args": {"units": 64, "use_bias": True, "activation": "relu"},
                "pooling_args": {"pooling_method": "segment_mean"}, "gather_args": {},
                "concat_args": {"axis": -1},
                "use_edge_features": True,
                "pooling_nodes_args": {"pooling_method": "sum"},
                "depth": 3, "verbose": 10,
                "output_embedding": "graph",
                "output_mlp": {"use_bias": [True, True, False], "units": [64, 32, 1],
                               "activation": ["relu", "relu", "linear"]},
            }
        },
        "training": {
            "fit": {"batch_size": 32, "epochs": 500, "validation_freq": 10, "verbose": 2,
                    "callbacks": [{"class_name": "kgcnn>LinearLearningRateScheduler",
                                   "config": {"learning_rate_start": 0.5e-3, "learning_rate_stop": 1e-5,
                                              "epo_min": 400, "epo": 500, "verbose": 0}}]
                    },
            "compile": {"optimizer": {"class_name": "Adam", "config": {"lr": 5e-3}},
                        "loss": "mean_absolute_error"
                        },
            "cross_validation": {"class_name": "KFold",
                                 "config": {"n_splits": 5, "random_state": None, "shuffle": True}},
            "scaler": {"class_name": "StandardScaler", "config": {"with_std": True, "with_mean": True, "copy": True}},
        },
        "data": {
            "dataset": {
                "class_name": "MoleculeNetDataset",
                "config": {},
                "methods": [
                ]
            },
            "data_unit": "mol/L"
        },
        "info": {
            "postfix": "",
            "kgcnn_version": "2.0.3"
        }
    }

if architecture_name == 'AttFP':
    hyper = {
        "model": {
            "class_name": "make_model",
            "module_name": "kgcnn.literature.AttentiveFP",
            "config": {
                "name": "AttentiveFP",
                "inputs": [{"shape": [None, 41], "name": "node_attributes", "dtype": "float32", "ragged": True},
                           {"shape": [None, 11], "name": "edge_attributes", "dtype": "float32", "ragged": True},
                           {"shape": [None, 2], "name": "edge_indices", "dtype": "int64", "ragged": True}],
                "input_embedding": {"node": {"input_dim": 95, "output_dim": 100},
                                    "edge": {"input_dim": 5, "output_dim": 100}},
                "attention_args": {"units": 200},
                "depth": 2,
                "dropout": 0.2,
                "verbose": 10,
                "output_embedding": "graph",
                "output_mlp": {"use_bias": [True, True, True], "units": [200, 100, 1],
                               "activation": ["kgcnn>leaky_relu", "selu", "linear"]}
            }
        },
        "training": {
            "fit": {
                "batch_size": 200, "epochs": 100, "validation_freq": 1, "verbose": 2,
                "callbacks": []
            },
            "compile": {
                "optimizer": {"class_name": "Addons>AdamW", "config": {"lr": 0.0031622776601683794,
                                                                       "weight_decay": 1e-05}
                              }
            },
            "cross_validation": {"class_name": "KFold",
                                 "config": {"n_splits": 5, "random_state": None, "shuffle": True}},
            "scaler": {"class_name": "StandardScaler", "config": {"with_std": True, "with_mean": True, "copy": True}},
            "execute_folds": 1
        },
        "data": {
            "dataset": {
                "class_name": "MoleculeNetDataset",
                "config": {},
                "methods": [
                ]
            },
            "data_unit": "mol/L"
        },
        "info": {
            "postfix": "",
            "kgcnn_version": "2.0.3"
        }
    }

# checked
if architecture_name == 'ChemProp':
    hyper = {
        "model": {
            "class_name": "make_model",
            "module_name": "kgcnn.literature.DMPNN",
            "config": {
                "name": "DMPNN",
                "inputs": [
                    {"shape": [None, 41], "name": "node_attributes", "dtype": "float32", "ragged": True},
                    {"shape": [None, 11], "name": "edge_attributes", "dtype": "float32", "ragged": True},
                    {"shape": [None, 2], "name": "edge_indices", "dtype": "int64", "ragged": True},
                    {"shape": [None, 1], "name": "edge_indices_reverse", "dtype": "int64", "ragged": True}
                ],
                "input_embedding": {
                    "node": {"input_dim": 95, "output_dim": 100},
                    "edge": {"input_dim": 5, "output_dim": 100}
                },
                "pooling_args": {"pooling_method": "sum"},
                "edge_initialize": {"units": 200, "use_bias": True, "activation": "relu"},
                "edge_dense": {"units": 200, "use_bias": True, "activation": "linear"},
                "edge_activation": {"activation": "relu"},
                "node_dense": {"units": 200, "use_bias": True, "activation": "relu"},
                "verbose": 10, "depth": 5,
                "dropout": {"rate": 0.2},
                "output_embedding": "graph",
                "output_mlp": {
                    "use_bias": [True, True, False], "units": [200, 100, 1],
                    "activation": ["kgcnn>leaky_relu", "selu", "linear"]
                }
            }
        },
        "training": {
            "fit": {"batch_size": 32, "epochs": 300, "validation_freq": 1, "verbose": 2, "callbacks": []
                    },
            "compile": {
                "optimizer": {"class_name": "Adam",
                              "config": {"lr": {
                                  "class_name": "ExponentialDecay",
                                  "config": {"initial_learning_rate": 0.001,
                                             "decay_steps": 1600,
                                             "decay_rate": 0.5, "staircase": False}
                              }
                              }
                              },
                "loss": "mean_absolute_error"
            },
            "cross_validation": {"class_name": "KFold",
                                 "config": {"n_splits": 5, "random_state": None, "shuffle": True}},
            "scaler": {"class_name": "StandardScaler", "config": {"with_std": True, "with_mean": True, "copy": True}},
            "execute_folds": 1
        },
        "data": {
            "dataset": {
                "class_name": "MoleculeNetDataset",
                "config": {},
                "methods": [
                    {"map_list": {"method": "set_edge_indices_reverse"}}
                ]
            },
            "data_unit": "mol/L"
        },
        "info": {
            "postfix": "",
            "kgcnn_version": "2.0.3"
        }
    }

if architecture_name == 'DimeNetPP':
    hyper = {
        "model": {
            "class_name": "make_model",
            "module_name": "kgcnn.literature.DimeNetPP",
            "config": {
                "name": "DimeNetPP",
                "inputs": [{"shape": [None], "name": "node_number", "dtype": "float32", "ragged": True},
                           {"shape": [None, 3], "name": "node_coordinates", "dtype": "float32", "ragged": True},
                           {"shape": [None, 2], "name": "range_indices", "dtype": "int64", "ragged": True},
                           {"shape": [None, 2], "name": "angle_indices", "dtype": "int64", "ragged": True}],
                "input_embedding": {"node": {"input_dim": 95, "output_dim": 128,
                                             "embeddings_initializer": {"class_name": "RandomUniform",
                                                                        "config": {"minval": -1.7320508075688772,
                                                                                   "maxval": 1.7320508075688772}}}},
                "emb_size": 128, "out_emb_size": 256, "int_emb_size": 64, "basis_emb_size": 8,
                "num_blocks": 4, "num_spherical": 7, "num_radial": 6,
                "cutoff": 5.0, "envelope_exponent": 5,
                "num_before_skip": 1, "num_after_skip": 2, "num_dense_output": 3,
                "num_targets": 128, "extensive": False, "output_init": "zeros",
                "activation": "swish", "verbose": 10,
                "output_embedding": "graph",
                "output_mlp": {"use_bias": [True, False], "units": [128, 1],
                               "activation": ["swish", "linear"]}
            }
        },
        "training": {
            "fit": {
                "batch_size": 10, "epochs": 872, "validation_freq": 10, "verbose": 2, "callbacks": []
            },
            "compile": {
                "optimizer": {
                    "class_name": "Addons>MovingAverage", "config": {
                        "optimizer": {
                            "class_name": "Adam", "config": {
                                "learning_rate": {
                                    "class_name": "kgcnn>LinearWarmupExponentialDecay", "config": {
                                        "learning_rate": 0.001, "warmup_steps": 30.0, "decay_steps": 40000.0,
                                        "decay_rate": 0.01
                                    }
                                }, "amsgrad": True
                            }
                        },
                        "average_decay": 0.999
                    }
                },
                "loss": "mean_absolute_error",
            },
            "cross_validation": {"class_name": "KFold",
                                 "config": {"n_splits": 5, "random_state": None, "shuffle": True}},
            "scaler": {"class_name": "StandardScaler", "config": {"with_std": True, "with_mean": True, "copy": True}},
        },
        "data": {
            "dataset": {
                "class_name": "MoleculeNetDataset",
                "config": {},
                "methods": [
                    {"map_list": {"method": "set_range", "max_distance": 4, "max_neighbours": 20}},
                    {"map_list": {"method": "set_angle"}}
                ]
            },
            "data_unit": "mol/L"
        },
        "info": {
            "postfix": "",
            "postfix_file": "",
            "kgcnn_version": "2.0.3"
        }
    }

if architecture_name == 'GIN':
    hyper = {
        "model": {
            "class_name": "make_model",
            "module_name": "kgcnn.literature.GIN",
            "config": {
                "name": "GIN",
                "inputs": [{"shape": [None, 41], "name": "node_attributes", "dtype": "float32", "ragged": True},
                           {"shape": [None, 2], "name": "edge_indices", "dtype": "int64", "ragged": True}],
                "input_embedding": {"node": {"input_dim": 96, "output_dim": 100}},
                "depth": 4,
                "dropout": 0.1,
                "gin_mlp": {"units": [100, 100], "use_bias": True, "activation": ["relu", "relu"],
                            "use_normalization": True, "normalization_technique": "batch"},
                "gin_args": {},
                "last_mlp": {"use_bias": [True, True, True], "units": [200, 100, 1],
                             "activation": ["kgcnn>leaky_relu", "selu", "linear"]},
                "output_embedding": "graph", "output_to_tensor": True,
                "output_mlp": {"use_bias": True, "units": 1,
                               "activation": "linear"}
            }
        },
        "training": {
            "fit": {"batch_size": 32, "epochs": 200, "validation_freq": 1, "verbose": 2, "callbacks": []
                    },
            "compile": {
                "optimizer": {"class_name": "Adam",
                              "config": {"lr": {
                                  "class_name": "ExponentialDecay",
                                  "config": {"initial_learning_rate": 0.001,
                                             "decay_steps": 1600,
                                             "decay_rate": 0.5, "staircase": False}
                              }
                              }
                              },
                "loss": "mean_absolute_error"
            },
            "cross_validation": {"class_name": "KFold",
                                 "config": {"n_splits": 5, "random_state": None, "shuffle": True}}
        },
        "data": {
            "dataset": {"class_name": "MoleculeNetDataset",
                        "config": {},
                        "methods": []
                        },
            "data_unit": "unit"
        },
        "info": {
            "postfix": "",
            "kgcnn_version": "2.0.3"
        }
    }

if architecture_name == 'GINE':
    hyper = {
        "model": {
            "class_name": "make_model_edge",
            "module_name": "kgcnn.literature.GIN",
            "config": {
                "name": "GIN",
                "inputs": [{"shape": [None, 41], "name": "node_attributes", "dtype": "float32", "ragged": True},
                           {"shape": [None, 11], "name": "edge_attributes", "dtype": "float32", "ragged": True},
                           {"shape": [None, 2], "name": "edge_indices", "dtype": "int64", "ragged": True}],
                "input_embedding": {"node": {"input_dim": 96, "output_dim": 100},
                                    "edge": {"input_dim": 5, "output_dim": 100}},
                "depth": 4,
                "dropout": 0.1,
                "gin_mlp": {"units": [100, 100], "use_bias": True, "activation": ["relu", "relu"],
                            "use_normalization": True, "normalization_technique": "batch"},
                "gin_args": {},
                "last_mlp": {"use_bias": [True, True, True], "units": [200, 100, 1],
                             "activation": ["kgcnn>leaky_relu", "selu", "linear"]},
                "output_embedding": "graph",
                "output_mlp": {"use_bias": True, "units": 1,
                               "activation": "linear"}
            }
        },
        "training": {
            "fit": {"batch_size": 32, "epochs": 200, "validation_freq": 1, "verbose": 2, "callbacks": []
                    },
            "compile": {
                "optimizer": {"class_name": "Adam",
                              "config": {"lr": {
                                  "class_name": "ExponentialDecay",
                                  "config": {"initial_learning_rate": 0.001,
                                             "decay_steps": 1600,
                                             "decay_rate": 0.5, "staircase": False}
                              }
                              }
                              },
                "loss": "mean_absolute_error"
            },
            "cross_validation": {"class_name": "KFold",
                                 "config": {"n_splits": 5, "random_state": None, "shuffle": True}}
        },
        "data": {
            "dataset": {
                "class_name": "MoleculeNetDataset",
                "config": {},
                "methods": []
            },
            "data_unit": "unit"
        },
        "info": {
            "postfix": "",
            "kgcnn_version": "2.0.3"
        }
    }

# Print to visually make sure we have parsed correctly the parameters
print("My parameters")
print("Loss", nn_loss)
print("LSTM", nn_lstmunits)
print("DENSE", nn_denseunits)
print("PROBA", nn_return_proba)
print("LR start", nn_lr_start)

# parameters of the model for OCHEM
random_seed = seed
np.random.seed(seed);

# Check against other keras models in OCHEM
log_filename = 'model.log';
modelname = "model.h5";

# +
if gpu >= 0:
    device = "gpu:%s/" % (gpu)
else:
    device = "cpu"

if gpu >= 0:
    config = tf.compat.v1.ConfigProto()
    config.gpu_options.allow_growth = True
    sess = tf.compat.v1.Session(config=config)

print('before', TRAIN)

if TRAIN == "True":
    # change the number of output (target) of the model
    # define columns names to grab from the input file
    cols = ['Result%s' % (i) for i in range(output_dim)]
    # I don't like this ... we need a proper way to assign output_dim!
    try:
        hyper['model']["config"]['output_mlp']['units'][-1] = output_dim
    except:
        hyper['model']["config"]['output_mlp']['units'] = output_dim
    try:
        hyper['model']["config"]['last_mlp']['units'][-1] = output_dim  # multi MLP normally so list of arrays!
    except:
        print('cannot change the last_mlp output dimension')

    # failed next line for GINE parameters
    hyperparams = HyperParameter(hyper,
                                 model_name=hyper["model"]["config"]["name"],
                                 model_module=hyper["model"]["module_name"],
                                 model_class=hyper["model"]["class_name"],
                                 dataset_name="MoleculeNetDataset")

    inputs = hyperparams["model"]["config"]["inputs"]

    print('training')
    dataset = MoleculeNetDataset(file_name=TRAIN_FILE, data_directory="", dataset_name="MoleculeNetDataset")
    dataset.prepare_data(overwrite=True, smiles_column_name="smiles",
                         make_conformers=True, add_hydrogen=True,
                         optimize_conformer=True, num_workers=None)
    dataset.read_in_memory(label_column_name=cols, add_hydrogen=False,
                           has_conformers=False)
    dataset.set_attributes(add_hydrogen=False,
                           has_conformers=True)
    print('*' * 100)

    # For geometric models find a way to deal with them...
    # need to identified which are 3D models properly ...
    # this will change the prepare_data function too 
    dataset.set_methods(hyper["data"]["dataset"]["methods"])

    # Test model input is in dataset.
    invalid_graphs = dataset.clean(hyperparams["model"]["config"]["inputs"])
    dataset.assert_valid_model_input(hyperparams["model"]["config"]["inputs"])  # failed for GCN code here

    # Model identification
    make_model = get_model_class(hyperparams["model"]["config"]["name"], hyperparams["model"]["class_name"])

    # check Dataset
    data_name = dataset.dataset_name
    data_unit = hyperparams["data"]["data_unit"]
    data_length = dataset.length
    k_fold_info = hyperparams["training"]["cross_validation"]["config"]

    labels = np.array(dataset.obtain_property("graph_labels"), dtype="float")

    print("dataset length:", len(labels))
    # is a classification or a regression problem
    isClass = isClassifer(labels)
    # change activation function ! we do use sigmoid or softmax in general
    if isClass:
        hyperparams['model']["config"]['output_mlp']['activation'][-1] = 'sigmoid'

    #kf = KFold(**k_fold_info)
    #split_indices = kf.split(X=np.arange(data_length - 1)[:, None])
    split_indices = np.uint8(np.floor(np.random.uniform(0,10,size=data_length))) 
    test = split_indices==0
    train = ~test
    
    # hyper_fit and epochs
    hyper_fit = hyperparams['training']['fit']
    epo = hyper_fit['epochs']
    epostep = hyper_fit['validation_freq']

    train_loss = []
    test_loss = []
    mae_5fold = []
    all_test_index = []
    model = None

    # how it works in Keras the Train / Val split
    filename = 'loss.csv'

    hl = CSVLogger(filename, separator=",", append=True)
    es = EarlyStopping(monitor='val_loss', mode='min', patience=earlystop_patience, restore_best_weights=True)
    reduce_lr = ReduceLROnPlateau(monitor='val_loss', factor=redlr_factor, verbose=0, patience=redlr_patience,
                                  min_lr=1e-5)

    # Make model
    model = make_model(**hyperparams['model']["config"])

    # data extraction is define here for a given Kfold => make them like in OCHEM classical call
    # TODO so just use the first split provided by KFold for the moment
    #for train_index, test_index in split_indices:
        # Select train and test data.
    #    xtrain, ytrain = dataset[train_index].tensor(inputs), labels[train_index]
    #    xtest, ytest = dataset[test_index].tensor(inputs), labels[test_index]
    #    break

     # data extraction is define here for a given Kfold => make them like in OCHEM classical call TODO so just use the first split provided by KFold for the moment
    testidx = np.transpose(np.where(test==True))
    trainidx = np.transpose(np.where(train==True))
    print("*"*10)
    print(data_length,np.sum(train), len(trainidx),np.sum(test), len(testidx))
    print("*"*10)
    xtrain, ytrain = dataset[trainidx].tensor(inputs), labels[trainidx]
    xtest, ytest = dataset[testidx].tensor(inputs), labels[testidx]
    print(len(xtrain),len(xtest))

    opt = tfa.optimizers.RectifiedAdam(learning_rate=0.001)
    # need to change this for class / regression 
    if isClass:
        model.compile(opt, loss=BCEmask if isClass else RMSEmask, metrics=['accuracy'])
    else:
        model.compile(opt, loss=BCEmask if isClass else RMSEmask, metrics=[tf.keras.metrics.RootMeanSquaredError()])

    print(model.summary())

    # Start and time training
    hyper_fit = hyperparams['training']['fit']
    start = time.process_time()

    # need to change that to have ragged not numpy or tensor error 
    hist = model.fit(xtrain, ytrain,
                     validation_data=(xtest, ytest),
                     batch_size=batch_size,
                     shuffle=True,
                     epochs=nbepochs,
                     verbose=2,
                     callbacks=[es, reduce_lr, hl]
                     )

    model.save_weights(modelname)
    print("Saved model to disk")
    pickle.dump(hyperparams, open("modelparameters.p", "wb"))

    # Probably this should become model.save_weights()
    tar = tarfile.open(MODEL_FILE, "w:gz");
    tar.add(modelname);
    tar.add("modelparameters.p")

    tar.close();

    try:
        os.remove(modelname);
        os.remove("modelparameters.p")
    except:
        pass;

    print("Relax!");

else:
    # Look in OCHEM how other Keras models are stored.
    tar = tarfile.open(MODEL_FILE);
    tar.extractall();
    tar.close();

    print("Loaded model from disk")
    hyper = pickle.load(open("modelparameters.p", "rb"))

    df = pd.read_csv(APPLY_FILE)

    validsmiles = df['smiles'].apply(lambda x: SmilesOK(x))
    # if 'Result0'  don't add this line
    if not 'Result0' in df.columns:
        df['Result0'] = 0

    df.to_csv(APPLY_FILE, index=False)
    cols = ['Result%s' % (i) for i in range(output_dim)]  # change this for MTL tasks not the case today
    print(cols)

    inputs = hyper["model"]["config"]["inputs"]

    print('evaluation')
    dataset = MoleculeNetDataset(file_name=APPLY_FILE, data_directory="", dataset_name="MoleculeNetDataset")
    # add Hydrogen ??? is it a good idea for 2D models ...
    dataset.prepare_data(overwrite=True, smiles_column_name="smiles",
                         make_conformers=True, add_hydrogen=True,
                         optimize_conformer=True, num_workers=None)
    
    dataset.read_in_memory(label_column_name="Result0", add_hydrogen=False,
                           has_conformers=False)
    dataset.set_attributes(add_hydrogen=False,
                           has_conformers=True)

    dataset.set_methods(hyper["data"]["dataset"]["methods"])



    invalid_graphs = dataset.clean(inputs)
    print('-' * 10)
    print('error graphs:', invalid_graphs)
    print('-' * 10)

    # Model creation    
    make_model = get_model_class(hyper["model"]["config"]["name"], hyper["model"]["class_name"])

    model = make_model(**hyper['model']["config"])
    # load stored best weights

    model.load_weights(modelname)

    x_pred = dataset.tensor(inputs)

    a_pred = model.predict(x_pred)

    dfres = pd.DataFrame(a_pred, columns=cols)

    # we need to check if the graph has not computed some cases and regenerate the full index 
    if len(invalid_graphs) > 0:
        dfresall = pd.DataFrame(index=df.index, columns=dfres.columns)
        idx = [i for i in range(len(df)) if i not in invalid_graphs]
        dfres.index = idx
        dfresall.loc[idx, :] = dfres.loc[idx, :]

    else:
        dfresall = dfres

    dfresall.to_csv(RESULT_FILE, index=False)

    try:
        os.remove(modelname);
        os.remove("modelparameters.p")
    except:
        pass;

    print("Relax!");
# -
