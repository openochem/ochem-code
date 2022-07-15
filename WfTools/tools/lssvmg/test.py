from model import model
import joblib
import numpy as np
import csv

descr = "ALOGPS"
m = model(center=False, whiten=False, reduce_dim=False, nfold=3)
m.load_data(path="runs/LSSVMG_"+descr+"/train.csv", target_columns=[-1])
m.fit_preproc()
m.train(kernel='rbf', gamma=2**9, lamb=2**0, rbf_gamma=2**-15, optimize=False, fast_optimize=True)

m.load_data(path="runs/LSSVMG_"+descr+"/apply.csv")
pred = m.apply("result.csv")
