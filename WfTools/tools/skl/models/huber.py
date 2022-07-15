from sklearn.linear_model import HuberRegressor as skmodel
from sklearn.multioutput import MultiOutputRegressor as multi

def builder (rng, njobs, mlt):
    m = skmodel()
    if mlt:
        m = multi(m, n_jobs=njobs)
    return m
