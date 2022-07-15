from sklearn.linear_model import LassoLarsCV as skmodel
from sklearn.multioutput import MultiOutputRegressor as multi

def builder (rng, njobs, mlt):
    m = skmodel(max_iter=500, normalize=False, n_jobs=njobs)
    if mlt:
        m = multi(m, n_jobs=njobs)
    return m
