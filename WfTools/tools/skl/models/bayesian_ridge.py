from sklearn.linear_model import BayesianRidge as skmodel
from sklearn.multioutput import MultiOutputRegressor as multi


def builder (rng, njobs, mlt):
    m = skmodel(n_iter=300, normalize=False)
    if mlt:
        m = multi(m, n_jobs=njobs)
    return m
