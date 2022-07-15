from sklearn.ensemble import BaggingRegressor as skmodel
from sklearn.multioutput import MultiOutputRegressor as multi
from sklearn.ensemble import ExtraTreesRegressor as base

def builder (rng, njobs, mlt):
    m = skmodel(base_estimator=base(n_jobs=njobs), n_estimators=30, n_jobs=njobs, random_state=rng)
    if mlt:
        m = multi(m, n_jobs=njobs)
    return m
