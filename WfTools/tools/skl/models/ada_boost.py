from sklearn.ensemble import AdaBoostRegressor as skmodel
from sklearn.multioutput import MultiOutputRegressor as multi

def builder (rng, njobs, mlt):
    m = skmodel(n_estimators=300, loss='linear', random_state=rng)
    if mlt:
        m = multi(m, n_jobs=njobs)
    return m
