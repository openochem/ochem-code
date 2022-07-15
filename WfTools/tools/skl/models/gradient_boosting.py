from sklearn.ensemble import GradientBoostingRegressor as skmodel
from sklearn.multioutput import MultiOutputRegressor as multi

def builder (rng, njobs, mlt):
    m = skmodel(n_estimators=500, random_state=rng, loss='huber')
    if mlt:
        m = multi(m, n_jobs=njobs)
    return m


