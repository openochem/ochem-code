from sklearn.linear_model import OrthogonalMatchingPursuitCV as skmodel
from sklearn.multioutput import MultiOutputRegressor as multi

def builder (rng, njobs, mlt):
    m = skmodel(n_jobs=njobs, normalize=True)
    if mlt:
        m = multi(m, n_jobs=njobs)
    return m
