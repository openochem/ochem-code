from sklearn.linear_model import RANSACRegressor as skmodel
from sklearn.ensemble import ExtraTreesRegressor as base

def builder (rng, njobs, mlt):
    return skmodel(base_estimator=base(), random_state=rng, loss='squared_loss')
