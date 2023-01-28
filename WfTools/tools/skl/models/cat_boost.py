from catboost import CatBoostRegressor as skmodel

def builder (rng, njobs, mlt):
    m = skmodel(loss_function='RMSE',depth=10,logging_level='Silent')
    return m
