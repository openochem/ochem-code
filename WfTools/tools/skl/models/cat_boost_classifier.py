from catboost import CatBoostClassifier as skmodel

def builder (rng, njobs, mlt):
    m = skmodel(depth=10,logging_level='Silent')
    return m
