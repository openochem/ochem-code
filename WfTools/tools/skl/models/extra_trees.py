from sklearn.ensemble import ExtraTreesRegressor as skmodel

def builder (rng, njobs, mlt):
    return skmodel(n_estimators=300, n_jobs=njobs, random_state=rng)


