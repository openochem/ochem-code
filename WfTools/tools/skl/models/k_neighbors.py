from sklearn.neighbors import KNeighborsRegressor as skmodel

def builder (rng, njobs, mlt):
    return skmodel(n_jobs=njobs, p=1)
