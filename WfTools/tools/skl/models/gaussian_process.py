from sklearn.gaussian_process import GaussianProcessRegressor as skmodel

def builder (rng, njobs, mlt):
    return skmodel(normalize_y=True, random_state=rng)
