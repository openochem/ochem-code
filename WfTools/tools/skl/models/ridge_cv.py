from sklearn.linear_model import RidgeCV as skmodel

def builder (rng, njobs, mlt):
    return skmodel(normalize=False)
