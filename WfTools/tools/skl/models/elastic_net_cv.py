from sklearn.linear_model import MultiTaskElasticNetCV as multi_skmodel
from sklearn.linear_model import ElasticNetCV as skmodel


def builder (rng, njobs, mlt):
    if mlt:
        return multi_skmodel(normalize=False, n_jobs=njobs, random_state=rng)
    else:
        return skmodel(normalize=False, n_jobs=njobs, random_state=rng)
