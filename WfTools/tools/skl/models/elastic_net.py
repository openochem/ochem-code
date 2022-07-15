from sklearn.linear_model import MultiTaskElasticNet as multi_skmodel
from sklearn.linear_model import ElasticNet as skmodel


def builder (rng, njobs, mlt):
    if mlt:
        return multi_skmodel(normalize=False, random_state=rng)
    else:
        return skmodel(normalize=False, random_state=rng)
