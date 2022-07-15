from sklearn.kernel_ridge import KernelRidge as skmodel

# tunable:
# alpha - smaller
# coef0: larger
#
# check numerical stability with aggressive values

def builder (rng, njobs, mlt):
     return skmodel(alpha=0.1, kernel='rbf', coef0=10)
