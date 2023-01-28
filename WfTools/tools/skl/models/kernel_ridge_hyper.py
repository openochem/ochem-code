from sklearn.kernel_ridge import KernelRidge
from sklearn.model_selection import GridSearchCV, train_test_split
import numpy as np

# tunable:
# alpha - smaller
# coef0: larger
#
# check numerical stability with aggressive values

def builder (rng, njobs, mlt):
    kernel_used = 'rbf'
    alpha = np.logspace(-10,-1,10)
    gamma = np.logspace(-10,-1,10)
    tuned_parameters = [{'kernel':[kernel_used],'alpha': alpha, 'gamma': gamma}]
    return GridSearchCV(KernelRidge(), tuned_parameters, cv=5, scoring='neg_mean_absolute_error', n_jobs=njobs, verbose=1)
