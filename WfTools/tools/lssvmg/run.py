from base import base
import joblib

if __name__ == "__main__":
    from load import read_task
    import sys
    import os
    task = read_task(sys.argv[1])
    if task['train_mode']:
        if task.get('gpu')==-1:
            from lssvm import LSSVM
        else:
            from lssvm_gpu import LSSVM_GPU as LSSVM
    
    if task['train_mode']:
        m = base(LSSVM, gpu = task['gpu'], center=task['center'], pca=task['pca'], seed=task['seed'],  metric=task['metric'])
        m.load_train_data(path=task['train_data_file'], ntargets=task['ntargets'])
        m.train(task['nfold'], glob_opt=task['glob_opt'])
        # workaround for running at no-gpu machine without cupy after gpu training
        if task.get('gpu'):
            from lssvm import LSSVM as LSSVM_cpu
            mcpu = LSSVM_cpu(gpu=-1)
            mcpu.set_params(gamma=m.lssvm.gamma, lamb=m.lssvm.lamb, rbf_gamma=m.lssvm.rbf_gamma)
            mcpu.fitted = True
            mcpu.outdim = m.lssvm.outdim
            mcpu.alpha = m.lssvm.alpha
            mcpu.b = m.lssvm.b
            mcpu.X_train = m.lssvm.X_train
            mcpu.Y_train = m.lssvm.Y_train
            m.lssvm = mcpu
        
        joblib.dump(m, task['model_file'])
        if os.path.exists(task['apply_data_file']):
            m.load_test_data(path=task['apply_data_file'])
            m.apply(task['result_file'])

    else:
        m = joblib.load(task['model_file'])
        m.load_test_data(path=task['apply_data_file'])
        m.apply(task['result_file'])
        
                              
