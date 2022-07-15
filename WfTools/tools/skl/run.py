from base import base
import joblib

if __name__ == "__main__":
    from load import read_task
    import sys
    import os
    task = read_task(sys.argv[1])

    if task['train_mode']:
        model = base(task)
        model.load_train_data()
        model.build_model()
        model.train()
        #print (">>> score: "+str(model.score(model.X, model.Y)))
        joblib.dump(model, task['model_file'])
        
        if os.path.exists(task['apply_data_file']):
            model.load_test_data()
            model.apply()
            #print(">>> validate score: "+str(model.validate(task['apply_data_file'])))
            
    else:
        model = joblib.load(task['model_file'])
        model.load_test_data()
        model.apply()
