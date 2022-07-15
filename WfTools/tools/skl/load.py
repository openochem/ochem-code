try:
    # python2
    import ConfigParser
except ImportError:
    # python3
    import configparser as ConfigParser

def check_option (config, section, name, default=None, type=str):
        if not config.has_section(section):
            raise Exception('No such config section: '+section) 
        
        if type==str:
            accessor = config.get
        elif type==int:
            accessor = config.getint
        elif type==float:
            accessor = config.getfloat
        elif type==bool:
            accessor = config.getboolean
        else:
            raise Exception('Internal exception #389423')
        
        if not config.has_option(section, name):
            if default is None:
                raise Exception('Mandatory config option '+name+' not found')
            else:
                return default
        else:
            try:
                return accessor(section, name)
            except:
                raise Exception('Bad value for config option: '+name+', expected type: '+str(type))

def read_task (path):
    config = ConfigParser.ConfigParser()
    with open(path, 'r') as fh:
        config.readfp(fh)
    task = {}

    task['model_file'] = check_option(config, 'Task', 'model_file')
    task['train_mode'] = check_option(config, 'Task', 'train_mode', type=bool, default=False)
    task['parallelize'] = check_option(config, 'Details', 'parallelize', type=bool, default=False)

    if task['train_mode']:
        task['seed'] = check_option(config, 'Details', 'seed', type=int, default=0)
       
        task['train_data_file'] = check_option(config, 'Task', 'train_data_file')
        task['apply_data_file'] = check_option(config, 'Task', 'apply_data_file', default="nonexistent")
        task['result_file'] = check_option(config, 'Task', 'result_file', default="nonexistent")

        task['base_model'] = check_option(config, 'Task', 'base_model', default="all")
        if task['base_model'] == 'all_cv':
            task['base_model'] = 'all'
            task['cv_only'] = True
        else:
            task['cv_only'] = False

    else:
        task['apply_data_file'] = check_option(config, 'Task', 'apply_data_file')
        task['result_file'] = check_option(config, 'Task', 'result_file')

    return task
