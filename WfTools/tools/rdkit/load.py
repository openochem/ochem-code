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

    task['input_file'] = check_option(config, 'Task', 'input_file')
    task['output_file'] = check_option(config, 'Task', 'output_file')
    
    task['autocorr2d']  = check_option(config, 'Task', 'autocorr2d', type=bool, default=False)
    
    task['autocorr3d']  = check_option(config, 'Task', 'autocorr3d', type=bool, default=False)
    
    task['getaway']     = check_option(config, 'Task', 'getaway', type=bool, default=False)
    
    task['morse']       = check_option(config, 'Task', 'morse', type=bool, default=False)
    
    task['rdf']         = check_option(config, 'Task', 'rdf', type=bool, default=False)
    
    task['whim']        = check_option(config, 'Task', 'whim', type=bool, default=False)
    task['whim_thresh'] = check_option(config, 'Task', 'whim_thresh', type=float, default=0.01)

    task['maccs']       = check_option(config, 'Task', 'maccs', type=bool, default=False)
    
    task['topological'] = check_option(config, 'Task', 'topological', type=bool, default=False)
    task['topological_nbits'] = check_option(config, 'Task', 'topological_nbits', type=int, default=1024)
    
    task['morgan_nbits'] = check_option(config, 'Task', 'morgan_nbits', type=int, default=1024)
    task['morgan_radius'] = check_option(config, 'Task', 'morgan_radius', type=int, default=2)
    task['morgan_fcfp'] = check_option(config, 'Task', 'morgan_fcfp', type=bool, default=False)
    task['morgan_counts'] = check_option(config, 'Task', 'morgan_counts', type=bool, default=False)
    
    task['morgan'] = check_option(config, 'Task', 'morgan', type=bool, default=False)
    if task['morgan']:
        if task['morgan_counts']:
            task['morgan'] = False
    else:
        task['morgan_counts'] = False

    task['rdkitdef']     = check_option(config, 'Task', 'rdkitdef', type=bool, default=False)

    task['avalon_nbits'] = check_option(config, 'Task', 'avalon_nbits', type=int, default=1024)
    task['avalon_counts'] = check_option(config, 'Task', 'avalon_counts', type=bool, default=False)
    task['avalon'] = check_option(config, 'Task', 'avalon', type=bool, default=False)
    if task['avalon']:
        if task['avalon_counts']:
            task['avalon'] = False
    else:
        task['avalon_counts'] = False

    task['atom_pairs'] = check_option(config, 'Task', 'atom_pairs', type=bool, default=False)
    task['btf'] = check_option(config, 'Task', 'btf', type=bool, default=False)
    task['bpf'] = check_option(config, 'Task', 'bpf', type=bool, default=False)
    task['torsions'] = check_option(config, 'Task', 'torsions', type=bool, default=False)
    
    task['sascore'] = check_option(config, 'Task', 'sascore', type=bool, default=False)
    task['scalars'] = check_option(config, 'Task', 'scalars', type=bool, default=False)
    task['scalars_secondary'] = check_option(config, 'Task', 'scalars_secondary', type=bool, default=False)
    task['scalars_fragments'] = check_option(config, 'Task', 'scalars_fragments', type=bool, default=False)
    
    return task
