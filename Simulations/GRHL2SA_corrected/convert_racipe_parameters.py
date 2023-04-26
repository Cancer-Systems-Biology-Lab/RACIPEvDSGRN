import numpy as np
import DSGRN
import os
import sys
from tqdm import tqdm
import multiprocessing as mp
from functools import partial


def dsgrn_FPs(param_file,param_name_file,network_file):
    """
    Create a file which gives the DSGRN FPs associated to a RACIPE parameter file.
    The file name is the same as param_file with _FPs appended to the name. 

    Input:
        param_file - (string) path to the RACIPE parameter file
        param_name_file - (string) path to a .prs file which describes which columns 
                          of the parameter file correspond to which parameter. 
        network_file - (string) path to a DSGRN network file.
    """
    racipe_params = np.loadtxt(param_file)
    num_params = racipe_params.shape[0]
    network = DSGRN.Network(network_file)
    pg = DSGRN.ParameterGraph(network)
    param_converter = Racipe2DSGRN(param_name_file,network)
    # racipe2dsgrn = make_racipe2dsgrn(param_name_file, network)
    print('Starting FP labeling of parameters in ' + param_file + '...\n', flush = True)
    # num_processors = mp.cpu_count()
    num_processors = 40
    pool = mp.Pool(num_processors)
    chunksize = max(2**12,num_params//num_processors)
    pbar = tqdm(total = num_params)
    work = partial(get_FPs,param_converter = param_converter,pg = pg,Network = network)
    with open(param_file.rstrip('.dat') + '_FPs.dat','w') as FP_writer, open(param_file.rstrip('.dat') + '_eqval.dat','w') as eq_writer:
        for k,vals in enumerate(pool.imap(work,racipe_params,chunksize = chunksize)):
        # for k in range(len(racipe_params)):
            # vals = get_FPs(racipe_params[k], pg = pg, param_converter=param_converter, Network=network)
            FPs = vals[0]
            eq_vals = vals[1]
            if FPs == -1: #RACIPE parameter does not give valid DSGRN parameter
                FP_writer.write('-1')
                eq_writer.write('-1')
            elif len(FPs) == 0: #no fixed points
                FP_writer.write('-2')
                eq_writer.write('-2')
            else:
                for (i,FP) in enumerate(FPs):
                    for coord in FP:
                        FP_writer.write(str(coord))
                    if i < len(FPs) - 1:
                        FP_writer.write('\t')
                for (i,eq) in enumerate(eq_vals):
                    for (j,coord) in enumerate(eq):
                        eq_writer.write('{:.3f}'.format(coord))
                        if j < len(eq) - 1:
                            eq_writer.write(',')
                    if i < len(eq_vals) - 1:
                        eq_writer.write('\t')
            if k < num_params - 1:
                FP_writer.write('\n')
                eq_writer.write('\n')
            pbar.update()
    print('done.\n',flush = True)

class Racipe2DSGRN:
    #pickleable version of racipe2DSGRN so that multiprocessing will work. 
    def __init__(self,param_name_file,network):
        self.init_attributes(param_name_file,network)
    def init_attributes(self,param_name_file,network):
        """
        Initialize attributes needed for convert.
        """
        N = network.size()
        node_names = [network.name(i) for i in range(N)]
        num_param_names = 0
        index_dict = dict()
        with open(param_name_file,'r') as reader:
            for i,line in enumerate(reader.readlines()):
                if i == 0:
                    #first line gives column names
                    continue
                num_param_names += 1
                param_type_end = line.find('_')
                param_type = line[:param_type_end]
                if param_type not in ['Prod','Deg','Trd','Act','Inh']:
                    continue
                node_start = param_type_end + 4
                if param_type in ['Prod','Deg']:
                    node_end = line.find('\t',node_start)
                    node = line[node_start:node_end]
                    node_index = node_names.index(node)
                    if param_type == 'Prod':
                        index_dict[i-1] = ('Prod',node_index)
                    else:
                        index_dict[i-1] = ('Deg',node_index)
                else:
                    source_end = line.find('To',node_start)
                    source = line[node_start:source_end]
                    target_start = source_end + 2
                    target_end = line.find('\t',target_start)
                    target = line[target_start:target_end]
                    source_index = node_names.index(source)
                    target_index = node_names.index(target)
                    if param_type == 'Trd':
                        index_dict[i-1] = ('Trd',source_index,target_index)
                    else:
                        index_dict[i-1] = ('lam',source_index,target_index) 
        num_inputs = np.array([[len(network.inputs(i))] if len(network.inputs(i))>0 else [1] for i in range(N)])
        self.N = N
        self.index_dict = index_dict
        self.num_inputs = num_inputs
        self.num_param_names = num_param_names
    
    def convert(self,param,return_deg = False):#, N=N, index_dict = index_dict, num_inputs = num_inputs):
        """
        Convert a RACIPE parameter to a DSGRN parameter. 

        Input: 
            param - (numpy vector) RACIPE parameter extracted from a row of a RACIPE parameter file
            return_deg - (optional,boolean) Returns the degradation rates Deg if True
        Output:
            L,U,theta - NxN numpy array indexed by [source,target]
            Deg - Nx1 numpy array. Only returned if return_deg == True
        """
        N = self.N
        index_dict = self.index_dict
        num_inputs = self.num_inputs
        num_param_names = self.num_param_names
        Prod = np.zeros([N,1])
        Trd = np.zeros([N,N])
        Deg = np.zeros([N,1])
        lam = np.zeros([N,N])
        K = len(param)
        for i, val in enumerate(reversed(param)):
            if num_param_names-1-i in index_dict:
                cur_entry = index_dict[num_param_names-1-i]
                if cur_entry[0] == 'Prod':
                    Prod[cur_entry[1],0] = val
                elif cur_entry[0] == 'Deg':
                    Deg[cur_entry[1],0] = val
                elif cur_entry[0] == 'Trd':
                    Trd[cur_entry[1],cur_entry[2]] = val
                else: #cur_entry[0] == 'lam'
                    if val > 1:
                        lam[cur_entry[1],cur_entry[2]] = 1/val
                    else: 
                        lam[cur_entry[1],cur_entry[2]] = val
        Prod = Prod**(1/num_inputs) #split production rate evenly across input nodes
        theta = Trd*Deg #multiply Trd[i,:] by d[i,0]
        L = lam*np.transpose(Prod) #multiply lam[:,i] by Prod[i,0]
        U = np.ones(lam.shape)*np.transpose(Prod) #U[j -> i] = U[j,i] = Prod[i,0]
        if return_deg:
            return L,U,theta,Deg
        else:
            return L,U,theta

def get_FPs(param,pg,param_converter,Network):
    L,U,theta,Deg = param_converter.convert(param,return_deg = True)
    # print(L)
    # print(theta)
    # L,U,theta = racipe2dsgrn(param)
    dsgrn_param_error = False
    try:
        dsgrn_param_index = DSGRN.par_index_from_sample(pg,L,U,theta)
    except TypeError:
        return -1,-1
    # print(dsgrn_param_index)
    if (dsgrn_param_index == -1):
        return -1,-1
    dsgrn_param = pg.parameter(dsgrn_param_index)
    eq_cells = DSGRN.EquilibriumCells(dsgrn_param,eqtype = 'topdim')
    eq_vals = get_eq_vals(L,U,theta,Deg,eq_cells,Network)
    return eq_cells, eq_vals

def get_eq_vals(L,U,theta,Deg,eq_cells,Network):
    N = Network.size()
    eq_vals = []
    theta = theta/Deg #undo multiplication of thresholds by the degradation rate
    sorted_thresholds = [sorted([theta[source,target] for target in Network.outputs(source)]) for source in range(N)]
    for cell in eq_cells:
        input_vec = np.zeros([N,1])
        for node in range(N):
            cur_prod = 1
            for source_set in Network.logic(node):
                cur_sum = 0
                for source in source_set:
                    activating = Network.interaction(source,node)
                    if cell[source] == 0:
                        reference_threshold = 0
                    elif cell[source] < len(sorted_thresholds[source]):
                        reference_threshold = sorted_thresholds[source][cell[source]-1]
                    else:
                        reference_threshold = np.inf
                    if theta[source,node] > reference_threshold:
                        if activating:
                            cur_sum += L[source,node]
                        else:
                            cur_sum += U[source,node]
                    else:
                        if activating:
                            cur_sum += U[source,node]
                        else:
                            cur_sum += L[source,node]
                cur_prod *= cur_sum
            input_vec[node,0] = cur_prod
        eq_vals.append(list((input_vec/Deg).flatten()))
    return eq_vals


def label_index(param_file,param_name_file,network_file):
    """
    Create a file which labels the parameters in a RACIPE parameter file with the 
    DSGRN parameter. The file name is the same as param_file with _index appended 
    to the name. 

    Input:
        param_file - (string) file path to the RACIPE parameter file
        param_name_file - (string) file path to a .prs file which describes which columns
                          of the parameter file correspond to which parameter
        network_file - (string) file path to a DSGRN network specification file
    Output:
        None (results are written to a .dat file)
    """
    racipe_params = np.loadtxt(param_file)
    num_params = racipe_params.shape[0]
    network = DSGRN.Network(network_file)
    pg = DSGRN.ParameterGraph(network)
    racipe2dsgrn = make_racipe2dsgrn(param_name_file,network)
    print('Starting index labeling of parameters in ' + param_file + '...', flush = True)
    with open(param_file.rstrip('.dat') + '_index.dat','w') as writer:
        for k in tqdm(range(num_params)):
            param = racipe_params[k,:]
            L,U,theta = racipe2dsgrn(param)
            dsgrn_param_error = False
            try:
                dsgrn_param_index = DSGRN.par_index_from_sample(pg,L,U,theta)
            except TypeError:
                dsgrn_param_error = True
                dsgrn_param_index = -1
            if k < num_params - 1:
                writer.write(str(dsgrn_param_index) + '\n')
            else:
                writer.write(str(dsgrn_param_index))
    print('done.',flush = True)

def make_racipe2dsgrn(param_name_file,network):
    """
    Returns a function (racipe2dsgrn) which converts a RACIPE parameter to its corresponding
    DSGRN parameter. Given that there are N parameters, assumes the last N columns of 
    param_name_file correspond to these parameters.

    Input:
        param_name_file - (string) file path to a .prs file which describes which columns
                          of the parameter file correspond to which parameter
        network - DSGRN network object. The RACIPE node names must agree with the 
                  DSGRN node names. 
    Output: 
        racipe2dsgrn (function)
            Input: 
                param - numpy vector with entries corresponding to a row of the RACIPE parameter file. 
                return_deg - (boolean,optional) If True, returns the vector of degradation rates. 
            Output:
                L,U,theta - (default) NxN numpy matrices. Note that entries are indexed by [source,target].
                Deg
                L,U,theta,Deg - (if return_deg == True) Deg is a Nx1 numpy array
    """
    N = network.size()
    node_names = [network.name(i) for i in range(N)]
    num_param_names = 0
    index_dict = dict()
    with open(param_name_file,'r') as reader:
        for i,line in enumerate(reader.readlines()):
            if i == 0:
                #first line gives column names
                continue
            num_param_names += 1
            param_type_end = line.find('_')
            param_type = line[:param_type_end]
            if param_type not in ['Prod','Deg','Trd','Act','Inh']:
                continue
            node_start = param_type_end + 4
            if param_type in ['Prod','Deg']:
                node_end = line.find('\t',node_start)
                node = line[node_start:node_end]
                node_index = node_names.index(node)
                if param_type == 'Prod':
                    index_dict[i-1] = ('Prod',node_index)
                else:
                    index_dict[i-1] = ('Deg',node_index)
            else:
                source_end = line.find('To',node_start)
                source = line[node_start:source_end]
                target_start = source_end + 2
                target_end = line.find('\t',target_start)
                target = line[target_start:target_end]
                source_index = node_names.index(source)
                target_index = node_names.index(target)
                if param_type == 'Trd':
                    index_dict[i-1] = ('Trd',source_index,target_index)
                else:
                    index_dict[i-1] = ('lam',source_index,target_index) 
    num_inputs = np.array([[len(network.inputs(i))] if len(network.inputs(i))>0 else [1] for i in range(N)])
    def racipe2dsgrn(param,return_deg = False):#, N=N, index_dict = index_dict, num_inputs = num_inputs):
        Prod = np.zeros([N,1])
        Trd = np.zeros([N,N])
        Deg = np.zeros([N,1])
        lam = np.zeros([N,N])
        K = len(param)
        for i, val in enumerate(reversed(param)):
            if num_param_names-1-i in index_dict:
                cur_entry = index_dict[num_param_names-1-i]
                if cur_entry[0] == 'Prod':
                    Prod[cur_entry[1],0] = val
                elif cur_entry[0] == 'Deg':
                    Deg[cur_entry[1],0] = val
                elif cur_entry[0] == 'Trd':
                    Trd[cur_entry[1],cur_entry[2]] = val
                else: #cur_entry[0] == 'lam'
                    if val > 1:
                        lam[cur_entry[1],cur_entry[2]] = 1/val
                    else: 
                        lam[cur_entry[1],cur_entry[2]] = val
        Prod = Prod**(1/num_inputs) #split production rate evenly across input nodes
        theta = Trd*Deg #multiply Trd[i,:] by d[i,0]
        L = lam*np.transpose(Prod) #multiply lam[:,i] by Prod[i,0]
        U = np.ones(lam.shape)*np.transpose(Prod) #U[j -> i] = U[j,i] = Prod[i,0]
        if return_deg:
            return L,U,theta,Deg
        else:
            return L,U,theta
    return racipe2dsgrn


network_file = "GRHL2 : (~ZEB)\nZEB : (~GRHL2)(~miR200)(ZEB)(SLUG)\nmiR200 : (~ZEB)(~SLUG)\nSLUG : (SLUG)"                                                                               
param_file = sys.argv[1]
param_name_file = sys.argv[2]
dsgrn_FPs(param_file, param_name_file, network_file)
