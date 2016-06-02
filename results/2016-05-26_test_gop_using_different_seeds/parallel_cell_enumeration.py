__author__ = 'Chauncey'

#import Queue
import interior_point as IP
import copy
import numpy as np
from time import time
import time
from multiprocessing import Pool
import itertools

def reflection(thetaB):
    ''' Get the thetaB for the reflected region: 0 ->1, 1->0
    thetaB: the sign of hyperplane '''
    
    ref_thetaB = []
    for i in xrange(len(thetaB)):
        ref_thetaB.append(thetaB[i])
        if thetaB[i] == 1:
            ref_thetaB[i] = 0
        elif thetaB[i] == 0:
            ref_thetaB[i] = 1
    return ref_thetaB

def initialize_region(coefficients,  nCuts, dim,  threshold):
    ''' Calculate the first feasible region
    coefficients: the coefficients of all hyperplanes.
    nCuts: the number of hyperplanes.
    Dim: the dimension of space.
    threshold: the tolerance of the distance from the point to hyperplane '''
    
    #test the region one by one unitl one feasible is found
    NBC = pow(2,  nCuts)
    for index in xrange(NBC):
        thetaBstr = "{0:b}".format(index).zfill(nCuts)
        marker = list(thetaBstr)
        result = [int(i) for i in marker]
        status,  xOpt,  obj = IP.interior_point(coefficients,  nCuts,  dim, result)
        if status == 2 and obj >= threshold:
            return result

def check_candidate_hyperplanes(coefficients,  candidate_hyperplanes, nCuts,  dim,  region,  threshold, pool = None):
    ''' Check if each candidate hyperplane is the strict hyperplane.
    coefficients: the coefficients of all hyperplanes.
    candidate_hyperplanes: check this hyperplane is strict or redundant.
    nCuts: the number of hyperplanes.
    Dim: the dimension of space.
    region: the sign of the hyperplanes. Adjacency regions of this region will be return.
    threshold: the tolerance of the distance from the point to a hyperplane '''
    
    #non-redundant hyperplane
    strict_hyperplanes = []
    
    if pool == None:
        #calcualte the norm 2 of the coefficients except the constant of all hyperplanes 
        """ norm_2 is the input of function interior_point_using_distance
        norm_2 = []
        for i in xrange(nCuts):
            weight = 0.0
            for j in xrange(dim):
                weight += coefficients[i][j]**2
            weight = math.sqrt(weight)
            norm_2.append(weight)
        """

        #check the candidate hyperplanes
        for can_hp in candidate_hyperplanes:
            new_region = flip(region, can_hp)
            status,  xOpt,  obj = IP.interior_point(coefficients,  nCuts,  dim,  new_region)
            #status,  xOpt,  obj = IP.interior_point_using_distance(coefficients,  nCuts,  dim,  new_region, norm_2)
            if status == 2 and obj >= threshold:
                strict_hyperplanes.append(can_hp)
    else:
        #Parallel
        results = [pool.apply_async(check_interior_point,  args = (coefficients,  nCuts,  dim,  region, can_hp,  threshold)) for can_hp in candidate_hyperplanes]
        for p in results:
                #result = [can_hp, True or false]
                result = p.get() 
                if result[1]:
                    strict_hyperplanes.append(result[0])
    return strict_hyperplanes

def check_interior_point(coefficients,  nCuts,  dim,  region, can_hp,  threshold):
    ''' In order to parallel, we use this function to check the new region exist or not'''
    new_region = flip(region,  can_hp)
    status,  xOpt,  obj = IP.interior_point(coefficients,  nCuts,  dim,  new_region)
    if status == 2 and obj >= threshold:
        return can_hp,  True
    else:
        return can_hp,  False
    

def cal_adjacency_regions(coefficients,  nCuts,  dim,  region,  certain_hyperplane, threshold,  pool =None,  closedset=[]):
    '''   Calculate the adjacency regions of the given region. This function icludes two parts:
    1) calculate the strict hyperplanes.
    2) get the adjacency regions based on the strict hyperplanes.
    
    coefficients: the coefficients of all hyperplanes.
    nCuts: the number of hyperplanes.
    Dim: the dimension of space.
    region: the sign of the hyperplanes. Adjacency regions of this region will be return.
    certain_hyperplane: shows which hyperplane we are considering.
    threshold: the tolerance of the distance from the point to a hyperplane'''
    
    #the calculated adjacency regions
    adj_regions = []
    
    #put all the hyperplanes into candidate_hyperplanes
    candidate_hyperplanes = list(xrange(nCuts))
    
    #test the other candidate hyperplanes.
    strict_hyperplanes = check_candidate_hyperplanes(coefficients,  candidate_hyperplanes, nCuts,  dim,  region,  threshold, pool)   
    #print "strict_hyperplanes:",  strict_hyperplanes
    
    #flip the strict hyperplane to get the adjacency regions.
    for i in strict_hyperplanes:
        if i != certain_hyperplane:
            flip_region = flip(region, i)
            if flip_region not in closedset:
                adj_regions.append(flip_region)
        
    #return the adjacency regions
    return adj_regions

def flip(region,  index):
    '''  flip the element of index. 1 -> 0,  0 -> 1
    region: the sign of hyperplanes
    index: the position should be fliped. '''
    
    flip_region = copy.copy(region)
    if flip_region[index] == 0:
        flip_region[index] = 1
    elif flip_region[index] == 1:
        flip_region[index] = 0
    else:
        assert(1 == 0),  'No 0 and 1 appears in the thetaB_list.'
    return flip_region

def parallel_cell_numeration(coefficients,  nCuts,  dim,  threshold, pool = None,  certain_hyperplane = 0):
    '''  Get the thetaB_list of all the regions.
    coefficients: the coefficients of all hyperplanes.
    nCuts: the number of hyperplanes.
    Dim: the dimension of space.
    threshold: the tolerance of the distance from the point to hyperplane
    certain_hyperplane: assume the sign of this hyperplane does not change. '''
    
    #sub_NUM_CORES = 10
    #Parallel
    #sub_pool = Pool(processes=sub_NUM_CORES)
    
    thetaB_list = []
    #the structure to manage the process of all the function.
    #openset = Queue.Queue()
    openset = []
    closedset = []

    Initialized_region = initialize_region(coefficients,  nCuts,  dim,  threshold)
    
    #openset.put(Initialized_region)
    openset.append(Initialized_region)
    print ('\nfinding adjacent regions...') 
	
    # non-parallel
    if pool == None:
        while len(openset) != 0:
            #get the region and delete it from the openset.
            region = openset.pop(0) 
            closedset.append(region)

            #calculate the adjacency regions
            adj_regions = cal_adjacency_regions(coefficients,  nCuts,  dim,  region, certain_hyperplane,  threshold)#, sub_pool)

            for region in adj_regions:
                #if region never appears.
                if region not in openset and region not in closedset:
                    openset.append(region)
    # parallel
    else:
        start = time.clock()
        while len(openset) != 0:
            closedset.extend(openset)
            results = [pool.apply_async(cal_adjacency_regions,  args = (coefficients,  nCuts,  dim,  region, certain_hyperplane,  threshold,  None, closedset)) for region in openset]

            #set openset is empty
            openset = []
            #get the result
            for p in results:
                openset.extend(p.get())
            #get the unique list
            openset.sort()
            openset = list(openset  for openset,  _ in itertools.groupby(openset))
    
        end = time.clock()
        time_adj_region = end - start
        print ('\nFinish finding adjacent regions using: %0.2f'  %(time_adj_region))  

    #Reflection 
    for region in closedset:
        thetaB_list.append(region)
        ref_thetaB = reflection(region)
        thetaB_list.append(ref_thetaB)

    #length = len(thetaB_list)
    return thetaB_list


if __name__ == "__main__":
    """
    #test 1
    coefficients = np.array([[1, 0, 0], [0, 1, 0], [2, -1, 0]])
    nCuts = 3
    dim = 2
    thetaB_list = parallel_cell_numeration(coefficients,  nCuts,  dim)
    print thetaB_list,  len(thetaB_list)
    """

    #test 2
    coefficients = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [1, -1, 0, 0]])
    nCuts = 4
    dim = 3

    """
    #test 3
    coefficients = np.array([[5, -1, 0], [0, 1, 0], [2, -1, 0]])
    nCuts = 3
    dim = 2
    """    
    #test4
    #coefficients = np.array([[1, 0, 0], [5, -1, 0], [2, -1, 0]])
    #nCuts = 3
    #dim = 2
    
    #test5
    coefficients = np.array([[ 0,           0,           1,          1,         -2,        ], 
                             [ 1,         0,           0,           0,           0,        ], 
                             [ 0,           1.0/3,   1,           1.0/3,  -2.0/3], 
                             [ 0,          0.5,          1,           0,          0,         ]])
    common_point = np.array([[ 0,   0], [ 0,   2]])
   
    nCuts = 4
    dim = 4
    thetaB_list = parallel_cell_numeration(coefficients,  nCuts,  dim,  np.spacing(1))
    print thetaB_list,  len(thetaB_list)
    pass
    
