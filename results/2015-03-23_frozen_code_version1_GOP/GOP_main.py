__author__ = 'fzhang'

import numpy as np
from multiprocessing import Pool
import h5py as h5
from time import time
import time
from tree import Tree

import parallel_cell_enumeration as p_cell

import parallel_masterprob as p_ma
import pre_process as pre
import subprob as sub


if __name__ == "__main__":
    
    # Set a random seed
    np.random.seed(19860522)
    
    # Define global optimization parameters
    #gb.setParam(gb.GRB.Param.Threads, NUM_THREADS)
    
    NUM_CORES = 10
    
    #Parallel
    pool = Pool(processes=NUM_CORES)
    pool = None
  
    ##################### Define the problem data. #################        
    '''M = 20
    K = 2
    N = 10
    # make theta_star
    x_star = np.random.normal(size=(M, K))  
    alpha = np.random.uniform(0,1,K)
    theta = np.random.dirichlet(alpha, N)
    theta_star = np.transpose(theta) 
    y = np.dot(x_star, theta_star)'''
    
    # illustrate data
    x_star = np.array([(0, -1)])
    theta_star = np.array([ (1, 0, 0.5),
                            (0, 1, 0.5)]) 
    y = np.dot(x_star, theta_star)
        
    # Define the sparsity of x
    # cons is used in the solve_master_s to constraint the sum of all the elements of |x_star|.  
    # It would be easier to set cons three times of the sum of all the absoute values for x.
    cons = np.sum(abs(x_star))    
    
    (M, N)  = np.shape(y)
    (M, K) = np.shape(x_star)
    
    #print M, K, N
    print "Optimal ||y-x*theta||_2^2 = %0.2f" % np.linalg.norm(y-np.dot(x_star, theta_star),2)**2
    
    # Define the problem parameters
    MAXITER = 200000000
    e = 0.01
    
    # Initialize the problem parameters   
    x_stor = []
    SUBD = np.inf
    MLBD = -np.inf

    # Initialize the problem decision variables
    #xBar = x_star + 0.1*np.random.normal(size=(M,K))

    # Randomly generate x between 0 and 1.
    xBar =  np.random.random_sample((M,K)) 
    print "Initial x is: %s" %xBar
        
    theta_U = np.zeros((K,N))
    theta_L = np.zeros((K,N))
    for i in xrange(N):
        for j in xrange(K):
            theta_U[j][i] = 1.0
            theta_L[j][i] = 0.0
    
    # record the optimal value 
    thetaBar = []
    lamBar = []
    muBar = []
    x_all = [xBar]
    MLBD_Bar = []
    SUBD_Bar = []
    # record all the MLBD of the generated nodes
    MLBD_all = []
    # record the chosen nodes with the lowest MLBD of each iteration
    nodes_all = []
    # record the time for each iteration
    time_iteration = []
	
    print "Start optimizing..."
    
    index=0
    global num_node
    num_node = 0
    current_node = 0
    
    #Claim a tree
    tree = Tree()
    node = tree.add_node(current_node, theta_L, theta_U)
    num_node = num_node + 1 

    start_all = time.clock()
    print x_all[-1]
    
    while index < MAXITER:
        start = time.clock()
        print "----------------------------iteration %d---------------------------" % index
        
        # Solve the subproblem
        objOpt, thetaOpt, lamOpt, muOpt = sub.solve_subproblem(y, xBar)   #(objOpt) upper bound

        thetaBar.append(thetaOpt)
        lamBar.append(lamOpt)
        muBar.append(muOpt)

        SUBD = np.amin([objOpt, SUBD])
                        
        print "THETA"
        print thetaOpt
        print "X"
        print xBar

        # preprocess: deal with duplicate hyperplanes, and remove all 0 coefficients hyperplanes.
        g_flag, replicated_marker, coefficients =  pre.pre_process(xBar, thetaOpt, lamOpt, muOpt, y)
        
        # print the flag and deplicate markers
        #print "g_flag", g_flag
        #print "replicated_marker", replicated_marker
        #print len(coefficients)
        #for co_index in xrange(len(coefficients)):
        #  print coefficients[co_index]

        # Get all the unique hyperplanes and save the coefficients of them.
        linker, unique_coefficients = pre.unique_coeff(g_flag,  replicated_marker,  coefficients,  M,  K,  N)

        # Set a threshold as the distance used in the cell enumeration
        distance = [np.spacing(1)]
        for i in xrange(len(unique_coefficients)):
            sum = 0.0
            sum += unique_coefficients[i][-1]
            for j in xrange(M*K):
                sum += xBar[j/K][j%K]* unique_coefficients[i][j]
            distance.append(np.fabs(sum))
        
        # Take the maximum of the 'distances' from the common point to all the hyperplanes
        # Make the threshold greater than or equal to np.spacing(1)
        threshold = max(distance)

        # Get the unique regions which are represented by thetaB_list (using parallel)
        pre_thetaB_list = p_cell.parallel_cell_numeration(unique_coefficients,  len(unique_coefficients),  M*K,  threshold,  pool)
        thetaB_list = pre.extend_back(pre_thetaB_list,  linker,  g_flag,  replicated_marker,  K,  N)
        print "\nthe length of thetaB_list is:",  len(thetaB_list)

        # Solve the master problems defined by the thetaB_list (using parallel)
        x_stor, Q_stor, next_node, num_node, MLBD_stor = p_ma.solve_master(tree, num_node, current_node, g_flag, thetaB_list, SUBD,  coefficients,  xBar, thetaOpt, lamOpt, muOpt, y,  cons, i, pool)  
        MLBD_all.extend(MLBD_stor)
        
        # Set the master problem lower bound and the next x value        
        current_node = tree.search_leaves(0, 0, np.inf)
        #print current_node, tree.nodes[current_node].parent, tree.nodes[current_node].thetaB

        nodes_all.append(current_node)
        xBar = tree.nodes[current_node].xOpt
        MLBD = tree.nodes[current_node].MLBD
        tree.nodes[current_node].MLBD = np.inf

        # Calculate the time for each iteration
        end = time.clock()

        # record the x_all, MLBD_Bar and SUBD_Bar of the chosen nodes with the lowest MLBD        
        x_all.append(xBar)
        MLBD_Bar.append(MLBD)
        SUBD_Bar.append(SUBD)
        
        time_iter = '%0.2f'  %(end - start)
        print ('\nTime used for this iteration: %0.2f'  %(end - start)) 
        time_iteration.append(time_iter)
  
        with h5.File('test.hdf5','w') as f:
            f.create_dataset('MLBD', data=MLBD_Bar)
            f.create_dataset('SUBD', data=SUBD_Bar)
            f.create_dataset('xOpt', data=x_all)
            f.create_dataset('thetaOpt', data=thetaBar)
            f.create_dataset('lamOpt', data=lamBar)
            f.create_dataset('muOpt', data=muBar)    
            f.create_dataset('MLBD_all_nodes', data=MLBD_all)
            f.create_dataset('selected_nodes', data=nodes_all)
            f.create_dataset('time', data=time_iteration)
                                                                    
        print('Current bounds: [%0.2f, %0.2f]' % (MLBD, SUBD))
        
        # add another convergence:
        #close = np.fabs(np.dot(x_all[-2], thetaOpt) - y) < 0.05
        
        # if SUBD - MLBD <= e or close.all():
        if SUBD - MLBD <= e: 
            print "\n==================Optimal x*theta================="
            print np.dot(x_all[-2], thetaOpt)
            print "======================Exact y====================="
            print y
            print "=====================Optimal x===================="
            print x_all[-2]
            print "======================Exact x====================="
            print x_star
            print "===================Optimal theta=================="
            print thetaOpt
            print "====================Exact theta==================="
            print theta_star

            index = MAXITER
            end_all = time.clock()
            print ('\nAll the iterations cost: %0.2f'  %(end_all - start_all)) 
        index += 1    
