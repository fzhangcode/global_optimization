import numpy as np
import gurobipy as gb
import itertools as it
import matplotlib.pyplot as plt
from multiprocessing import Pool

import h5py as h5

from time import time

def solve_subproblem(y, x):
	'''Solve the GOP subproblem (dual) for a fixed value of x '''

	(M, N) = np.shape(y)
	K = np.shape(x)[1]

	# Create a new model
	m = gb.Model("sub")
	m.setParam('OutputFlag', False)
	
	# Create variables
	theta = [[0 for i in xrange(N)] for k in xrange(K)]
	for (k,i) in it.product(range(K), range(N)):
		theta[k][i] = m.addVar(lb = -gb.GRB.INFINITY, ub = gb.GRB.INFINITY,
									vtype = gb.GRB.CONTINUOUS,
									name = "theta_%d_%d" % (k,i) )
									
	#Integrate new variables
	m.update()
	
	#Construct the objective min_theta sum_i (y_i - x*theta_i)^2
	X2 = np.dot(x.T, x)
	
	obj = gb.QuadExpr()
	for i in xrange(N):
		yx = np.dot(y[:,i], x)
		
		# convert numpy.float64(0) to a native Python type.
		obj.addConstant(np.asscalar(np.dot(y[:,i].T, y[:,i]))) 
		
		for k1 in xrange(K):
			obj.addTerms(-2*yx[k1], theta[k1][i])
			for k2 in xrange(K):
				obj.addTerms(X2[k1][k2], theta[k1][i], theta[k2][i])
				
	m.setObjective(obj, gb.GRB.MINIMIZE)
	
	# Add constraint: sum_k theta_{ki} = 1 for all i
	for i in xrange(N):
		m.addConstr(gb.quicksum(theta[k][i] for k in xrange(K)) == 1)
	
	# Add constraint: theta_{ki} >= 0 for all (k,i)
	for (k,i) in it.product(range(K), range(N)):
		m.addConstr(theta[k][i] >= 0)
        
	# Optimize the model
	m.optimize()
	assert(m.getAttr('Status') == 2)
    
	# Package the optimizing value of theta
	thetaOpt = np.empty(np.shape(theta))
	for (k,i) in it.product(range(K), range(N)):
		thetaOpt[k,i] = theta[k][i].x
			
	# Package the Lagrange dual variables values
	lamOpt = np.empty(N)
	muOpt = np.empty(shape=(K,N))
	constrs = m.getConstrs()

	for i in xrange(N):
		lamOpt[i] = constrs[i].getAttr("pi")
		
        for (k,i) in it.product(range(K), range(N)):
            muOpt[k,i] = constrs[N+k*N+i].getAttr("pi")
            
	# Return the objective(Upper bound), theta, lambda, mu
	return m.objVal, thetaOpt, lamOpt, muOpt
	


if __name__ == '__main__':
	np.random.seed(1)
	
	#Define global optimization parameters
	#gb.setParam(gb.GRB.Param.Threads, NUM_THREADS)
	#pool = pool(processes = NUM_CORES)
	
	#Define the problem data
	y = np.array( [ (0.05, 0.2, -0.4, -0.85, -0.25),
					(-0.35, -0.4, -0.2, -0.05, -0.25),
					(0.6, 0.4, 1.2, 1.8, 1) ] );
	x_star = np.array( [ (-1, 0.5), 
						(0, -0.5),
						(2, 0)]);
	theta_star = np.array( [ (0.3, 0.2, 0.6, 0.9, 0.5), 
							(0.7, 0.8, 0.4, 0.1, 0.5) ] )
						 
	(M, N) = np.shape(y) # Dimesions of y are (M,N)
	
	#Define the problem parameters
	K = 2  # x(M*K) theta(K*N)
	MAXITER = 500
	e = 1.0
	
	#Initialize the problem parameters
	NBC = 2**(K*N)
	
	Qstor = []
	xstor = []
	
	SUBD = np.array([np.inf]) # supper boundary
	MLBD = np.array([-np.inf]) # lower boundary
	
	#Initialize the problem decision variables
	xBar = [x_star + np.random.normal(size = (M,K))]
	
	objOpt, thetaOpt, lamOpt, muOpt = solve_subproblem(y, xBar[-1])
	thetaBar = []
	lamBar = []
	print "Start Optimizing..."
	print objOpt
	print thetaOpt
	print "-"*50
	print lamOpt
	print "-"*50
	print muOpt
'''	
	i = 0
	while i < MAXITER:
	
	sttime = time()
	
		# Solve the subproblem
		objOpt, thetaOpt, lamOpt = solve_subproblem(y, xBar[-1])
		thetaBar.append(thetaOpt)
		lamBar.append(lamOpt)
		SUBD = np.append(SUBD, np.amin([objOpt, SUBD[-1]]))# the mini  
		
	sutime = time()	
		
		# Get the qualifying Lagrangian constraints from the previous iterations.
		thetaBK = []
		for k in xrange(i):
			lagrange_args = lagrange
'''	
	
"""	
if __name__ == '__main__':
	
    np.random.seed(1)
	
    #Define global optimization parameters
    #gb.setParam(gb.GRB.Param.Threads, NUM_THREADS)
    #pool = pool(processes = NUM_CORES)
	
    #Define the problem data						 
    M = 3
    N = 5
    S = 1000
            
    #Define the problem parameters
    K = 2  
	
    SUBD = np.array([np.inf]) # subproblem upper boundary
    MLBD = np.array([-np.inf]) # master problem lower boundary

    '''fig, ax = plt.subplots()
    y =  0.4*np.random.randint(10, size=(M, N)) - 0.4*np.random.randint(10, size=(M, N)) 
    print "y:\n", y
    
    # sample different x
    x = make_x(y, S, K)     
    
    # plot SUBD of different x
    for index in xrange(S):
        objOpt, thetaOpt, lamOpt, muOpt = solve_subproblem(y, x[index])
        plt.plot(index, objOpt, 'ro')        

    ylb = plt.ylabel('SUBD', fontsize = 24)        
    xlb = plt.xlabel('Random Samples', fontsize = 24)
    plt.setp(plt.gca().get_xticklabels(), fontsize=22)
    plt.setp(plt.gca().get_yticklabels(), fontsize=22)

    plt.show()'''
    
    fig, ax = plt.subplots()
    y = - 0.5*np.random.randint(10, size=(M, N)) #- 0.5*np.random.randint(10, size=(M, N)) 
    print "y:\n", y

    # sample different x
    x = np.zeros(shape = (M, K))
    index = 0
    for i1 in np.arange(-4, 5, 2):       # array([-4, -2,  0,  2,  4])  
        #print "i1", i1
        for i2 in np.arange(-4, 5, 2):
            #print "i2", i2
            for i3 in np.arange(-4, 5, 2):         
                #print "i3", i3
                for i4 in np.arange(-4, 5, 2): 
                    #print "i4", i4
                    for i5 in np.arange(-4, 5, 2):  
                        #print "i5", i5
                        for i6 in np.arange(-4, 5, 2):
                                #print "i6", i6
                                x[0][0] = i1
                                x[0][1] = i2
                                x[1][0] = i3
                                x[1][1] = i4
                                x[2][0] = i5
                                x[2][1] = i6
                             
                                objOpt, thetaOpt, lamOpt, muOpt = solve_subproblem(y, x)
                                #index = i1* 5^5 + i2* 5^4 + i3* 5^3 + i4* 5^2 + i5* 5 + i6
                                index += 1
                                
                                if index%1000 == 0:
                                    print index
                                    
                                plt.plot(index, objOpt, 'ro')    

                                                                                    
    ylb = plt.ylabel('SUBD', fontsize = 24)        
    xlb = plt.xlabel('Random Samples', fontsize = 24)
    plt.setp(plt.gca().get_xticklabels(), fontsize=22)
    plt.setp(plt.gca().get_yticklabels(), fontsize=22)

    plt.show()
    
"""    

