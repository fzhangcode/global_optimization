__author__ = 'fzhang'

import numpy as np
import gurobipy as gb
import itertools as it
from itertools import repeat
import time
#from time import time

def add_qualifying_constraint(m,  coefficients, M,  K,  N,  thetaB,  g_flag,  t):
    ''' Add the qualifying constraints to model m.  The qualifying 
        constraint is formed by linearizing the Lagrangian with respect 
        to x about x0. Then linearizing that with respect to theta 
        about thetat.  '''

    qualifying_constraint = []
    x =  [[m.getVarByName("x_%d_%d" % (j,k)) for k in xrange(K)] for j in xrange(M)]

    
    for i in xrange(len(coefficients)):
        #calcuate the qualifying constraint formula
        g_expr = gb.LinExpr()
        g_expr.addConstant(coefficients[i][-1])
        for j in xrange(M*K):
            g_expr.add(x[j/K][j%K]* coefficients[i][j])
        qualifying_constraint.append(g_expr)
        #add the qualifying constraints
        if g_flag[i%K][i/K] != 0:
            ##################### Add constraints: g_expr #########################
            if thetaB[i%K][i/K] == 1:
                m.addConstr(g_expr <= np.spacing(1), name="qc%d_%d_%d" % (t,k,i))
            elif thetaB[i%K][i/K] == 0:
                m.addConstr(g_expr >= -np.spacing(1), name="qc%d_%d_%d" % (t,k,i))
            m.update() 
            

    #return qualifying constraints to calculate lagrange constraint
    return qualifying_constraint

def add_previous_lagrangian_constraint(m,  lagrangian_coefficient,  M,  K,  N,  t):
    '''Add the previous linearized lagrangian constraints to model m.'''
    
    #Extract the variables from the model.
    Q = m.getVarByName("Q")
    x = [[m.getVarByName("x_%d_%d" % (j,k)) for k in xrange(K)] for j in xrange(M)]
    
    #Initialize the constraint.
    L = gb.LinExpr()
    for j in xrange(M*K):
        L.add(x[j/K][j%K] * lagrangian_coefficient[j])
    L.addConstant(lagrangian_coefficient[-1])
    m.addConstr(Q >= L,  name = "Lagrangian_%d"%t)
    m.update()
    

def add_lagrangian_constraint(m,  qualifying_constraint,  xt,  thetat,  thetaB,  lam,  mu,  y,  t):
    ''' Add the linearized Lagrangian constraint to model m.
        The lagrangian is linearized with respect to theta0 and calculated at theta_i   '''

    (M,N) = np.shape(y)
    K = np.shape(thetat)[0]
    
    # Extract the variables from the model
    Q = m.getVarByName("Q")
    x = [[m.getVarByName("x_%d_%d" % (j,k)) for k in xrange(K)] for j in xrange(M)]
    
    # Initialize the constraint. The lagrangian constraint is linear in theta
    L = gb.LinExpr()
    x0_x = np.dot(xt.T, x)
    x0_x0 = np.dot(xt.T, xt)

    for i in xrange(N):
        L0 = np.dot(y[:,i].T, y[:,i])
        L0 += -np.dot(thetat[:,i], np.dot(x0_x0, thetat[:,i]))
        L0 += -2 * np.dot( np.dot(y[:,i].T, x), thetat[:,i])
        L0 += 2 * np.dot(thetat[:,i], np.dot( x0_x, thetat[:,i]))
        L.add(L0)
        
        L.addConstant(-lam[i] * (np.sum(thetat[:,i]) - 1))
        L.addConstant( -np.dot(mu[:,i].T, thetat[:,i]) )

        for j in xrange(K):
            L.add( qualifying_constraint[i*K + j]* (thetaB[j, i] - thetat[j, i]))

    m.addConstr(Q >= L, name = "Lagrangian_%d" % t)
    m.update()

    
    #Get the coefficients of lagrangian constraints
    num = L.size()
    coeff_constant = L.getConstant()
    coeff = np.zeros((M,  K))
    for index in xrange(num):
        var_name = L.getVar(index).VarName.encode('ascii',  'ignore')
        indexes = [ind for ind,  ltr in enumerate(var_name) if ltr == '_']
        i = int(var_name[indexes[0] + 1 : indexes[1]])
        j = int(var_name[indexes[1] + 1 : ])
        coeff[i][j] += L.getCoeff(index)
    
    lagrangian_coef = []
    for i in xrange(M):
        for j in xrange(K):
            lagrangian_coef.append(float(coeff[i][j]))
    lagrangian_coef.append(float(coeff_constant))
    return lagrangian_coef
    

def solve_master(tree, num_node, Current_node, g_flag, thetaB_list, SUBD,  coefficients, xBar, thetaOpt, lamOpt, muOpt, y,  cons, iteration,  pool=None):
    '''we solve the relaxed master problems based on thetaB_list, then select the infimum of all minimum values.
       Parameters: About the tree : tree, Current_node
                   About the subproblem: SUBD, xBar, thetaOpt, lamOpt, muOpt, y
                   About the boundary: theta_L, theta_U''' 
                   
    (M, N) = np.shape(y)
    K = np.shape(xBar[-1])[0]
    
    x_stor = None
    Q_stor = np.inf
    next_node = -1
    
    #store all the MLBD
    MLBD_stor = [] 

    #store all the MLBD
    if pool == None:
        tree.nodes[Current_node].set_parameters_qualifying_constraint(lamOpt,  thetaOpt, muOpt,  xBar, SUBD,  g_flag,  coefficients)
        #check whether the coefficients are already stored into the parents or not.
        
        print ('\n%d master problems are solving...'  %len(thetaB_list))

        for index in xrange(len(thetaB_list)):
            thetaB = thetaB_list[index].copy()
            status, objVal, xOpt,  thetaB, lagrangian_coefficient= solve_master_s(tree, Current_node, coefficients, thetaOpt, xBar, lamOpt, muOpt, thetaB.copy(), y, g_flag, cons)
            #print objVal, xOpt
            
            if status == 2 and objVal < SUBD - np.spacing(1):
                node = tree.add_node(num_node, 0, 1, Current_node)
                node.set_parameters_thetaB(thetaB,  xOpt, objVal, lagrangian_coefficient)
                MLBD_stor.append(objVal)
                if objVal < Q_stor:
                    Q_stor = objVal
                    next_node = num_node
                    x_stor = xOpt
                num_node = num_node + 1

    else:
        tree.nodes[Current_node].set_parameters_qualifying_constraint(lamOpt,  thetaOpt, muOpt,  xBar, SUBD,  g_flag,  coefficients)
        len_thetaB = len(thetaB_list)
        print ('\n%d master problems are solving...'  %len_thetaB)
        results = [pool.apply_async(solve_master_s,  args = (tree, Current_node, coefficients, thetaOpt, xBar, lamOpt, muOpt, thetaB.copy(), y, g_flag, cons)) for thetaB in thetaB_list]

        #put all the result into the tree.
        for p in results:
            #result = [status, objVal, xOpt, thetaB]
            result = p.get() 
            if result[0] == 2 and result[1] < SUBD - np.spacing(1):
                node = tree.add_node(num_node,  0,  1,  Current_node)
                node.set_parameters_thetaB(result[3],  result[2],  result[1], result[4])
                #node.set_parameter(lamOpt,  thetaOpt, result[3],  muOpt,  xBar,  result[2],  SUBD,  result[1], g_flag,  coefficients)
                MLBD_stor.append(result[1])
                if result[1] < Q_stor:
                    Q_stor = result[1]
                    next_node = num_node
                    x_stor =  result[2]
                num_node += 1

    return x_stor, Q_stor, next_node, num_node, MLBD_stor

def solve_master_s(tree, Current_node, coefficients, thetaOpt, xBar, lamOpt, muOpt, thetaB, y, g_flag, cons):
    '''Solve one master problem using Gurobipy'''
    
    (M,  K) = np.shape(xBar)
    (M,  N) = np.shape(y)

    ######################## Create a new model ##########################
    m = gb.Model("master")
    
    ##################### Set optimization parameters ####################   
    m.setParam('OutputFlag', False)

    ################### Create decision variables ########################
    x = [[0 for k in xrange(K)] for j in xrange(M)]
    for (j,k) in it.product(range(M), range(K)):
        x[j][k] = m.addVar(lb=-100, ub=100,
                           vtype=gb.GRB.CONTINUOUS, name="x_%d_%d" % (j,k) )
        
    # Create the slack variable for the objective    
    Q = m.addVar(lb=-gb.GRB.INFINITY, ub=gb.GRB.INFINITY, name="Q" )
    
    # Integrate decision variables
    m.update()
 
       ###################### Create the objective: min_x Q ##################
    m.setObjective(Q, gb.GRB.MINIMIZE)

    ####################### Add constraints ###############################
    # Add the sparsity constraint
    s = [[2 for k in xrange(K)] for j in xrange(M)]
    for (j,k) in it.product(range(M), range(K)):
        s[j][k] = m.addVar(lb=0.0, ub = gb.GRB.INFINITY, name="s_%d_%d" % (j,k) )
    m.update()
    
    for (j,k) in it.product(range(M), range(K)):
        m.addConstr( x[j][k] <= s[j][k], name="L1ub_s_%d_%d" % (j,k) )
        m.addConstr( x[j][k] >= -s[j][k], name="L1_lb_s_%d_%d" %(j,k) )

    obj = gb.LinExpr()
    for (j,k) in it.product(range(M), range(K)):
        obj.add(s[j][k])
    m.addConstr( obj <= cons, name="L1_sum" )    
    m.update()  
    

    #Compose the Lagrangian and qualifying constraints
    identifier = Current_node
    t = 0

    while identifier != 0:
        parent = tree.return_parent(identifier)
        qualifying_constraint= add_qualifying_constraint(m,  tree.nodes[parent].coefficients,  M,  K,  N, tree.nodes[identifier].thetaB,  tree.nodes[parent].g_flag,  t)   
        add_previous_lagrangian_constraint(m,  tree.nodes[identifier].lagrangian_coef,  M,  K,  N,  t)
        #lagrangian_coefficient = add_lagrangian_constraint(m, qualifying_constraint,   tree.nodes[parent].xBar, tree.nodes[parent].thetaBar,  tree.nodes[identifier].thetaB, tree.nodes[parent].lamBar, tree.nodes[parent].muBar, y, t)
        identifier = parent
        t += 1
    
    qualifying_constraint = add_qualifying_constraint(m,  coefficients, M,  K,  N, thetaB,  g_flag,  t)
    lagrangian_coefficient = add_lagrangian_constraint(m, qualifying_constraint,  xBar, thetaOpt, thetaB, lamOpt, muOpt, y, t)


    ############################ Optimize the master problem #####################    
    try:
        m.optimize()
    except gb.GurobiError as e:
        print e.message
        
    ############################ Check optimization results ######################
    if m.Status == gb.GRB.OPTIMAL:
        #Extract to optimal value and the optimum
        xOpt = np.empty((M, K))
        for (j, k) in it.product(range(M), range(K)):
            xOpt[j, k] = round(x[j][k].x,  4)

        return (m.Status, m.objVal, xOpt,  thetaB, lagrangian_coefficient)
    else:
        return(m.Status,  np.inf,  np.nan,  np.nan, lagrangian_coefficient)
