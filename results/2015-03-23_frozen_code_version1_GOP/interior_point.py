__author__ = 'Chauncey'

import numpy as np
import gurobipy as gb
import math

def interior_point(coefficients, nCuts,  dim,  marker):
    '''  Get the interior point from the certain region described by marker. The returned point should be far from all the hyperplanes.
    coefficients: the coefficients of hyperplanes with the same order. For example, (a1, a2, a3, c) in the hyperplane (a1 * x1 + a2 * x2 + a3 * x3 + c = 0)
    nCuts: the number of the hyerplanes.
    Dim: the dimension of the space.
    marker: the sign of the hyperplane.
    index: the hyperplane we want to trim. Usually index = -1, which means no hyperplane will be trimed.'''
    
    #Create a new model
    m = gb.Model('interior_point')
    
    #Set parameters
    m.setParam('OutputFlag',  False)
    
    x = [0 for i in xrange(dim)]
    for i in xrange(dim):
        x[i] = m.addVar(lb = -1e7,  ub = 1e7,  vtype=gb.GRB.CONTINUOUS,  name = 'x_%d'%(i))
    
    obj = m.addVar(lb = 0.0,  vtype = gb.GRB.CONTINUOUS,  name = 'objective')
    m.update()
    
    for i in xrange(nCuts):
        g_expr = gb.LinExpr()
        g_expr.addConstant(coefficients[i][-1])
        for j in xrange(dim):
            g_expr.add(x[j] * coefficients[i][j])
        if marker[i] == 0:
            m.addConstr(-1 * g_expr + obj <= np.spacing(0),  name = 'qc_%d' %(i))
        elif marker[i] == 1:
            m.addConstr(g_expr + obj <= np.spacing(0),  name = 'qc_%d' %(i))
    m.update()
    
    #Create the objective : maximize obj
    m.setObjective(obj,  gb.GRB.MAXIMIZE)
    m.update()
    
    #Optimize the test problem.
    try:
        m.optimize()
    except gb.GurobiError as e:
        print e.message

    if m.Status == gb.GRB.OPTIMAL:
        xOpt = np.empty(dim)
        for i in xrange(dim):
            xOpt[i] = x[i].x
        return m.Status,  xOpt,  obj.x
    else:
        return m.Status,  np.nan,  np.nan

def interior_point_using_distance(coefficients, nCuts,  dim,  marker, weight, index = -1):
    ''' Get the interior point from the certain region described by marker. The returned point should be far from all the hyperplanes. 
    coefficients: the coefficients of hyperplanes with the same order. For example, (a1, a2, a3, c) in the hyperplane (a1 * x1 + a2 * x2 + a3 * x3 + c = 0)
    nCuts: the number of the hyerplanes.
    Dim: the dimension of the space.
    marker: the sign of the hyperplane.
    weight: calcualte the norm 2 of the coefficients except the constant of all hyperplanes  
    index: the hyperplane we want to trim. Usually index = -1, which means no hyperplane will be trimed.
    
    Here, we add the different weights for different hyperplanes when we calculate the distance.

    maximize z
     subject to  A*x + b + weight*z <= np.spacing(1), for all sign >=0
                       -1*(A*x + b) + weight*z <= np.spacing(1), for all sign <= 0
    '''
    
    #Create a new model
    m = gb.Model('Interior_point')
    
    #Set parameters
    m.setParam('OutputFlag',  False)
    
    x = [0 for i in xrange(dim)]
    for i in xrange(dim):
        x[i] = m.addVar(lb = -1e7,  ub = 1e7,  vtype=gb.GRB.CONTINUOUS,  name = 'x_%d'%(i))
    
    obj = m.addVar(lb = 0.0,  vtype = gb.GRB.CONTINUOUS,  name = 'objective')
    m.update()
    
    for i in xrange(nCuts):
        #hyperplane index is trimed, so there is no constraint for index hyperplane.
        if index != i:
            g_expr = gb.LinExpr()
            g_expr.addConstant(coefficients[i][-1])
            for j in xrange(dim):
                g_expr.add(x[j] * coefficients[i][j])
            if marker[i] == 0:
                m.addConstr(-1 * g_expr + weight[i]*obj <= np.spacing(0),  name = 'qc_%d' %(i))
            elif marker[i] == 1:
                m.addConstr(g_expr + weight[i]*obj <= np.spacing(0),  name = 'qc_%d' %(i))
    m.update()
    
    #Create the objective : maximize obj
    m.setObjective(obj,  gb.GRB.MAXIMIZE)
    m.update()
    
    #Optimize the test problem.
    try:
        m.optimize()
    except gb.GurobiError as e:
        print e.message

    if m.Status == gb.GRB.OPTIMAL:
        xOpt = np.empty(dim)
        for i in xrange(dim):
            xOpt[i] = x[i].x
        return m.Status,  xOpt,  obj.x
    else:
        return m.Status,  np.nan,  np.nan

if __name__ == "__main__":
    #test 1
    coefficients = np.array([[1, 0, 0], [0, 1, 0], [2, -1, 0]])
    nCuts = 3
    dim = 2
    marker = [0, 0, 0]
    
    #test 2
    #coefficients = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [1, -1, 0, 0]])
    #nCuts = 4
    #dim = 3
    #marker = [0, 0, 0, 1]
    
    #test 3
    coefficients = np.array([[5, -1, 0], [0, 1, 0], [2, -1, 0]])
    nCuts = 3
    dim = 2
    marker = [0, 0, 0]   
    
    #Status,  Xopt,  obj = interior_point(coefficients, nCuts,  dim,  marker)
    Status,  Xopt,  obj= interior_point_using_distance(coefficients, nCuts,  dim,  marker)
    print Status,  Xopt,  obj
    
    for i in xrange(nCuts):
        sum = 0.0
        norminator = 0.0
        sum += coefficients[i][-1]
        norminator = coefficients[i][-1]**2
        for j in xrange(dim):
            sum += Xopt[j] * coefficients[i][j]
            norminator += coefficients[i][j]**2
        #print sum
        print sum/math.sqrt(norminator)
    
