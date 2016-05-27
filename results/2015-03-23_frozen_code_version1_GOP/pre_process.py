__author__ = 'fzhang'

import numpy as np
import gurobipy as gb
import itertools as it
import sympy

def get_the_coefficient(g_expr,  M,  K,  m):
    '''  Get the coefficients of one qualifying constraint.
    More efficiently! We observe that VarName is named by "x_i_j". 
    Therefore, we can extract i and j from the name.'''
    num = g_expr.size()
    x = [[m.getVarByName("x_%d_%d" % (j,k)) for k in xrange(K)] for j in xrange(M)]
    
    coeff = np.zeros((M,  K))
    coeff_constant = g_expr.getConstant()
    for index in xrange(num):
        var_name = g_expr.getVar(index).VarName.encode('ascii',  'ignore')
        indexes = [ind for ind,  ltr in enumerate(var_name) if ltr == '_']
        i = int(var_name[indexes[0] + 1 : indexes[1]])
        j = int(var_name[indexes[1] + 1 : ])
        coeff[i][j] += g_expr.getCoeff(index)
    
    coefficient = []
    for i in xrange(M):
        for j in xrange(K):
            coefficient.append(float(coeff[i][j]))#set the precision 
    coefficient.append(float(coeff_constant))#set the precision 
    return coefficient


def normalize_coefficient(coef):
    '''Normalization is not used in the current version.'''
    '''Normalize the coefficients of Ab by the largest one of A.'''
    #get the largest magnitude of A
    max_coef = max(coef[:-1])
    min_coef = min(coef[:-1])
    normalization = 0.0
    if np.fabs(max_coef) > np.fabs(min_coef):
        normalization = np.fabs(max_coef)
    else:
        normalization = np.fabs(min_coef)
    assert(normalization != 0.0), "normalization should not equal to zero"
    #divide the largest magnitude 
    normalized_coef = []
    for i in xrange(len(coef)):
        normalized_coef.append(float(coef[i]/normalization))
    return normalized_coef
    

def check_line(coeff,  M,  K):
    #Since we used the truncation, here we compare the coeff with 0.0001 directly.
    for index in xrange(M*K + 1):
        if np.fabs(coeff[index]) >= 0.0001:
            return False
    return True


def replicate_line_with_threshold(g_expr,  g_expr_c,  M,  K,  threshold):
    '''Compare whether two lines are the same lines.'''
    coeffi_1 = np.zeros((M,  K))
    coeffi_2 = np.zeros((M,  K))
    coeff1_constant = g_expr[-1]
    coeff2_constant = g_expr_c[-1]
    #calculate the sum of all coefficients
    sum_coeff1 = coeff1_constant
    sum_coeff2 = coeff2_constant
    for i in xrange(M):
        for j in xrange(K):
            coeffi_1[i][j] = g_expr[i*K + j]
            coeffi_2[i][j] = g_expr_c[i*K + j]
            sum_coeff1 += coeffi_1[i][j]
            sum_coeff2 += coeffi_2[i][j]
    
    #check constant
    if np.fabs(sum_coeff1 * coeff2_constant - sum_coeff2 * coeff1_constant) > threshold:
        return 0
    #check all other coefficients
    for i in xrange(M):
        for j in xrange(K):
            if np.fabs(sum_coeff1 * coeffi_2[i][j] - sum_coeff2 * coeffi_1[i][j]) > threshold:
                #sum_coeff1 * coeffi_2[i][j] != sum_coeff2 * coeffi_1[i][j]
                return 0 
    #1 represents the two lines are the same equation
    return 1
    
    
def pre_process(xBar,  thetaOpt,  lamOpt,  muOpt,  y):
    '''Preprocessing step'''
    
    (M,  K) = np.shape(xBar)
    (K,  N) = np.shape(thetaOpt)
    #qualify_flag: mark whether the coefficients of the qualifying constraint is all zero(0) or not(1)
    qualify_flag = np.zeros((K, N))
    #Coefficients: store the coefficients for KN lines. The size if KN * (MK + 1)
    coefficients = []
    
    #Create a new model without optimization
    m = gb.Model("Just_get_the_coefficients")
    
    #Add variables
    x = [[0 for k in xrange(K)] for j in xrange(M)]
    for (j,k) in it.product(range(M), range(K)):
        x[j][k] = m.addVar(lb=-gb.GRB.INFINITY, ub=gb.GRB.INFINITY,
                           vtype=gb.GRB.CONTINUOUS, name="x_%d_%d" % (j,k) )
    m.update()
    
    #get the coefficients
    x0_x = np.dot(xBar.T, x)
    x0_x0 = np.dot(xBar.T, xBar)
    x0_x_x = x0_x + x0_x.T
    for i in xrange(N):
        for k in xrange(K):
            g_expr = gb.LinExpr()
            g_expr.addConstant(-2* np.dot(x0_x0[:, k], thetaOpt[:,i]) - lamOpt[i] - muOpt[k,i])
            S = 0
            for j in xrange(M):
                S += x[j][k] * float(y[j,i])

            g_expr.add(-2*S)
            g_expr.add(2 * np.dot(x0_x_x[:,k], thetaOpt[:,i]))
            coef = get_the_coefficient(g_expr, M, K, m)
            # Check and mark if a line has all the 0 coefficients 
            Flag = check_line(coef, M, K)
            if Flag:
                qualify_flag[k][i] = 0 # zero line
                coefficients.append(coef)
            else:
                qualify_flag[k][i] = 1 #not zero line
                coefficients.append(coef)
                
    ########################## Check a pair of two lines are the same or not. ###########################
    # Save the index of the duplicate qualifying constraints in the replicated_marker.
    replicated_marker = [[] for x in xrange(N*K)]
    for i in xrange(len(coefficients)):
        if qualify_flag[i%K][i/K] != 0:
            for j in xrange(i+1, len(coefficients)):
                if qualify_flag[j%K][j/K] != 0:
                    if replicate_line_with_threshold(coefficients[i], coefficients[j], M, K, 1e-4):
                        replicated_marker[i].extend([(j%K, j/K)])
                        replicated_marker[j].extend([(i%K, i/K)])
    return qualify_flag, replicated_marker, coefficients


def unique_coeff(g_flag, replicated_marker, coefficients, M, K, N):
    '''Get the list of unique coefficients'''
    
    #Get the marker of the unqiue coefficient
    marker = [0 for i in xrange(len(coefficients))]
    for i in xrange(len(coefficients)):
        if marker[i] == 0:
            if g_flag[i%K][i/K] == 0:
                marker[i] = -1
            else:
                ## replicated_marker[i] == empty
                if not replicated_marker[i]:
                    marker[i] = 0
                else:
                    for j in xrange(len(replicated_marker[i])):
                        marker[replicated_marker[i][j][0] + replicated_marker[i][j][1]*K] = -1
    #Get the unique coefficient based on the marker
    linker = []
    Unique_coefficient = []
    for i in xrange(len(marker)):
        if marker[i] == 0:
            linker.append(i)
            Unique_coefficient.append(coefficients[i])
    return linker, np.array(Unique_coefficient)


def extend_back(thetaB_list, linker, g_flag, replicated_marker, K, N):
    '''change the formula of thetaB_list to the matrix.'''
    
    extend_list = []
    for i in xrange(len(thetaB_list)):
        thetaB = np.zeros( (K, N))
        #set zero lines     
        for n in xrange(N):
            for k in xrange(K):
                thetaB[k][n] = -2
                if g_flag[k][n] == 0:
                    thetaB[k][n] = -1
        #set value
        for j in xrange(len(linker)):
            row = linker[j]/K
            col = linker[j]%K
            thetaB[col][row] = int(thetaB_list[i][j])
        for j in xrange(len(replicated_marker)):
            if replicated_marker[j] != []:
                if thetaB[j%K][j/K] == -2:
                    for n in xrange(len(replicated_marker[j])):
                        value = thetaB[replicated_marker[j][n][0]][replicated_marker[j][n][1]]
                        if value != -2:
                           thetaB[j%K][j/K] = value
        extend_list.append(thetaB)
    return extend_list

if __name__ == "__main__":
    pass
