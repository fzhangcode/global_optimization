__author__ = 'fzhang'

import numpy as np

class Node:
    '''def __init__(self,  identifier):
        self.__identifier = identifier
        self.__children = []
        self.__thetaBar = None
        self.__lamBar = None
        self.__muBar = None
        self.__xBar = None
        self.__SUBD = np.inf
        self.__MLBD = -np.inf
        
        self.__theta_L = None
        self.__theta_U = None
       # self.thetaBk =[]'''

    def __init__(self, identifier, theta_L, theta_U):
        self.__identifier = identifier
        self.__children = []
        self.__parent = -1
        self.__thetaB = None
        self.__thetaBar = None
        self.__lamBar = None
        self.__muBar = None
        self.__xBar = None
        self.__xOpt = None
        self.__SUBD = -np.inf
        self.__MLBD = np.inf
        self.__theta_L = theta_L
        self.__theta_U = theta_U
        self.__Depth = 0
        self.__g_flag = None
        self.__coefficients = None
        self.__lagrangian_coef = None
        
    
    def set_parameter(self,  lam,  theta, thetaB, mu,  xBar, xOpt, SUBD,  MLBD, g_flag,  coefficients,  lagrangian_coef):
        self.__lamBar = lam
        self.__thetaBar = theta
        self.__thetaB = thetaB
        self.__muBar = mu
        self.__xBar = xBar
        self.__xOpt = xOpt
        self.__SUBD = SUBD
        self.__MLBD = MLBD
        self.__g_flag = g_flag
        self.__coefficients = coefficients
        self.__lagrangian_coef = lagrangian_coef
        
    def set_parameters_qualifying_constraint(self, lam,  theta,  mu,  xBar,  SUBD,  g_flag,  coefficients):
        self.__lamBar = lam
        self.__thetaBar = theta
        self.__muBar = mu
        self.__xBar = xBar
        self.__SUBD = SUBD
        self.__g_flag = g_flag
        self.__coefficients = coefficients


    def set_parameters_thetaB(self,  thetaB,  xOpt,  MLBD,  lagrangian_coef):
        self.__thetaB = thetaB
        self.__xOpt = xOpt
        self.__MLBD = MLBD
        self.__lagrangian_coef = lagrangian_coef


    @property
    def lagrangian_coef(self):
        return self.__lagrangian_coef
    
    @property
    def Depth(self):
        return self.__Depth

    @property
    def thetaB(self):
        return self.__thetaB

    @property
    def theta_L(self):
        return self.__theta_L
	
    @property
    def theta_U(self):
        return self.__theta_U
        
    @property
    def parent(self):
        return self.__parent
        
    @property
    def thetaBar(self):
        return self.__thetaBar
        
    @property
    def lamBar(self):
        return self.__lamBar
        
    @property
    def muBar(self):
        return self.__muBar
        
    @property
    def xBar(self):
        return self.__xBar

    @property
    def xOpt(self):
        return self.__xOpt 
		
    @property
    def SUBD(self):
        return self.__SUBD
    
    @property
    def MLBD(self):
        return self.__MLBD
    
    @property
    def identifier(self):
        return self.__identifier
        
    @property
    def g_flag(self):
        return self.__g_flag
        
    @property
    def children(self):
        return self.__children
        
    @property
    def coefficients(self):
        return self.__coefficients

    
    def add_child(self,  identifier):
        self.__children.append(identifier)
        
    def set_parent(self,  parent):
        self.__parent = parent
        

    def set_theta_Upper_and_Lower(self,  theta_L,  theta_U):
         self.__theta_L = theta_L
         self.__theta_U = theta_U
