__author__ = 'fzhang'

from node import Node

(_ROOT,  _DEPTH,  _BREADTH) = range(3)

class Tree:
    
    def __init__(self):
        self.__nodes = {}
        
    @property
    def nodes(self):
        return self.__nodes
    '''''
    def add_node(self,  identifier,  parent = None):
        node = Node(identifier)
        self[identifier] = node
        if parent is not None:
            self[parent].add_child(identifier)
        return node'''
        
    def add_node(self, identifier, theta_L, theta_U, parent = None):
        #print "-"*40, "add a new node"
        node = Node(identifier, theta_L, theta_U);
        self[identifier] = node
        if parent is not None:
            self[parent].add_child(identifier)
            
        self[identifier].set_parent(parent)
        return node
        
    def return_parent(self,  identifier):
        return self[identifier].parent
    
    
    def search_leaves(self, identifier, Current_node, infimum, depth=_ROOT):
        
        children = self[identifier].children
        #print identifier, Current_node
        if len(children) == 0:
            #print self[identifier].MLBD, identifier
            if self[identifier].MLBD < infimum:
                Current_node = identifier
                infimum = self[identifier].MLBD
        depth += 1
        
        for child in children:
            node = self.search_leaves(child, Current_node, infimum, depth)
            if self[node].MLBD < self[Current_node].MLBD:
                Current_node = node
             
        return Current_node
    
    
    def display(self, identifier, depth=_ROOT):
        children = self[identifier].children
        if depth == _ROOT:
            print("{0}".format(identifier))
        else:
            print("t"*depth, "{0}".format(identifier)), self[identifier].MLBD, self[identifier].thetaB

        depth += 1
        for child in children:
            self.display(child, depth)  # recursive call
    
    def traverse(self,  identifier,  mode = _DEPTH):
        #Python generator
        yield identifier
        queue = self[identifier].children
        while queue:
            yield queue[0]
            expansion = self[queue[0]].children
            if mode is _DEPTH:
                queue = expansion + queue[1:] # Depth first
            elif mode is _BREADTH:
                queue = queue[1:] + expansion # Width first
    
    def __getitem__(self,  key):
        return self.__nodes[key]
        
    def __setitem__(self,  key,  item):
        self.__nodes[key] = item
        

    
