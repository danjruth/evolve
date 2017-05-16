# -*- coding: utf-8 -*-
"""
Created on Sun May 14 17:02:02 2017

@author: danjr
"""

import random
import pandas as pd

class member_rules:
    '''
    Class to define the rules specifying what values genes can take.
    '''
    
    def __init__(self,gene_dict):
        '''
        gene_dict: a dict in which each gene is one entry. Each entry is
            another dict, with min and max specified            
        '''
        
        self.gene_dict = gene_dict
        
        self.genes = list(gene_dict.keys())
        
    def check_adherence(self,member):
        '''
        Ensure that a member's genes adhere to the rules.
        '''
        for gene in self.genes:
            rule = self.gene_dict[gene]
            if (rule['max'] > member.genes[gene]) or (rule['min'] < member.genes[gene]):
                return False
        else:
            return True
        
    def random_values(self):
        gene_vals = {}
        for gene in self.genes:
            rule = self.gene_dict[gene]
            gene_vals[gene] = random.uniform(rule['min'],rule['max'])
        return gene_vals
    
class member:
    '''
    Class for each individual member of a population.
    '''
    
    def __init__(self):
        self.genes = None
        self.rules = None
        self.score = None
        
    def rand_init(self,rules):
        '''
        Randomly give genes to this member, following the rules as specified in
        an instance of a member_rules class.
        '''
        self.rules = rules
        self.genes = self.rules.random_values()
        
    def evaluate(self,opt_fun):
        self.score = opt_fun(self.genes)
        return self.score

class population:
    '''
    Class to contain population members and evaluate a population.
    '''
    
    def __init__(self,rules,pop_size=20):
        members = [member() for m in range(pop_size)]
        [members[m].rand_init(rules) for m in range(pop_size)]
        self.members = members
        
    def evaluate(self,opt_fun):
        [mem.evaluate(opt_fun) for mem in self.members]
        
    def sort_via_df(self):
        
        genes_list = list()
        [genes_list.append(mem.genes) for mem in self.members]
        
        # initialize a dataframe to store the genes and score of 
        df = pd.DataFrame(genes_list)
        for m,mem in enumerate(self.members):
            df.loc[m,'score'] = mem.score
        
        df.sort_values('score',ascending=False,inplace=True)
        df.index = range(len(df.index))
            
        self.df = df
        
    def df_to_memberlist(self):
            

class genetic_algorithm:
    '''
    Class which defines and runs the genetic algorithm.
    '''    
    
    def __init__(self,opt_fun):
        self.opt_fun = opt_fun
        self.max_iters = 200
        self.converge_criteria = (0.000001,4) # opt solution does not change by x1 over x2 iterations
        
    def run(self,initial_population=None,viz_update=None):
        None