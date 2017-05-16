# -*- coding: utf-8 -*-
"""
Created on Sun May 14 17:02:02 2017

@author: danjr
"""

import random
import pandas as pd
import numpy as np

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
        self.rules = rules
        
    def evaluate(self,opt_fun):
        [mem.evaluate(opt_fun) for mem in self.members]
        # now that it's evaluated, make it a tuple
        self.members = tuple(self.members)
        
    def sort_via_df(self):
        
        genes_list = list()
        [genes_list.append(mem.genes) for mem in self.members]
        
        # initialize a dataframe to store the genes and score of 
        df = pd.DataFrame(genes_list)
        for m,mem in enumerate(self.members):
            df.loc[m,'score'] = mem.score
        
        df['has_changed'] = False
        df.sort_values('score',ascending=False,inplace=True)
        df.index = range(len(df.index))
            
        self.df = df
        
    def reinit_some(self,reinit_frac=0.5):
        '''
        Delete and randomly reinitialize the worse-performing members.
        '''
        
        df = self.df.copy()
        
        for ix in np.arange(round(len(self.members)*(1-reinit_frac)),len(self.members)):
            new_member = member()
            new_member.rand_init(self.rules)
            new_genes = new_member.genes
            df.loc[ix,:] = None
            for gene in list(new_genes.keys()):
                df.loc[ix,gene] = new_genes[gene]
            df.loc[ix,'has_changed'] = True
        
        self.df = df        
        
    def mutate_pop(self,mut_frac=0.8,mut_prob=0.3,mut_num=0.2):
        '''
        Perform mutation on members of the population.
        '''
        
        new_df = self.df.copy()
        
        for ix_float in np.arange(1,round(len(self.members)*mut_frac)):
            
            ix = int(ix_float)
            mem = self.members[ix]
            
            if random.uniform(0,1) > mut_prob:
                
                # calculate the new value a gene will mutate to
                gene_to_mut = random.choice(mem.rules.genes)
                rel_mut_amount = (0.5*(float(ix) / float(len(self.members))))**1.5
                mut_amount = rel_mut_amount * (mem.rules.gene_dict[gene_to_mut]['max'] - mem.rules.gene_dict[gene_to_mut]['min']) * np.random.normal()
                new_val = mem.genes[gene_to_mut] + mut_amount
                
                # assign this value to both the member object and the dataframe
                #print('   mutating '+str(self.members[ix].genes[gene_to_mut])+' to '+str(new_val))
                self.members[ix].genes[gene_to_mut] = new_val
                new_df.loc[ix,gene_to_mut] = new_val
                
        self.df = new_df
        
    def cross_pop(self,cross_frac=0.8,cross_prob=0.3,cross_num=0.2):
        '''
        Perform crossover on members of the population.
        '''
        
        new_df = self.df.copy()
        
        for ix_float in np.arange(1,round(len(self.members)*cross_frac)):
            
            ix = int(ix_float)
            mem = self.members[ix]
            
            if random.uniform(0,1) > cross_prob:
                
                parent_ix = random.choice(range(len(self.members)))
                
                # calculate the new value a gene will mutate to
                gene_to_cross = random.choice(mem.rules.genes)
                new_val = self.members[parent_ix].genes[gene_to_cross]
                
                # assign this value to both the member object and the dataframe
                #print('   mutating '+str(self.members[ix].genes[gene_to_mut])+' to '+str(new_val))
                self.members[ix].genes[gene_to_cross] = new_val
                new_df.loc[ix,gene_to_cross] = new_val
                
        self.df = new_df
                
    def df_to_memberlist(self):
        
        new_list = list()
        
        for m,mem in enumerate(self.members):
            
            # update the member's genes if they have changed
            #if (self.df.loc[m,'has_changed'] == True):
            if True:
                #print(mem.genes)
                for gene in mem.rules.genes:                    
                    mem.genes[gene] = self.df.loc[m,gene]
                #print(mem.genes)
                #print('-------------')
                    
            new_list.append(mem)
            
        self.members = new_list         
            

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