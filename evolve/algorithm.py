# -*- coding: utf-8 -*-
"""
Created on Sun May 14 17:02:02 2017

@author: danjr
"""

import random
import pandas as pd
import numpy as np
import copy
import matplotlib.pyplot as plt

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
        
    def random_values(self):
        '''
        Generate a set of genes given the rules.
        '''
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
        
    def force_adherence(self):
        '''
        Clamp gene values between the min/max (in case they were mutated 
        outside of this range).
        '''
        for gene in self.genes:
            rule = self.rules.gene_dict[gene]
            self.genes[gene] = max(self.genes[gene],rule['min'])
            self.genes[gene] = min(self.genes[gene],rule['max'])
        
    def evaluate(self,opt_fun):
        '''
        Score this member with the function to optimize.
        '''
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
        
    def plot_state(self,fig=None,all_configs_df=None):
        '''
        Make a matrix plot of all the members and encode the scores with color.
        
        If all_configs_df (pandas dataframe) is given, those members will be 
        plotted semi-transparently.
        '''
        if fig is None:
            fig = plt.figure(figsize=(9,8))
        
        n_genes = len(self.rules.genes)
        print(n_genes)
        
        fig.clear()
        
        ai = 0
        for yi,yg in enumerate(self.rules.genes):
            for xi,xg in enumerate(self.rules.genes):
                ai = ai+1
                ax = fig.add_subplot(int(n_genes),int(n_genes),int(ai))
                if all_configs_df is not None:
                    ax.scatter(all_configs_df[xg],all_configs_df[yg],c=all_configs_df['score'],alpha=0.3,s=20)
                ax.scatter(self.df[xg],self.df[yg],c=self.df['score'],s=20)
                ax.set_ylim([self.rules.gene_dict[yg]['min'],self.rules.gene_dict[yg]['max']])
                ax.set_xlim([self.rules.gene_dict[xg]['min'],self.rules.gene_dict[xg]['max']])
                
                if xi==0:
                    ax.set_ylabel(yg)
                else:
                    ax.yaxis.set_ticklabels([])
                    
                if yi==n_genes-1:
                    ax.set_xlabel(xg)
                else:
                    ax.xaxis.set_ticklabels([])
        
        plt.show()
        plt.pause(.1)
        plt.show()
        
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
        
    def reinit_some(self,reinit_frac=0.3):
        '''
        Delete and randomly reinitialize the worse-performing members.
        '''
        
        df = self.df.copy()
        
        for ix in np.arange(round(len(self.members)*(1-reinit_frac)),len(self.members)):
            
            dup_best = np.random.uniform() < 0.3
            if dup_best:
                new_member = copy.deepcopy(self.members[0])
            else:
                new_member = member()
                new_member.rand_init(self.rules)
            new_genes = new_member.genes
            df.loc[ix,:] = None
            for gene in list(new_genes.keys()):
                df.loc[ix,gene] = new_genes[gene]
            df.loc[ix,'has_changed'] = True
        
        self.df = df        
        
    '''
    For mutation and crossover, first just update the df storing the genes.
    Then, use df_to_memberlist to update the list of members with their new
    genes.
    '''
        
    def mutate_pop(self,mut_frac=0.8,mut_prob=0.6):
        '''
        Perform mutation on members of the population.
        '''
        
        new_df = self.df.copy()
        
        for ix_float in np.arange(1,round(len(self.members)*mut_frac)):
            
            ix = int(ix_float)
            mem = self.members[ix]
            
            if random.uniform(0,1) < mut_prob:
                
                num_mut = random.choice(range(len(mem.rules.genes))) + 1
                for gi in range(num_mut):
                
                    new_df.loc[ix,'has_changed'] = True
                    
                    # calculate the new value a gene will mutate to
                    gene_to_mut = random.choice(mem.rules.genes)
                    rel_mut_amount = ((float(ix) / float(len(self.members))))**1.0
                    mut_amount = rel_mut_amount * (mem.rules.gene_dict[gene_to_mut]['max'] - mem.rules.gene_dict[gene_to_mut]['min']) * np.random.normal()
                    new_val = new_df.loc[ix,gene_to_mut] + mut_amount
                    
                    #print(str(ix)+': '+str(mem.genes[gene_to_mut])+' ('+str(new_df.loc[ix,gene_to_mut])+') to '+str(new_val))
                    
                    # assign this value to the dataframe
                    new_df.loc[ix,gene_to_mut] = new_val
                
        self.df = new_df
        
    def cross_pop(self,cross_frac=0.8,cross_prob=0.6,cross_from_frac=0.4):
        '''
        Perform crossover on members of the population.
        '''
        
        new_df = self.df.copy()
        
        for ix_float in np.arange(1,round(len(self.members)*cross_frac)):
            
            ix = int(ix_float)
            mem = self.members[ix]
            
            if random.uniform(0,1) < cross_prob:
                
                new_df.loc[ix,'has_changed'] = True
                
                parent_ix = random.choice(range(int(len(self.members)*cross_from_frac)))
                
                # get the gene from a parent
                gene_to_cross = random.choice(mem.rules.genes)
                new_val = self.df.loc[parent_ix,gene_to_cross]
                
                # assign this value to the dataframe
                new_df.loc[ix,gene_to_cross] = new_val
                
        self.df = new_df
                
    def df_to_memberlist(self):
        
        new_list = list()        
        for m,mem in enumerate(self.members):

            for gene in mem.rules.genes:                    
                mem.genes[gene] = self.df.loc[m,gene]

            mem.force_adherence()
            new_list.append(mem)
            
        self.members = new_list
            

class genetic_algorithm:
    '''
    Class which defines and runs the genetic algorithm.
    '''    
    
    def __init__(self,opt_fun,max_iters = 20):
        self.opt_fun = opt_fun
        self.max_iters = max_iters
        self.converge_criteria = (0.000001,4) # opt solution does not change by x1 over x2 iterations
        
    def run(self,initial_population,ax=None):
        
        self.pop_list = [initial_population]
        
        pop = copy.deepcopy(initial_population)
        
        all_configs_df = pd.DataFrame()
        
        generation = 0
        has_converged = False
        fig = plt.figure(figsize=(12,8))
        while has_converged==False:
            
            generation = generation+1
            
            pop = copy.deepcopy(pop)
                                    
            pop.evaluate(self.opt_fun)
            pop.sort_via_df()
            all_configs_df = all_configs_df.append(pop.df.copy())
            
            pop.plot_state(fig=fig,all_configs_df=all_configs_df)
            
            print('Reinitializing...')
            pop.reinit_some()
            
            print('Mutation and crossover...')
            pop.mutate_pop()
            pop.cross_pop()
            
            print('Back to the list of members...')
            pop.df_to_memberlist()
            
            # show the top score--the actual top score is that stored at the
            # top of the dataframe, as the list entries have not been updated
            print('TOP SCORE:')
            print(pop.df.loc[0,:])
                        
            self.pop_list.append(copy.deepcopy(pop))
            
            has_converged = (generation>=self.max_iters)                
            
        self.optimal_result = pop.members[0]
        self.all_configs_df = all_configs_df
        
        return self.optimal_result.genes
        
def get_variable_importance(df):
    
    from sklearn.tree import DecisionTreeRegressor
    
    x = df.drop(['score','has_changed'],axis=1)
    print(x)
    y = df['score']
    
    clf = DecisionTreeRegressor(min_samples_leaf=5).fit(x.values,y)
    
    imp = {}
    for g,gene in enumerate(x.columns):
        imp[gene] = clf.feature_importances_[g]
    
    return imp