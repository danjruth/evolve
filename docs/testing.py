# -*- coding: utf-8 -*-
"""
Created on Sun May 14 17:30:14 2017

@author: danjr
"""

import algorithm as ga

rules_dict = {'a':{'min':2,'max':3},
              'b':{'min':0,'max':1},
              'c':{'min':-1,'max':1},
              'd':{'min':-5,'max':1}}
my_rules = ga.member_rules(rules_dict)
my_member = ga.member()
my_member.rand_init(my_rules)

def fun_to_opt(genes):
    return genes['a']**2 - genes['b']*genes['a'] + genes['c']*genes['d'] +  genes['b'] - genes['b']**2

my_pop = ga.population(my_rules,pop_size=100)

for i in range(50):

    my_pop.evaluate(fun_to_opt)
    my_pop.sort_via_df()
    #print(my_pop.df)
    
    print('Reinitializing...')
    my_pop.reinit_some()
    #print(my_pop.df)
    
    print('Mutating and crossover...')
    my_pop.mutate_pop()
    my_pop.cross_pop()
    #print(my_pop.df)
    
    print('Back to list...')
    my_pop.df_to_memberlist()
    print(my_pop.df.head())
    #for mem in my_pop.members:
        #print(mem.genes)
    
    print('TOP SCORE:')
    print(my_pop.df.loc[0,'score'])