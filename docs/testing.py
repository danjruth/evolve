# -*- coding: utf-8 -*-
"""
Created on Sun May 14 17:30:14 2017

@author: danjr
"""

import evolve.algorithm as ga
import numpy as np

# define some complicated function to optimize
def fun_to_opt(genes):
    
    a = genes['a']
    b = genes['b']
    c = genes['c']
    d = genes['d']    
    
    score = np.sin(a/d) + 3*np.sin(a*c) - (b-3)**2 * (d-a)
    
    return score

# define values our parameters to optimize can take
rules_dict = {'a':{'min':1,'max':2},
              'b':{'min':0,'max':6},
              'c':{'min':-5,'max':2},
              'd':{'min':0.5,'max':1}}
my_rules = ga.member_rules(rules_dict)

# come up with an initial population
my_pop = ga.population(my_rules,pop_size=50)

# create and run the genetic algorithm
alg = ga.genetic_algorithm(fun_to_opt)
best_genes = alg.run(my_pop)