# -*- coding: utf-8 -*-
"""
Created on Sun May 14 17:30:14 2017

@author: danjr
"""

import evolve.algorithm as ga
import numpy as np

def fun_to_opt(genes):
    
    a = float(genes['a'])
    b = float(genes['b'])
    c = float(genes['c'])
    d = float(genes['d'])
    
    score = 2*np.sin((1.5*a+b)*3) + 3*np.sin((a+c**2)/2) + (np.sin(c+d+b+a))**2 - ((a+c-b-d)/4)**2 - 0.5*(d-4)**2
    
    return score

# define values our parameters to optimize can take
rules_dict = {'a':{'min':.6,'max':10},
              'b':{'min':5,'max':15},
              'c':{'min':1,'max':5},
              'd':{'min':2,'max':5}}
my_rules = ga.member_rules(rules_dict)

# come up with an initial population
my_pop = ga.population(my_rules,pop_size=60)

# create and run the genetic algorithm
alg = ga.genetic_algorithm(fun_to_opt,max_iters=20)
best_genes = alg.run(my_pop)