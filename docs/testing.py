# -*- coding: utf-8 -*-
"""
Created on Sun May 14 17:30:14 2017

@author: danjr
"""

import algorithm as ga

rules_dict = {'a':{'min':0,'max':1},
              'b':{'min':5,'max':12}}
my_rules = ga.member_rules(rules_dict)
my_member = ga.member()
my_member.rand_init(my_rules)

print(my_rules.genes)
print(my_member.genes)

def fun_to_opt(genes):
    return genes['a'] + genes['b']

my_member.evaluate(fun_to_opt)
print(my_member.score)

my_pop = ga.population(my_rules)
my_pop.evaluate(fun_to_opt)
my_pop.sort_via_df()

print(my_pop.df)