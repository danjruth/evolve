# evolve
Run and visualize genetic algorithms in Python.

## Usage
First define the function that you want to maximize. It should take just a `dict` as an input, where each key/value corresponds to one of the parameters to optimize. For example:

```python
import numpy as np

def fun_to_opt(genes):
    
    a = float(genes['a'])
    b = float(genes['b'])
    c = float(genes['c'])
    d = float(genes['d'])
    
    score = 2*np.sin((1.5*a+b)*3) + 3*np.sin((a+c**2)/2) + (np.sin(c+d+b+a))**2 - ((a+c-b-d)/4)**2
    
    return score    
```

Then, create an instance of a `member_rules` class specifying the range of values each gene can take:

```python
import evolve.algorithm as ga

rules_dict = {'a':{'min':.6,'max':10},
              'b':{'min':5,'max':15},
              'c':{'min':1,'max':5},
              'd':{'min':2,'max':5}}
my_rules = ga.member_rules(rules_dict)
```

Finally, come up with an initial population of members, then initialize and run the algorithm:

```python
initial_pop = ga.population(my_rules,pop_size=60)
alg = ga.genetic_algorithm(fun_to_opt,max_iters=50)
best_genes = alg.run(initial_pop)
```

As the algorithm runs, it'll update a figure showing the results of each combination of genes it's tried over all the iterations (blue is a low score; yellow is a high score):
![Example visualization](/docs/example_output.PNG "Example visualization")

Once the algorithm's done, you can get the relative importance of each variable via a `scikit-learn` regression tree:

```python
imp = ga.get_variable_importance(alg.all_configs_df)
```
