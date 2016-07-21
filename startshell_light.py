import matplotlib.pyplot as plt
import numpy as np
import sys
from importlib import reload
home_dir='C:\\users\\sabrina\\documents\\research\\code for real user\\'
sys.path.append(home_dir+'tests')
import utilities as ut
import copula_analysis as ca
import copulaModels as mod
import vines
from itertools import combinations
import results

parameter={'offsets':[11,14,17],'date_range':('2010-07-01 18:00:00','2016-07-05 01:00:00'),'first_hour':(0,0)}
#titles=[{'type':'Wind','location':'total','kind':'forecast'},{'type':'Wind','location':'NP','kind':'error'},{'type':'Solar','location':'NP','kind':'error'}]
titles=[{'type':'Solar','location':'total','kind':'error'}]
copula=ca.create_copulaManager(titles,parameter)
list_models=[vines.cop2d_gumbel, vines.cop2d_clayton, vines.cop2d_frank, vines.cop2d_gaussian, vines.cop2d_student, vines.cop2d_uniform]
combs2 = [lambda unifs, comb=comb: vines.WeightedCopula(unifs, comb, precise=True) for comb in combinations(list_models, 2)]
combs2b = [lambda unifs, comb=comb: vines.WeightedCopula(unifs, comb, precise=False) for comb in combinations(list_models, 2)]
combs3 = [lambda unifs, comb=comb: vines.WeightedCopula(unifs, comb, precise=True) for comb in combinations(list_models, 3)]
combs3b = [lambda unifs, comb=comb: vines.WeightedCopula(unifs, comb, precise=False) for comb in combinations(list_models, 3)]
#models = [mod.cop_gaussian, mod.cop_student, lambda unifs: vines.D_vine(unifs, list_models=list_models),
#          lambda unifs: vines.D_vine(unifs, list_models=[list_models + combs2, list_models]),
#          lambda unifs: vines.D_vine(unifs, list_models=[list_models + combs2 + combs3, list_models])]
models = list_models + combs2 + combs2b + combs3 + combs3b
res=ca.test_models(copula,list_models=models)
#results.compile_results(res)
