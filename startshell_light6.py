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
import results

parameter = {'offsets': [11, 14, 16, 18, 20],'date_range': ('2010-07-01 18:00:00', '2016-07-05 01:00:00'),
             'first_hour': (0,0)}
#titles=[{'type':'Wind','location':'total','kind':'forecast'},{'type':'Wind','location':'NP','kind':'error'},{'type':'Solar','location':'NP','kind':'error'}]
parameter2 = parameter.copy()
parameter2['offsets'] = [17, 14]
param_list = [parameter, parameter2]
titles = [{'type': 'Wind','location': 'total','kind': 'error'}]
#titles.append({'type':'Solar','location':'SP','kind':'error'})

copula = ca.create_copulaManager(titles, parameter)

gu = vines.cop2d_gumbel
ga = vines.cop2d_gaussian
st = vines.cop2d_student
fr = vines.cop2d_frank
un = vines.cop2d_uniform
cl = vines.cop2d_clayton

list1 = [gu, ga, st, fr, un, cl]

list2 = [gu, ga, st, fr, un, cl,
         lambda unifs: vines.WeightedCopula(unifs, [gu, ga, un], precise=False),
         lambda unifs: vines.WeightedCopula(unifs, [gu, st], precise=True),
         lambda unifs: vines.WeightedCopula(unifs, [gu, ga, st], precise=True),
         lambda unifs: vines.WeightedCopula(unifs, [gu, ga, st], precise=False),
         lambda unifs: vines.WeightedCopula(unifs, [gu, fr, st], precise=True),
         lambda unifs: vines.WeightedCopula(unifs, [gu, st, un], precise=False),
         lambda unifs: vines.WeightedCopula(unifs, [gu, fr, st], precise=False),
         lambda unifs: vines.WeightedCopula(unifs, [gu, ga], precise=False),
         lambda unifs: vines.WeightedCopula(unifs, [gu, un], precise=False),
         lambda unifs: vines.WeightedCopula(unifs, [gu, cl, ga], precise=False)]

list_models = [mod.cop_gaussian, mod.cop_student]
list_models.extend([lambda x, i=i: vines.D_vine(x, list_models=[[i], list1, list1]) for i in list2])
list_models.extend([lambda x, i=i: vines.D_vine(x, list_models=[[i], [un], [un]]) for i in list2])
list_models.append(lambda x: vines.D_vine(x, list_models=[list2, list1, list1]))

res = ca.test_models(copula, list_models=list_models)
results.compile_results(res, 'horse_races/more_dims/Solar_4D_NP_SP')
