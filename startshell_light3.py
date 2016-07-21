import matplotlib.pyplot as plt
import numpy as np
import sys
from importlib import reload
home_dir='/home/ambroiseidoine/UCD/'
sys.path.append(home_dir+'tests')
import utilities as ut
import copula_analysis as ca
import copulaModels as mod
import vines


parameter={'offsets':[14,17],'date_range':('2010-07-01 18:00:00','2016-07-05 01:00:00'),'first_hour':(0,0)}
#titles=[{'type':'Wind','location':'total','kind':'forecast'},{'type':'Wind','location':'NP','kind':'error'},{'type':'Solar','location':'NP','kind':'error'}]
titles=[{'type':'Wind','location':'total','kind':'error'}]
copula=ca.create_copulaManager(titles,parameter)

def create_list_models(nb_max=3):
    mod_tp=[vines.cop2d_frank,vines.cop2d_gumbel,vines.cop2d_clayton,vines.cop2d_gaussian,vines.cop2d_student]
    length=len(mod_tp)
    res=mod_tp.copy()
    l=[0 for i in range(length)]
    while l is not None:
        if 1<sum(l)<=nb_max:
            temp=[mod_tp[i] for i in range(length) if l[i]==1]
            def f(x,temp=temp):
                return vines.WeightedCopula(x,temp)
            res.append( f )
        l=ut.table_increment(l)
    return res

list_models=[vines.cop2d_uniform,vines.cop2d_frank,vines.cop2d_gumbel,vines.cop2d_clayton,vines.cop2d_gaussian,vines.cop2d_student]

def vine3(unifs): 
	return vines.D_vine(unifs,list_models=create_list_models(nb_max=3))

def vine2(unifs): 
	return vines.D_vine(unifs,list_models=create_list_models(nb_max=2))

def vine1(unifs): 
	return vines.D_vine(unifs,list_models=create_list_models(nb_max=1))

final_list=[mod.cop_gaussian,mod.cop_student,
lambda x: vines.D_vine(x,list_models=[create_list_models(nb_max=3),list_models,list_models]),
lambda x: vines.D_vine(x,list_models=[create_list_models(nb_max=2),list_models,list_models]),
lambda x: vines.D_vine(x,list_models=[create_list_models(nb_max=1),list_models,list_models]),
lambda x: vines.D_vine(x,list_models=[create_list_models(nb_max=3),list_models,[vines.cop2d_uniform]]),
lambda x: vines.D_vine(x,list_models=[create_list_models(nb_max=2),list_models,[vines.cop2d_uniform]]),
lambda x: vines.D_vine(x,list_models=[create_list_models(nb_max=1),list_models,[vines.cop2d_uniform]]),
lambda x: mod.cop_independent_combination(x,[vine3,vine3],separators=[2]),
lambda x: mod.cop_independent_combination(x,[vine2,vine2],separators=[2]),
lambda x: mod.cop_independent_combination(x,[vine1,vine1],separators=[2])]

res=ca.test_models(copula,list_models=final_list,end_incr=900)



