### start shell

import matplotlib.pyplot as plt
import numpy as np
import sys
from importlib import reload
home_dir='C:\\Users\\Sabrina\\Documents\\Research\\code for real user\\'
sys.path.append(home_dir+'tests')
import utilities as ut
import copula_analysis as ca
import copulaModels as mod
import vines


parameter={'offsets':[14,17],'date_range':('2010-07-01 18:00:00','2016-07-05 01:00:00'),'first_hour':(0,0)}
#titles=[{'type':'Wind','location':'total','kind':'forecast'},{'type':'Wind','location':'NP','kind':'error'},{'type':'Solar','location':'NP','kind':'error'}]
titles=[{'type':'Wind','location':'total','kind':'error'}]
copula=ca.create_copulaManager(titles,parameter)

list_models=[vines.cop2d_uniform,vines.cop2d_frank,vines.cop2d_gumbel,vines.cop2d_clayton,vines.cop2d_gaussian,vines.cop2d_student]


def create_list_models(nb_max=3):
    list_models=[vines.cop2d_uniform,vines.cop2d_frank,vines.cop2d_gumbel,vines.cop2d_clayton,vines.cop2d_gaussian,vines.cop2d_student]

    length=len(list_models)
    res=list_models.copy()
    l=[0 for i in range(length)]
    while l is not None:
        if 1<sum(l)<=nb_max:
            temp=[list_models[i] for i in range(length) if l[i]==1]
            def f(x,temp=temp):
                return vines.WeightedCopula(x,temp)
            res.append( f )
        l=ut.table_increment(l)
    return res



res=ca.test_models(copula,list_models=create_list_models(nb_max=5),end_incr=10000)
