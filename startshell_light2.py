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
titles=[{'type':'Solar','location':'total','kind':'error'}]
copula=ca.create_copulaManager(titles,parameter)

gu=vines.cop2d_gumbel
ga=vines.cop2d_gaussian
st=vines.cop2d_student
fr=vines.cop2d_frank
un=vines.cop2d_uniform
cl=vines.cop2d_clayton

list1=[gu,ga,st,fr,un,cl]

list_models=[gu,ga,st,fr,un,cl,
lambda unifs: vines.WeightedCopula(unifs,[gu,ga,st]),
lambda unifs: vines.WeightedCopula(unifs,[gu,ga,st,un]),
lambda unifs: vines.WeightedCopula(unifs,[gu,st,fr]),
lambda unifs: vines.WeightedCopula(unifs,[fr,gu,ga,st]),
lambda unifs: vines.WeightedCopula(unifs,[gu,ga,st,cl]),
lambda unifs: vines.WeightedCopula(unifs,list1,precise=True,max_models=2),
lambda unifs: vines.WeightedCopula(unifs,list1,precise=True,max_models=3),
lambda unifs: vines.WeightedCopula(unifs,list1,precise=True,max_models=4),
lambda unifs: vines.WeightedCopula(unifs,list1,precise=True,max_models=10),
lambda unifs: vines.WeightedCopula(unifs,list1,precise=False,max_models=3)]

#res=ca.test_models(copula,list_models=list_models,end_incr=92)


