import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pyemd import emd
import math
from pylab import rcParams
import numpy as np
import scipy as sc
from scipy import stats
from scipy.integrate import quad as integ
import dateutil.parser as dt
from dateutil.relativedelta import relativedelta
import sys
import ephem
from importlib import reload
home_dir='C:/Users/Sabrina/Documents/Research/Code for Real User/'
sys.path.append(home_dir+'prescient/release/Prescient_1.0/')
sys.path.append(home_dir+'prescient/release/Prescient_1.0/exec')
sys.path.append(home_dir+'tests')
import ScenGen.workflow.EpiModel as EpiModel
import exec.MasterOptions as MasterOptions
import utilities as ut
import copula_analysis as ca
import test as test
import copulaModels as mod
import subsetModel as sub
import copulaMS as cms
import random



parameter1={'offsets':[8],'date_range':('2013-07-01 18:00:00','2015-07-05 01:00:00'),'first_hour':(-5,5)}
parameter2={'offsets':[8,11,15],'date_range':('2013-07-01 18:00:00','2015-07-05 01:00:00'),'first_hour':(-5,5)}
parameter2bis={'offsets':[8,11,15],'date_range':('2013-07-01 18:00:00','2015-07-05 01:00:00'),'first_hour':(-2,2)}
parameter3={'offsets':[15],'date_range':('2013-07-01 18:00:00','2015-07-05 01:00:00'),'first_hour':(-2,2)}
#parameter4={'offsets':[8,11],'date_range':('2013-07-01 18:00:00','2015-07-05 01:00:00'),'first_hour':(-2,2)}
parameter4={'offsets':[14,17],'date_range':('2013-07-01 18:00:00','2015-07-05 01:00:00'),'first_hour':(-2,2)}
titles=[{'type':'Wind','location':'NP','kind':'forecast'},{'type':'Wind','location':'NP','kind':'error'},{'type':'Solar','location':'NP','kind':'error'}]
parameter_list=[parameter1,parameter2,parameter3]
#copula=ca.create_copulaManager(titles,parameter1,list_parameters=parameter_list)

titles=[{'type':'Wind','location':'NP','kind':'error'}]
#copula1=ca.create_copulaManager(titles,parameter4)
#copula2=ca.create_copulaManager(titles,parameter2)
copula3=ca.create_copulaManager(titles,parameter2bis)

#reload(ut);reload(mod);reload(sub);reload(ca);reload(test)

# parameter4={'offsets':[14,17],'date_range':('2013-07-01 18:00:00','2015-07-05 01:00:00'),'first_hour':(-2,2)}
# copula1.update(parameter4)
