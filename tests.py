import matplotlib.pyplot as plt
import math
import random
import numpy as np
from scipy import stats
import scipy.special as spe
from scipy.optimize import root
from scipy.optimize import minimize_scalar as minim
import subsetModel as sub
import utilities as ut
from scipy.integrate import quad as integ

def density1(x):
    return 1/(2*math.sqrt(2*3.1415)*3)*math.exp(-(x-3)**2/18) + 1/(2*math.sqrt(2*3.1415))*math.exp(-(x-3.5)**2/2)

def density2(x):
    if 1<x<5:
        return 0.25
    else:
        return 0

# coeffs: 2 , 5
def cost1(x,a=2,b=5):
    return -a*x+(a+b)*max(0,x)

def cost2(x,a=2,b=5):
    return a*x**2+(b-a)*(max(0,x))**2


def expectancy(cost,density):

    def res(y):
        def g(x,y=y):
            return cost(y-x)*density(x)
        # ut.curve(g,inter=[-10,10])
        return integ(g,-10,10)[0]
    return res

def pdf_to_quantile(f,precision=0.00001):
    def tmp(x):
        return x*f(x)
    mn=integ(tmp,-10000,10000)[0]
    val=integ(f,-10000,mn)[0]
    max_iter=100
    def quantile(x,precision=precision):
        cur=mn
        quant=val
        incr=0
        while (abs(quant-x)>precision)&(incr<max_iter):
            past_cur=cur
            cur=(x-quant)/f(past_cur)+cur
            quant=quant+(integ(f,past_cur,cur)[0])
            incr+=1

        quant=integ(f,-10000,cur)[0]

        while (abs(quant-x)>precision)&(incr<max_iter):
            past_cur=cur
            cur=(x-quant)/f(past_cur)+cur
            quant=quant+integ(f,past_cur,cur)[0]
            incr+=1

        print('nb iteration; %d, quantile: %f'%(incr,quant))
        return float(cur)
    return quantile

def f(x):
    if 0<=x<0.5: return 0.75
    if 0.5<=x<1: return 0.375+0.5*x+0.5*((1-x)**2+(0.5-x)**2)
    if 1<=x: return x
    else: return f(-x)

def potential():
    mean=[2.2,3.6,5.9,5.7,5.1,6.2,6.45,6.1,4.32,2.3]
    var=[(10-x)*(x+1)/10 for x in range(10)]
    f=[]
    for i in range(10):
        def temp(x,i=i):
            return var[i]*(x-mean[i])**2
        f.append(temp)

    def func(l):
        x,y=l[0],l[1]
        if x<0:
            return f[0](y)
        elif x>=9:
            return f[9](y)
        else:
            ind=int(math.floor(x))
            return (1-(x-ind))*f[ind](y) + (x-ind)*f[ind+1](y)

    return func