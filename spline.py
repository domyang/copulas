import numpy as np
import scipy as sc
from scipy.special import binom
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
import sys

home_dir='/home/ambroiseidoine/UCD/'

sys.path.append(home_dir+'code/prescient/release/Prescient_1.0/')
sys.path.append(home_dir+'code/prescient/release/Prescient_1.0/exec')
import ScenGen.workflow.EpiModel as EpiModel
import exec.MasterOptions as MasterOptions

def draw(self,x_diag,y_diag,x_anti,y_anti,degree):
    s=UnivariateSpline(x_diag,y_diag,k=degree)
    coeffs_diag=s.get_coeffs()

#bernstein polynoms
def BPoly(coeffs,knots):
    def f(x):
        print(coeffs)
        print(knots)
        res=0
        length=len(coeffs)-1
        for i in range(length+1):
            res=res+binom(length,i)*np.power((x-knots[0]),i)*np.power((knots[1]-x),length-i)/np.power((knots[1]-knots[0]),length)*coeffs[i]
        return res
    return f

### returns the spline interpolation coefficients
# returns
#   () dico coefficients for the interpolation
# arguments:
#   () x,y points to be interpolated
#   () visualize: if you wish to plot an example
def find_spline(x,y,options,visualize=True):

    print('seg_N: %d, seg_s: %.2f, seg_kappa: %.2f ,epifit_error_norm: %s, L1Linf_solver: %s, L2Norm_solver: %s'%
          (options.seg_N, options.seg_s, options.seg_kappa, options.epifit_error_norm, options.L1Linf_solver, options.L2Norm_solver))

    model=EpiModel.FitEpispline(x,y,options)

    approx=[]
    for i in model.s.keys():
        approx.append((model.s.get_values())[i])

    # getting coefficients for the estimating function
    N=options.seg_N
    min_x,max_x=(min(x),max(x))
    s0,v0=(model.s0.get_values()[None],model.v0.get_values()[None])

    a=[]
    keys=model.a.keys()
    for i in keys:
        a.append((model.a.get_values())[i])

    dico={'a':a,'s0':s0,'v0':v0,'min_x':min_x,'max_x':max_x,'N':N}

    if visualize:
        f=create_spline(dico)
        plt.plot(x,y,color='blue')
        plt.plot(x,approx,color='red')
        index=np.arange(min_x,max_x,0.1)
        val=[f(i) for i in index]
        plt.plot(index,val,color='green')
        plt.show()

    return dico


### returns a spline function
# returns
#   () f : spline function
# arguments
#   () dico coefficients for the interpolation
def create_spline(dico):

    a=dico['a'];s0=dico['s0'];v0=dico['v0'];min_x=dico['min_x'];max_x=dico['max_x'];N=dico['N']

        #estimating function
    def f(t):
        delta=(max_x-min_x)/N
        len_a=len(a)
        t-=min_x
        k=max(0,min(int(np.floor(t/delta)),len_a-1))
        res=0
        for i in range(k):
            res+=(t+(0.5-i-1)*delta)*a[i]
        res*=delta
        res+=0.5*(t-k*delta)**2*a[k]
        res+=s0+t*v0
        return (res)

    return(f)





def main(args=None):
    # Parse command-line options.
    print ("hello")
    try:
        options_parser, guiOverride = MasterOptions.construct_options_parser()
        (options, args) = options_parser.parse_args(args=args)
    except SystemExit:
        # the parser throws a system exit if "-h" is specified - catch
        # it to exit gracefully.
        return
    c=[0,1,2,3,4,5,6,7,8]
    d=[1,0.1,0.9,2,3.6,4,4.2,4.3,4.1]
    find_spline(c,d,options)

# MAIN ROUTINE STARTS NOW #

# main(sys.argv)