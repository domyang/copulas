import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import math
import random
import utilities as ut
import copulaMS as cms


class SubsetModel(object):


    unifs=[[]]
    model=None
    dim=0
    length=0
    name=''
    precision=30
    quad_points=0
    acr=''

    def is_in_subset(self,x,d):
        return False

    def inverse(self,x):
        return (-1,0)

    def uniform(self,subset_number,interval,count):
        def f(x):
            return 0
        return {'size':1,'density':f,'points':[]}

    def area(self,interval):
        return interval[1]-interval[0]

    def density(self):
        def f(x):
            return 1
        return f

    def pprint(self):
        s='\n### Subset Model: '+self.name+' ###\n\n'
        s+='unifs: (%d,%d)\n'%(len(self.unifs),len(self.unifs[0]))
        for i in range(self.dim):
            s+='    unifs[%d]: %r\n'%(i,self.unifs[i][:6])
        print(s)

#           ############
            ### TAIL ###
            ############


class tail(SubsetModel):


    name='tail'
    acr='TL'

    def __init__(self,unifs,model=None,precision=None):
        self.unifs=unifs
        self.model=model
        self.dim=len(unifs)
        self.length=len(unifs[0])
        if (precision is not None)&(type(precision)==int):
            self.precision=precision
        else:
            self.precision=int(1+np.sqrt(self.length))

        quad_points=[0 for i in range(2**self.dim)]
        for pt in list(zip(*unifs)):
            quad_points[self.inverse(pt)[0]]+=1
        self.quad_points=quad_points

    def is_in_subset(self,x,d):
        res=False
        p=1
        for i in range(self.dim):
            if 0.5*(2-d)<x[i]<=1:
                res+=p
                p*=2
            if (x[i]>1)|(x[i]<0):
                return -1
        return res

    def inverse(self,x):
        sub=0
        norm=0
        p=1
        for i in x:
            norm=max(norm,1-2*abs(0.5-i))
            if 0.5<i<=1:
                sub+=p
            p*=2
            if (i>1)|(i<0):
                return -1

        return(sub,norm)

    def density(self,precision=None):


        if (precision is None)|(type(precision)!=int):
            precision=self.precision
        f1=cms.tails_ut(self.unifs,precision=precision,interpolation='linear', cumulative=False)
        # tail_func=cms.tails_ut(self.unifs,precision=precision,interpolation='epi_spline', cumulative=False)
        # def res(x):
        #     inv=self.inverse(x)
        #     return tail_func[inv[0]](inv[1])
        # return res

    def area(self,interval):
        return interval[1]**2-interval[0]**2

    def uniform(self,corner,interval,count):
        def aux(dim):
            a,b=interval[0],interval[1]
            if dim==1:
                return [random.uniform(a,b)]
            else:
                x=random.random()
                if x<((b-a)*b**(dim-1))/(b**dim-a**dim):
                    res=[]
                    for i in range(dim-1):
                        res.append(random.uniform(0,b))
                    res.append(random.uniform(a,b))
                else:
                    res=aux(dim-1)
                    res.append(random.uniform(0,a))
                return res

        result=[]
        for i in range(count):
            result.append(aux(self.dim))

        cor=[]
        for i in range(self.dim):
            cor.append(corner%2)
            corner//=2
        cor.reverse()

        for i in range(self.dim):
            if cor[i]==1:
                for k in result:
                    k[i]=1-0.5*k[i]
            else:
                for k in result:
                    k[i]=0.5*k[i]
        return result

# t=sub.tail(copula2.unifM)
# T=t.uniform(6,[0.7,0.8],1000)
# res=[[],[],[]]
# for i in T:
#     for k in range(3):
#         res[k].append(i[k])
# fig=plt.figure()
# ax=fig.add_subplot(111, projection='3d')
# ax.scatter3D(res[0],res[1],res[2])

#           ################
            ### DIAGONAL ###
            ################


class diagonal(SubsetModel):

    name='diagonal'
    acr='DG'

    def __init__(self,unifs,concerned=[0,1],model=None,precision=None):
        self.unifs=unifs
        self.model=model
        self.dim=len(unifs)
        self.length=len(unifs[0])
        b=( self.dim>=len(concerned)>1)
        for i in concerned:
            if not b:
                break
            b&=self.dim>i>=0
            b&=type(i)==int
        b=len(concerned)==len(set(concerned)) # check for doublons

        if not b:
            raise(RuntimeError('at lest two dimensions must be considered'))
        self.concerned=concerned
        self.dim_c=len(concerned)
        if (precision is not None)&(type(precision)==int):
            self.precision=precision
        else:
            self.precision=int(1+np.sqrt(self.length))
        self.quad_points=[self.length]


    def is_in_subset(self,x,d):
        y=[x[i] for i in self.concerned]
        prod=0
        for i in range(self.dim_c):
            prod+=y[i]
        prod/=self.dim_c

        dist=0
        for i in y:
            dist+=(i-prod)**2
        dist=math.sqrt(dist)

        return dist<d

    def inverse(self,x):
        try:
            y=[x[i] for i in self.concerned]
        except IndexError:
            print('x: %r' %x)
            print('concerned: %r' %self.concerned)
        prod=0
        for i in range(self.dim_c):
            prod+=y[i]
        prod/=self.dim_c
        dist=0
        for i in y:
            dist+=(i-prod)**2
        dist=math.sqrt(dist)
        return(0,dist*math.sqrt(self.dim_c))

    def density(self,precision=None):
        if (precision is None)|(type(precision)!=int):
            precision=self.precision
        abs=[(i+0.5)/precision for i in range(precision)]
        bins=[0 for i in range(precision)]
        incr=precision/self.length
        list_inverse=[]
        for x in zip(*self.unifs):
            d=self.inverse(x)[1]
            list_inverse.append(d)
            bins[min(precision-1,int(d*precision))]+=incr

        # return abs,bins
        general_approximate=ut.create_spline(ut.find_spline(abs,bins,positiveness_constraint=True))
        tail_approximate,transition_start=estimate_tail(list_inverse,precision=precision)

        def f(x,tail_approximate=tail_approximate,general_approximate=general_approximate,transition_start=transition_start):
            transition_end=(2+transition_start)/3
            coeff = max(0,min(1,(transition_end-x)/(transition_end-transition_start)))
            return coeff*general_approximate(x)+(1-coeff)*tail_approximate(x)

        return([f])

    def uniform(self,subset_number,interval,count,returnArea=False):
        dim_c=self.dim_c
        a,b=interval[0],interval[1]

        def basis(dim):
            proj=[]
            for i in range(dim-1):
                u=np.array([-1/(dim-1) for j in range(dim)])
                u[i]=1
                for j in proj:
                    u=u-sum(j*u)*j
                u=u/math.sqrt(sum(u*u))
                proj.append(u)
            proj=np.matrix(proj)*math.sqrt(dim_c/(dim_c-1))
            return proj

        proj=basis(dim_c)
        l=np.matrix([1 for i in range(dim_c)])

        def aux(dim,aa,bb):
            if dim<=0:
                return []
            else:
                test=aa**(dim-1)*(bb-aa)/(bb**dim-aa**dim)
                r=random.random()
                if r<test/2:
                    # print(1)
                    res=[random.uniform(-aa/2,aa/2) for i in range(dim-1)]
                    res.append(-random.uniform(aa/2,bb/2))
                    return res
                elif r>1-test/2:
                    # print(2)
                    res=[random.uniform(-aa/2,aa/2) for i in range(dim-1)]
                    res.append(random.uniform(aa/2,bb/2))
                    return res
                else:
                    # print(3)
                    res=aux(dim-1,aa,bb)
                    res.append(random.uniform(-bb/2,bb/2))
                    return res

        incr=0
        res=[]
        stop=0

        while (incr<count)&(stop<100000):
            x=np.matrix(aux( dim_c-1,a/math.sqrt(dim_c-1),b ))*proj
            x+=l*random.random()
            x=x.tolist()[0]
            temp=[random.random() for i in range(self.dim)]
            for i in range(dim_c):
                temp[self.concerned[i]]=x[i]
            x=temp

            inside=a<=self.inverse(x)[1]<=b
            for i in x:
                inside&=0<=i<=1

            if inside:
                incr+=1
                res.append(tuple(x))
            stop+=1

        if returnArea:
            return incr/stop

        return res

    def area(self,interval):
        def sphere_volume(dim,r):
            if dim%2==0:
                return math.pi**(dim//2)/(math.factorial(dim//2))*r**dim
            else:
                return (4*math.pi)**(dim//2)*2*math.factorial(dim//2)/math.factorial(dim)*r**dim

        frac=self.uniform(0,interval,2000,returnArea=True)
        vol=sphere_volume(self.dim_c-1,math.sqrt(self.dim_c)*0.5)

        return vol*(interval[1]**(self.dim_c-1)-interval[0]**(self.dim_c-1))*frac+0.0001

#           #################
            ### MAIN AXIS ###
            #################


class main_axis(SubsetModel):

    name='main_axis'
    acr='MA'

    def __init__(self,unifs,model=None,precision=None):
        self.unifs=unifs
        self.model=model
        self.dim=len(unifs)
        self.length=len(unifs[0])
        if (precision is not None)&(type(precision)==int):
            self.precision=precision
        else:
            self.precision=int(1+np.sqrt(self.length))
        self.axis=np.empty(self.dim)
        self.axis[:]=1
        self.quad_points=[self.length]

    def is_in_subset(self,x,d):
        return sum(np.array(x)*self.axis)<=d*self.dim

    def inverse(self,x):
        return(0,sum(np.array(x)*self.axis/self.dim))

    def density(self,precision=None):
        if (precision is None)|(type(precision)!=int):
            precision=self.precision
        abs=[(i+0.5)/precision for i in range(precision)]
        bins=[0 for i in range(precision)]
        incr=precision/self.length
        list_inverse=[]
        for x in zip(*self.unifs):
            d=self.inverse(x)[1]
            list_inverse.append(d)
            bins[min(precision-1,int(d*precision))]+=incr

        general_estimate=ut.create_spline(ut.find_spline(abs,bins,positiveness_constraint=True))
        # tail_low,start_transition_low=estimate_tail([1-i for i in list_inverse])
        # tail_high,start_transition_high=estimate_tail(list_inverse)
        #
        # start_transition_low=min(start_transition_high,start_transition_low)
        # start_transition_high=max(start_transition_high,start_transition_low)
        #
        # def f(x,general_estimate=general_estimate,tail_low=tail_low,tail_high=tail_high,
        #       start_transition_low=start_transition_low,start_transition_high=start_transition_high):
        #
        #     end_transition_low=1/3*start_transition_low
        #     end_transition_high=(2+start_transition_high)/3
        #
        #     if x<end_transition_low:
        #         return 1-tail_low(1-x)
        #     elif x<start_transition_low:
        #         coeff=(x-end_transition_low)/(start_transition_low - end_transition_low)
        #         return coeff*general_estimate(x) + (1-coeff)*(1-tail_low(1-x))
        #     elif x<start_transition_high:
        #         return general_estimate(x)
        #     elif x<end_transition_high:
        #         coeff=(end_transition_high-x)/(end_transition_high-start_transition_high)
        #         return coeff*general_estimate(x) +(1-coeff)*tail_high(x)
        #     else:
        #         return tail_high(x)

        return([general_estimate])

    def uniform(self,subsetNumber,interval,count,returnArea=False):
        radius=min(min(interval[1],1-interval[0])*self.dim*math.sqrt(self.dim-1),math.sqrt(self.dim/2))

        def basis(dim):
            proj=[]
            for i in range(dim-1):
                u=np.array([-1/(dim-1) for j in range(dim)])
                u[i]=1
                for j in proj:
                    u=u-sum(j*u)*j
                u=u/math.sqrt(sum(u*u))
                proj.append(u)
            proj=np.matrix(proj)
            return proj

        incr=0
        stop=0
        res=[]
        # laux=[]
        while (incr<count)&(stop<10000):
            pt=(random.uniform(interval[0],interval[1])*self.axis).tolist()
            # print(pt)
            aux=np.matrix([[random.uniform(-radius,radius) for i in range(self.dim-1)]])
            # print(aux)
            aux= aux*basis(self.dim)
            # print(aux)
            # laux.append(aux.tolist()[0])

            for i in range(0,self.dim):
                pt[i]+=aux[0,i]

            b=True
            for i in pt:
                if not 0<=i<=1:
                    b=False
                    break
            if b:
                res.append(pt)
                incr+=1
            stop+=1

        # ut.plot(laux)
        if returnArea:
            return(interval[1]-interval[0])*math.sqrt(self.dim)*(2*radius)**(self.dim-1)*incr/stop
        else:
            return res

    def area(self,interval):
        return(self.uniform(0,interval,10000,returnArea=True))






#-------------------------------------------------------

def create_diagonal(concerned):
    def aux(unifs,model=None):
        return diagonal(unifs,concerned=concerned,model=model)
    return aux

# returns a list of subsets to be used in copulaModel.cop_customized
def all_subsets(dim):
    # subsets=[tail]
    subsets=[]
    for i in range(dim):
        for j in range(i+1,dim):
            subsets.append(create_diagonal([i,j]))
    subsets.append(main_axis)
    return subsets


### returns an estimate of the tail density with the form y= a(1-x)^b = exp( log(a) + b*log(1-x) )
### it will compute a and b for 1-(empirical_cdf) and then derive it
def estimate_tail(obs,precision=30):

    temp=[]
    for i in obs:
        if 0<i<1:
            temp.append(1-i)
        elif i>=1:
            temp.append(0.000000001)
        else:
            temp.append(0.999999999)
    obs=temp

    def default(s):
        return 1.5*math.sqrt(1-s)

    min_pt=10
    obs_kept=[]
    for i in obs:
        if i<2/precision:
            obs_kept.append(i)

    if len(obs_kept)<min_pt:
        obs_kept=sorted(obs)[:min_pt]
    else:
        obs_kept.sort()

    length0=len(obs)
    if length0==0:
        print('tail estimation impossible with only one observation')
        return default

    length=len(obs_kept)
    yy=list(map(math.log,[(1+i)/(length0+1) for i in range(length)]))
    xx=list(map(math.log,obs_kept))

    mean_x=np.mean(xx)
    mean_y=np.mean(yy)

    squares_x=np.sum([(x-mean_x)**2 for x in xx])
    if squares_x==0:
        print('tail estimation impossible with only one observation')
        return default

    b=np.sum([(xx[i]-mean_x)*(yy[i]-mean_y) for i in range(length)])/squares_x
    log_a=mean_y-b*mean_x
    a=math.exp(log_a)

    b=b-1
    a*=(b+1)

    def f(s):
        return a*(1-s)**b

    return f,min(obs_kept)

#
# def makef(a,b):
#     def f(x):
#         return (x*(b+1)/a)**(1/(b+1))
#     return f
#
# f_inv=makef(2,3)
# x0=[(i+1)/1000 for i in range(999)]
# y0=x0
# x0=[f_inv(i) for i in y0]
# plt.figure(),plt.plot(x0,y0)