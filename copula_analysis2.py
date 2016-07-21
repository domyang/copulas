import time
import csv
import sys
import matplotlib.pyplot as plt
import math
import numpy as np
from scipy import stats
import dateutil.parser as dt
from dateutil.relativedelta import relativedelta
import utilities as ut
import copulaMS as cms
import copulaModels as mod
import subsetModel as sub
import vines

home_dir = 'C:\\Users\\Sabrina\\Documents\\Research\\Code for Real User\\'

### this function fetches the data to create a MS (multi series) copula manager.
### for solar data, it additionally takes into account the solar hour factor
### (normalize the variance to 1, and keep the variance of errors as a function of solar hour in memory)
# returns a MS copula manager
# arguments:
#   () titles: {}
def create_copulaManager(titles,parameters_def,list_parameters=None,just_parameters=None):
    length=len(titles)
    series=[]

    good_format=False
    if (type(titles)==list):
        if len(titles)>0:
            good_format=True
    if not good_format:
        raise(RuntimeError('titles must be a list of dictionaries, with length>0'))
    if not list_parameters is None:
        if len(list_parameters)!=length:
            raise(RuntimeError('list_parameters must be either None or of same length than titles'))
    if not just_parameters is None:
        if len(just_parameters)!=length:
            raise(RuntimeError('just_parameters must be either None or of same length than titles'))

    opts={'type':{'Wind','Solar'},'location':{'NP','SP','total'},'kind':{'error','forecast','actual'}}
    mandatory_keys={'type','location','kind'}
    for t in titles:
        good_format=False
        if type(t)==dict:
            if set(t.keys()).issubset(mandatory_keys):
                good_format2=True
                for key in mandatory_keys:
                    if not t[key] in opts[key]:
                        good_format2=False
                good_format=good_format2
    if not good_format:
        raise(RuntimeError('please check the arguments specifying the kind of data'))

    for i in range(length):
        t=titles[i]
        index=titles.index(t)
        if(index<i):
            temp={}
            for key in series[index].keys():
                temp[key]=series[index][key].copy()
            series.append(temp)
        else:
            temp={}
            res=get_data(type=t['type'],location=t['location'])
            if t['type']=='Solar':
                # call to take into account the importance of solar hour in solar power forecast errors
                res2=ut.prepare_solar(res,visualize=False)
                temp['date']=res2['date']
                temp['vect']=res2[t['kind']]
                temp['data']={'forecast':res2['for'],'forecast_d':ut.list_derivative(res2['for']), 'hour_sol':res2['hour_sol']}
                temp['var_function']=res2['var_function']
            else:
                temp['date']=res['date']
                temp['var_function']=None
                if t['kind']=='error':
                    err=[]
                    a,b=res['act'].copy(),res['for'].copy()
                    while a!=[]:
                        err.append(a.pop()-b.pop())
                    err.reverse()
                    temp['vect']=err
                elif t['kind']=='forecast':
                    temp['vect']=res['for']
                else:
                    temp['vect']=res['act']
                temp['data']={'forecast':res['for'],'forecast_d':ut.list_derivative(res['for'])}
            temp['title']=t
            series.append(temp)

    return cms.CopulaManagerMS(series,parameters_def,list_parameters=list_parameters,just_parameters=just_parameters)




#------------------------------------ comparing functions ------------------------------------------------------




### computes the L2 distance between two copula
def compare_copula_l2(to_compare,visualize=True):
    res=[]
    for i in to_compare['unifs']:
        res.append(ut.compute_distance_l(to_compare['unifs'][0],i))

    if visualize:
        ind=sorted(range(len(res)),key=res.__getitem__,reverse=False)
        nm=list(map(to_compare['names'].__getitem__,ind))
        scores=list(map(res.__getitem__,ind))
        wid=0.5
        abs=np.array(range(len(scores)))
        plt.figure()
        plt.title('L2 comparison of the underlying copulae')
        plt.bar(abs,scores,width=0.5)
        plt.xticks(abs + wid/2., nm,rotation=40)

    return res


### compares models fitted to a copula, using emd distance.
# returns:
#   () a matrix featuring the distance between the various models
#arguments:
#   () copula: a copulaMS manager
#   () sample_size: the size of the sample to be used to compute emd distance.
#       The final distance will be an average of the EMD computed on these samples
def compare_distributions_emd(copula,to_compare=None,sample_size=50,visualize=True):
    length=copula.lengthM
    dim=copula.dim
    vects=[]

    if (to_compare is None):
        vects.append(copula.vectM)
        # creating a fitted gaussian
        covariance=np.cov(vects[0])
        means=np.mean(vects[0],axis=1)
        gaussian=np.transpose(np.random.multivariate_normal(means, covariance, length))
        temp=[]
        for i in range(dim):
            temp.append(gaussian[i].tolist())
        vects.append(temp)

        # computing the inverse of marginals' CDF
        F_inv=ut.marginals_cdf_inv(vects[0])

        # creating copulae and the corresponding distribution
        x=np.arange(0,1,0.02)

        unifs=copula.unifM
        subsets=[]
        for i in range(dim):
            for j in range(i+1,dim):
                subsets.append(sub.create_diagonal([i,j]))
        subsets.append(sub.main_axis)


        for cop in [mod.cop_gaussian(unifs),mod.cop_student(unifs)]: #,mod.cop_student_custom(unifs),
                    #mod.cop_customized(mod.cop_gaussian,unifs,subsets),mod.cop_customized(mod.cop_student,unifs,subsets)]:
            f=cop.f
            unifs=f(length)
            vect_temp=[]
            if visualize:
                plt.figure()
                plt.title('empirical_copula')
                plt.plot(unifs[0],unifs[1],'.')
            for i in range(dim):
                vect_temp.append(list(map(F_inv[i],unifs[i])))

            vects.append(vect_temp)
    else:
        for i in to_compare['vects']:
            vects.append(i.copy())

    # reshuffling the vectors
    points=[]
    for i in vects:
        temp=list(zip(*i))
        np.random.shuffle(temp)
        points.append(temp)

    # Computing the mean of emd distance between samples of distribution
    nb_dist=len(vects)

    for i in range(nb_dist):
        l_temp=len(points[i])
        if length>l_temp:
            print("distribution %d is length %d instead of %d"%(i,l_temp,length))
            length=l_temp

    nb_sample=int(length/sample_size)
    res=np.identity(nb_dist)
    for i in range(nb_dist):
        for j in range(i,nb_dist):
            p1=points[i]
            p2=points[j]
            print('len1 %d len2 %d'%(len(p1),len(p2)))
            print(p1[:10])
            print(p2[:10])
            res[i,j]=0;res[j,i]=0
            for k in range(nb_sample):
                a=k*sample_size
                b=((k+1)%nb_sample)*sample_size
                print(nb_sample)
                print('dist %d %d --- k: %.2f a: %.2f b: %.2f'%(i,j,k,a,b))
                res[i,j]+=ut.compute_emd(p1[a:(a+sample_size)],p2[b:(b+sample_size)])
            res[j,i]=res[i,j]

    # if visualize:
    #     for v in vects:
    #         plt.figure()
    #         plt.plot(v[0],v[1],'.')
    #         print(res)

    for i in range(nb_dist-1):
        for j in range(i+1,nb_dist):
            res[i,j]=round(res[i,j]/math.sqrt(res[i,i]*res[j,j]),4)
            res[j,i]=res[i,j]
    for i in range(nb_dist):
        res[i,i]=1
        res.tolist()

    print(res)
    if visualize:
        scores=res[0].copy()
        ind=sorted(range(len(scores)),key=scores.__getitem__,reverse=False)
        if to_compare is None:
            nm=['' for i in range(len(scores))]
        else:
            nm=list(map(to_compare['names'].__getitem__,ind))
        scores=list(map(scores.__getitem__,ind))
        wid=0.5
        ab=np.array(range(len(scores)))
        plt.figure()
        plt.title('emd comparison of the distributions')
        print(ab,scores)
        plt.bar(ab,scores,width=0.5)
        plt.xticks(ab + wid/2., nm,rotation=40)


    return(res)


### compares models fitted to a c
# opula,using log-likelihood:
def compare_distributions_log(to_compare,visualize=True):
    vects=to_compare['vects']
    length=len(vects)
    density=to_compare['dis_density']

    res=[]
    for i in range(length):
        print(i)
        res.append(ut.emp_log_likelihood(vects[0],vects[i],density1=density[0],density2=density[i]))

    if visualize:
        for result in list(zip(*res)):
            ind=sorted(range(len(result)),key=result.__getitem__,reverse=True)
            nm=list(map(to_compare['names'].__getitem__,ind))
            scores=list(map(result.__getitem__,ind))
            scores=[scores[0]-i for i in scores]
            wid=0.5
            abs=np.array(range(len(scores)))
            plt.figure()
            plt.title('log-likelihood comparison \n (-L(model)+L(original))')
            plt.bar(abs,scores,width=0.5)
            plt.xticks(abs + wid/2., nm,rotation=40)
    return res


### creates a 'to_compare' dictionary and compares different models using the 3 previous distances 
def compare_distributions(copula,list_models=None,visualize=True):
    unifs=copula.unifM
    vects=[copula.vectM]
    names=['real']
    cop_density=[None]
    length=copula.lengthM
    dim=copula.dim

    # creating a fitted gaussian
    covariance=np.cov(vects[0])
    means=np.mean(vects[0],axis=1)
    gaussian=np.transpose(np.random.multivariate_normal(means, covariance, length))
    temp=[]
    for i in range(dim):
        temp.append(gaussian[i].tolist())
    vects.append(temp)
    names.append('gaussian (w/o copula)')
    cov_inv=np.matrix(covariance)**(-1)
    fact_gau=np.sqrt(np.linalg.det(cov_inv)/(2*math.pi)**dim)
    def den_gau(vec):
        if type(vec[0]) in {int,np.float,float,np.int}:
            vec=[[i] for i in vec]
        res=[]
        for i in zip(*vec):
            i=np.matrix(i)
            res.append(fact_gau*math.exp(-i*cov_inv*np.transpose(i))/2)
        return res

    cop_density.append(den_gau)
    del(temp)

    # print(vects)
    # raise(RuntimeError())
    if list_models is None:
        # different models
        list_models=[mod.cop_emp_redistributed(unifs),mod.cop_gaussian(unifs),mod.cop_student(unifs),
            mod.cop_customized(mod.cop_gaussian,unifs),mod.cop_customized(mod.cop_student,unifs)]
        # renormalizing the customized copulae
        for i in list_models[-2:]:
            list_models.append(mod.cop_emp_redistributed(i.val,model=i))
        # return list_models

    unifs=[]
    for i in list_models:
        unifs.append(i.val)
        names.append(i.name)
        cop_density.append(i.pdf)

    cop_density[2]=None
    dis_density=ut.copula_to_densities(vects[0],cop_density)


    to_dist=ut.copula_to_distribution(vects[0],visualize=False)
    for i in unifs:
        vects.append(to_dist(i))

    unifs.reverse()
    unifs.append(mod.uniforms(vects[1],rand=False))
    unifs.append(mod.uniforms(vects[0],rand=False))
    unifs.reverse()

    to_compare={'vects':vects,'unifs':unifs,'names':names,'dis_density':dis_density}


    #comparing the distributions using various measures
    res=[]
    # res_tp=[ut.emp_log_likelihood(vects[0],vects[i],density1=dis_density[0],density2=dis_density[i]) for i in range(len(vects))]
    # if visualize:
    #     temp=[i[0] for i in res_tp]
    #     ind=sorted(range(len(temp)),key=temp.__getitem__,reverse=True)
    #     nm=list(map(to_compare['names'].__getitem__,ind))
    #     scores=list(map(temp.__getitem__,ind))
    #     scores=[scores[0]-i for i in scores]
    #     wid=0.5
    #     abs=np.array(range(len(scores)))
    #     plt.figure()
    #     plt.title('log-likelihood comparison \n (-L(model)+L(original))')
    #     plt.bar(abs,scores,width=0.5)
    #     plt.xticks(abs + wid/2., nm,rotation=40)
    #
    #
    # res.append(res_tp)

    res.append(compare_distributions_log(to_compare,visualize=True))
    res.append(compare_copula_l2(to_compare))
    res.append((compare_distributions_emd(copula,to_compare=to_compare).tolist())[0])
    # res.append(compare_distributions_log(to_compare))

    return res,to_compare


### This function compares the predicting accuracy of the various models
# arguments:
#   () copula is a copulaManagerMS
#   () win_days is a window parameter specifying the number of days before and after the current time of the year 
#      that should be taken into account. (if=45, it translates into 91 days in the past years and 45 in the current year)
#   () win_forecast is a window parameter specifying the width of the forecast window: if q=CDF_forecast(current_forecast), 
#      we will consider dates at which the forecast was in [CDF^-1(q-'win_forecast'),CDF^-1(q+'win_forecast')]
def test_models(copula,win_days=45,repeat_win_days=False,win_forecast=0.2,visualize=False,start_incr=None,end_incr=None,keep_vec=False):
    start_time=time.time()

    # keeping the old window parameters in param_fixed
    param_fixed=[]
    for par in copula.parameters:
        dic={}
        dic['date_range']=par['date_range']
        dic['offsets']=par['offsets'].copy()
        dic['first_hour']=par['first_hour']
        param_fixed.append(dic)

    # initializing variables: parameters, dim (dimension of the copula), forecast, errors, dates, 
    #                         forecastCDF and forecastQuantile (inversse function)
    
    parameters=copula.parameters.copy()
    dim=copula.copulae[0].dim

    if not 'forecast' in copula.copulae[0].dataM.keys():
        parameters[0]['forecast']=[(-100000,100000) for i in parameters[0]['offsets']]
        copula.update(parameters[0],list_parameters=parameters)

    forecasts=[]
    for i in copula.copulae[0].dataM['forecast']:
        forecasts.append(i.copy())
    errors=[]
    for i in copula.vectM:
        errors.append(i.copy())
    dates=copula.dateM.copy()
    print('forecasts (%d) %r\nerrors (%d) %r\ndates (%d) %r'%(len(forecasts),forecasts[:8],
                                                              len(errors),errors[:8],len(dates),dates[:8]))

    forecastCDF=ut.empirical_CDF_scalar(forecasts[0])
    forecastQuantile=ut.empirical_CDF_inv(forecasts[0])


    ### BIG LOOP ###


    incr=0
    first=True
    last_time=time.time()
    res={'len':[],'log':[],'rank':[],'sum_pdf':[],'vec':[],'problem':[],'past_log':[]}

    # loop over each considered hour: each time, 
    #       The copula manager is updated to fit the window (using 'win_days' and 'win_forecast')
    #       Models are created to fit the updated copula
    #       The log_likelihood of the observation is computed for all different models
    for obs in zip(*[forecasts[0],list(zip(*errors)),dates]):

        # selecting the observation range
        incr+=1
        if start_incr is None:
            if incr<400:
                continue
        else:
            if incr<start_incr:
                continue
        if end_incr is not None:
            if incr>=end_incr:
                continue


        t_print=[]
        for i in (time.time()-last_time,time.time()-start_time):
            i=int(i)
            t_print.append((i//3600,i//60%60,i%60))
        last_time=time.time()
        print('\n\n   ####################   \n\niter %d: forecast %r, error %r, date %s\n\n'
              'time elapsed in the last loop: %d:%d:%d, time since start: %d:%d:%d'
              '\n\n   ####################   \n\n'
              %(incr,obs[0],obs[1],obs[2],t_print[0][0],t_print[0][1],t_print[0][2],t_print[1][0],t_print[1][1],t_print[1][2]))

        ### selecting past observations using the window ###

        print('selecting observations using the window')
        mid=forecastCDF(obs[0])
        temp=[(-np.inf,np.inf) for obs in range(dim)]
        temp[0]=(float(forecastQuantile(max(0.0001,mid-win_forecast))),float(forecastQuantile(min(0.9999,mid+win_forecast))))
        parameters[0]['forecast']=temp
        del temp
        parameters[0]['predicted_day']=obs[2]
        if not repeat_win_days:
            parameters[0]['date_range']=(str(dt.parse(obs[2])-relativedelta(days=win_days)),str(dt.parse(obs[2])-relativedelta(hours=1)))
        else:
            parameters[0]['date_range']=ut.intersect_dates(param_fixed[0]['date_range'],obs[2],win_days)


        copula.update(parameters[0],list_parameters=parameters)


        ### fitting models to the distribution ###

        print('fitting models to the distribution')
        length=copula.lengthM
        dim=copula.dim

        if length<50:
            res['problem'].append('length <50 at iteration %d'%incr)
            continue


        try:
            # creating the density of a fitted gaussian
            def create_gaussian_density():
                covariance=np.cov(copula.vectM)
                means=np.mean(copula.vectM,axis=1)

                cov_inv=np.linalg.inv(np.matrix(covariance))
                fact_gau=np.sqrt(np.linalg.det(cov_inv)/(2*math.pi)**dim)
                def den_gau(vec):
                    if type(vec[0]) in {int,np.float,float,np.int}:
                        vec=[[j] for j in vec]
                    res=[]
                    for j in zip(*vec):
                        j=np.matrix(j)
                        res.append(fact_gau*math.exp(-(j-means)*cov_inv*np.transpose(j-means)/2))
                    return res
                return den_gau

            densities=[create_gaussian_density()]

            # creating a list of copula models
            list_models=[mod.cop_gaussian(copula.unifM),mod.cop_student(copula.unifM),
                         mod.cop_customized(mod.cop_student,copula.unifM),
                         vines.D_vine(copula.unifM)]#,vines.C_vine(copula.unifM)]

            # computing the densities of the models, their log-likelihood, and selecting the 'best candidate'

            names=['gaussian']
            cop_densities=[]
            best_model_past=mod.cop_gaussian
            best_log_past=0
            log_past=[]
            for j in list_models:
                if first:
                    names.append(j.name)
                cop_densities.append(j.pdf)
                lld=sum([math.log(k) for k in j.pdf(copula.unifM)])
                log_past.append(lld)
                if lld>best_log_past:
                    best_log_past=lld
                    best_model_past=j.pdf

            names.append('selected model')
            cop_densities.append(best_model_past)
            log_past.append(best_log_past)

            res['past_log'].append(log_past)

            # computing the rank of 'obs' among the window points
            CDFs=ut.marginals_cdf(copula.vectM)
            rank=[CDFs[i](obs[1][i]) for i in range(dim)]


            try:
                ### computing the log likelihood ###
                res_log=[math.log(den(obs[1])[0]) for den in densities]
                res_log.extend([den(obs[1])[0] for den in   ut.copula_to_densities(copula.vectM,cop_densities,log_return=True)])

                res['log'].append(res_log)
            except:
                res['log'].append(None)
                res['problem'].append('incr: %d, problem in the log: %r'%(incr,sys.exc_info()[0]))
            res['len'].append(length)
            res['rank'].append(rank)

            if keep_vec:
                res['vec'].append([copula.unifM,list_models[0].val,list_models[2].val])
        except:
            res['problem'].append('incr %d general problem: %r'%(incr,sys.exc_info()[0]))

        if first:
            res['names']=names
            first=False

    copula.update(param_fixed[0],list_parameters=param_fixed)
    return res


### This function looks at whether, knowing error at t1, errors at t0 and t2 are independent
def conditional_independance(copula,win_days=45,repeat_win_days=False,win_forecast=0.2,start_incr=None,end_incr=None):

    start_time=time.time()
    # keeping the old window parameters in param_fixed
    param_fixed=[]
    for i in copula.parameters.copy():
        for key in i.keys():
            if type(i[key]) in {set,list,dict}:
                i[key]=i[key].copy()
        param_fixed.append(i)

    # initializing variables: parameters, dim (dimension of the copula), forecast, errors, dates,
    #                         forecastCDF and forecastQuantile (inversse function)

    parameters=copula.parameters.copy()
    dim=copula.copulae[0].dim

    if not 'forecast' in copula.copulae[0].dataM.keys():
        parameters[0]['forecast']=[(-100000,100000) for i in parameters[0]['offsets']]
        copula.update(parameters[0],list_parameters=parameters)

    forecasts=[]
    for i in copula.copulae[0].dataM['forecast']:
        forecasts.append(i.copy())
    errors=[]
    for i in copula.vectM:
        errors.append(i.copy())
    dates=copula.dateM.copy()

    incr=0
    res=[]

    forecastCDF=ut.empirical_CDF_scalar(forecasts[0])
    forecastQuantile=ut.empirical_CDF_inv(forecasts[0])

    # loop over each considered hour: each time,
    #       The copula manager is updated to fit the window (using 'win_days' and 'win_forecast')
    #       Models are created to fit the updated copula
    #       The log_likelihood of the observation is computed for all different models
    # try:
    last_time=time.time()
    for obs in zip(*[forecasts[0],list(zip(*errors)),dates]):

        incr+=1
        if start_incr is None:
            if incr<400:
                continue
        else:
            if incr<start_incr:
                continue
        if end_incr is not None:
            if incr>=end_incr:
                continue

        try:
            t_print=[]
            for i in (time.time()-last_time,time.time()-start_time):
                i=int(i)
                t_print.append((i//3600,i//60%60,i%60))
            last_time=time.time()
            print('\n\n   ####################   \n\niter %d: forecast %r, error %r, date %s\n\n'
                  'time elapsed in the last loop: %d:%d:%d, time since start: %d:%d:%d'
                  '\n\n   ####################   \n\n'
                  %(incr,obs[0],obs[1],obs[2],t_print[0][0],t_print[0][1],t_print[0][2],t_print[1][0],t_print[1][1],t_print[1][2]))

            ### selecting observations using the window ###

            print('selecting observations using the window')
            mid=forecastCDF(obs[0])
            temp=[(-math.inf,math.inf) for obs in range(dim)]
            temp[0]=(float(forecastQuantile(max(0.0001,mid-win_forecast))),float(forecastQuantile(min(0.9999,mid+win_forecast))))
            parameters[0]['forecast']=temp
            del temp
            parameters[0]['predicted_day']=obs[2]
            if not repeat_win_days:
                parameters[0]['date_range']=(str(dt.parse(obs[2])-relativedelta(days=win_days)),str(dt.parse(obs[2])-relativedelta(hours=1)))
            else:
                parameters[0]['date_range']=ut.intersect_dates(param_fixed[0]['date_range'],obs[2],win_days)


            copula.update(parameters[0],list_parameters=parameters)

            gau=mod.cop_gaussian(copula.unifM)
            res.append([gau.par['cov'][0,1],gau.par['cov'][1,2],gau.par['cov'][0,2]])
        except:
            pass
    return res






#-------------------------------------------- other functions --------------------------------------------------


### This functions takes a list of (same length) vectors in arguments, and return a copula based on these vectors
def fake_copulaManager(vects):
    length=len(vects[0])
    d=dt.parse('2012/01/01 00:00')
    hour=relativedelta(hours=1)
    date=[]
    for i in range(length):
        date.append(str(d))
        d+=hour
    series=[]
    for v in vects:
        ser={'date':date,'vect':v,'data':{},'title':{'type':'Wind', 'location':'NP', 'kind':'forecast'}, 'var_function':None}
        series.append(ser)
    parameter={'offsets':[0]}
    return cms.CopulaManagerMS(series,parameter)


### This function fetches the data from external files
# returns a dictionary, 'data' with keys:   - 'act' (actuals)
#                                           - 'for' (forecasts)
#                                           - 'date'
def get_data(type='Wind',location='total',filename=''):
    if filename=='':
        csv_read=csv.reader(open(home_dir+'tests/data/sources.csv'))
        index=-1
        index_dir=-1
        break1=False
        for line in csv_read:
            if (len(line)<1) or (line[1]=='#'):
                continue
            elif line[0]=='title':
                for i in range(len(line)):
                    if line[i]==type+'_'+location:
                        index=i
                    if line[i]=='home_dir':
                        index_dir=i
            if(break1):
                break
        if (index==-1)|(index_dir==-1):
            print('data specification is not valid 1')
            return -1
        else:
            csv_read=csv.reader(open(home_dir+'tests/data/sources.csv'))
            for line in csv_read:
                if (len(line)<1) or (line[1]=='#'):
                    continue
                if line[0]=='data':
                    filename='%s%s' % (str(line[index_dir]),str(line[index]))
                    break
        if filename=='':
            print('data specification is not valid 2')
            return -1
        else:
            print('retrieving data from %s \n' % filename)

    data={'date':[],'for':[],'act':[]}
    csv_read=csv.reader(open(filename))
    for line in csv_read:
        if (len(line)<1) or (line[1]=='#'):
            continue
        data['act'].append(float(line.pop()))
        data['for'].append(float(line.pop()))
        data['date'].append(str(dt.parse(line.pop())))

    print('successfully retrieved %s data \n' % (type+'_'+location) )
    return data