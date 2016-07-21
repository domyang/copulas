## to analyse the correlation between errors, both in time and space
import csv
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pyemd import emd
import math
from pylab import rcParams
import numpy as np
import scipy as sc
import dateutil.parser as dt
from dateutil.relativedelta import relativedelta
import spline


home_dir='/home/ambroiseidoine/UCD/'

#------------------------------ introductory analysis -----------------------------------------------------------------


# arguments:
#   () type: either Solar or Wind
#   () date_range: optional date range to be considered (list of two dates)
# Calls to:
#   () get_error (returns prediction errors)
#   () restrain_to_date_range (to treat a specified date range)
# Returns:
#   () correlation between north and south error forecast
def space_correlation(type,date_range=[]):

    #getting the data
    date1,errors1=get_data(type,'NP')
    print('date len=%d,errors len=%d'%(len(date1),len(errors1)))
    date2,errors2=get_data(type,'SP')
    print('date len=%d,errors len=%d'%(len(date2),len(errors2)))

    if len(date_range)==2:
        restrain_to_date_range(date1,errors1,date_range)
        restrain_to_date_range(date1,errors1,date_range)

    #estimate correlation, checking that dates correspond
    correlation=0
    count=0
    len1=len(date1)
    len2=len(date2)
    diff=0
    for i in range(len1):
        i=len1-i-1
        d1=date1.pop()
        date_match=d1==date2[i+diff]
        if not date_match:
            for j in range(max(0,i+diff-50),min(len2,i+diff+50)):
                if d1==date2[j]:
                    diff=j-i
                    date_match=True
                    break
        if date_match:
            correlation=correlation+errors1[i]*errors2[i+diff]
            count=count+1

    #final result
    correlation=correlation/(count*np.sqrt(sc.var(errors1)*sc.var(errors2)))
    print('space correlation between North and South for %s is %.3f \n'%(type,correlation))
    return correlation


### This functions draws the time correlation of prediction errors
# arguments:
#   () type, location: specify the errors wanted
#   () date_and_errors: if you prefer to put in directly the dates and errors, a tuple containing these 2 lists
#   () range_max : it will analyse the time correlated errors for time intervals ranging from 0 hour to range_max hour
#   () date_range: 2-elements-list to specify a datetime range for the data. by default, all data is considered
#   () fig_nb: figure number for plotting
# Calls to:
#   () get_error (returns prediction errors)
#   () restrain_to_date_range (to treat a specified date range)
# Returns nothing
def time_correlation(type,location,date_and_errors=(),range_max=120, date_range=[], fig_nb=1):

    # getting the date and errors
    try:
        if len(date_and_errors)==2 & len(date_and_errors[0])>0:
            date,errors=date_and_errors
        else:
            raise IndexError("no good arguments")
    except IndexError:
        date,errors=get_data(type,location)
        print('date len=%d,errors len=%d'%(len(date),len(errors)))
        if len(date)<1 | len(errors)<1:
            print('Problem with the provided date & error. Taking a default file instead')
            file_name_wind='/home/ambroiseidoine/UCD/code/prescient/release/Prescient_1.0/exec/confcomp_data/Wind/Wind_actual_forecast_SP15_070114_063015.csv'
            date,errors=get_data(file_name_wind)

    # restraining data to 'date_range' if provided
    if (len(date_range)==2):
        date,errors=restrain_to_date_range(date,errors,date_range)

    # basic statistical analysis
    mean_error=sc.mean(errors);
    var_error=sc.var(errors);
    correl=[]

    nb_days=int(np.ceil(range_max/24))
    for i in range (range_max):
        correl.append((sc.mean([errors[j]*errors[i+j] for j in range(len(errors)-range_max) ])- np.power(mean_error,2))/var_error)

    # when is the correlation under 0.5, 0.2, 0.1?
    temp=correl.copy()
    temp.reverse()
    incr=0
    levels=[0.1,0.2,0.5]
    temp_levels=levels.copy()
    temp_levels.sort(reverse=True)
    crossing_times=[]
    for lev in temp_levels:
        while temp.pop()>lev:
            incr=incr+1
        crossing_times.append(incr)
        print('correlation under %.2f after %d hours'% (lev,incr))
    del(temp,incr)


    # additional information
    location_names={'total':'California','NP':'North California', 'SP':'South California'}
    print("mean error : %d \n" % mean_error)


    # plotting
    plt.figure(fig_nb)
    ax=plt.subplot()
    for i in range(3):
        i=2-i
        ax.fill_between([0,range_max],[-levels[i],-levels[i]],[levels[i],levels[i]], color=(0.7+i*0.1,0.7+i*0.1,0.7+i*0.1))
    plt.plot([0,range_max],[0,0],color='black')
    plt.plot(correl);
    plt.xticks(range(0,24*nb_days,24))
    plt.title("time correlation in %s forecast errors (%s)" % (type,location_names[location]))
    plt.xlabel("time interval")
    plt.ylabel("correlation")
    #plt.show()


#------------------------------ adjacent functions for introductory analysis ------------------------------------------


### This functions reads data as specified by the type (Solar,Wind) and location (total,NP,SP0
# It returns 2 lists:
#   () the date & time (formated string)
#   () the difference between forecast and actual (at this time)
def get_data(type='Solar',location='total', kind='errors'):
    if (type in {'Solar','Wind'}) & (location in {'NP','SP','total'}):
        if location!= 'total':
            location=location+'15' # file names compatibility
#
        data_dir=home_dir+'code/prescient/release/Prescient_1.0/exec/confcomp_data/'
        data_files=['_actual_forecast_','_070114_063015.csv']
        date=[]
        values={'for':[],'act':[]}
#
        if (type=='Wind') & (location=='total'):
            # we add NP15 and SP15 (beware of the missing data: we check that the dates are the same)
            values_temp={'for':[],'act':[]}
            date_temp=[]
            count=0
            for loc in ['SP15','NP15']:
                file_name=data_dir+type+'/'+type+data_files[0]+loc+data_files[1]
                csv_read=csv.reader(open(file_name,'r'))
                if loc=='SP15':
                    for line in csv_read:
                        if (len(line)<1) or (line[1]=='#'):
                            continue
                        date_temp.append(str(dt.parse(line[0])))
                        values_temp['for'].append(float(line[1]))
                        values_temp['act'].append(float(line[2]))
                else:
                    len_value_temp= len(values_temp['act'])
                    for line in csv_read:
                        if (len(line)<1) or (line[1]=='#'):
                            continue
                        if count>=len_value_temp:
                            break
                        if dt.parse(date_temp[count])!=dt.parse(line[0]): #check that the dates correspond
                            print('missing data at date %s for %s' % (dt.parse(line[0]),type))
                            go_on=True
                            for i in range(max(-10,-count),min(10,len_value_temp-1-count)): # Searching neighbouring dates
                                if dt.parse(date_temp[count+i])==dt.parse(line[0]):
                                    count=count+i
                                    go_on=False
                                    continue
                            if go_on:
                                continue
                            else:
                                date.append(str(dt.parse(line[0])))
                        else:
                            date.append(str(dt.parse(line[0])))
#
                        values['for'].append(float(line[1])+values_temp['for'][count])
                        values['act'].append(float(line[2])+values_temp['act'][count])
                        count=count+1
#
        elif (type=='Solar') & (location=='NP15'):
            # we subtract total and SP15 (beware of the missing data: we check that dates correspond)
            values_temp={'for':[],'act':[]}
            date_temp=[]
            count=0
            for loc in ['SP15','total']:
                file_name=data_dir+type+'/'+type+data_files[0]+loc+data_files[1]
                csv_read=csv.reader(open(file_name,'r'))
                if loc=='SP15':
                    for line in csv_read:
                        if (len(line)<1) or (line[1]=='#'):
                            continue
                        date_temp.append(str(dt.parse(line[0])))
                        values_temp['for'].append(float(line[1]))
                        values_temp['act'].append(float(line[2]))
                else:
                    len_value_temp= len(values_temp['act'])
                    for line in csv_read:
                        if (len(line)<1) or (line[1]=='#'):
                            continue
                        if count>=len_value_temp:
                            break
                        if dt.parse(date_temp[count])!=dt.parse(line[0]): #check that the dates correspond
                            print('missing data at date %s for %s' % (dt.parse(line[0]),type))
                            go_on=True
                            for i in range(max(-10,-count),min(10,len_value_temp-1-count)): # Searching neighbouring dates
                                if dt.parse(date_temp[count+i])==dt.parse(line[0]):
                                    count=count+i
                                    go_on=False
                                    continue
                            if go_on:
                                continue
                            else:
                                date.append(str(dt.parse(line[0])))
                        else:
                            date.append(str(dt.parse(line[0])))
#
                        values['for'].append(float(line[1])-values_temp['for'][count])
                        values['act'].append(float(line[2])-values_temp['act'][count])
                        count=count+1
#
        else:
            file_name=data_dir+type+'/'+type+data_files[0]+location+data_files[1]
            csv_read=csv.reader(open(file_name,'r'))
            for line in csv_read:
                if (len(line)<1) or (line[1]=='#'):
                    continue
                date.append(str(dt.parse(line[0])))
                values['for'].append(float(line[1]))
                values['act'].append(float(line[2]))
        print(file_name)
#
        print('kind ; %s' % kind)
        if kind=='errors':
            print('successfully retrieved data for %s in region: %s \n' % (type,location))
            return date,map_list(values['act'],values['for'])
        elif kind=='forecasts':
            print('successfully retrieved data for %s in region: %s \n' % (type,location))
            return date,values['for']
        elif kind=='actuals':
            print('successfully retrieved data for %s in region: %s \n' % (type,location))
            return date,values['act']
        else:
            print('kind is not correct, here is the retrieved data for %s in region: %s \n' % (type,location))
        return date,values
    else:
        print('foo')
        return ()

### This function maps a function 'f' on two lists, returning a list
# Returns a list:
#   () res : [f(vect1[0],vect2[0]), ... ]
# Arguments:
#   () vect1,vect2, 2 lists to map f on
#   () f: a function
def map_list(vect1,vect2,f=None,action='subtract'):
    if (action=='subtract')&(f is None):
        def f(x,y):
            return x-y;
    elif f is None:
        print('no valid action specified')
        return -1
    diff=len(vect1)-len(vect2)
    if(diff):
        print('possible error: vectors do not have same length')
#
    v1,v2=vect1.copy(),vect2.copy()
    if(diff>0):
        v1=v1[0:-diff]
    elif(diff<0):
        v2=v2[0:diff]
#
    v1.reverse(),v2.reverse()
#
    if type(f(0,0))==list:
        length=len(f(0,0))
        res=[[] for i in f(0,0)]
        while (v1!=[])&(v2!=0):
            temp=(f(v1.pop(),v2.pop()))
            for i in range(length):
                res[i].append(temp[i])
    else:
        res=[]
        while (v1!=[])&(v2!=[]):
            res.append(f(v1.pop(),v2.pop()))
    return res

# arguments:
#   () date: list of ordered dates
#   () errors: list of corresponding values
#   () date_range: 2-elements-list to specify a datetime range for the data.
# returns:
#   () a tuple: (date,errors) shortened to fit the date range
def restrain_to_date_range (date,errors, date_range):
    date=date.copy()
    errors=errors.copy()

    if len(date)!=len(errors):
        return ([],[])
    else:
        date_min,date_max=(max(dt.parse(date[0]),dt.parse(date_range[0])),min(dt.parse(date[len(date)-1]),dt.parse(date_range[1])))
        date_temp=[]
        error_temp=[]
        while(dt.parse(date.pop())>date_max):
            errors.pop()
        last_date=dt.parse(date.pop())
        while(last_date>date_min):
            date_temp.append(str(last_date))
            error_temp.append(errors.pop())
            last_date=dt.parse(date.pop())
        if(last_date==date_min):
            date_temp.append(str(last_date))
            error_temp.append(errors.pop())

        error_temp.reverse()
        date_temp.reverse()
        print('restrained to specified date range')
        return(date_temp,error_temp)



#------------------------------ copula analysis -----------------------------------------------------------------------


### This function plots and return the empirical copula of 2 error forecasts with possible hour offsets
# arguments:
#   () type1,location1,type2,location2: specifies the 2 data range to be compared with
#   () hour_offset: it will compare errors1 at date1 with errors2 at date2+offset (in hours)
#   () redistribute: if different from Identity, it will call a function to redistribute the data according to a certain law (gaussian for example)
#   () parameter: a dictionary specifying the parameters: {'main': ... ,'loc': ... ,'kin': ... } with main in {'date','Wind','Solar'}, loc in {'total','NP','SP'} and kin in {'errors','forecasts','actuals'}
#   () which_data
# Calls to:
#   () prepare_data (returns two lists of errors)
#   () draw_copula (to create and plot the empirical copula points and density)
# Returns a dict as {'unif':[[unif1,unif2],...], 'characteristics': [(mean, 10th percentile , 90th percentile) , ...] ] with:
#   () (unif1,unif2): the error distribution, with uniform marginals
def plot_copula_2d(type1,location1,type2,location2,options,hour_offset=0,date_range=[], fig_nb=1, redistribute='Identity',parameter={},kind='errors',first_hour=None):
    result={'unif':[],'characteristics':[[],[]]}
    returned=prepare_data(type1,location1,type2,location2,options,hour_offset=hour_offset,date_range=date_range,kind=kind,first_hour=first_hour)
    for ret in returned:
        errors1,errors2,date,count,offset=ret

        title='Copula of forecast error \n between (%s in %s) and (%s in %s), \n \n hour offset: %d'% (type1,location1,type2,location2,offset)
        xlabels='%s in %s'% (type1,location1)
        ylabels='%s in %s'% (type2,location2)

        if (type(redistribute)==str)&(redistribute!='Identity'):
            if redistribute=='Gaussian':
                print('this is what a gaussian fitted to the data would give')
                errors1,errors2,count=redistribute_gaussian(errors1,errors2)
                title+=('\n redistributed: gaussian')
            elif redistribute=='Diagonal':
                print('we look at the correlation between the sum and subtraction of the two vectors')
                errors1,errors2,count=redistribute_diagonal(errors1,errors2)
                title+=('\n redistributed: diagonals')
            elif redistribute=='Spline':
                print('this is what a spline copula would give')
                errors1,errors2,count=redistribute_spline(errors1,errors2,options)
                title+=('\n redistributed: spline copula')
            else:
                print('redistribute error: no function for %s; plotting default copula'% redistribute)
                title+=('\n redistributed: NO')
        else:
            title+=('\n redistributed: NO')

        # calls draw_copula to plot the copula and get the results to return
        if parameter=={}:
            res=draw_copula(errors1,errors2,title=title,xlabels=xlabels,ylabels=ylabels,fig_nb=fig_nb)
            result['unif'].append(res)
        else:
            res=draw_copula_parameter(errors1,errors2,date,type1,location1,type2,location2,hour_offset=hour_offset,date_range=date_range,title=title,xlabels=xlabels,ylabels=ylabels,fig_nb=fig_nb, visualize=True, parameter=parameter)
            result['unif'].append(res)

        # what is the upper and lower quantile of the error distribution?
        l=np.array(errors1)
        temp=(np.mean(l),np.percentile(l,10),np.percentile(l,90))
        print('%s mean: %d quantiles 1: (10,%d) (90,%d)' % (kind,temp[0],temp[1],temp[2]))
        result['characteristics'][0].append(temp)
        l=np.array(errors2)
        temp=(np.mean(l),np.percentile(l,10),np.percentile(l,90))
        print('%s mean: %d quantiles 2: (10,%d) (90,%d)' % (kind,temp[0],temp[1],temp[2]))
        result['characteristics'][1].append(temp)
        fig_nb+=1

        print(result['unif'])

    return result


### This function plots and return the empirical copula of 2 error forecasts, but also
### the distribution along the diagonal and anti-diagonal axis
# Returns a list as [ [diag,anti_diag],... ] with:
#   () diag,anti_diag: the distribution along the diagonal and anti_diagonal axis
# arguments:
#   () type1,location1,type2,location2: specifies the 2 data range to compare
#   () hour_offset: it will compare errors1 at date1 with errors2 at date2+offset (in hours)
# Calls to:
#   () plot_copula_2d (returns the empirical copula)
#   () compute_diag (computes the distribution of 2D points along the diagonal and anti-diagonal axis)
def copula_analysis(type1,location1,type2,location2,options,hour_offset=0,date_range=[], fig_nb=1):
    if type(hour_offset)!=list:
        hour_offset=[hour_offset]

    returned=plot_copula_2d(type1,location1,type2,location2,options,hour_offset=hour_offset,date_range=date_range, fig_nb=fig_nb)['unif']
    len_offset=len(hour_offset)
    step=300
    incr=0
    result=[]
    for ret in returned:
        unif1,unif2=ret
        diag,anti_diag=compute_diag(unif1,unif2,step)

        index=[i for i in range(step)]
        f_anti_diag=spline.create_spline(spline.find_spline(index,anti_diag,options,visualize=False))
        temp=options.seg_N
        options.seg_N//=2
        f_diag=spline.create_spline(spline.find_spline(index[0:(step//2)],diag[0:(step//2)],options,visualize=False))
        options.seg_N=temp
        del(temp)

        print(len(hour_offset)-incr-1)

        plt.figure(fig_nb+2*incr+len_offset-1+1)
        plt.title('diagonal distribution, offset: %d' % (hour_offset[incr]))
        val=[f_diag(step//2-abs(step//2-i)) for i in index]
        plt.plot(index,val,color='red')
        plt.plot(diag,'.')

        plt.figure(fig_nb+2*incr+len_offset-1+2)
        plt.title('anti-diagonal distribution, offset: %d' % (hour_offset[incr]))
        val=[f_anti_diag(i) for i in index]
        plt.plot(index,val,color='red')
        plt.plot(anti_diag,'.')

        incr=incr+1
        result.append([diag,anti_diag])


# this function works with copula points:
# it computes quadrants with an equal number of points, and the barycenter of the extreme points in each quadrant
# returns as a list [[x0,x1,..],[y0,y1,..]]:
#   () coordinates of the barycenters
#arguments:
#   () type1,location1,type2,location2: specifies the 2 data range to compare
#   () quantile= specifies the extreme points to be taken into account
#   () nbPoints= number of quadrants
#   () distance function= a norm used to determine extreme points ('l1','l2',linf')
def representative_points(type1,location1,type2,location2,options,quantile=0.1,nbPoints=4,distance_func='l1',hour_offset=0,date_range=[], fig_nb=1):

    # retrieving copula data
    returned=(plot_copula_2d(type1,location1,type2,location2,options,hour_offset=hour_offset,date_range=date_range, fig_nb=fig_nb))['unif']
    if returned!=[]:
        ret=returned.pop()
        [unif1,unif2]=ret
        length=len(unif1)
        if(len(unif2)!=length):
            print('coordinate vectors must be same length')
            return -1

    # computing coor (coordinate along the edge of a 1 length square)
    #           and quant (distance to the center of the square)
        temp1,temp2=unif1.copy(),unif2.copy()
        coor=[]
        quant=[]

        def comp_angle(x,y):
            valAbs=abs(x)+abs(y)
            if valAbs==0:
                return-1
            if y-x<0:
                if y+x<0:
                    res=(0+(1-x/y)/2)
                else:
                    res=(1+(y/x+1)/2)
            else:
                if y+x>0:
                    res=(2+(1-x/y)/2)
                else:
                    res=(3+(y/x+1)/2)
            return res

        def comp_angle_inv(ang):
            num=ang//1
            if num==0:
                return (ang,0)
            elif num==1:
                return (1,ang-num)
            elif num==2:
                return (1-ang+num,1)
            else:
                return (0,max(0,1-ang+num))

        while temp1!=[]:
            x,y=2*temp1.pop()-1,2*temp2.pop()-1
            if(distance_func=='l1'):
                valAbs=l1(x,y)
            elif(distance_func=='l2'):
                valAbs=l2(x,y)
            elif(distance_func=='linf'):
                valAbs=linf(x,y)
            else:
                print('distance function %s is not defined, using l1 norm as default'% distance_func)
                valAbs=l1(x,y)

            quant.append(valAbs)
            coor.append(comp_angle(x,y))

        # sorting the copula points according to coor
        coor.reverse()
        quant.reverse()
        unif1_temp=unif1.copy()
        unif2_temp=unif2.copy()

        all=[]
        for i in range(length):
            all.append([unif1_temp.pop(),unif2_temp.pop(),quant.pop(),coor.pop()])
        def temp_sort(a):
            return a[3]
        all.sort(key=temp_sort)

        # defining the separators of the quadrants
        separators=[]
        for i in range(nbPoints):
            index=int(np.floor((i+1/2)*length/nbPoints))
            separators.append((all[index][3],index))


        def findSep(i):
            res=0
            if i>=separators[nbPoints-1][1]:
                res=0
            else:
                for j in range(nbPoints):
                    if i<separators[j][1]:
                        res=j
                        break
            return res

        points=[[0 for i in range(nbPoints)],[0 for i in range(nbPoints)]]
        quadrants_pts=[[] for i in range(nbPoints)]
        for i in range(length):
            x,y,radius,angle=all.pop()
            index = findSep(length-i-1)
            quadrants_pts[index].append((x,y,radius))

        def temp_sort(a):
            return a[2]

        # computing the barycenters and plotting
        for i in range (nbPoints):
            quadrants_pts[i].sort(key=temp_sort)
            length_temp=int(np.floor(len(quadrants_pts[i])*quantile))
            list_temp=[[],[]]
            for j in range(length_temp):
                x,y,radius=quadrants_pts[i].pop()
                points[0][i]+=x
                points[1][i]+=y
                list_temp[0].append(x)
                list_temp[1].append(y)
            points[0][i]/=length_temp
            points[1][i]/=length_temp
            plt.figure(fig_nb)
            # plt.subplot(211)
            plt.plot(list_temp[0],list_temp[1],'.',color='green')
        plt.plot(points[0],points[1],'+',color='red',markersize=10,markeredgewidth=2)
        for i in range (nbPoints):
            line=comp_angle_inv(separators[i][0])
            plt.plot([line[0],0.5],[line[1],0.5],color='red',linewidth=2)

        print('\n %d quadrants, quantile= %0.3f' %(nbPoints,quantile))
        for i in range(nbPoints):
            print('point of the %d th quadrant: %0.3f , %0.3f'%(i,points[0][i],points[1][i]) )
        print('\n');

        return points
    return []


def best_copula(type1,location1,type2,location2,options,nb_parts=4,nb_points_parts=0,hour_offset=0,date_range=[], fig_nb=1, precision=20):

    # retrieving and preparing the data
    result=prepare_data(type1,location1,type2,location2,options,hour_offset=hour_offset,date_range=date_range)
    if(result==[]):
        print('No data to deal with')
        return -1;
    errors1,errors2,date,count,offset=result.pop()
    copula=['Identity','Gaussian','Gaussian no fit','Spline']

    unif1,unif2=(draw_copula(errors1,errors2,visualize=False))
    length=len(unif1)
    parts=[int(length*i//nb_parts) for i in range(nb_parts)]
    parts.append(length)

    values=[[[],[]]for i in range(len(copula))]

    for i in range(nb_parts-1):
        i+=1
        print(i)
        err1=errors1[parts[i-1]:parts[i]]
        err2=errors2[parts[i-1]:parts[i]]

        if (nb_points_parts<10)|(nb_points_parts>50000):
            l=min(parts[i+1],2*parts[i]-parts[i-1])
        else:
            l=parts[i]+nb_points_parts

        # generating the copula to compare to
        for j in range(len(copula)):
            redistribute=copula[j]
            if (type(redistribute)==str)&(redistribute!='Identity'):
                if redistribute=='Gaussian':
                    e1,e2,count=redistribute_gaussian(err1,err2,length=nb_points_parts)
                if redistribute=='Gaussian no fit':
                    e1,e2,count=redistribute_gaussian(err1,err2,length=nb_points_parts,fit=False)
                elif redistribute=='Diagonal':
                    e1,e2,count=redistribute_diagonal(err1,err2,length=nb_points_parts)
                elif redistribute=='Spline':
                    e1,e2,count=redistribute_spline(err1,err2,options,length=nb_points_parts)
            else:
                lbis=min(parts[i+1],2*parts[i]-parts[i-1])
                e1,e2=errors1[parts[i]:lbis],errors2[parts[i]:lbis]

            unif1,unif2=(draw_copula(e1,e2,visualize=False))
            print('\n \n \n %d: length: %d %d'%(i,sum(unif1),sum(unif2)))


            (values[j][0]).extend(unif1)
            (values[j][1]).extend(unif2)

            # title='Copula of forecast error \n between (%s in %s) and (%s in %s), \n \n hour offset: %d'% (type1,location1,type2,location2,offset)
            # xlabels='%s in %s'% (type1,location1)
            # ylabels='%s in %s'% (type2,location2)
            #
            # if (type(redistribute)==str)&(redistribute!='Identity'):
            #     if redistribute=='Gaussian':
            #         print('this is what a gaussian fitted to the data would give')
            #         title+=('\n redistributed: gaussian')
            #     elif redistribute=='Gaussian no fit':
            #         print('this is what a random gaussian would give')
            #         title+=('\n redistributed: gaussian no fit')
            #     elif redistribute=='Diagonal':
            #         print('we look at the correlation between the sum and subtraction of the two vectors')
            #         title+=('\n redistributed: diagonals')
            #     elif redistribute=='Spline':
            #         print('this is what a spline copula would give')
            #         title+=('\n redistributed: spline copula')
            #     else:
            #         print('redistribute error: no function for %s; plotting default copula'% redistribute)
            #         title+=('\n redistributed: NO')
            # else:

            # if fig_nb<=9:
            #     plt.figure(fig_nb)
            #     plt.plot(unif1,unif2,'.')
            #     plt.xlabel(xlabels)
            #     plt.title(title)
            #     plt.ylabel(ylabels)
            #     fig_nb+=1


    for j in range(len(copula)):
        redistribute=copula[j]
        title='Copula of forecast error \n between (%s in %s) and (%s in %s), \n \n hour offset: %d'% (type1,location1,type2,location2,offset)
        xlabels='%s in %s'% (type1,location1)
        ylabels='%s in %s'% (type2,location2)

        if (type(redistribute)==str)&(redistribute!='Identity'):
            if redistribute=='Gaussian':
                print('this is what a gaussian fitted to the data would give')
                title+=('\n redistributed: gaussian')
            elif redistribute=='Gaussian no fit':
                print('this is what a random gaussian would give')
                title+=('\n redistributed: gaussian no fit')
            elif redistribute=='Diagonal':
                print('we look at the correlation between the sum and subtraction of the two vectors')
                title+=('\n redistributed: diagonals')
            elif redistribute=='Spline':
                print('this is what a spline copula would give')
                title+=('\n redistributed: spline copula')
            else:
                print('redistribute error: no function for %s; plotting default copula'% redistribute)
                title+=('\n redistributed: NO')
        else:
            title+=('\n redistributed: NO')

        plt.figure(fig_nb)
        plt.plot(values[j][0][0:(min(parts[2],2*parts[1])-parts[1])],values[j][1][0:(min(parts[2],2*parts[1])-parts[1])],'.')
        plt.xlabel(xlabels)
        plt.title(title)
        plt.ylabel(ylabels)

        fig_nb+=1

    print('\n########')
    for j in range(len(copula)-1):
        j+=1
        redistribute=copula[j]
        distance=0
        for i in range(nb_parts-1):
            i+=1
            begin=parts[i]-parts[1]
            if (nb_points_parts<10)|(nb_points_parts>50000):
                end=min(parts[i+1],2*parts[i]-parts[i-1])-parts[1]
            else:
                end=begin+nb_points_parts
            distance+=compute_distance_emd(values[0][0][begin:end],values[0][1][begin:end],values[j][0][begin:end],values[j][1][begin:end],precision=precision)

        print('distance between empirical and %s copula: %f' % (redistribute,distance))
    print('########\n')




#------------------------------ adjacent functions for copula analysis ------------------------------------------------


### This function prepares the data for a later copula analysis
# returns:
#   () errors1,errors2 : the errors between forecasts for the two data files
#   () date: the corresponding date (for the first data range, if there is an offset)
#   () count: the length of these 4 vectors
#   () offset: the hour offset between the two data ranges
# calls to:
#   () restrain_to_date_range (to treat a specified date range)
def prepare_data(type1,location1,type2,location2,options,hour_offset=0,date_range=[],kind='errors',first_hour=None):
    date_original1,errors_original1=get_data(type1,location1,kind=kind)
    if (type1==type2)&(location1==location2):
        date_original2,errors_original2=(date_original1.copy(),errors_original1.copy())
    else:
        date_original2,errors_original2=get_data(type2,location2,kind=kind)


    # Iterating on hour_offset
    if type(hour_offset)!=list:
        hour_offset=[hour_offset]
    result=[]
    for offset in hour_offset:
        date_temp1,errors_temp1=(date_original1.copy(),errors_original1.copy())
        date_temp2,errors_temp2=(date_original2.copy(),errors_original2.copy())

        # setting an hour offset to data
        if (offset!=0) & (type(offset)==int):
            date_temp2=[str(dt.parse(date_temp2[i])+relativedelta(hours=offset)) for i in range(len(date_temp2))]

        if len(date_range)==2:
            date_temp1,errors_temp1=restrain_to_date_range(date_temp1,errors_temp1,date_range)
            date_temp2,errors_temp2=restrain_to_date_range(date_temp2,errors_temp2,date_range)

        # checking for corresponding dates
        len1=len(date_temp1)
        len2=len(date_temp2)
        errors1=[]
        errors2=[]
        date=[]
        diff=(dt.parse(date_temp1[0])-dt.parse(date_temp2[0])).seconds//3600 + (dt.parse(date_temp1[0])-dt.parse(date_temp2[0])).days*24
        print('diff %d' % diff)
        count=0
        for i in range(len1):
            if (i+diff>len2-1)|(i+diff<0):
                continue

            d1=date_temp1[i]
            date_match= (d1==date_temp2[i+diff])
            if not date_match:
                for j in range(max(0,i+diff-50),min(len2,i+diff+50)):
                    if d1==date_temp2[j]:
                        diff=j-i
                        date_match=True
                        break
            if date_match:
                date.append(d1)
                errors1.append(errors_temp1[i])
                errors2.append(errors_temp2[i+diff])
                count=count+1


        isfh,first_hour=is_list_of(first_hour,child=int)
        if isfh:
            b=[]
            for i in date:
                b.append(dt.parse(i).hour in first_hour)
            [errors1,errors2,date]=remove_in_list([errors1,errors2,date],b,when=[False])
            count=len(errors1)

        result.append((errors1,errors2,date,count,offset))
    return result


### This function plots and return the empirical copula of 2 error forecasts
# arguments:
#   () vect1,vect2: specifies the 2 data range to be compared with
# Returns a list whose elements are:
#   () (unif1,unif2): the error distribution, with uniform marginals
def draw_copula(vect1,vect2,title='',xlabels='',ylabels='',fig_nb=1, visualize=True, parameter={}):
    count=len(vect1)
    print('count : %d '% count)
    if(len(vect2)!=count):
        print("vect1 and vect2 must have the same length")
        return []
    else:
        index1=sorted(range(count), key=vect1.__getitem__,reverse=False)
        index2=sorted(range(count), key=vect2.__getitem__,reverse=False)

        unif1=[0 for i in range(count)]
        unif2=[0 for i in range(count)]

        print('len pour 1: %d len pour 2: %d count: %d' % (len(index1),len(index2),count))

        for i in range(count):
            unif1[index1[i]]=i/count
            unif2[index2[i]]=i/count

        # estimating correlation
        correlation=(sum([unif1[i]*unif2[i] for i in range(count)])-0.25*count)/(((count+1)*(2*count+1)/(count*6))-0.25*count) # divided by the sum of squares formula
        title=title+('\n estimated correlation: %.4f' % correlation)


        #estimating copula density
        #
        # dim=int(np.floor(np.sqrt(count/10)))
        # density=[[0 for i in range(dim)] for j in range(dim)]
        # temp1=unif1.copy()
        # temp2=unif2.copy()
        # for i in range(count):
        #     x=int(np.floor(temp1.pop()*dim))
        #     y=int(np.floor(temp2.pop()*dim))
        #     density[x][y]=density[x][y]+1
        # del(temp1,temp2)


        # plotting copula points
        if visualize:
            rcParams['figure.figsize'] = 7,7


            plt.figure(fig_nb)
            # plt.subplot(211)
            plt.title(title)
            plt.plot(unif1,unif2,'.')
            plt.xlabel(xlabels)
            plt.ylabel(ylabels)
            axes = plt.gca()
            axes.set_xlim([0,1])
            axes.set_ylim([0,1])


            # plotting copula density
            # plt.subplot(212)
            # ax=plt.subplot(212)
            # patches=[]
            # inv=1.0/dim
            # inv_max_density=1/np.max(density)
            # i=0
            # for den in density:
            #     j=0
            #     for val in den:
            #         col=1-val*inv_max_density*0.8
            #         rect = mpatches.Rectangle([i,j],inv,inv, color=(col,col,col))
            #         ax.add_patch(rect)
            #         patches.append(rect)
            #         j=j+inv
            #
            #     i=i+inv
            #
            # plt.xlabel(xlabels)
            # plt.ylabel(ylabels)
            # axes = plt.gca()
            # axes.set_xlim([0,1])
            # axes.set_ylim([0,1])

            rcParams['figure.figsize'] = 7,6

        return (unif1,unif2)


### This function plots and return the empirical copula of 2 error forecasts, coloured according to some parameter (date,forecasts,...)
# arguments:
#   () vect1,vect2: specifies the 2 data range to be compared with
#   () parameter: a dictionary specifying the parameters: {'main': ... ,'loc': ... ,'kin': ... } with main in {'date','Wind','Solar'}, loc in {'total','NP','SP'} and kin in {'errors','forecasts','actuals'}
# Returns:
#   () (unif1,unif2): the error distribution, with uniform marginals
def draw_copula_parameter(vect1,vect2,date,type1,location1,type2,location2,hour_offset=0,date_range=[],title='',xlabels='',ylabels='',fig_nb=1, visualize=True, parameter={}):
    vect1,vect2=vect1.copy(),vect2.copy()
    count=len(vect1)
    print('count : %d '% count)
    if(len(vect2)!=count):
        print("vect1 and vect2 must have the same length")
        return []
    else:
        index1=sorted(range(count), key=vect1.__getitem__,reverse=False)
        index2=sorted(range(count), key=vect2.__getitem__,reverse=False)

        unif1=[0 for i in range(count)]
        unif2=[0 for i in range(count)]

        print('len pour 1: %d len pour 2: %d count: %d' % (len(index1),len(index2),count))

        for i in range(count):
            unif1[index1[i]]=i/count
            unif2[index2[i]]=i/count

        # estimating correlation
        correlation=(sum([unif1[i]*unif2[i] for i in range(count)])-0.25*count)/(((count+1)*(2*count+1)/(count*6))-0.25*count) # divided by the sum of squares formula
        title=title+('\n estimated correlation: %.4f' % correlation)


        # initializing param
        param=[]
        dic_main={'date','Wind','Solar'}
        dic_loc={'total','NP','SP'}
        dic_kin={'errors','forecasts','actuals'}
        keys=parameter.keys()

        # getting parameters from 'parameter'
        if (parameter!={})&('main' in keys):
            if parameter['main']=='date':
                temp=date
                temp.reverse()
                t0=dt.parse('01/01/2010 00:00')
                while temp!=[]:
                    t=temp.pop()
                    param.append((dt.parse(t)-t0).days%365)
                del(t)
                title=title+('\n parameter: date')
            else:
                typ=type1
                loc=location1
                kin='forecasts'

                if 'loc' in keys:
                    if parameter['loc'] in dic_loc:
                        typ=parameter['loc']
                    else:
                        print('Wrong \'loc\' parameter, using %s' % loc)
                else:
                    print('No \'loc\' parameter, using %s' % loc)

                if parameter['main'] in dic_main:
                    if parameter['main']=='Solar':
                        loc='Solar'
                    else:
                        loc='Wind'

                if 'kind' in keys:
                    if parameter['kind']=='errors':
                        kin='errors'
                    elif parameter['kind']=='actuals':
                        kin='actuals'
                    elif parameter['kind']=='forecasts':
                        kin='forecasts'
                    else:
                        print('wrong \'kind\' parameter, using %s' % kin)
                else:
                    print('No \'kind\' parameter, using %s' % kin)
                title=title+('\n parameter: %s %s in %s'%(typ,kin,loc))

                # retrieving the data for param, and checking that the dates correspond to those in 'date'
                if (typ==type1)&(loc==location1)&(kin=='errors'):
                    param=vect1
                elif (typ==type2)&(loc==location2)&(kin=='errors'):
                    param=vect2
                else:
                    date_temp,param_temp=get_data(type=typ,location=loc,kind=kin)
                    # date checking
                    len_temp=len(date_temp)
                    print(date_temp)
                    difference=dt.parse(date[0])-dt.parse(date_temp[0])
                    diff=difference.seconds//3600+difference.days*24
                    param_temp_bis=[]
                    for i in range(len(date)):
                        current=date[i]
                        if dt.parse(current)==dt.parse(date_temp[i+diff]):
                            param_temp_bis.append(param_temp[i+diff])
                        else:
                            found=False
                            for j in range(max(0,i+diff-50),min(len_temp,i+diff+50)):
                                if dt.parse(current)==dt.parse(date_temp[j]):
                                    param_temp_bis.append(param_temp[i+diff])
                                    diff=j-i
                                    found=True
                                    break
                            if not found:
                                param_temp_bis.append(None)
                    param_temp_bis.reverse()
                    unif1_temp=unif1.copy()
                    unif2_temp=unif2.copy()
                    unif1_temp.reverse()
                    unif2_temp.reverse()
                    unif1=[]
                    unif2=[]
                    while(param_temp_bis!=[]):
                        t=param_temp_bis.pop()
                        if t is None:
                            unif1_temp.pop()
                            unif2_temp.pop()
                        else:
                            param.append(t)
                            unif1.append(unif1_temp.pop())
                            unif2.append(unif2_temp.pop())
        elif parameter!={}:
            print('wrong parameter, missing \'main\' \n (try using this format: {\'main\': ... ,\'loc\': ... ,\'kin\': ... } )')
        else:
            print('empty parameter')
            return -1
        countBis=len(param)
        par_max=max(param)
        par_min = min(param)
        # param is initialized

        print(param)
        # plotting copula points
        if visualize:

            rcParams['figure.figsize'] = 7,7
            plt.figure(fig_nb)
            # plt.subplot(211)
            plt.title(title)
            col=[]
            typ_col='default'
            intervalls_col=[]
            if (parameter!={})&('main' in keys)&(parameter['main']=='date'):
                col=color_value(param,typ=typ_col,date=True,intervalls=intervalls_col)
            else:
                col=color_value(param,typ=typ_col,date=False,intervalls=intervalls_col)

            [unif1,unif2,col]=remove_in_list([unif1,unif2,col],col,when=['transparent'])
            plt.scatter(unif1,unif2,color=col)
            plt.xlabel(xlabels)
            plt.ylabel(ylabels)
            axes = plt.gca()
            axes.set_xlim([0,1])
            axes.set_ylim([0,1])


            # plotting correlation as a function of param
            # plt.subplot(212)
            # precision=500
            # smoothing=5
            # correlation=[0 for i in range(precision)]
            # number=[0 for i in range(precision)]
            # cor_final=[0 for i in range(precision)]
            # for i in range(countBis):
            #     fact=int(min(max((param[i]-par_min)/(par_max-par_min)*precision//1,0),precision-1))
            #     correlation[fact]+=unif1[fact]*unif2[fact]
            #     number[fact]+=1
            # for i in range(precision):
            #     if number[i]!=0:
            #         correlation[i]/=number[i]
            # for i in range(precision):
            #     min_temp,max_temp=max(0,i-smoothing),min(precision-1,i+smoothing)
            #     if sum(number[min_temp:max_temp])!=0:
            #         cor_final[i]=np.average(correlation[min_temp:max_temp],weights=number[min_temp:max_temp])
            #
            # plt.bar([i*(par_max-par_min)+par_min for i in range(precision)],cor_final)
            # plt.xlabel(xlabels)
            # plt.ylabel(ylabels)
            # # axes = plt.gca()
            # # axes.set_xlim([0,1])
            # # axes.set_ylim([0,1])

            # rcParams['figure.figsize'] = 7,6

        return (unif1,unif2)


### removes values in 'l' when the value of 'k' is in 'when'
# returns a list of list:
#   () res : same structure as k if k is a list of list, a list of list anyway
# arguments:
#   () l: a list or a list of lists (in this case children of same length as k)
#   () k: a list to specify when to remove elements
def remove_in_list(l,k,when=[],reverse=False):
    temp=k.copy()
    templ=l.copy()
    if (type(l)!= list ):
        print('first argument must be of type list')
        return-1
    elif type(l[0])!=list:
        templ=[templ]
    length=len(k)
    nblists=len(templ)
    for i in range(nblists):
        if len(l[i])!=length:
            print('first argument\'s childs must of same length than second argument (num %i: %d, 2nd arg: %d)'%(i,len(l[i]),length) )
            return -1
    res=[[] for i in range(nblists)]
    while True:
        b= temp.pop() in when
        if reverse:
            b=not b
        if b:
            for i in range(nblists):
                templ[i].pop()
        else:
            for i in range(nblists):
                res[i].append(templ[i].pop())
        if temp==[]:
            break

    for i in res:
        i.reverse()
    return res


def is_list_of(l,child=int):
    res=False
    if l is not None:
        if type(l)==child:
            l=[l]
        if type(l)==list:
            if len(l)>0:
                res=True
                for i in l:
                    if type(i)!=child:
                        res=False
                        break
    return (res,l)

### returns a color list corresponding to 'vect' values
# returns:
#   () col : list of 3-elements tuples (->rgb colors)
# arguments:
#   () vect : values to match the color to
#   () typ : type of color pattern (if int <0.5 and>=0, exclude the values in the interval 10% - 90%)
#   () intervalls : to get constant colors within the intervalls
def color_value(vect,typ='default',date=False,intervalls=[]):
    vect=vect.copy()
    length=len(vect)
    col=[]
    ma,mi=max(vect),min(vect)

    if intervalls!=[]:
        b=True
        t=0
        for i in intervalls:
            print(b)
            print(i)
            b&=(type(i)==float)&(i>=t)
            t=i
        b&=t<=1
        if b:
            def norm(x):
                x=(x-mi)/(ma-mi)
                t=0
                inter=[0]
                inter.extend(intervalls)
                for i in range(len(inter)):
                    if inter[i]<x:
                        t=i
                    else:
                        break
                inter.append(1)
                t=(inter[t]+inter[t+1])/2
                return t
        else:
            def norm(x):
                return (x-mi)/(ma-mi)
    else:
        def norm(x):
            return (x-mi)/(ma-mi)

    if date:
        vect=[2*min(i-mi,ma-i) for i in vect]


    if (type(typ)==float):
        if (typ<0.5)&(typ>0):
            for i in vect:
                if norm(i)<typ:
                    col.append('blue')
                elif norm(i)>1-typ:
                    col.append('red')
                else:
                    col.append('transparent')
        else:
            print('error in argument \'typ\'')
            return -1

    elif typ=='old':
        for i in vect:
            fact=norm(i)
            col.append((0.5*fact+ 2*fact*(1-fact), 2*fact*(1-fact) , 0.5-0.5*fact+2*fact*(1-fact)) )
    elif typ=='blueshade':
        for i in vect:
            fact=norm(i)
            col.append((1-fact,1-fact,1))
    else:
        for i in vect:
            fact=norm(i)
            col.append((fact**2,2*fact*(1-fact),(1-fact)**2))

    return col

### This function computes the distribution of points along the diagonal and anti-diagonal axis
# Returns a tuple whose elements are:
#   () diag: a vector containing the number of points of pi/4 oriented bands
#   () anti_diag: a vector containing the number of points of -pi/4 oriented bands
# arguments
#   () (unif1,unif2): two vectors of same length (-> 2d points)
def compute_diag(unif1,unif2,step,min=0,max=1):
    unif1=unif1.copy()
    unif2=unif2.copy()
    length=len(unif1)
    diag=[0 for i in range(step)]
    anti_diag=[0 for i in range(step)]
    for i in range(length):
        x,y=(unif1.pop(),unif2.pop())


        fact=float(1-np.abs(1-(x+y))) # normalizing factor along the anti diagonal axis
        if (x+y>2*min)&(x+y<2*max):
            temp=int(np.floor(((x-y)/(2*fact)+0.5)*step))
            if temp==step:
                temp=step-1
                print('care: unif1=1 and unif2=0 at index %d' % (step-i-1))
            diag[temp]=diag[temp]+1

            temp=int(np.floor((x+y)/2*step))
            if temp==step:
                temp=step-1
                print('care: unif1=1 and unif2=1 at index %d' % (step-i-1))
            anti_diag[temp]=anti_diag[temp]+1

    return (diag,anti_diag)


# Returns a function [(x,y)->int]
#   () f: the spline 2D-repartition function
# arguments
#   () (unif1,unif2): two vectors of same length (-> 2d points)
def compute_distribution(unif1,unif2):
    count=len(unif1)
    if count==len(unif1):
        dim=max(1,int(np.floor(np.sqrt(count/16))))
        distr=[[0 for i in range(dim)] for j in range(dim)]
        points=[[[] for i in range(dim)] for j in range(dim)]
        temp1=unif1.copy()
        temp2=unif2.copy()
        for i in range(count):
            a=temp1.pop()
            b=temp2.pop()
            x=min(max(0,int(np.floor(a*dim))),dim-1)
            y=min(max(0,int(np.floor(b*dim))),dim-1)
            distr[x][y]=distr[x][y]+1
            points[x][y].append((a,b))
        del(temp1,temp2)
        for i in range(dim-1):
            distr[1][i+1]+=distr[1][i]
            distr[i+1][1]+=distr[i][1]
        for i in range(dim-1):
            for j in range(dim-1):
                distr[i+1][j+1]+=distr[i][j+1]+distr[i+1][j]-distr[i][j]

        def f(a,b):
            x=min(max(0,int(np.floor(a*dim))),dim-1)
            y=min(max(0,int(np.floor(b*dim))),dim-1)
            res=distr[x][y]
            for t in points[x][y]:
                if(a>=t[0])&(b>=t[1]):
                    res+=1
            return res/count
        return f
    else:
        print('unif1 and unif2 must be same length')
        return -1

# Returns a float:
#   () res: the L1-distance between the two repartition function, as estimated from the points
# arguments
#   () (unif_a1,unif_a2): two vectors of same length (-> 2d points)
#   () (unif_b1,unif_b2): two other vectors of same length (-> 2d points)
def compute_distance(unif_a1,unif_a2,unif_b1, unif_b2, precision=300):
    fa=compute_distribution(unif_a1,unif_a2)
    fb=compute_distribution(unif_b1,unif_b2)
    res=0
    for i in range(precision):
        for j in range(precision):
            res+=abs(fa(i/precision,j/precision)-fb(i/precision,j/precision))
    return res/precision**2

    # len1=len(unif_a1)
    # len2=len(unif_b1)
    # if(len1!=len(unif_a2))|(len2!=len(unif_b2)):
    #     print('X and Y coordinates must be same length')
    #     return -1
    #
    # if len1>len2:
    #     unif_a1,unif_a2=unif_a1[0:len2],unif_a2[0:len2]
    # elif len1<len2:
    #     unif_b1,unif_b2=unif_b1[0:len2],unif_b2[0:len2]
    #
    # def aux(prec,range):
    #     s1=sum([(unif_a1[i]>=range[0][0])&(unif_a1[i]<range[0][1])&
    #                 (unif_a2[i]>=range[1][0])&(unif_a2[i]<range[1][1]) for i in range(len1)])
    #     s2=sum([(unif_b1[i]>=range[0][0])&(unif_b1[i]<range[0][1])&
    #                 (unif_b2[i]>=range[1][0])&(unif_b2[i]<range[1][1]) for i in range(len2)])
    #     res=abs(s1-s2)
    #     if prec==0:
    #         return res
    #     else:
    #         wid=(range[0][1]-range[0][0])/2
    #         res+=0.5*aux(prec-1,[[range[0][0],range[0][0]+wid],[range[1][0],range[1][0]+wid]])
    #         res+=0.5*aux(prec-1,[[range[0][0]+wid,range[0][1]],[range[1][0],range[1][0]+wid]])
    #         res+=0.5*aux(prec-1,[[range[0][0],range[0][0]+wid],[range[1][0]+wid,range[1][1]]])
    #         res+=0.5*aux(prec-1,[[range[0][0]+wid,range[0][1]],[range[1][0]+wid,range[1][1]]])
    #         return res
    # return aux(precision,[[0,1],[0,1]])


# Returns a float:
#   () res: the L1-distance between the two repartition function, as estimated from the points
# arguments
#   () (unif_a1,unif_a2): two vectors of same length (-> 2d points)
#   () (unif_b1,unif_b2): two other vectors of same length (-> 2d points)
def compute_distance_emd(a1,a2,b1, b2, precision=10):
    # with the earth mover distance
    counta=len(a1)
    if (len(a2)!=counta):
        print('given vectors must be same length, here: %d and %d' % (len(a1),len(a2)))
        return -1
    countb=len(b1)
    if (len(b2)!=countb):
        print('given vectors must be same length, here: %d and %d' % (len(b1),len(b2)))
        return -1
    prec2=precision**2
    distance=np.array([[float(np.sqrt((i%precision-j%precision)**2+(i//precision-j//precision)**2))/precision for j in range(precision**2)] for i in range(precision**2)])
    sa=np.array([0. for i in range(prec2)])
    sb=np.array([0. for i in range(prec2)])
    inva=float(1/counta)
    for i in range(counta):
        index=min(int(np.floor(a1[i]*precision)),precision-1)+precision*min(int(np.floor(a2[i]*precision)),precision-1)
        sa[index]=sa[index]+inva
#
    invb=float(1/countb)
    for i in range(countb):
        index=min(int(np.floor(b1[i]*precision)),precision-1)+precision*min(int(np.floor(b2[i]*precision)),precision-1)
        sb[index]= sb[index]+invb
    return emd(sa,sb,distance)


#------------------------------ redistributing functions for copula analysis -------------------------------------------


# arguments:
#   () vect1,vect2: specifies the 2 data range of same length to fit the normal law
# returns:
#   () points distributed according to a normal law, fitted to the input data
#   () the length of the output vectors
def redistribute_gaussian(vect1,vect2,length=0,fit=True):

    if (length>50000)|(length<10):
        count=len(vect1)
    else:
        count=length
    if fit:
        cov=np.cov(np.matrix([vect1,vect2]))
        print(cov)
    else:
        cov=np.array([[1,0],[0,1]])
        print(cov)
    print('in gaussian, count=%d' % count)
    redistributed=np.transpose(np.random.multivariate_normal([0,0], cov, count))

    return (redistributed[0],redistributed[1],count)


# arguments:
#   () vect1,vect2: specifies the 2 data range of same length
# returns:
#   () (vect1+vect2,vect1-vect2)
#   () the length of the output vectors
def redistribute_diagonal(vect1,vect2,length=0):
    vect1=vect1.copy()
    vect2=vect2.copy()

    count=len(vect1)
    if(len(vect2)!=count):
        raise ValueError("vect1 and vect2 must have same first dimension")
    else:
        temp1=[]
        temp2=[]
        for i in range(count):
            x=vect1.pop()
            y=vect2.pop()
            temp1.append(x+y)
            temp2.append(x-y)
        if (length>50000)|(length<10):
            return(temp1,temp2,count)
        else:
            t1,t2=[],[]
            for i in range(length//count):
                t1.extend(temp1)
                t2.extend(temp2)
            t1.extend(temp1[0:(length%count)])
            t1.extend(temp1[0:(length%count)])
            print('length: %d , count: %d'% (length,count))
            return(t1,t2,length)


# arguments:
#   () vect1,vect2: specifies the 2 data range of same length
# returns:
#   () x,y (lists): points distributed according to an empiricial copula
#   () the length of the output vectors
def redistribute_spline_old(vect1,vect2,options):

    step=300
    step+=step%2 #make sure step is even
    length=len(vect1)
    unif1,unif2=draw_copula(vect1,vect2,visualize=False)
    diag,anti_diag=compute_diag(unif1,unif2,step)
    length_d,length_ad=(len(diag),len(anti_diag))

    # sum
    diag[0]+=diag[length_d-1]
    for i in range(length_d//2-1):
        diag[i+1]+=diag[i]+diag[length_d-i-2]
    for i in range(length_ad-1):
        anti_diag[i+1]+=anti_diag[i]

    # renormalize
    max=diag[length_d//2-1]
    for i in range(length_d//2):
        diag[i]/=max
    max=anti_diag[length_ad-1]
    for i in range(length_ad):
        anti_diag[i]/=max

    plt.figure(9)
    plt.plot(diag[1:150])
    plt.figure(8)
    plt.plot(anti_diag)


    index=[i/step+0.5/step for i in range(step)]
    x=[0];x.extend(anti_diag)
    y=[0];y.extend(index)
    f_anti_diag=spline.create_spline(spline.find_spline(x,y,options,visualize=False))
    temp=options.seg_N
    options.seg_N//=2

    x=[0];y=[0]
    for i in range(step//2):
        x.append(diag[i])
        y.append(index[i]*2)

    f_diag=spline.create_spline(spline.find_spline(x,y,options,visualize=False))
    options.seg_N=temp
    del(temp)

    temp_x=[]
    temp_y=[]
    count=len(vect1)
    AD,D=[],[]
    u1,u2,b=(np.random.uniform(0,1,count),np.random.uniform(0,1,count),2*np.random.binomial(1,0.5,count)-1)
    for i in range(count):
        ad=f_anti_diag(u1[i])
        d=f_diag(u2[i])
        temp_x.append(ad+b[i]*(1-d)/2*(1-2*abs(0.5-ad)))
        temp_y.append(ad-b[i]*(1-d)/2*(1-2*abs(0.5-ad)))

        AD.append(ad)
        D.append(d)

    plt.figure(10)
    plt.plot(AD,'.')
    plt.figure(11)
    plt.plot(D,'.')

    return(temp_x,temp_y,count)

# idem old
def redistribute_spline(vect1,vect2,options,length=0):

    step=300
    step+=step%4 #make sure step is even
    if (length>50000)|(length<10):
        length1=len(vect1)
    else:
        length1=length
    unif1,unif2=draw_copula(vect1,vect2,visualize=False)
    length_d,length_ad=(step//4,step)

    # compute anti diagonal density
    anti_diag=compute_diag(unif1,unif2,step)[1]
    max_ad=max(anti_diag)
    anti_diag=[i/max_ad for i in anti_diag]
    index=[(i+0.5)/length_ad for i in range(length_ad)]

    f_anti_diag=(spline.create_spline(spline.find_spline(index,anti_diag,options,visualize=False)))

    # plt.figure(9)
    # plt.plot(range(100),[f_anti_diag(i/100) for i in range(100)])
    # plt.figure(8)
    # plt.hist(anti_diag,50)

    # compute diagonal density coefficients for different 'anti-diagonal bands'
    index=[(i+0.5)/length_d for i in range(length_d)]
    a=[];s0=[];v0=[];dico={}
    nb_bands=5
    bands=[i/nb_bands for i in range(nb_bands)]
    for i in bands:
        diag_temp=(compute_diag(unif1,unif2,step//2,min=i,max=i+1/nb_bands))[0]

        #assuming that diag is symetrical
        diag_bis=diag_temp.copy()
        diag_bis.reverse()
        diag=[]
        for i in range(length_d):
            diag.append(diag_temp.pop()+diag_bis.pop())
        del(diag_temp,diag_bis)

        #normalize
        max_d=max(diag)
        diag=[i/max_d for i in diag]

        dico=spline.find_spline(index,diag,options,visualize=False)
        a.append(dico['a'])
        s0.append(dico['s0'])
        v0.append(dico['v0'])


    count=0
    x,y=[],[]
    D,AD=[],[]
    d,ad,b=0,0,0
    go_on=True
    while(go_on):
        u1,u2=list(np.random.uniform(0,1,100)),list(np.random.uniform(0,1,100))
        v1,v2=list(np.random.uniform(0,1,100)),list(np.random.uniform(0,1,100))
        binom=list(2*np.random.binomial(1,0.5,100)-1)
        go_on_bis=True
        for i in range(100):

            if(count>=length1):
                go_on=False
                break

            # reject method to generate the diagonal and anti_diagonal coordinates (d and ad),
            # b is to determine on which side of the first diagonal the point will be
            # x,y are the coordinates of the generated points to be returned

            #generation of anti-diagonal coordinate
            U1,V1=(u1.pop(),v1.pop())
            if(V1<f_anti_diag(U1)):
                ad=U1

                # estimation of the diagonal density
                def create_dico(dico,i,j,wi,wj):
                    dico['a']=[(wi*a[i][s]+wj*a[j][s])/(wi+wj) for s in range(len(a[0]))]
                    dico['s0']=(wi*s0[i]+wj*s0[j])/(wi+wj)
                    dico['v0']=(wi*v0[i]+wj*v0[j])/(wi+wj)
                    return dico

                if ad<=0.5/nb_bands:
                    dico=create_dico(dico,0,0,1,1)
                elif ad>=1-0.5/nb_bands:
                    dico=create_dico(dico,0,0,1,1)
                else:
                    k=min(int(np.floor(ad*nb_bands-0.5)),nb_bands-2)
                    mod=(ad*nb_bands-0.5)%1
                    dico=create_dico(dico,k,k+1,1-mod,mod)
                f_diag=spline.create_spline(dico)

                # generation of diagonal coordinate
                point_added=False
                while(not point_added):
                    if(u2==[])|(v2==[])|(binom==[]):
                        go_on_bis=False
                        break
                    U2,V2=(u2.pop(),v2.pop())

                    b=binom.pop()

                    if(V2<f_diag(U2)):
                        d=U2
                        point_added=True
                        count+=1

                # now we compute x and y
                if point_added:
                    D.append((1-d)*b/2+0.5)
                    AD.append(ad)

                    x.append(ad+b*(1-d)*(0.5-abs(0.5-ad)))
                    y.append(ad-b*(1-d)*(0.5-abs(0.5-ad)))

            if(not go_on_bis):
                break

    # plt.figure(10)
    # plt.hist(AD,50)
    # plt.figure(11)
    # plt.hist(D,50)
    return (x,y,count)


#------------------------------ distance functions for representative points -------------------------------------------

def linf(x,y):
    return max(abs(x),abs(y))

def l1(x,y):
    return abs(x)+abs(y)

def l2(x,y):
    return x**2+y**2
