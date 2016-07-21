import error_correlation as ec
import matplotlib.pyplot as plt
import sys
import dateutil.parser as dt
from dateutil.relativedelta import relativedelta

home_dir='/home/ambroiseidoine/UCD/'
sys.path.append(home_dir+'code/prescient/release/Prescient_1.0/')

import exec.MasterOptions as MasterOptions

##############################
##  This section is to run  ##
##############################

### time correlation

def test(options,action):
    fn=0;

    if 'solar correlation' in action:
        ec.time_correlation('Solar','total',fig_nb=fn); fn=fn+1
        ec.time_correlation('Solar','SP',fig_nb=fn); fn=fn+1
        ec.time_correlation('Solar','NP',fig_nb=fn); fn=fn+1

    if 'seasons comparison' in action:
        # ec.time_correlation('Wind','total',date_range=['07/01/14 00:00','10/01/14 00:00'],fig_nb=fn); fn=fn+1
        # ec.time_correlation('Wind','total',date_range=['01/01/15 00:00','03/01/15 00:00'],range_max=24*10,fig_nb=fn); fn=fn+1

        nb_period=6
        date1=dt.parse('07/01/14 00:00');date2=dt.parse('06/30/15 23:00')
        delta=(date2-date1).days//nb_period
        dates_list=[str(date1+relativedelta(days=i*delta)) for i in range(nb_period)]
        dates_list.append(str(date2))
        print(dates_list)
        offsets=[5]
        for i in range(nb_period):
            ec.plot_copula_2d('Wind','NP','Wind','NP',options,date_range=[dates_list[i],dates_list[i+1]],hour_offset=offsets,fig_nb=fn);fn=fn+len(offsets);

    if 'time correlation' in action:
        ec.time_correlation('Wind','total',fig_nb=fn); fn=fn+1
        ec.time_correlation('Wind','SP',fig_nb=fn); fn=fn+1
        ec.time_correlation('Wind','NP',fig_nb=fn); fn=fn+1

    if 'space correlation' in action:
        ec.space_correlation('Wind')
        ec.space_correlation('Solar')

    if 'copula wind' in action:
        offsets=[11,23,24,35]
        ec.plot_copula_2d('Wind','NP','Wind','NP',options,hour_offset=offsets,fig_nb=fn);fn=fn+len(offsets)

    if 'copula approximation' in action:
        offsets=[5]
        # unif1,unif2=(ec.plot_copula_2d('Wind','NP','Wind','NP',options,hour_offset=offsets,fig_nb=fn))[0][0];fn=fn+len(offsets);
        # unif1_emp,unif2_emp=(ec.plot_copula_2d('Wind','NP','Wind','NP',options,hour_offset=offsets,fig_nb=fn, redistribute='Spline'))[0][0];fn=fn+len(offsets)
        # unif1_gau,unif2_gau=(ec.plot_copula_2d('Wind','NP','Wind','NP',options,hour_offset=offsets,fig_nb=fn,redistribute='Gaussian'))[0][0];fn=fn+len(offsets)
        #
        # print('##### \n distance between gaussian and empirical: %f \n' %ec.compute_distance_emd(unif1,unif2,unif1_gau,unif2_gau,precision=30))
        # print('##### \n distance between spline and empirical: %f \n' %ec.compute_distance_emd(unif1,unif2,unif1_emp,unif2_emp, precision=30))

        ec.best_copula('Wind','NP','Wind','NP',options,nb_parts=6 ,nb_points_parts=2000,hour_offset=offsets,fig_nb=fn,precision=20);fn=fn+3*len(offsets)

    if 'copula spline estimation' in action:
        offsets=[2]
        ec.copula_analysis('Wind','NP','Wind','NP',options,hour_offset=offsets,fig_nb=fn);fn=fn+3*len(offsets);
        ec.plot_copula_2d('Wind','NP','Wind','NP',options,hour_offset=offsets,fig_nb=fn, redistribute='Spline');fn=fn+len(offsets)
        ec.plot_copula_2d('Wind','NP','Wind','NP',options,hour_offset=offsets,fig_nb=fn,redistribute='Gaussian');fn=fn+len(offsets)

    if 'representative points' in action:
        offsets=[11,17]
        list=ec.representative_points('Solar','NP','Solar','NP',options,nbPoints=8,distance_func='l2',hour_offset=offsets);fn+=len(offsets)

    if 'parameter dependence' in action:
        offsets=[1,2,4,8,16,24]
        # ec.plot_copula_2d('Wind','NP','Wind','NP',options,hour_offset=offsets,fig_nb=fn,parameter={'main':'date'});fn=fn+len(offsets);
        # ec.plot_copula_2d('Wind','NP','Wind','NP',options,hour_offset=offsets,fig_nb=fn,parameter={'main':'wind','kind':'forecasts'});fn=fn+len(offsets);
        # ec.plot_copula_2d('Wind','NP','Wind','NP',options,hour_offset=offsets,fig_nb=fn,parameter={'main':'wind','kind':'actuals'});fn=fn+len(offsets);
        # offsets=[11,23,24,35]
        # offsets=[11]
        res=[]
        res.append(ec.plot_copula_2d('Wind','NP','Wind','NP',options,hour_offset=offsets,fig_nb=fn,kind='errors')['characteristics']);fn+=len(offsets)

        print(res)
        for i in res:
            print(i)

    plt.show()


##############################
## This section is to check ##
##############################

def foo(options):
    fn=0;
    # u1=[0.1,0.11,0.13,0.2,0.35,0.36,0.4,0.41,0.43,0.7,0.81,0.95]
    # u1=[i for i in range(20)]
    # u2=[1 for i in u1]
    # col=ec.color_value([i for i in range(len(u1))],typ='blueshade',intervalls=[0.,0.5,0.9])
    # [u1,u2,col]=ec.remove_in_list([u1,u2,col],col,when=['transparent'])
    # u2.reverse()
    # f=ec.compute_distribution(u1,u2)
    # print(ec.remove_in_list([[2,-2,34,2,-1],[1,2,3,4,5]],[1,2,3,4,5],when={3,4}))
    offsets=[11,23,24,35]
    ec.plot_copula_2d('Wind','NP','Wind','NP',options,hour_offset=offsets,fig_nb=fn,kind='forecasts');fn=fn+len(offsets)
    plt.show()


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


    diff_actions={'solar correlation', 'seasons comparison', 'time correlation', 'space correlation', 'copula wind',
        'copula approximation','copula spline estimation', 'representative points'}


    # test(options,{'seasons comparison'})
#     test(options,{'copula approximation'})
#     test(options,{'copula spline estimation'})
#     test(options,{'representative points'})
#     foo(options)
    print(options)


# MAIN ROUTINE STARTS NOW #

main(sys.argv)

