import csv
import math
import sys
import time
from datetime import timedelta
from copy import deepcopy
import os

import copulaMS as cms
import copulaModels as mod
import dateutil.parser as dt
import matplotlib.pyplot as plt
import numpy as np
import subsetModel as sub
import utilities as ut
import vines
from dateutil.relativedelta import relativedelta

home_dir='c:\\users\\sabrina\\documents\\research\\code for real user\\'

### this function fetches the data to create a MS (multi series) copula manager.
### for solar data, it additionally takes into account the solar hour factor
### (normalize the variance to 1, and keep the variance of errors as a function of solar hour in memory)
# returns a MS copula manager
# arguments:
#   () titles: {}
def create_copulaManager(titles, parameters_def, list_parameters=None,
                         just_parameters=None, filename=None, sourcefile=None):
    length = len(titles)
    series = []

    if (not isinstance(titles, list)) and (len(titles) > 0):
        raise (RuntimeError('titles must be a list of dictionaries, with length>0'))
    if (list_parameters is not None) and (len(list_parameters) != length):
        raise (RuntimeError('list_parameters must be either None or of same length than titles'))
    if (just_parameters is not None) and (len(just_parameters) != length):
        raise (RuntimeError('just_parameters must be either None or of same length than titles'))

    opts = {'type': {'Wind', 'Solar'}, 'location': {'NP', 'SP', 'total'}, 'kind': {'error', 'forecast', 'actual'}}
    mandatory_keys = {'type', 'location', 'kind'}
    for title in titles:
        if isinstance(title, dict) and set(title.keys()).issubset(mandatory_keys):
            for key in mandatory_keys:
                if title[key] not in opts[key]:
                    raise RuntimeError('please check the arguments specifying the kind of data')
        else:
            raise RuntimeError('please check the arguments specifying the kind of data')

    for i, title in enumerate(titles):
        index = titles.index(title)
        if index < i: # if repeat of title?
            serie = {}
            for key in series[index]:
                serie[key] = series[index][key].copy()
            series.append(serie)
        else:
            serie = {}
            res = get_data(type=title['type'], location=title['location'], filename=filename, sourcefile=sourcefile)
            if title['type'] == 'Solar':
                # call to take into account the importance of solar hour in solar power forecast errors
                res2 = ut.prepare_solar(res, visualize=False)
                serie['date'] = res2['date']
                serie['vect'] = res2[title['kind']]
                serie['data'] = {'forecast': res2['for'], 'forecast_d': ut.list_derivative(res2['for']),
                                'hour_sol': res2['hour_sol']}
                serie['var_function'] = res2['var_function']
            else:
                serie['date'] = res['date']
                serie['var_function'] = None
                if title['kind'] == 'error':
                    serie['vect'] = [actual - forecast for (actual, forecast) in zip(res['act'], res['for'])]
                elif title['kind'] == 'forecast':
                    serie['vect'] = res['for']
                else:
                    serie['vect'] = res['act']
                serie['data'] = {'forecast': res['for'], 'forecast_d': ut.list_derivative(res['for'])}
            serie['title'] = title
            series.append(serie)

    return cms.CopulaManagerMS(series, parameters_def, list_parameters=list_parameters, just_parameters=just_parameters)


### This function compares the predicting accuracy of the various models
# arguments:
#   () copula is a copulaManagerMS
#   () win_days is a window parameter specifying the number of days before and after the current time of the year
#      that should be taken into account. (if=45, it translates into 91 days in the past years and 45 in the current year)
#   () win_forecast is a window parameter specifying the width of the forecast window: if q=CDF_forecast(current_forecast),
#      we will consider dates at which the forecast was in [CDF^-1(q-'win_forecast'),CDF^-1(q+'win_forecast')]
def test_models(copula, win_days=90, list_models=None, repeat_win_days=False, win_forecast=0.2,
                start_incr=None, end_incr=None, compare_dist=False):
    min_points = 30
    win_days = max(min_points, win_days)

    if start_incr is None:
        start_incr = win_days

    start_time = time.time()

    # keeping the old window parameters in param_fixed
    param_fixed = []
    for par in copula.parameters:
        dic = {'date_range': par['date_range'], 'offsets': par['offsets'], 'first_hour': par['first_hour']}
        param_fixed.append(dic)

    # initializing variables
    parameters = copula.parameters
    nb_series = len(copula.copulae)
    forecasts = [copula.dataM[i]['forecast'][0] for i in range(nb_series)]
    errors = copula.vectM
    dates = copula.dateM

    print('forecasts (%d) %r\nerrors (%d) %r\ndates (%d) %r' % (len(forecasts[0]), forecasts[0][:8],
                                                                len(errors), errors[:8], len(dates), dates[:8]))

    incr = 0
    first = True
    last_time = time.time()
    res = {'len': [], 'log': [], 'rank': [], 'problem': [], 'past_log': [], 'selected_model': [],
           'proj_emd': [], 'proj_quantile': [], 'wei_par': [], 'dates': []}

    res['parameters'] = deepcopy(parameters[0])
    res['parameters'].update({'win_days': win_days, 'win_forecast': win_forecast, 'repeat_win_days': repeat_win_days})

    # loop over each considered hour: each time,
    #       The copula manager is updated to fit the window (using 'win_days' and 'win_forecast')
    #       Models are created to fit the updated copula
    #       The log_likelihood of the observation is computed for all different models
    for forecast, error, date in zip(zip(*forecasts), zip(*errors), dates):

        # selecting the observation range
        incr += 1
        if incr < start_incr:
            continue
        if end_incr is not None:
            if incr >= end_incr:
                continue

        # printing information
        def print_info(last_time=last_time):
            t_print = []
            for i in (time.time() - last_time, time.time() - start_time):
                i = int(i)
                t_print.append((i // 3600, i // 60 % 60, i % 60))
            last_time = time.time()
            print('\n\n   ####################   \n\niter %d: forecast %r, error %r, date %s\n\n'
                  'time elapsed in the last loop: %d:%d:%d, time since start: %d:%d:%d'
                  '\n\n   ####################   \n\n'
                  % (incr, forecast, error, date, t_print[0][0], t_print[0][1], t_print[0][2], t_print[1][0],
                     t_print[1][1], t_print[1][2]))
            return last_time

        last_time = print_info(last_time=last_time)

        observations = select_observations(copula, date, repeat_win_days=repeat_win_days,
                                           win_days=win_days, win_forecast=win_forecast, param_fixed=param_fixed)

        if copula.lengthM < min_points:
            print('Error updating copula at incr %d, too few points' % incr, file=sys.stderr)
            continue

        if observations is None:
            continue
        else:
            vects, unifs = observations

        # fitting models to the distribution
        # print('fitting models to the distribution')
        length = len(vects[0])
        dim = len(vects)

        # creating the density of a fitted gaussian
        densities = [ut.create_gaussian_density(vects)]

        # creating a list of copula models
        if list_models is None:
            list_models = [mod.cop_gaussian, mod.cop_student]
        models = []
        for i, model in enumerate(list_models):
            # print('\n### %d th model###' % i)
            created = model(unifs)
            models.append(created)
            if created.ACR == 'WE':
                res['wei_par'].append(created.get_names())

        # computing the densities of the models, their log-likelihood, and selecting the 'best candidate'
        names = ['gaussian distribution']
        cop_densities = []
        best_model_past = models[0]
        best_log_past = 0
        log_past = []
        for model in models:
            if first:
                names.append(model.name)
            cop_densities.append(model.pdf)
            lld = sum([np.log(k) for k in model.pdf(unifs)])
            log_past.append(lld)
            if lld > best_log_past:
                best_log_past = lld
                best_model_past = model

        names.append('selected model')
        cop_densities.append(best_model_past.pdf)
        log_past.append(best_log_past)

        res['past_log'].append(log_past)
        res['selected_model'].append(best_model_past.name)
        if first:
            res['names'] = names
            first = False

        # computing the rank of 'obs' among the window points
        CDFs = ut.marginals_cdf(vects)
        rank = [float(CDFs[i](error[i])) for i in range(dim)]

        # computing the tail metrics:
        C_to_D = ut.copula_to_distribution(vects)
        simulations = [cop.simulate(10000) for cop in models]

        """
        for i in range(len(models)):
            if not 0 < np.min(simulations[i]) < np.max(simulations[i]) < 1:
                print(i)
                print(simulations[i])
                print(copula.lengthM)
                print(models[i].print_names())
                print(models[i].print_par())
                return models[i]
        """

        simulations = [C_to_D(sim) for sim in simulations]
        tail_metrics = ut.compare_tails(simulations, vects, error, quantile=0.1)
        res['proj_emd'].append(tail_metrics[0])
        res['proj_quantile'].append(tail_metrics[1])

        # computing the log likelihood
        try:
            if compare_dist:
                res_log = [np.log(den(error)[0]) for den in densities]
                res_log.extend(
                    [den(error)[0] for den in ut.copula_to_densities(vects, cop_densities, log_return=True)])
            else:
                res_log = [den(rank)[0] for den in
                           ut.distribution_to_copula_densities(vects, densities, log_return=True)]
                res_log.extend([np.log(den([[r] for r in rank])[0]) for den in cop_densities])

            res['log'].append(res_log)
        except:
            res['log'].append(None)
            res['problem'].append('incr: %d, problem in the log: %r' % (incr, sys.exc_info()[0]))

        res['len'].append(length)
        res['rank'].append(rank)
        res['dates'].append(date)

        # res['problem'].append('incr %d general problem: %r' % (incr, sys.exc_info()[0]))
        # print(incr, sys.exc_info()[0], file=sys.stderr)

    copula.update(param_fixed[0], list_parameters=param_fixed)
    return res


# Visualization and printing of the results
def visualize_result(res, add_title='', save=True, quantile=0.1):
    nb_models = len(res['names']) - 2
    length = len(res['log'])
    dim = int(np.log(np.shape(res['proj_quantile'])[2]) / np.log(2)) + 1

    titles = []
    designations = []
    comparisons = []

    # general log_likelihood comparison
    title = 'Mean log-likelihood of day-ahead forecast errors'
    titles.append(title)
    title = '%s \n%s' % (title, add_title)
    temp = np.mean(res['log'], 0)
    comparisons.append(temp[1:-1])
    designations.append('log-likelihood')
    print(title)
    print(temp)
    ut.sorted_barplot(temp, xlabels=res['names'], title=title)
    if save:
        plt.savefig('%s.png' % title.replace(' ', '_').replace('\n', '_'))

    # emd comparison of the projections (FITTING)
    title = 'EMD comparison of the tails (all diagonals - FITTING)'
    temp = np.mean(np.mean(np.mean(res['proj_emd'], 0), 1), 1)
    titles.append(title)
    comparisons.append(temp)
    designations.append('EMD fit all')

    # emd comparison of the main diagonal projection (FITTING)
    title = 'EMD comparison of the tails (main diagonal - FITTING)'
    temp = [i[-1] for i in np.mean(np.mean(res['proj_emd'], 0), 2)]
    titles.append(title)
    comparisons.append(temp)
    designations.append('EMD fit main')

    # emd comparison of the projections
    vec = [[[] for _ in range(nb_models)] for _ in range(2 ** (dim - 1))]
    for i in res['proj_quantile']:
        for j in range(nb_models):
            for k in range(2 ** (dim - 1)):
                vec[k][j].append(i[j][k])

    control = [(i + 1) / (length + 1) for i in range(length)]
    emd_vec = [[ut.univariate_EMD_in_tails(j, control, quantile=quantile) for j in i] for i in vec]

    temp = np.mean(np.mean(emd_vec, 0), 1)
    title = 'EMD comparison of the tails (all diagonals)'
    titles.append(title)
    comparisons.append(temp)
    designations.append('EMD all')

    temp = [i[0] for i in emd_vec[-1]]
    title = 'EMD comparison of the lower tail'
    titles.append(title)
    comparisons.append(temp)
    designations.append('EMD low')

    temp = [i[1] for i in emd_vec[-1]]
    title = 'EMD comparison of the upper tail'
    titles.append(title)
    comparisons.append(temp)
    designations.append('EMD up')

    incr = 0
    for (title, temp) in zip(*[titles, comparisons]):
        if incr > 0:
            title = '%s\n%s' % (title, add_title)
            print(title)
            print(temp)
            ut.sorted_barplot(temp, xlabels=res['names'][1:-1], title=title)
            if save:
                plt.savefig('%s.png' % title.replace(' ', '_').replace('\n', '_'))
        incr += 1

    ut.table_latex(comparisons, xlabels=res['names'][1:-1], ylabels=designations,
                   title='comparison of the models: ' + add_title)

    winner_list = []
    incr = 0
    for temp in comparisons:
        if incr == 0:
            index = sorted(range(nb_models), key=temp.__getitem__, reverse=True)
        else:
            index = sorted(range(nb_models), key=temp.__getitem__)
        winner_list.append(list(map(res['names'][1:-1].__getitem__, index)))
        incr += 1

    ut.table_latex(winner_list, xlabels=list(range(nb_models)), ylabels=designations,
                   title='comparison of the models: ' + add_title)
    return titles, comparisons, designations, winner_list


def select_observations(copula, date, repeat_win_days=False, win_days=90, win_forecast=0.2, param_fixed=None):
    if isinstance(date, str):
        date = dt.parse(date)


    nb_series = len(copula.copulae)
    parameters = copula.parameters
    original_parameters = deepcopy(parameters)
    min_points = 30
    index = copula.dateM.index(date)
    forecast = [copula.dataM[i]['forecast'][0][index] for i in range(nb_series)]
    for ser in range(nb_series):
        parameters[ser]['win_forecast'] = [(-np.inf, np.inf) for _ in range(copula.copulae[ser].dim)]
        parameters[ser]['predicted_day'] = date
        if not repeat_win_days:
            parameters[ser]['date_range'] = (str(date - timedelta(days=win_days)),
                                             str(date - timedelta(hours=1)))
        else:
            parameters[ser]['date_range'] = ut.intersect_dates(param_fixed[0]['date_range'], date, win_days)

    copula.update(parameters[0], list_parameters=parameters)

    forecast_tp = [copula.dataM[i]['forecast'][0] for i in range(nb_series)]
    for i in range(nb_series):
        forecast_tp[i].append(forecast[i])
    sorted_indexes = [sorted(range(copula.lengthM + 1), key=forecast_tp[i].__getitem__) for i in
                      range(nb_series)]

    # gets the indexes of the nearest points (enough of them)
    def get_indexes():
        nb = int(np.ceil(win_forecast * copula.lengthM)) * 2
        indexes_tp = []

        while len(indexes_tp) < min_points:

            indexes_tp = around_obs(sorted_indexes[0], nb)
            for i in range(1, nb_series):
                indexes_tp = np.intersect1d(indexes_tp, around_obs(sorted_indexes[i], nb))
            nb += 2
        return indexes_tp

    indexes = get_indexes()
    vects = [[vect[j] for j in indexes] for vect in copula.vectM]
    unifs = ut.uniforms(vects, rand=False)
    copula.update(original_parameters[0])
    return vects, unifs


# gets the nb nearest neighbours of ths observation (whose value is len(x)-1) in x
def around_obs(x, nb):
    half = nb // 2
    index = x.index(len(x) - 1)
    start = max(0, min(index - half, len(x) - nb - 1))
    res_tp = x[start:index]
    res_tp.extend(x[index + 1:start + nb + 1])
    return res_tp

# -------------------------------------------- other functions --------------------------------------------------


### This functions takes a list of (same length) vectors in arguments, and return a copula based on these vectors
def fake_copulaManager(vects):
    length = len(vects[0])
    d = dt.parse('2012/01/01 00:00')
    hour = relativedelta(hours=1)
    date = []
    for i in range(length):
        date.append(str(d))
        d += hour
    series = []
    for v in vects:
        forecast = list(np.random.uniform(0, 4000, length))
        ser = {'date': date, 'vect': v, 'data': {}, 'title': {'type': 'Wind', 'location': 'NP', 'kind': 'error'},
               'var_function': None}
        ser['data']['forecast'] = forecast
        ser['data']['forecast_d'] = ut.list_derivative(forecast)
        series.append(ser)
    parameter = {'offsets': [14], 'date_range': ('2010-07-01 18:00:00', '2016-07-05 01:00:00'), 'first_hour': (0, 0)}
    return cms.CopulaManagerMS(series, parameter)


### This function looks at the evolution of the copula through the time
def copula_evolution(copula, win_days=45, day_interval=60, nb_max=5, method='color', uniforms=True, box=True):
    date_range = copula.date_range
    start = dt.parse(date_range[0])
    end = dt.parse(date_range[-1])
    par = copula.parameters
    original_date_range = par[0]['date_range']
    box_data = []

    if uniforms:
        low_all, high_all = 0, 1
    else:
        vector = copula.vectM
        low_all, high_all = min(min(vector[0]), min(vector[1])), max(max(vector[0]), max(vector[1]))
        low_all, high_all = 1.1 * low_all - 0.1 * high_all, 1.1 * high_all - 0.1 * low_all

    for i in range(nb_max):
        if (end - start).days < 20:
            break

        bound = start + relativedelta(days=win_days)
        par[0]['date_range'] = (str(start), str(bound))
        copula.update(par[0], list_parameters=par)

        if uniforms:
            vector = copula.unifM
            low, high = 0, 1
        else:
            vector = copula.vectM
            low, high = min(min(vector[0]), min(vector[1])), max(max(vector[0]), max(vector[1]))
            low, high = 1.1 * low - 0.1 * high, 1.1 * high - 0.1 * low

        if box:
            box_data.append(vector[1])
        print('###')
        print((low, high))
        print(vector)

        title = ''
        fig_num = 0
        if method == 'manichean':
            nb_pt = int(copula.lengthM / 3)
            vec1 = [vector[0][:nb_pt], vector[1][:nb_pt]]
            vec2 = [vector[0][-nb_pt:], vector[1][-nb_pt:]]
            fig = plt.figure(figsize=(8, 4))
            fig_num = fig.number
            plt.subplot(121)
            plt.scatter(vec1[0], vec1[1])
            plt.xlim(low, high)
            plt.ylim(low, high)
            plt.subplot(122)
            plt.scatter(vec2[0], vec2[1])
            plt.xlim(low, high)
            plt.ylim(low, high)
            title += 'first and last third of the copula points'
            fig.tight_layout()
            fig.subplots_adjust(top=0.85)
        elif method == 'None':
            pass
        elif method == 'same_color':
            fig_num = plt.figure().number
            plt.scatter(vector[0], vector[1])
            plt.xlim(low_all, high_all)
            plt.ylim(low_all, high_all)
            title += 'copula points in chronological order'


        elif not method == 'color':
            print('Wrong method argument, using default: \'color\'')
            method = 'color'

        if method == 'color':
            date_start = copula.dateM[0]
            date = [(date.timestamp() - date_start.timestamp()) / 86400 for date in copula.dateM]
            fig_num = plt.figure().number
            plot = plt.scatter(vector[0], vector[1], c=date, cmap='bwr')
            plt.xlim(low, high)
            plt.ylim(low, high)
            cbar = plt.colorbar(plot)
            title += 'copula points in chronological order'
            cbar.set_label('days from start', rotation=270)

        title += '\ntime range: (%s,%s)' % (start.strftime("%B"), bound.strftime("%B"))
        plt.figure(fig_num)
        plt.suptitle(title)

        start = start + relativedelta(days=day_interval)

    par[0]['date_range'] = original_date_range
    copula.update(par[0], list_parameters=par)

    if box:
        plt.figure()
        plt.subplot(211)
        plt.boxplot(box_data)
        plt.title('Evolution of the data')
        plt.ylabel('Power (MW)')
        plt.xlabel('Quarters')
        plt.subplot(212)
        plt.scatter(copula.dateM, copula.vectM[1])
        plt.ylabel('Power (MW)')
        plt.xlabel('Days')

    return None


### This function fetches the data from external files
# returns a dictionary, 'data' with keys:   - 'act' (actuals)
#                                           - 'for' (forecasts)
#                                           - 'date'
# If you have a file which dictates which files to draw from, you can specify with the sourcefile parameter
# When doing this, you also must specify the type and location of the data
def get_data(filename=None, type='Wind', location='total', sourcefile=None):
    if filename is None:
        if sourcefile is None:
            sourcefile = home_dir + 'tests/data/sources.csv'

        csv_read = csv.reader(open(sourcefile))
        index = -1
        index_dir = -1
        break1 = False
        for line in csv_read:
            if (len(line) < 1) or (line[1] == '#'):
                continue
            elif line[0] == 'title':
                for i in range(len(line)):
                    if line[i] == type + '_' + location:
                        index = i
                    if line[i] == 'home_dir':
                        index_dir = i
            if (break1):
                break
        if (index == -1) or (index_dir == -1):
            print('data specification is not valid 1')
            return -1
        else:
            csv_read = csv.reader(open(sourcefile))
            for line in csv_read:
                if (len(line) < 1) or (line[1] == '#'):
                    continue
                if line[0] == 'data':
                    filename = '%s%s' % (str(line[index_dir]), str(line[index]))
                    break
        if filename == '':
            print('data specification is not valid 2')
            return -1
        else:
            print('retrieving data from %s \n' % filename)

    data = {'date': [], 'for': [], 'act': []}
    csv_read = csv.reader(open(filename))
    for line in csv_read:
        if (len(line) < 1) or (line[1] == '#'):
            continue
        date, forecast, actual = line
        data['date'].append(dt.parse(date))
        data['for'].append(float(forecast))
        data['act'].append(float(actual))

    print('successfully retrieved %s data \n' % (type + '_' + location))
    return data


def create_list_models(nb_max=3):
    list_models = [vines.cop2d_uniform, vines.cop2d_frank, vines.cop2d_gumbel, vines.cop2d_clayton,
                   vines.cop2d_gaussian, vines.cop2d_student]

    length = len(list_models)
    res = list_models.copy()
    l = [0 for i in range(length)]
    while l is not None:
        if 1 < sum(l) <= nb_max:
            temp = [list_models[i] for i in range(length) if l[i] == 1]

            def f(x, temp=temp):
                return vines.WeightedCopula(x, temp)

            res.append(f)
        l = ut.table_increment(l)
    return res


# ------------------------------------ old comparing functions ------------------------------------------------------


### This function compares the predicting accuracy of the various models
# arguments:
#   () copula is a copulaManagerMS
#   () win_days is a window parameter specifying the number of days before and after the current time of the year
#      that should be taken into account. (if=45, it translates into 91 days in the past years and 45 in the current year)
#   () win_forecast is a window parameter specifying the width of the forecast window: if q=CDF_forecast(current_forecast),
#      we will consider dates at which the forecast was in [CDF^-1(q-'win_forecast'),CDF^-1(q+'win_forecast')]
def test_models_old(copula, win_days=45, repeat_win_days=False, win_forecast=0.2, visualize=False, start_incr=None,
                    end_incr=None, keep_vec=False, compare_dist=True):
    start_time = time.time()

    # keeping the old window parameters in param_fixed
    param_fixed = []
    for par in copula.parameters:
        dic = {'date_range': par['date_range'], 'offsets': par['offsets'].copy(), 'first_hour': par['first_hour']}
        param_fixed.append(dic)

    # initializing variables: parameters, dim (dimension of the copula), forecast, errors, dates,
    #                         forecastCDF and forecastQuantile (inversse function)

    parameters = copula.parameters.copy()
    dim = copula.copulae[0].dim

    if not 'forecast' in copula.copulae[0].dataM.keys():
        parameters[0]['forecast'] = [(-100000, 100000) for i in parameters[0]['offsets']]
        copula.update(parameters[0], list_parameters=parameters)

    forecasts = []
    for i in copula.copulae[0].dataM['forecast']:
        forecasts.append(i.copy())
    errors = []
    for i in copula.vectM:
        errors.append(i.copy())
    dates = copula.dateM.copy()
    print('forecasts (%d) %r\nerrors (%d) %r\ndates (%d) %r' % (len(forecasts), forecasts[:8],
                                                                len(errors), errors[:8], len(dates), dates[:8]))

    forecastCDF = ut.empirical_CDF_scalar(forecasts[0])
    forecastQuantile = ut.empirical_CDF_inv(forecasts[0])

    ### BIG LOOP ###


    incr = 0
    first = True
    last_time = time.time()
    res = {'len': [], 'log': [], 'rank': [], 'sum_pdf': [], 'vec': [], 'problem': [], 'past_log': [],
           'selected_model': [],
           'proj_emd': [], 'proj_quantile': []}

    # loop over each considered hour: each time,
    #       The copula manager is updated to fit the window (using 'win_days' and 'win_forecast')
    #       Models are created to fit the updated copula
    #       The log_likelihood of the observation is computed for all different models
    for obs in zip(*[forecasts[0], list(zip(*errors)), dates]):

        obs = list(obs)
        dt_object = dt.parse('200-1-1 00:00')
        obs[2] = dt_object.fromtimestamp(obs[2]).__str__()

        # selecting the observation range
        incr += 1
        if start_incr is None:
            if incr < 400:
                continue
        else:
            if incr < start_incr:
                continue
        if end_incr is not None:
            if incr >= end_incr:
                continue

        t_print = []
        for i in (time.time() - last_time, time.time() - start_time):
            i = int(i)
            t_print.append((i // 3600, i // 60 % 60, i % 60))
        last_time = time.time()
        print('\n\n   ####################   \n\niter %d: forecast %r, error %r, date %s\n\n'
              'time elapsed in the last loop: %d:%d:%d, time since start: %d:%d:%d'
              '\n\n   ####################   \n\n'
              % (
              incr, obs[0], obs[1], obs[2], t_print[0][0], t_print[0][1], t_print[0][2], t_print[1][0], t_print[1][1],
              t_print[1][2]))

        ### selecting past observations using the window ###

        print('selecting observations using the window')
        mid = forecastCDF(obs[0])
        temp = [(-np.inf, np.inf) for obs in range(dim)]
        temp[0] = (float(forecastQuantile(max(0.0001, mid - win_forecast))),
                   float(forecastQuantile(min(0.9999, mid + win_forecast))))
        parameters[0]['forecast'] = temp
        del temp
        parameters[0]['predicted_day'] = obs[2]
        if not repeat_win_days:
            parameters[0]['date_range'] = (
            str(dt.parse(obs[2]) - relativedelta(days=win_days)), str(dt.parse(obs[2]) - relativedelta(hours=1)))
        else:
            parameters[0]['date_range'] = ut.intersect_dates(param_fixed[0]['date_range'], obs[2], win_days)

        try:
            copula.update(parameters[0], list_parameters=parameters)
        except:
            parameters[0]['forecast'] = [(-np.inf, np.inf) for obs in range(dim)]
            copula.update(parameters[0], list_parameters=parameters)

        if copula.lengthM < 30:
            parameters[0]['forecast'] = [(-np.inf, np.inf) for obs in range(dim)]
            res['problem'].append('not applying win_forecast at iteration %d' % incr)
            copula.update(parameters[0], list_parameters=parameters)

        ### fitting models to the distribution ###

        print('fitting models to the distribution')
        length = copula.lengthM
        dim = copula.dim

        if length < 30:
            res['problem'].append('length <30 at iteration %d' % incr)
            continue

        # try:
        # creating the density of a fitted gaussian
        def create_gaussian_density():
            covariance = np.cov(copula.vectM)
            means = np.mean(copula.vectM, axis=1)

            cov_inv = np.linalg.inv(np.matrix(covariance))
            fact_gau = np.sqrt(np.linalg.det(cov_inv) / (2 * math.pi) ** dim)

            def den_gau(vec):
                if type(vec[0]) in {int, np.float, float, np.int}:
                    vec = [[j] for j in vec]
                res = []
                for j in zip(*vec):
                    j = np.matrix(j)
                    res.append(fact_gau * math.exp(-(j - means) * cov_inv * np.transpose(j - means) / 2))
                return res

            return den_gau

        densities = [create_gaussian_density()]

        # creating a list of copula models
        list_models = [mod.cop_gaussian(copula.unifM), mod.cop_student(copula.unifM),
                       # mod.cop_customized(mod.cop_student,copula.unifM,redistribute=True),
                       vines.D_vine(copula.unifM), mod.cop_uniform(copula.unifM)]  # ,vines.C_vine(copula.unifM)]

        # list_models=[vines.D_vine(copula.unifM,rearrange=False),vines.D_vine(copula.unifM,rearrange=True)]

        # list_models=[vines.cop2d_gaussian(copula.unifM),vines.cop2d_student(copula.unifM)]

        # computing the densities of the models, their log-likelihood, and selecting the 'best candidate'

        names = ['gaussian']
        cop_densities = []
        best_model_past = list_models[0]
        best_log_past = 0
        log_past = []
        for j in list_models:
            if first:
                names.append(j.name)
            cop_densities.append(j.pdf)
            lld = sum([math.log(k) for k in j.pdf(copula.unifM)])
            log_past.append(lld)
            if lld > best_log_past:
                best_log_past = lld
                best_model_past = j

        names.append('selected model')
        cop_densities.append(best_model_past.pdf)
        log_past.append(best_log_past)

        res['past_log'].append(log_past)
        res['selected_model'].append(best_model_past.name)
        if first:
            res['names'] = names
            first = False

        # computing the rank of 'obs' among the window points
        CDFs = ut.marginals_cdf(copula.vectM)
        rank = [float(CDFs[i](obs[1][i])) for i in range(dim)]

        # computing the tail metrics:
        C_to_D = ut.copula_to_distribution(copula.vectM)
        simulations = [C_to_D(cop.simulate(1000)) for cop in list_models]
        tail_metrics = ut.compare_tails(simulations, copula.vectM, obs[1], quantile=0.1)
        res['proj_emd'].append(tail_metrics[0])
        res['proj_quantile'].append(tail_metrics[1])

        try:
            ### computing the log likelihood ###
            if compare_dist:
                res_log = [math.log(den(obs[1])[0]) for den in densities]
                res_log.extend(
                    [den(obs[1])[0] for den in ut.copula_to_densities(copula.vectM, cop_densities, log_return=True)])
            else:
                res_log = [den(rank)[0] for den in
                           ut.distribution_to_copula_densities(copula.vectM, densities, log_return=True)]
                res_log.extend([math.log(den([[r] for r in rank])[0]) for den in cop_densities])

            res['log'].append(res_log)
        except:
            res['log'].append(None)
            res['problem'].append('incr: %d, problem in the log: %r' % (incr, sys.exc_info()[0]))
        res['len'].append(length)
        res['rank'].append(rank)

        if keep_vec:
            res['vec'].append([copula.unifM, list_models[0].val, list_models[2].val])
            # except:
            #     res['problem'].append('incr %d general problem: %r'%(incr,sys.exc_info()[0]))

    copula.update(param_fixed[0], list_parameters=param_fixed)
    return res


### computes the L2 distance between two copula
def compare_copula_l2(to_compare, visualize=True):
    res = []
    for i in to_compare['unifs']:
        res.append(ut.compute_distance_l(to_compare['unifs'][0], i))

    if visualize:
        ind = sorted(range(len(res)), key=res.__getitem__, reverse=False)
        nm = list(map(to_compare['names'].__getitem__, ind))
        scores = list(map(res.__getitem__, ind))
        wid = 0.5
        abs = np.array(range(len(scores)))
        plt.figure()
        plt.title('L2 comparison of the underlying copulae')
        plt.bar(abs, scores, width=0.5)
        plt.xticks(abs + wid / 2., nm, rotation=40)

    return res


### compares models fitted to a copula, using emd distance.
# returns:
#   () a matrix featuring the distance between the various models
# arguments:
#   () copula: a copulaMS manager
#   () sample_size: the size of the sample to be used to compute emd distance.
#       The final distance will be an average of the EMD computed on these samples
def compare_distributions_emd(copula, to_compare=None, sample_size=50, visualize=True):
    length = copula.lengthM
    dim = copula.dim
    vects = []

    if (to_compare is None):
        vects.append(copula.vectM)
        # creating a fitted gaussian
        covariance = np.cov(vects[0])
        means = np.mean(vects[0], axis=1)
        gaussian = np.transpose(np.random.multivariate_normal(means, covariance, length))
        temp = []
        for i in range(dim):
            temp.append(gaussian[i].tolist())
        vects.append(temp)

        # computing the inverse of marginals' CDF
        F_inv = ut.marginals_cdf_inv(vects[0])

        # creating copulae and the corresponding distribution
        x = np.arange(0, 1, 0.02)

        unifs = copula.unifM
        subsets = []
        for i in range(dim):
            for j in range(i + 1, dim):
                subsets.append(sub.create_diagonal([i, j]))
        subsets.append(sub.main_axis)

        for cop in [mod.cop_gaussian(unifs), mod.cop_student(unifs)]:  # ,mod.cop_student_custom(unifs),
            # mod.cop_customized(mod.cop_gaussian,unifs,subsets),mod.cop_customized(mod.cop_student,unifs,subsets)]:
            f = cop.f
            unifs = f(length)
            vect_temp = []
            if visualize:
                plt.figure()
                plt.title('empirical_copula')
                plt.plot(unifs[0], unifs[1], '.')
            for i in range(dim):
                vect_temp.append(list(map(F_inv[i], unifs[i])))

            vects.append(vect_temp)
    else:
        for i in to_compare['vects']:
            vects.append(i.copy())

    # reshuffling the vectors
    points = []
    for i in vects:
        temp = list(zip(*i))
        np.random.shuffle(temp)
        points.append(temp)

    # Computing the mean of emd distance between samples of distribution
    nb_dist = len(vects)

    for i in range(nb_dist):
        l_temp = len(points[i])
        if length > l_temp:
            print("distribution %d is length %d instead of %d" % (i, l_temp, length))
            length = l_temp

    nb_sample = int(length / sample_size)
    res = np.identity(nb_dist)
    for i in range(nb_dist):
        for j in range(i, nb_dist):
            p1 = points[i]
            p2 = points[j]
            print('len1 %d len2 %d' % (len(p1), len(p2)))
            print(p1[:10])
            print(p2[:10])
            res[i, j] = 0;
            res[j, i] = 0
            for k in range(nb_sample):
                a = k * sample_size
                b = ((k + 1) % nb_sample) * sample_size
                print(nb_sample)
                print('dist %d %d --- k: %.2f a: %.2f b: %.2f' % (i, j, k, a, b))
                res[i, j] += ut.compute_emd(p1[a:(a + sample_size)], p2[b:(b + sample_size)])
            res[j, i] = res[i, j]

    # if visualize:
    #     for v in vects:
    #         plt.figure()
    #         plt.plot(v[0],v[1],'.')
    #         print(res)

    for i in range(nb_dist - 1):
        for j in range(i + 1, nb_dist):
            res[i, j] = round(res[i, j] / math.sqrt(res[i, i] * res[j, j]), 4)
            res[j, i] = res[i, j]
    for i in range(nb_dist):
        res[i, i] = 1
        res.tolist()

    print(res)
    if visualize:
        scores = res[0].copy()
        ind = sorted(range(len(scores)), key=scores.__getitem__, reverse=False)
        if to_compare is None:
            nm = ['' for i in range(len(scores))]
        else:
            nm = list(map(to_compare['names'].__getitem__, ind))
        scores = list(map(scores.__getitem__, ind))
        wid = 0.5
        ab = np.array(range(len(scores)))
        plt.figure()
        plt.title('emd comparison of the distributions')
        print(ab, scores)
        plt.bar(ab, scores, width=0.5)
        plt.xticks(ab + wid / 2., nm, rotation=40)

    return (res)


### compares models fitted to a c
# opula,using log-likelihood:
def compare_distributions_log(to_compare, visualize=True):
    vects = to_compare['vects']
    length = len(vects)
    density = to_compare['dis_density']

    res = []
    for i in range(length):
        print(i)
        res.append(ut.emp_log_likelihood(vects[0], vects[i], density1=density[0], density2=density[i]))

    if visualize:
        for result in list(zip(*res)):
            ind = sorted(range(len(result)), key=result.__getitem__, reverse=True)
            nm = list(map(to_compare['names'].__getitem__, ind))
            scores = list(map(result.__getitem__, ind))
            scores = [scores[0] - i for i in scores]
            wid = 0.5
            abs = np.array(range(len(scores)))
            plt.figure()
            plt.title('log-likelihood comparison \n (-L(model)+L(original))')
            plt.bar(abs, scores, width=0.5)
            plt.xticks(abs + wid / 2., nm, rotation=40)
    return res


### creates a 'to_compare' dictionary and compares different models using the 3 previous distances
def compare_distributions(copula, list_models=None, visualize=True):
    unifs = copula.unifM
    vects = [copula.vectM]
    names = ['real']
    cop_density = [None]
    length = copula.lengthM
    dim = copula.dim

    # creating a fitted gaussian
    covariance = np.cov(vects[0])
    means = np.mean(vects[0], axis=1)
    gaussian = np.transpose(np.random.multivariate_normal(means, covariance, length))
    temp = []
    for i in range(dim):
        temp.append(gaussian[i].tolist())
    vects.append(temp)
    names.append('gaussian (w/o copula)')
    cov_inv = np.matrix(covariance) ** (-1)
    fact_gau = np.sqrt(np.linalg.det(cov_inv) / (2 * math.pi) ** dim)

    def den_gau(vec):
        if type(vec[0]) in {int, np.float, float, np.int}:
            vec = [[i] for i in vec]
        res = []
        for i in zip(*vec):
            i = np.matrix(i)
            res.append(fact_gau * math.exp(-i * cov_inv * np.transpose(i)) / 2)
        return res

    cop_density.append(den_gau)
    del (temp)

    # print(vects)
    # raise(RuntimeError())
    if list_models is None:
        # different models
        list_models = [mod.cop_emp_redistributed(unifs), mod.cop_gaussian(unifs), mod.cop_student(unifs),
                       mod.cop_customized(mod.cop_gaussian, unifs), mod.cop_customized(mod.cop_student, unifs)]
        # renormalizing the customized copulae
        for i in list_models[-2:]:
            list_models.append(mod.cop_emp_redistributed(i.val, model=i))
            # return list_models

    unifs = []
    for i in list_models:
        unifs.append(i.val)
        names.append(i.name)
        cop_density.append(i.pdf)

    cop_density[2] = None
    dis_density = ut.copula_to_densities(vects[0], cop_density)

    to_dist = ut.copula_to_distribution(vects[0], visualize=False)
    for i in unifs:
        vects.append(to_dist(i))

    unifs.reverse()
    unifs.append(mod.uniforms(vects[1], rand=False))
    unifs.append(mod.uniforms(vects[0], rand=False))
    unifs.reverse()

    to_compare = {'vects': vects, 'unifs': unifs, 'names': names, 'dis_density': dis_density}

    # comparing the distributions using various measures
    res = []
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

    res.append(compare_distributions_log(to_compare, visualize=True))
    res.append(compare_copula_l2(to_compare))
    res.append((compare_distributions_emd(copula, to_compare=to_compare).tolist())[0])
    # res.append(compare_distributions_log(to_compare))

    return res, to_compare
