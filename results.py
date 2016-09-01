import csv
import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
import utilities as ut
from scipy import stats
import dateutil.parser as dt

from copula_analysis import select_observations
from vines import cop2d_gumbel, cop2d_frank, cop2d_uniform, cop2d_gaussian, cop2d_clayton, cop2d_student

simple_copulas = ['frank', 'gumbel', 'uniform', 'gaussian', 'student', 'clayton']
simple_copulas.sort()

copula_constructors = [cop2d_clayton, cop2d_frank, cop2d_gaussian, cop2d_gumbel, cop2d_student]

class Event:
    """
    This is an object which contains all the pertinent details stored after each cycle of the test.
    It stores the log-likelihood, quantile,
    """
    def __init__(self, names, date, logs, quantile, emd, rank):
        self.names = names
        self.date = date
        self.logs = logs
        self.quantile = quantile
        self.emd = emd
        self.rank = rank

    def log_by_name(self, name):
        name_index = self.names.index(name)
        return self.logs[name_index]

    def __str__(self):
        return "Event: {}".format(self.date)

    __repr__ = __str__

class Results:
    def __init__(self, results_dict, manager=None):
        self.manager = manager
        self.nb_models = nb_models = len(results_dict['names']) - 2
        self.length = length = len(results_dict['log'])
        self.dim = dim = int(np.log(np.shape(results_dict['proj_quantile'])[2]) / np.log(2)) + 1
        self.proj_quantile = results_dict['proj_quantile']
        self.proj_emd = results_dict['proj_emd']
        self.ranks = results_dict['rank']
        self.log = results_dict['log']
        self.names = results_dict['names']
        self.dates = results_dict['dates']
        self.parameters = results_dict['parameters']
        self.events = [Event(self.names, date, logs, quantile, emd, rank) for date, logs, quantile, emd, rank in zip(self.dates,
                                                                                                         self.log,
                                                                                                         self.proj_quantile,
                                                                                                         self.proj_emd,
                                                                                                         self.ranks)]
        for key in results_dict:
            setattr(self, key, results_dict[key])

        vec = [[[] for _ in range(nb_models)] for _ in range(2 ** (dim - 1))]
        for i in results_dict['proj_quantile']:
            for j in range(nb_models):
                for k in range(2 ** (dim - 1)):
                    vec[k][j].append(i[j][k])

        control = [(i + 1) / (length + 1) for i in range(length)]
        self.emd_vec = [[ut.univariate_EMD_in_tails(j, control, quantile=0.1) for j in i] for i in vec]

        self.mean_log_likelihood = np.mean(self.log, 0)
        self.emd_fit_all = np.mean(np.mean(np.mean(self.proj_emd, 0), 1), 1)
        self.emd_fit_main = [i[-1] for i in np.mean(np.mean(self.proj_emd, 0), 2)]
        self.emd_all = np.mean(np.mean(self.emd_vec, 0), 1)
        self.emd_low = [low for (low, _) in self.emd_vec[-1]]
        self.emd_up = [up for (_, up) in self.emd_vec[-1]]

        self.attrs = ['log-likelihood', 'EMD fit all', 'EMD fit main', 'EMD all', 'EMD low', 'EMD up']

    def ranks_by_name(self, name):
        """
        Returns a dictionary of all the ranks of the experiment by name of distribution, usually use
        for the main simple distributions as the other distributions are more complicated.
        Returns dict {'log-likelihood':, 'EMD Fit all': , 'EMD Fit Main': , 'EMD UP': , 'EMD low': 'EMD all'
        """
        name_index = self.names[1:-1].index(name)
        return {'log-likelihood': self.nb_models - stats.rankdata(self.mean_log_likelihood[1:-1])[name_index],
                # lower has higher rank
                'EMD fit all': stats.rankdata(self.emd_fit_all)[name_index],  # higher has higher rank
                'EMD fit main': stats.rankdata(self.emd_fit_main)[name_index],
                'EMD all': stats.rankdata(self.emd_all)[name_index],
                'EMD low': stats.rankdata(self.emd_low)[name_index],
                'EMD up': stats.rankdata(self.emd_up)[name_index]}

    def stats_by_name(self, name):
        """
        Returns a dictionary of all the stats of the experiment by name of distribution, usually use
        for the main simple distributions as the other distribution names are more complicated.
        Returns dict {'log-likelihood':, 'EMD Fit all': , 'EMD Fit Main': , 'EMD UP': , 'EMD low': 'EMD all'
        """
        name_index = self.names[1:-1].index(name)
        return {'log-likelihood': self.mean_log_likelihood[1:-1][name_index],
                'EMD fit all': self.emd_fit_all[name_index],
                'EMD fit main': self.emd_fit_main[name_index],
                'EMD all': self.emd_all[name_index],
                'EMD low': self.emd_low[name_index],
                'EMD up': self.emd_up[name_index]}

    def get_event(self, date):
        """
        This method accepts a date as either a string or a datetime object and returns the corresponding
        Event with the details for that day
        """
        if isinstance(date, str):
            date = dt.parse(date)

        date_index = self.dates.index(str(date))
        return self.events[date_index]

    def get_range(self, start_date, end_date):
        if isinstance(start_date, str):
            start_date = dt.parse(start_date)
        if isinstance(end_date, str):
            end_date = dt.parse(end_date)

        start_index = self.dates.index(str(start_date))
        end_index = self.dates.index(str(end_date))
        return self.events[start_index:end_index+1]


    def plot_day(self, date, simple=False):
        """
        Creates a scatterplot for the points in the window used to construct the copulas for any given date.
        If simple is True, it also gives 3D plots of the simple copulas
        """
        if self.manager is None:
            raise RuntimeError("Results object must have associated Copula Manager to retrieve data")

        if isinstance(date, str):
            date = dt.parse(date)

        _, unifs = select_observations(self.manager, date, win_days=self.parameters['win_days'],
                                       win_forecast=self.parameters['win_forecast'])

        ut.hist3D(*unifs, bins=10)

        date_index = self.dates.index(str(date))
        point = self.ranks[date_index]

        plt.figure()
        plt.scatter(*unifs)
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.title("{} win_days: {} win_forecast: {}".format(date, self.parameters['win_days'],
                                                            self.parameters['win_forecast']))
        plt.plot(*point, color='r', marker='x')

        if simple:
            for name, copula in zip(simple_copulas, copula_constructors):
                model = copula(unifs)
                ut.curve3d(model.pdf, points=False, title=name)

    def log_plot(self, start_date=None, end_date=None):
        if start_date is None:
            start_date = self.dates[0]
        if end_date is None:
            end_date = self.dates[-1]

        if isinstance(start_date, str):
            start_date = dt.parse(start_date)
        if isinstance(end_date, str):
            end_date = dt.parse(end_date)

        start_index = self.dates.index(str(start_date))
        end_index = self.dates.index(str(end_date))
        plt.figure()
        lines = []
        for i, copula in enumerate(self.names[1:7]):
            logs = list(zip(*self.log))[i+1][start_index:end_index+1]
            lines.extend(plt.plot(list(map(dt.parse, self.dates))[start_index:end_index+1], logs, label=copula))
        plt.ylabel("Log-Likelihood")
        plt.xlabel("Date")
        plt.legend(handles=lines, loc=0)
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.title("{} - {} Log Likelihood".format(start_date, end_date))

    def general_table(self, simple=True, title=None):
        rows = []
        if simple:
            names = simple_copulas
        else:
            names = self.names

        name_stats = {name: self.stats_by_name(name) for name in names}

        for attr in self.attrs:
            row = []
            for name in names:
                row.append("{:7.4f}".format(name_stats[name][attr]))
            rows.append(row)

        return ut.table_latex(rows, xlabels=names, ylabels=self.attrs, title=title)

    def winner_table(self, attr, simple=False, upto=10):
        """
        Returns a Latex Table as a string ranking each of the copulas by one of the statistics.
        Possible values for attr are:
        - 'EMD fit all'
        - 'EMD fit main'
        - 'EMD all'
        - 'EMD low'
        - 'EMD up'
        - 'log-likelihood'
        """
        rows = []
        reverse = False
        if attr == 'log-likelihood':
            stat = self.mean_log_likelihood[1:-1]
            reverse = True
        elif attr == 'EMD fit all':
            stat = self.emd_fit_all
        elif attr == 'EMD fit main':
            stat = self.emd_fit_main
        elif attr == 'EMD all':
            stat = self.emd_all
        elif attr == 'EMD low':
            stat = self.emd_low
        else:
            stat = self.emd_up

        sorted_by_stat = sorted(zip(stat, self.names[1:-1]), reverse=reverse)

        for score, name in sorted_by_stat[:upto]:
            rows.append([name, "{:7.4f}".format(score)])

        ylabels = list(range(upto))
        if simple:
            for score, name in sorted_by_stat[upto + 1:]:
                if name in simple_copulas:
                    rows.append([name, "{:7.4f}".format(score)])
                    ylabels.append(int(self.ranks_by_name(name)[attr]))
        return ut.table_latex(rows, ylabels=ylabels, xlabels=["Model", attr])

    def write_csv(self, filename):
        """Writes the main statistics to a specified csv file"""
        with open(filename, 'w') as f:
            writer = csv.writer(f, lineterminator='\n')
            rows = []
            for param in sorted(self.parameters):
                if param == 'date_range':
                    rows.append(('date_range',)+self.date_range)
                    continue
                rows.append([param, self.parameters[param]])
            rows.append(['Copula', 'Mean Log-likelihood', 'EMD fit all', 'EMD fit main', 'EMD all',
                             'EMD low', 'EMD up'])
            for name, log, fit_all, fit_main, all, low, up, in zip(self.names[1:-1], self.mean_log_likelihood[1:-1],
                                                                   self.emd_fit_all, self.emd_fit_main, self.emd_all,
                                                                   self.emd_low, self.emd_up):
                rows.append(map(str, [name, log, fit_all, fit_main, all, low, up]))

            writer.writerows(rows)

    @property
    def date_range(self):
        return self.dates[0], self.dates[-1]

    def __iter__(self):
        return iter(self.events)

def compile_results(results_dict, dir_name=None, manager=None):
    if dir_name is None:  # Only used if the parameters key exists
        params = results_dict['parameters']
        dir_name = "{}_{}_{}_win_days_{}_win_forecast_{}".format(params['type'], params['location'], params['kind'],
                                                                 params['win_days'], params['win_forecast'])
    if not(os.path.isdir(dir_name)):
        os.mkdir(dir_name)

    ## PICKLE PICKLE PICKLE PICKLE PICKLE
    pickle.dump(results_dict, open(dir_name + os.sep + 'pickle', 'wb'))

    means = [np.mean(copula) for copula in zip(*results_dict['log'])]
    pairs = list(sorted(zip(results_dict['names'], means), key=lambda x: x[1]))

    # bar plot of means
    indices = np.arange(len(means))
    fig, ax = plt.subplots()
    ax.bar(indices, [pair[1] for pair in pairs])
    ax.set_ylabel('Average Log Likelihood')
    ax.set_title('Average Log Likelihood of Various Copulas Over a Series of Observations')
    ax.set_xticks(indices + 0.4)
    ax.set_xticklabels([pair[0] for pair in pairs], rotation=45)
    fig.tight_layout()
    plt.savefig(dir_name + os.sep + 'Average_log_likelihoods.png')

    res = Results(results_dict, manager=manager)

    res.write_csv(dir_name + os.sep + 'log_likelihoods_and_emds.csv')

"""
def write_csv(results_dict, filename):
    Writes the results of the copula analysis experiment to a csv file
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        if results_dict.get('parameters', None):
            for param in sorted(results_dict['parameters']):
                writer.writerow([param, results_dict['parameters'][param]])
        means = [np.mean(copula) for copula in zip(*results_dict['log'])]
        row0 = ['Average Log Likelihoods'] + means
        writer.writerow(row0)
        row1 = [''] + results_dict['names'] + ['selected model'] + ['rank1', 'rank2']
        row1 = row1 + sum([[name + ' emd'] * 4 for name in results_dict['names'][1:8]], [])
        row1 = row1 + sum([[name + ' quantile'] * 2 for name in results_dict['names'][1:8]], [])
        writer.writerow(row1)
        rows = []
        for i, logs in enumerate(results_dict['log']):
            if logs is None:
                continue
            row = [i] + logs + [results_dict['selected_model'][i],
                                results_dict['rank'][i][0], results_dict['rank'][i][1]]
            row = row + sum([[results_dict['proj_emd'][i][j][0][0]] +
                         [results_dict['proj_emd'][i][j][0][1]] +
                         [results_dict['proj_emd'][i][j][1][0]] +
                         [results_dict['proj_emd'][i][j][1][1]]for j in range(len(results_dict['names'])-2)], [])
            row = row + sum([[results_dict['proj_quantile'][i][j][0]] +
                             [results_dict['proj_quantile'][i][j][1]] for j in range(len(results_dict['names'])-2)], [])
            rows.append(row)
        writer.writerows(rows)
"""

def log_table(list_of_results):
    table = 'Trial #: '
    for name in list_of_results[0]['names'][1:-1]:
        table += '{:>7} '.format(name[:7])
    table += '\n'
    for i, trial in enumerate(list_of_results):
        table += '{:<9}'.format(i)
        for j, _ in enumerate(trial['names'][1:-1]):
            table += '{:7.4f} '.format(np.mean(list(zip(*trial['log']))[j+1]))
        table += '\n'

    return table


def averages(list_of_results):
    average_results = []
    for trial in list_of_results:
        average_results.append(np.mean(trial['log'], 0).tolist())
    return list(zip(*average_results))


def rank_table(names, values):
    table = 'Trial #: '
    for name in names:
        table += '{:>7} '.format(name[:7])
    table += '\n'
    for i, trial in enumerate(values[0]):
        table += '{:<9}'.format(i)
        row = [-col[i] for col in values]  # negative to reverse order
        ranks = stats.rankdata(row)
        for rank in ranks:
            table += '{:>7d} '.format(int(rank))
        table += '\n'

    return table


def join_results(list_of_results):
    new_results = list_of_results[0].copy()
    for result in list_of_results[1:]:
        for name in ['past_log', 'rank', 'log', 'len', 'proj_emd', 'selected_model', 'proj_quantile']:
            new_results[name].extend(result[name])

    return new_results

def mean_difference(data):
    """Returns the mean of the difference between consecutive values in a list"""
    return np.mean([x - y for (x,y) in zip(data, data[1:])])


def ranks(winner_lists):
    winner_ranks = {}
    for name in winner_lists[0]:
        winner_ranks[name] = []

    for winner_list in winner_lists:
        for i, winner in enumerate(winner_list):
            winner_ranks[winner].append(i)

    return winner_ranks

