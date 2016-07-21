import csv
import os
import matplotlib.pyplot as plt
import numpy as np
import pickle
from scipy import stats

def compile_results(results_dict, dir_name=None):
    if dir_name is None:  # Only used if the parameters key exists
        params = results_dict['parameters']
        dir_name = "{}_{}_{}_win_days_{}_win_forecast_{}".format(params['type'], params['location'], params['kind'],
                                                                 params['win_days'], params['win_forecast'])
    if not(os.path.isdir(dir_name)):
        os.mkdir(dir_name)

    #PICKLE PICKLE PICKLE PICKLE PICKLE
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

    fig2, ax2 = plt.subplots()
    m = min(means)
    adjusted_pairs = [(x, y - m) for (x, y) in pairs]
    ax2.bar(indices, [pair[1] for pair in adjusted_pairs])
    ax2.set_ylabel('Average Log Likelihood Minus Least')
    ax2.set_title('Average Log Likelihood of Various Copulas Over a Series of Observations')
    ax2.set_xticks(indices + 0.4)
    ax2.set_xticklabels([pair[0] for pair in pairs], rotation=45)
    fig2.tight_layout()
    plt.savefig(dir_name + os.sep + 'Adjusted_log_likelihoods.png')

    write_csv(results_dict, dir_name + os.sep + 'data.csv')


def write_csv(results_dict, filename):
    """Writes the results of the copula analysis experiment to a csv file"""
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
def merge_results(*results):
    new_results = {'log':[], 'proj_emd':[], 'proj_quantile':[], 'names':[]}
    for result in results:
        for i, log in enumerate(result['log']):
            new_results[]
"""

def get_emds(res):
    emds = res['proj_emd']
    new_emds = [[] for _ in emds[0]]
    for i, copula in enumerate(emds[0]):
        new_emds[i].append([emds[j][i][0][0] for j in range(len(emds))])
        new_emds[i].append([emds[j][i][0][1] for j in range(len(emds))])
        new_emds[i].append([emds[j][i][1][0] for j in range(len(emds))])
        new_emds[i].append([emds[j][i][1][1] for j in range(len(emds))])
    return new_emds

def z_interval(data, variance, confidence=.95):
    xbar = np.mean(data)
    moe = stats.norm.ppf((1+confidence)/2) * np.sqrt(variance/len(data))
    return xbar - moe, xbar + moe


def t_interval(data, confidence=.95):
    xbar = np.mean(data)
    var = np.var(data)
    N = len(data)
    moe = stats.t.ppf((1+confidence)/2, N - 1) * np.sqrt(var/N)
    return xbar - moe, xbar + moe

def fix_name(name):
    if name.startswith('weighted_ST'):
        return name.replace('ST', 'UN', 1)

def mean_difference(data):
    return np.mean([x - y for (x,y) in zip(data, data[1:])])

def ranks(winner_lists):
    winner_ranks = {}
    for name in winner_lists[0]:
        winner_ranks[name] = []

    for winner_list in winner_lists:
        for i, winner in enumerate(winner_list):
            winner_ranks[winner].append(i)

    return winner_ranks

