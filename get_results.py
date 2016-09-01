import sys

import copula_analysis as ca
import numpy as np
import results

home = 'C:\\users\\sabrina\\desktop\\horse_races\\Solar_Wind_14'



def format_results(list_of_results, filename=None, ranking_method=0):
    """
    This function accepts a list of multiple results experiments from test_models.
    It outputs a list of useful statistics into a specified file.
    These include the mean difference over consecutive values for each statistic, the standard deviations
    for each copula model of their rankings over each statistic and a final combined ranking of copulas by a method
    which can be specified.
    0 indicates ranking by max(rank(log-likelihood), rank(emd all))
    1 indicates ranking by mean(all statistics ranks)
    """

    if filename is None:
        file = sys.stdout
    else:
        file = open(filename, 'w')

    exp_results = [ca.visualize_result(res, save=False, quantile=0.1) for res in list_of_results]
    # viz results returns [titles, comparisons, designations, winners] all lists

    names = exp_results[0][3][0]

    mean_diffs = {}
    ranks = {}
    rank_stds = {}
    ranks_by_name = {name: {designation: [] for designation in exp_results[0][2]}
                       for name in names}
    results_by_name = {name: {designation: [] for designation in exp_results[0][2]}
                     for name in names}

    for i, designation in enumerate(exp_results[0][2]):
        for experiment in exp_results:
            for j, name in enumerate(list_of_results[0]['names'][1:-1]):
                results_by_name[name][designation].append(experiment[1][i][j])

    for i, designation in enumerate(exp_results[0][2]):
        for experiment in exp_results:
            for j, name in enumerate(experiment[3][i]):
                ranks_by_name[name][designation].append(j)  # Compiles the ranks for each experiment by name of distro
                ranks_by_name[name]['averages'] = []
                ranks_by_name[name]['worsts'] = []

    # Calculates mean difference for each measured statistic
    # Also compiles the ranks into a dictionary
    for test in exp_results:
        for i, designation in enumerate(test[2]):
            mean_diff = results.mean_difference(sorted(test[1][i]))
            if mean_diffs.get(designation):
                mean_diffs[designation].append(mean_diff)
            else:
                mean_diffs[designation] = [mean_diff]

            if ranks.get(designation):
                ranks[designation].append(test[3][i])
            else:
                ranks[designation] = [test[3][i]]

    rankings = {}

    for statistic in ranks:
        rankings[statistic] = results.ranks(ranks[statistic])
        rank_stds[statistic] = {name: np.std(rankings[statistic][name]) for name in rankings[statistic]}

    # Output
    # for statistic in sorted(mean_diffs):
    #     print('-----------------{}-------------------'.format(statistic), file=file)
    #     print('Mean Differences: ' + '  '.join(map(str, mean_diffs[statistic])) +
    #           ' Standard Deviation: ' + str(np.std(mean_diffs[statistic])), file=file)
    #     print('Standard Deviation of Ranks:', file=file)
    #     for name in sorted(rank_stds[statistic], key=lambda x:np.mean(rankings[statistic][x])):
    #         print('{:40} {:60}: {}'.format(name, str(rankings[statistic][name]), rank_stds[statistic][name]), file=file)
    #     print("Mean of STDS: {}".format(np.mean(list(rank_stds[statistic].values()))), file=file)
    #     print('Gaussian Ranks: {}\n'.format(ranks_by_name['gaussian'][statistic]), file=file)
    #
    # print('Combined Result: ', file=file)
    #
    # for i, experiment in enumerate(exp_results):
    #     print("\nExperiment {}:\n".format(i), file=file)
    #     if ranking_method == 0:
    #         print('-------------Ranked by max of EMD all and mean-log-likelihood-------------------', file=file)
    #         for j, name in enumerate(sorted(simple,
    #                                         key=lambda x: max(ranks_by_name[x]['EMD all'][i],
    #                                                           ranks_by_name[x]['log-likelihood'][i]))):
    #             if j < 20:
    #                 print('{:2}: {:40} {:>3} {:30} {:30} {} {:>3} {:30}'.format(j, name, 'EMD (all)',
    #                                                                    ranks_by_name[name]['EMD all'][i],
    #                                                                    results_by_name[name]['EMD all'][i],
    #                                                                    'mean log-likelihood',
    #                                                                    ranks_by_name[name]['log-likelihood'][i],
    #                                                                    results_by_name[name]['log-likelihood'][i]),
    #                       file=file, end=' ')
    #                 print("Difference:", abs(ranks_by_name[name]['EMD all'][i] - ranks_by_name[name]['log-likelihood'][i]), file=file)
    #             ranks_by_name[name]['worsts'].append(j)
    #     elif ranking_method == 1:
    #         print('--------------Ranked by average of EMD all and mean-log-likelihood----------------', file=file)
    #         for j, name in enumerate(sorted(ranks_by_name,
    #                                         key=lambda x: np.mean([ranks_by_name[x]['EMD all'][i],
    #                                                                ranks_by_name[x]['log-likelihood'][i]]))):
    #             if j < 20:
    #                 print('{:2}: {:40} {:30} {:>3} {:30} {:>3}'.format(j, name, 'EMD (all)',
    #                                                                    ranks_by_name[name]['EMD all'][i],
    #                                                                    'mean log-likelihood',
    #                                                                    ranks_by_name[name]['log-likelihood'][i]),
    #                       file=file, end='')
    #                 print("Difference:", abs(ranks_by_name[name]['EMD all'][i] - ranks_by_name[name]['log-likelihood'][i]),
    #                       file=file)
    #             ranks_by_name[name]['averages'].append(j)
    #     print("Average difference between log rank and EMD all rank is {}".format(np.mean([abs(ranks_by_name[name]['EMD all'][i] -
    #                                                                                        ranks_by_name[name]['log-likelihood'][i])
    #                                                                                        for name in ranks_by_name])), file=file)
    #     print("Average simple copula distance is {}".format(np.mean([abs(ranks_by_name[name]['EMD all'][i] -
    #                                                                      ranks_by_name[name]['log-likelihood'][i])
    #                                                                      for name in simple])), file=file)
    #
    # table_labels = ['Mean Log-Likelihood'] + ['EMD All ' + str(i+1) for i in range(len(exp_results))]
    # rows = [[round(results_by_name[name]['log-likelihood'][0], 3) for name in simple]]
    # for i, exp in enumerate(exp_results):
    #     rows.append([int(round(results_by_name[name]['EMD all'][i] * 10**5)) for name in simple])
    # print(rows)
    # ut.table_latex(rows, xlabels=simple, ylabels=table_labels, title="Solar Experiment Results with Quantile 0.5")
    #
    # """
    # table_labels = ['Mean Log-Likelihood'] + ['EMD All ' + str(i + 1) for i in range(len(exp_results))]
    # rows = [[results_by_name[name]['log-likelihood'][0] for name in simple[3:]]]
    # for i, exp in enumerate(exp_results):
    #     rows.append([results_by_name[name]['EMD all'][i] for name in simple[3:]])
    # print(rows)
    # ut.table_latex(rows, xlabels=simple[3:], ylabels=table_labels, title="Solar Experiment Results with Quantile 0.1")
    # """
    # """
    # print('-----------------FINAL RANKINGS----------------', file=file)
    # if ranking_method == 0:
    #     print('Average rank over the experiments:\n', file=file)
    #     for j, name in enumerate(sorted(ranks_by_name,
    #                                     key=lambda x: np.mean(ranks_by_name[x]['worsts']))):
    #         print('{:2}: {:40} {:.2f}'.format(j, name, np.mean(ranks_by_name[name]['worsts'])), file=file)
    #
    #     print('Worst rank over the experiments:\n', file=file)
    #     for j, name in enumerate(sorted(ranks_by_name,
    #                                     key=lambda x: max(ranks_by_name[x]['worsts']))):
    #         print('{:2}: {:40} {:.2f}'.format(j, name, max(ranks_by_name[name]['worsts'])), file=file)
    #
    # if ranking_method == 1:
    #     print('Average rank over the experiments:\n', file=file)
    #     for j, name in enumerate(sorted(ranks_by_name,
    #                                     key=lambda x: np.mean(ranks_by_name[x]['averages']))):
    #         print('{:2}: {:40} {:.2f}'.format(j, name, np.mean(ranks_by_name[name]['averages'])), file=file)
    #
    #     print('Worst rank over the experiments:\n', file=file)
    #     for j, name in enumerate(sorted(ranks_by_name,
    #                                     key=lambda x: max(ranks_by_name[x]['averages']))):
    #         print('{:2}: {:40} {:.2f}'.format(j, name, max(ranks_by_name[name]['averages'])), file=file)
    # """

