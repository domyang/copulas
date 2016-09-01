import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import math
import utilities as ut
import copula as cop
from datetime import timedelta

# READ ME:
#
# In this class, we define a copula class for MULTIPLE SERIES of data and multiple offsets
# The goal is to take into account the dependencies between what happened at t0 and at t0+offset, AS WELL AS SPACE DEPENDENCIES
# Say: wind power forecasts errors, at t=0 and t=5 hours, for South and North California
#
# The window is the data that we consider at one point: it can be specified by a date range (for the first date),
# the value of the forecasts, of the forecasts derivative, an hour range in the day (-> 'first_hour'), etc...

# In this Multi series managers, we use the intersection of the single series windows (dates that were considered everywhere)
#
# Therefore you can add an exotic constraint, like one on the derivative of Solar power forecast in Antarctic, by:
#       - adding a fictive serie 'Solar errors in Antarctic'
#       - setting the offsets to [] (this will make it virtual)
#       - setting the parameters to ...'forecast_d':(min,max)... for this series
# This can also be done by creating a boolean list, 'just_parameter', with True at the same index as the Antarctic serie in series



class CopulaManagerMS:
    nbSeries = 0
    dim = 0
    parameters = []
    # list telling which copulae are only to create the window
    just_parameters = None

    copulae = []
    vectM = [[]]
    unifM = [[]]
    lengthM = 0
    date = []
    indexes = [[]]

    # initialisation: this will create single-series copula managers and aggregate their data
    # 1 - a few necessary checks and initialisations
    # 2 - creation of the single series copula-managers
    # 3 - taking only the dates that correspond
    # 4 - retrieving the values (saved in vectM)
    # 5 - computing the copula points again
    # arguments:
    #   () series: specifies the data - FORMAT: [ {'date':[],'vect':[],'data':{...},'title':{} }, ...] WITH:
    #                                   ... 'data' : {'forecast':[],'forecast_d':[]}
    #                                   ... 'title' : {'type':'', 'location':'', 'kind':''}
    #   () parameters_def : default parameters which specifies the window
    #   () list_parameters : a list of dictionaries (->parameters for each element of series) overriding the default
    #   () just_parameters : boolean list; if true, the corresponding element of series is only used to build the window
    def __init__(self, series, parameters_def, list_parameters=None, just_parameters=None):
        self.parameters = [ser['title'] for ser in series]

        if series is None:
            raise RuntimeError("you need to give a 'series' argument to create a copulaManager")

        self.nbSeries = nbSeries = len(series)

        if nbSeries < 1:
            raise RuntimeError('Series should contain at least 1 element')

        for ser in series:
            if not {'date', 'vect', 'data', 'title'}.issubset(ser.keys()):
                raise (RuntimeError(
                    "each element of series should have at least these entry key: {'date','vect','data','title'}"))

        self.parameters = [parameters_def.copy() for _ in range(nbSeries)]

        for ser, par in zip(series, self.parameters):
            par.update(ser['title'])

        copulae = []
        for ser, par in zip(series, self.parameters):
            # print('%r\n%r\n%r\n%r'% (ser['date'],ser['vect'],ser['data'],par[i]))
            copulae.append(
                cop.CopulaManager(ser['date'], ser['vect'], ser['data'], par,
                                  var_function=ser['var_function']))
        self.copulae = copulae

        self.update(parameters_def, list_parameters=list_parameters, just_parameters=just_parameters)

    def update(self, parameters_def, list_parameters=None, just_parameters=None):
        # 1 - a few necessary checks and initialisations
        if just_parameters is None:
            just_parameters = self.just_parameters

        nbSeries = self.nbSeries

        if just_parameters is None:
            just_parameters = [False] * nbSeries

        par = [parameters_def.copy() for _ in range(nbSeries)]

        if list_parameters is not None:
            if (not isinstance(list_parameters, list)) or (len(list_parameters) != self.nbSeries):
                print(
                    "WARNING: list_parameters should be either None or a list same length as series")
                list_parameters = None

        if list_parameters is not None:
            for j in range(nbSeries):
                i = list_parameters[j]
                if i is not None and isinstance(i, dict):
                    par[j] = i

        for i in range(nbSeries):
            for key in {'type', 'location', 'kind'}:
                par[i][key] = self.parameters[i][key]
            # defines sub_parameters in case the serie's offsets are an empty list
            if (not isinstance(par[i]['offsets'], list)) or (len(par[i]['offsets']) == 0):
                par[i]['offsets'] = [0]
                just_parameters[i] = True

        # 2 - creation of the single series copula-managers

        copulae = self.copulae
        for i in range(nbSeries):
            copulae[i].update(par[i])

        # 3 - taking only the dates that correspond

        date_max = min([copula.dateM[-1] for copula in copulae])
        date_min = max([copula.dateM[0] for copula in copulae])

        # print('date range: %s %s' % (str(date_min), str(date_max)))

        date_iters = [iter(copula.dateM) for copula in copulae]

        current_indexes = [0] * nbSeries
        date_indexes = [[] for _ in range(nbSeries)]
        dates_element = [next(date_iter) for date_iter in date_iters]
        dateM = []

        cur = date_min
        while cur < date_max:

            b = True
            for i in range(nbSeries):
                while dates_element[i] < cur:
                    dates_element[i] = next(date_iters[i])
                    current_indexes[i] += 1

                if dates_element[i] != cur:
                    b = False

            if b:
                dateM.append(cur)
                for i in range(nbSeries):
                    date_indexes[i].append(current_indexes[i])
                cur += timedelta(hours=1)

            else:
                cur_tp = max(dates_element)
                if cur < cur_tp:
                    cur = cur_tp
                else:
                    cur += timedelta(hours=1)

        # 4 - retrieving the values (saved in vectM)
        series_to_consider = [i for i in range(nbSeries) if not just_parameters[i]]

        vectM = []
        for i in series_to_consider:
            for vec in self.copulae[i].vectM:
                vectM.append([vec[j] for j in date_indexes[i]])

        dataM = []
        for ser in series_to_consider:
            dataM_tp = {}
            for key in copulae[ser].dataM:
                dataM_tp[key] = [[copulae[ser].dataM[key][dim][i] for i in date_indexes[ser]] for dim in
                                 range(copulae[ser].dim)]
            dataM.append(dataM_tp)

        # 5 - computing the copula points again

        self.dim = len(vectM)
        self.dateM = dateM
        self.vectM = vectM
        self.dataM = dataM
        self.unifM = ut.uniforms(vectM, rand=False)
        self.lengthM = len(dateM)
        self.indexes = date_indexes
        self.parameters = par
        self.just_parameters = just_parameters
        self.copulae = copulae

    @property
    def date_range(self):
        a = str(self.dateM[0])
        b = str(self.dateM[-1])
        return (a, b)

    def get_corner_points(self, cornerFraction=0.5, endFraction=0.3, visualize=True):

        length = self.lengthM
        dim = self.dim
        nb_axes = 2 ** (dim - 1)
        unifs = [[] for i in range(length)]
        vects = [[] for i in range(length)]
        visualize &= (dim == 2) | (dim == 3)

        if not 0.05 < cornerFraction <= 1:
            raise (RuntimeError('cornerFraction must be between 0.05 and 1'))
        else:
            cornerFraction = int(length * cornerFraction // (nb_axes * 2) + 1)

        if not 0.05 < endFraction <= 1:
            raise (RuntimeError('cornerFraction must be between 0.05 and 1'))
        else:
            endFraction = int(math.ceil(endFraction * cornerFraction))

        unifM = []
        vectM = []
        for v in self.unifM:
            unifM.append(v.copy())

        for v in self.vectM:
            vectM.append(v.copy())

        # print('len: %d, unifM: %r'% (len(unifM),unifM))


        for i in range(length):
            for v in unifM:
                unifs[i].append(v.pop())

            for v in vectM:
                vects[i].append(v.pop())

        # axis=[[]for i in range(nb_axes)]
        axis = []
        cosinus = [[] for _ in range(nb_axes)]
        norms = []
        corners = []

        # only for visualization
        corners_unif = []
        considered_pts = []
        considered_pts_unif = []

        def norm(l):
            temp = 0
            for i in l:
                temp += i * i
            temp = math.sqrt(temp)
            if temp == 0:
                print('null vector: %r' % l)
                return ([0 for _ in l], 0)
            return ([i / temp for i in l], temp)

        print('creation of axes')
        norm_axis = math.sqrt(dim)
        for i in range(nb_axes):
            str = bin(i)[2:]
            temp = [norm_axis for i in range(dim - len(str))]
            for i in str:
                temp.append(-(2 * int(i) - 1) * norm_axis)
            axis.append(temp)

        print('projection of points')
        for pt in unifs:
            pt = [2 * i - 1 for i in pt]
            pt, n = norm(pt)
            norms.append(n)
            for j in range(nb_axes):
                temp = 0
                a = axis[j]
                for k in range(dim):
                    temp += pt[k] * a[k]
                cosinus[j].append(temp)

        print('selecting the data & creating representative points')
        for i in cosinus:
            indexes = sorted(range(length), key=i.__getitem__)
            norms_temp = list(map(norms.__getitem__, indexes))
            vects_temp = list(map(vects.__getitem__, indexes))
            print(norms_temp)
            print(vects_temp)

            v_pos = vects_temp[-cornerFraction:]
            ind_pos = sorted(range(cornerFraction), key=norms_temp[-cornerFraction:].__getitem__)
            pos = list(map(v_pos.__getitem__, ind_pos))
            pos = pos[-endFraction:]
            v_neg = vects_temp[:cornerFraction]
            ind_neg = sorted(range(cornerFraction), key=norms_temp[:cornerFraction].__getitem__)
            neg = list(map(v_neg.__getitem__, ind_neg))
            neg = neg[-endFraction:]

            pt = [0] * dim
            for j in pos:
                for k in range(dim):
                    pt[k] += j[k]
            corners.append([j / endFraction for j in pt])

            pt = [0] * dim
            for j in neg:
                for k in range(dim):
                    pt[k] += j[k]
            corners.append([j / endFraction for j in pt])

            if visualize:
                considered_pts.extend(pos)
                considered_pts.extend(neg)

                unifs_temp = list(map(unifs.__getitem__, indexes))

                u_pos = unifs_temp[-cornerFraction:]
                pos = list(map(u_pos.__getitem__, ind_pos))
                pos = pos[-endFraction:]
                u_neg = unifs_temp[:cornerFraction]
                neg = list(map(u_neg.__getitem__, ind_neg))
                neg = neg[-endFraction:]

                pt = [0] * dim
                for j in pos:
                    for k in range(dim):
                        pt[k] += j[k]
                corners_unif.append([j / endFraction for j in pt])

                pt = [0] * dim
                for j in neg:
                    for k in range(dim):
                        pt[k] += j[k]
                corners_unif.append([j / endFraction for j in pt])

                considered_pts_unif.extend(pos)
                considered_pts_unif.extend(neg)

        print('      number of points per quadrants: %d' % cornerFraction)
        print('      number of points considered per quadrants:%d' % endFraction)

        if visualize and (dim == 3):
            print('plotting')
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.set_title('representative points on the real distribution')
            ax.set_xlabel('X axis')
            ax.set_ylabel('Y axis')
            ax.set_zlabel('Z axis')
            temp = list(zip(*tuple(vects)))
            ax.scatter(temp[0], temp[1], temp[2], c='grey', color='grey')
            temp = list(zip(*tuple(considered_pts)))
            ax.scatter(temp[0], temp[1], temp[2], c='blue', color='blue')
            temp = list(zip(*tuple(corners)))
            ax.scatter(temp[0], temp[1], temp[2], c='red', color='red')

            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.set_title('representative points on the renormalized distribution')
            ax.set_xlabel('X axis')
            ax.set_ylabel('Y axis')
            ax.set_zlabel('Z axis')
            ax.set_xlim3d(0, 1)
            ax.set_ylim3d(0, 1)
            ax.set_zlim3d(0, 1)
            temp = list(zip(*tuple(unifs)))
            ax.scatter(temp[0], temp[1], temp[2], c='grey', color='grey')
            temp = list(zip(*tuple(considered_pts_unif)))
            ax.scatter(temp[0], temp[1], temp[2], c='blue', color='blue')
            temp = list(zip(*tuple(corners_unif)))
            ax.scatter(temp[0], temp[1], temp[2], c='red', color='red')

        if visualize and (dim == 2):
            print('plotting')
            plt.figure()
            plt.title('representative points on the real distribution')
            temp = list(zip(*tuple(vects)))
            plt.plot(list(temp[0]), list(temp[1]), '.', color='grey')
            temp = list(zip(*tuple(considered_pts)))
            plt.plot(temp[0], temp[1], '.', color='blue')
            temp = list(zip(*tuple(corners)))
            plt.plot(temp[0], temp[1], '+', color='red', ms=10, mew=3)

            plt.figure()
            plt.title('representative points on the renormalized distribution')
            temp = list(zip(*tuple(unifs)))
            plt.plot(temp[0], temp[1], '.', color='grey')
            temp = list(zip(*tuple(considered_pts_unif)))
            plt.plot(temp[0], temp[1], '.', color='blue')
            temp = list(zip(*tuple(corners_unif)))
            plt.plot(temp[0], temp[1], '+', color='red', ms=10, mew=3)

        print(unifs[:20])
        print(considered_pts_unif)
        for c in corners:
            print('corner: %r' % c)

        return corners

    def plot2D(self):
        if self.dim >= 2:
            v = []
            incr = 0
            for i in self.vectM:
                v.append(i)
                incr += 1
                if incr > 1:
                    break
            plt.figure()
            plt.plot(v[0], v[1], '.')

    def plot(self, proj=[[]], which='unif', date=True):
        if proj == [[]]:
            dim = self.dim
            for i in range(dim - 1):
                x = [0 for _ in range(dim)]
                y = [0 for _ in range(dim)]
                x[i] = 1
                y[i + 1] = 1
                print([x, y])
                self.plot(proj=[x, y], which=which)
        else:
            argument_check = False
            if type(proj) == np.matrixlib.defmatrix.matrix:
                if proj.shape == (2, self.dim):
                    argument_check = True
            elif type(proj == list):
                if len(proj) == 2:
                    if type(proj[0]) == list:
                        bool_temp = True
                        for i in proj:
                            if len(i) != self.dim:
                                bool_temp = False
                                break
                        argument_check = bool_temp
                if argument_check:
                    proj = np.matrix(proj)
            if not argument_check:
                raise (RuntimeError(
                    'projection matrix must be either a numpy matrix or a list of list, with dimensions (2,%d)' % self.dim))

            # normalisation of the line vectors
            for i in range(2):
                div = float(proj[i] * (proj[i].transpose()))
                proj[i] = proj[i] / div

            # projection des points sur les axes
            if which == 'unif':
                points = (proj * np.matrix(self.unifM)).tolist()
            else:
                points = (proj * np.matrix(self.vectM)).tolist()

            # plotting
            plt.figure()
            if date:
                plt.colorbar(plt.scatter(points[0], points[1], c=self.dateM, cmap='bwr'))
            else:
                plt.plot(points[0], points[1], '.')
            if which == 'unif':
                for i in range(self.dim):
                    x = proj[0, i]
                    y = proj[1, i]
                    plt.plot([0, x], [0, y], color='grey')
                    plt.text(x, y, 'axe %d' % i)

    def tails(self, precision=30, restrict=None, interpolation='epi_spline'):
        return tails_ut(self.unifM, precision=precision, restrict=restrict, interpolation=interpolation)

    def __str__(self):
        maxN = 6
        # print('\n### SCALAR PARAMETERS ### \n \n')
        # print('-> nbSeries= %d \n-> dim= %d' % (self.nbSeries, self.dim))
        # print('\n\n### PARAMETERS AND TERTIARY DATA ###\n')
        # print('copulae (len=%d) ...' % len(self.copulae))
        string = ''
        for i in range(self.nbSeries):
            par = self.parameters[i]
            string = '%s\n   -> %s %s in %s' % (string, par['type'], par['kind'], par['location'])
            if self.just_parameters[i]:
                string = '%s (just parameter)' % string
            temp = par['offsets']
            string = '%s\n      offsets (len=%d): [' % (string, len(temp))
            for j in temp:
                string += '%d, ' % j
            string = string[:-2]
            string += ']\n'
            if 'date_range' in par:
                start_date, end_date = par['date_range']
                string += '      date_range: (' + str(start_date) + ', ' + str(end_date) + ') \n'
            if 'first_hour' in par:
                start_hour, end_hour = par['first_hour']
                string += '      first_hour: (' + str(start_hour) + ', ' + str(end_hour) + ') \n'
            if 'forecast' in par:
                temp = par['forecast']
                if temp != []:
                    string += '      forecast (len=%d): [' % len(temp)
                    for j in temp:
                        if (len(j) == 2) & (type(j) == tuple):
                            string += '(%0.2f, %0.2f), ' % j
                        else:
                            string += '--ERROR--, '
                    string = string[:-2] + ']\n'
            if 'forecast_d' in par:
                temp = par['forecast_d']
                if temp != []:
                    string += '      forecast_d (len=%d): [' % len(temp)
                    for j in temp:
                        if (len(j) == 2) & (type(j) == tuple):
                            string += '(%0.2f, %0.2f), ' % j
                        else:
                            string = '%s--ERROR--, ' % string
                    string = string[:-2] + ']\n'

            string = '%s \n --- values --- \n\n' % string

            length = len(self.dateM)
            string = '%s\n      dateM[%d] (len=%d): [' % (string, i, length)
            for j in range(min(maxN, length)):
                string = '%s \'%s\' ,' % (string, self.dateM[j])
            if length <= maxN:
                string = '%s]\n' % string[:-1]
            else:
                string = '%s...]\n' % string

            ind = 0
            for v in self.vectM:
                ind += 1
                length = len(v)
                string = '%s\n      vectM[%d] (len=%d): [' % (string, ind, length)
                for j in range(min(maxN, length)):
                    string = '%s %0.2f ,' % (string, v[j])
                if length <= maxN:
                    string = '%s]\n' % string[:-1]
                else:
                    string = '%s...]\n' % string
            ind = 0
            for v in self.unifM:
                ind += 1
                length = len(v)
                string = '%s\n      unifM[%d] (len=%d): [' % (string, ind, length)
                for j in range(min(maxN, length)):
                    string = '%s %0.2f ,' % (string, v[j])
                if length <= maxN:
                    string = '%s]\n' % string[:-1]
                else:
                    string = '%s...]\n' % string

        return string

    def pprint(self):
        print(self)


### This adjacent function computes the tail distribution for each corner of the copula specified by unifs
def tails_ut(unifs, precision=30, restrict=None, interpolation='epi_spline', cumulative=True, visualize=False):
    dim = len(unifs)

    ### use restrict parameter to consider only a projection of the copula
    if restrict is not None:
        if type(restrict) == list:
            dim = len(restrict)
    if dim <= 0:
        return None
    nb = 2 ** dim
    quadrants = [[0 for i in range(precision)] for i in range(nb)]

    cop_temp = unifs
    if restrict is not None:
        if type(restrict) == list:
            cop_temp = list(map(unifs.__getitem__, restrict))

    # computes the density of the bands B_ab = { x : 0<=a< dist(x,corner) <b<=0.5 }
    for i in zip(*cop_temp):
        p = 1
        quad = 0
        dist = 0
        for j in i:
            if j > 0.5:
                quad += p
                dist = max(dist, 1 - j)
            else:
                dist = max(dist, j)
            p *= 2
        quadrants[quad][min(int(dist * 2 * precision), precision - 1)] += 1

    # normalizes the density
    quadrants_temp = []
    for quad in quadrants:
        s = sum(quad)
        quadrants_temp.append([i / s * precision for i in quad])
    quadrants = quadrants_temp
    del (quadrants_temp)

    x = [1 - (0.5 + i) / precision for i in range(precision)]
    x_rev = [1 - j for j in x]

    def transform(a, b):
        return math.log(a)

    if visualize:
        fig1 = plt.figure()
        fig2 = plt.figure()

    # computes the tail functions
    tail_functions = [None for j in range(nb)]
    for i in range(nb // 2):
        for ii in [i, nb - 1 - i]:

            # values for cumulative/density function
            if cumulative:
                temp = [0 for j in range(precision)]
                temp[0] = quadrants[ii][0]
                for j in range(precision - 1):
                    temp[j + 1] = temp[j] + quadrants[ii][j + 1]
            else:
                temp = quadrants[ii]

            # log transformed values

            log_x = []
            log_temp = []
            for k in zip(x_rev, temp):
                if (k[1] > 0) & (k[0] > 0):
                    log_x.append(math.log(k[0]))
                    log_temp.append(math.log(k[1]))

            # creates approximate cumulative(f) or density(df) functions for the tails
            if interpolation == 'linear':
                slope, intercept = stats.linregress(log_x, log_temp)[:2]
                print('%r , %r ' % (slope, intercept))
                print('f(x)= %r x^%r' % (math.exp(intercept), slope))

                def f(x):
                    if x > 0:
                        return math.exp(intercept) * x ** slope
                    else:
                        return (0)
            else:
                title = 'corner %d' % ii
                dic = ut.find_spline(x_rev, temp, visualize=visualize, title=title, increasingness_constraint=True,
                                     positiveness_constraint=True)
                f = ut.create_spline(dic)

            tail_functions[ii] = f

            if visualize:
                plt.figure(fig1.number)
                plt.plot(x_rev, temp)
                plt.figure(fig2.number)
                interpol = [f(j) for j in x_rev]
                plt.plot(x_rev, interpol)

    return tail_functions
