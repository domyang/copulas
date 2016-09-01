import math
import random
import sys

import dateutil.parser as dt
import ephem
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from dateutil.relativedelta import relativedelta
from pylab import rcParams
from scipy import optimize as opt
from scipy import stats
from scipy.interpolate import interp1d

home_dir = 'C:\\users\\sabrina\\documents\\research\\code for real user\\'

sys.path.append(home_dir + 'prescient/release/Prescient_1.0/')
sys.path.append(home_dir + 'prescient/release/Prescient_1.0/exec')
import ScenGen.workflow.EpiModel as EpiModel


# ------------------------------ index ---------------------------------------------------------------

# 1 # interpolation
# 2 # data handling
# 3 # plots
# 4 # distributions
# 5 # marginals
# 6 # matrix
# 7 # solar hours
# 8 # other


# ------------------------------ utils for interpolation ------------------------------------------------

class override_options:
    epifit_error_norm = 'L2'
    seg_N = 5
    seg_s = 0.01
    seg_kappa = 50
    L1Linf_solver = 'gurobi'
    L2Norm_solver = 'gurobi'


### returns the spline interpolation coefficients
# returns
#   () dico coefficients for the interpolation
# arguments:
#   () x,y points to be interpolated
#   () visualize: if you wish to plot an example
def find_spline(x, y, options=None, visualize=False, seg_N=5, positiveness_constraint=False,
                increasingness_constraint=False, title=''):
    if options is None:
        options = override_options()
        options.segN = seg_N

    if visualize:
        print('seg_N: %d, seg_s: %.2f, seg_kappa: %.2f ,epifit_error_norm: %s, L1Linf_solver: %s, L2Norm_solver: %s' %
              (options.seg_N, options.seg_s, options.seg_kappa, options.epifit_error_norm, options.L1Linf_solver,
               options.L2Norm_solver))

    # offset_x=min(x)
    # offset_y=min(y)
    # factor_x=(max(x)-offset_x+0.000001)/5
    # factor_y=(max(y)-offset_y+0.000001)/5
    # print('factor %.3f %.3f' %(factor_x,factor_y))
    # x_temp=[(i-offset_x)/factor_x for i in x]
    # y_temp=[(i-offset_y)/factor_y for i in y]
    #
    # print(x_temp);print(y_temp)

    offset = min(x)
    factor_x = (max(x) - offset + 0.000001)
    factor_y = (max(y) - min(y) + 0.000001)
    x_t = [float((i - offset) / factor_x) for i in x]
    y_t = [float(i / factor_y) for i in y]
    x_temp = []
    y_temp = []
    for i in range(len(x)):
        if not x_t[i] in x_temp:
            x_temp.append(x_t[i])
            y_temp.append(y_t[i])

    model = EpiModel.FitEpispline(x_temp, y_temp,
                                  options)  # ,positiveness_constraint=positiveness_constraint,increasingness_constraint=increasingness_constraint)

    approx = []
    for i in model.s.keys():
        approx.append(factor_y * (model.s.get_values())[i])

    # getting coefficients for the estimating function
    N = options.seg_N
    min_x, max_x = (0, max(x_temp))
    s0, v0 = (model.s0.get_values()[None], model.v0.get_values()[None])

    a = []
    keys = model.a.keys()
    for i in keys:
        a.append((model.a.get_values())[i])

    dico = {'a': a, 's0': s0, 'v0': v0, 'min_x': min_x, 'max_x': max_x, 'N': N, 'factor_x': factor_x,
            'factor_y': factor_y, 'offset': offset, 'positiveness_constraint': positiveness_constraint}

    if visualize:
        print(dico)
        plt.figure()
        plt.title(title)
        f = create_spline(dico)
        plt.plot(x, y, color='blue')  # True points
        plt.plot(x, approx, color='red')  #
        index = np.arange(offset, max_x * factor_x + offset, max_x * factor_x / 100)
        val = [f(i) for i in index]
        # print(val)
        plt.plot(index, val, color='green')

    return dico


### returns a spline function
# returns
#   () f : spline function
# arguments
#   () dico coefficients for the interpolation
def create_spline(dico):
    a = dico['a']
    s0 = dico['s0']
    v0 = dico['v0']
    min_x = dico['min_x']
    max_x = dico['max_x']
    N = dico['N']
    offset = dico['offset']
    factor_x = dico['factor_x']
    factor_y = dico['factor_y']

    # estimating function
    def f(t):
        t_temp = (t - offset) / factor_x
        delta = (max_x - min_x) / N
        len_a = len(a)
        k = max(0, min(int(np.floor(t_temp / delta)), len_a - 1))
        res = 0
        for i in range(k):
            res += (t_temp + (0.5 - i - 1) * delta) * a[i]
        res *= delta
        res += 0.5 * (t_temp - k * delta) ** 2 * a[k]
        res += s0 + t_temp * v0
        res *= factor_y
        if dico['positiveness_constraint']:
            res = max(res, 0)
        return res

    return (f)


### returns the derivative corresponding to a spline function
def create_spline_d(dico):
    a = dico['a']
    s0 = dico['s0']
    v0 = dico['v0']
    min_x = dico['min_x']
    max_x = dico['max_x']
    N = dico['N']
    offset = dico['offset']
    factor_x = dico['factor_x']
    factor_y = dico['factor_y']

    # estimating function
    def f(t):
        t = (t - offset) / factor_x
        delta = (max_x - min_x) / N
        len_a = len(a)
        k = max(0, min(int(np.floor(t / delta)), len_a - 1))
        res = 0
        for i in range(k):
            res += a[i]
        res *= delta
        res += (t - k * delta) * a[k]
        res += v0
        res *= factor_y / factor_x
        return (res)

    return f


### if x is sorted, and y=f(x), return a linear interpolation function
def interpol_simple(x, y, positive=True, exp_tails=True):
    length = len(x)
    if (len(y) != length):
        print('WARNING: x and y should be same length')
        return None

    last = x[0]
    x_temp, y_temp, length_temp = [], [], 0
    for i in range(1, length):
        xi = x[i]
        if xi > last:
            x_temp.append(xi)
            y_temp.append(y[i])
            length_temp += 1
            last = xi
    x, y, length = x_temp, y_temp, length_temp

    if length < 2:
        print('WARNING: not enough data to compute interpolation. Is the x-coordinate sorted?')
        return None

    mi, ma = x[0], x[-1]

    def f(a):
        if mi < a < ma:
            ind = np.searchsorted(x, a, side='left')
            s = a - x[ind - 1]
            t = x[ind] - a
            return (t * y[ind - 1] + s * y[ind]) / (s + t)
        elif a >= ma:
            s, t = a - x[-2], x[-1] - a
            return (t * y[-2] + s * y[-1]) / (s + t)
        else:
            s, t = a - x[1], x[2] - a
            return (t * y[1] + s * y[2]) / (s + t)

    if not positive:
        return f
    else:
        def g(a):
            return max(f(a), 0)

        return g


### use interpolation splines to approximate a scipy.stats object fast
class fast_stats(object):
    def __init__(self, obj, prec=5000, method='linear'):

        # for the ppf
        max_derivative = 20
        min_ppf = opt.brentq(lambda x: 1 / obj.pdf(obj.ppf(x)) - max_derivative, 0.01, 0.5)
        max_ppf = opt.brentq(lambda x: 1 / obj.pdf(obj.ppf(x)) - max_derivative, 0.5, 0.99)

        # for the cdf
        min_cdf = obj.ppf(min_ppf)
        max_cdf = obj.ppf(max_ppf)
        xx = np.linspace(min_cdf, max_cdf, prec)
        bounds_cdf = (xx[6], xx[-7])

        yy = list(map(obj.cdf, xx))
        bounds_ppf = (yy[6], yy[-7])

        ppf_estimate = interp1d(yy, xx, kind=method, assume_sorted=True)
        cdf_estimate = interp1d(xx, yy, kind=method, assume_sorted=True)

        def ppf(x):
            is_list = np.ndim(x) > 0
            if is_list:
                if len(x) == 1:
                    x = x[0]
                    is_list = False

            if not is_list:
                if bounds_ppf[0] < x < bounds_ppf[1]:
                    return [float(ppf_estimate(x))]
                else:
                    return [float(obj.ppf(x))]
            else:
                return obj.ppf(x)

        def cdf(x):
            is_list = np.ndim(x) > 0
            if is_list:
                if len(x) == 1:
                    x = x[0]
                    is_list = False

            if not is_list:
                if bounds_cdf[0] < x < bounds_cdf[1]:
                    return [float(cdf_estimate(x))]
                else:
                    return [float(obj.cdf(x))]
            else:
                return obj.cdf(x)

        zz = list(map(obj.pdf, xx))
        pdf_estimate = interp1d(xx, zz, kind=method, assume_sorted=True)

        def pdf(x):
            is_list = np.ndim(x) > 0
            if is_list:
                if len(x) == 1:
                    x = x[0]
                    is_list = False

            if not is_list:
                if bounds_cdf[0] < x < bounds_cdf[1]:
                    return [float(pdf_estimate(x))]
                else:
                    return [float(obj.pdf(x))]
            else:
                return obj.pdf(x)

        self.obj = obj
        self.ppf = ppf
        self.cdf = cdf
        self.pdf = pdf


### finds the quantiles of a function that returns the CDF value of different points
def simultaneous_quantile(f, obs, epsilon=1e-7, prec=1e-7):
    nb_iter = int(np.ceil(-np.log(prec) / np.log(2)))
    length = len(obs)
    floor_obs = np.array([epsilon] * length)
    width = 1 - 2 * epsilon

    def zerone(x):
        if x > 0:
            return 1
        else:
            return 0

    for i in range(nb_iter):
        width = width / 2
        values = f(floor_obs + width)
        floor_obs = np.array([zerone(j[2] - j[0]) * width + j[1] for j in zip(*[values, floor_obs, obs])])

    return floor_obs + width / 2


# ------------------------------ utils for data handling ------------------------------------------------


### removes values in 'l' when the value of 'k' is in 'when'
# returns a list of list:
#   () res : same structure as k if k is a list of list, a list of list anyway
# arguments:
#   () l: a list or a list of lists (in this case children of same length as k)
#   () k: a list to specify when to remove elements
def remove_in_list(l, k, when=[], reverse=False):
    temp = k.copy()
    templ = l.copy()
    if (type(l) != list):
        print('first argument must be of type list')
        return -1
    elif type(l[0]) != list:
        templ = [templ]
    length = len(k)
    nblists = len(templ)
    for i in range(nblists):
        if len(l[i]) != length:
            print('first argument\'s childs must of same length than second argument (num %i: %d, 2nd arg: %d)' % (
                i, len(l[i]), length))
            return -1
    res = [[] for i in range(nblists)]
    while True:
        b = temp.pop() in when
        if reverse:
            b = not b
        if b:
            for i in range(nblists):
                templ[i].pop()
        else:
            for i in range(nblists):
                res[i].append(templ[i].pop())
        if temp == []:
            break

    for i in res:
        i.reverse()
    return res


### checks the type (list of child)
# returns a tuple:
#   () res: boolean indicating if the returned list has the good format
#   () l: list cast to the right format if possible
def is_list_of(l, child=int):
    res = False
    if l is not None:
        if isinstance(l, child):
            l = [l]
        if isinstance(l, list):
            if len(l) > 0:
                res = True
                for i in l:
                    if type(i) != child:
                        res = False
                        break
    return (res, l)


### check that vects is a list of non empty lists of the same size
# returns a tuple:
#   () length (length of the elements of vects)
#   () format (boolean value)
def good_format(vects):
    length = -1
    format = np.ndim(vects) == 2
    if format:
        format &= len(vects) > 0
    if format:
        length = len(vects[0])
    return length, format


### recursively maps a function on the leaves of an n-dimensional array
def map_list(f, l):
    if type(l) in {list, np.ndarray}:
        for i in range(len(l)):
            l[i] = map_list(f, l[i])
        return (l)
    else:
        return f(l)


### assuming the list represents values of a function at regular intervals (of width 1),
### returns approximate values of the derivative
def list_derivative(l, width=5):
    l = l.copy()
    smooth = [math.exp(-(i * 0.5 / width) ** 2) for i in range(width)]
    s = sum(smooth)
    for i in range(width):
        smooth[i] = smooth[i] / s
    res = []
    wid = width * 2 + 1
    buff = []
    for i in range(wid):
        buff.append(l.pop())
    l.append(buff[-1])
    index = width - 1
    pof = []
    for i in range(len(l) + wid - 1 - 2 * (width)):
        buff[(index - width) % wid] = l.pop()
        index += 1
        index %= wid
        temp = 0
        for j in range(1, 1 + width):
            temp -= (buff[(index + j) % wid] - buff[(index - j) % wid]) * smooth[j - 1] / (2 * j)
        res.append(temp)
        pof.append(buff[index])
    temp = res[-1]
    res.extend([temp for i in range(width)])
    res.reverse()
    temp = res[-1]
    res.extend([temp for i in range(width)])
    return res


# ------------------------------ utils for plots ------------------------------------------------


# returns a color list corresponding to 'vect' values
# returns:
#   () col : list of 3-elements tuples (->rgb colors)
# arguments:
#   () vect : values to match the color to
#   () typ : type of color pattern (if int <0.5 and>=0, exclude the values in the interval 10% - 90%)
#   () intervalls : to get constant colors within the intervalls
def color_value(vect, typ='default', date=False, intervalls=[]):
    vect = vect.copy()
    length = len(vect)
    col = []
    ma, mi = max(vect), min(vect)

    if intervalls != []:
        b = True
        t = 0
        for i in intervalls:
            print(b)
            print(i)
            b &= (type(i) == float) & (i >= t)
            t = i
        b &= t <= 1
        if b:
            def norm(x):
                x = (x - mi) / (ma - mi)
                t = 0
                inter = [0]
                inter.extend(intervalls)
                for i in range(len(inter)):
                    if inter[i] < x:
                        t = i
                    else:
                        break
                inter.append(1)
                t = (inter[t] + inter[t + 1]) / 2
                return t
        else:
            def norm(x):
                return (x - mi) / (ma - mi)
    else:
        def norm(x):
            return (x - mi) / (ma - mi)

    if date:
        vect = [2 * min(i - mi, ma - i) for i in vect]

    if (type(typ) == float):
        if (typ < 0.5) & (typ > 0):
            for i in vect:
                if norm(i) < typ:
                    col.append('blue')
                elif norm(i) > 1 - typ:
                    col.append('red')
                else:
                    col.append('transparent')
        else:
            print('error in argument \'typ\'')
            return -1

    elif typ == 'old':
        for i in vect:
            fact = norm(i)
            col.append(
                (0.5 * fact + 2 * fact * (1 - fact), 2 * fact * (1 - fact), 0.5 - 0.5 * fact + 2 * fact * (1 - fact)))
    elif typ == 'blueshade':
        for i in vect:
            fact = norm(i)
            col.append((1 - fact, 1 - fact, 1))
    else:
        for i in vect:
            fact = norm(i)
            col.append((fact ** 2, 2 * fact * (1 - fact), (1 - fact) ** 2))

    return col


### plots a list of points:
# if dimension is 1,    plotted on a line
# if dimension is 3,    3D plot
# else                  2D plot of the first two coordinates
# use noTranspose=True if you have a list of coordinate lists.
def plot(points, noTranspose=False):
    if noTranspose:
        vec = points
        dim = len(points)
        length = len(points[0])
    else:
        vec = list(zip(*points))
        dim = len(points[0])
        length = len(points)

    if dim == 3:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.scatter3D(vec[0], vec[1], vec[2])
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.scatter3D([0], [0], [0], c='red')
    elif dim > 1:
        plt.figure()
        plt.plot(vec[0], vec[1], '.')
    elif dim == 1:
        plt.plot(vec[0], [1 for i in range(length)], '.')
    else:
        print('points has not the good format')


### plots the values of a function on a given interval
def curve3d(f, inter1=(0, 1), inter2=(0, 1), points=True, nb=25, title=''):
    xx = [(inter1[1] - inter1[0]) * (i % nb + 0.5) / nb + inter1[0] for i in range(nb ** 2)]
    yy = [(inter2[1] - inter2[0]) * (i // nb + 0.5) / nb + inter2[0] for i in range(nb ** 2)]
    if points:
        zz = list(map(f, list(zip(*[xx, yy]))))
    else:
        zz = f([xx, yy])

    print("%r , %r, %r" % (xx, yy, zz))
    print("%r , %r, %r" % (len(xx), len(yy), len(zz)))
    # plot 3d
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_trisurf(xx, yy, zz)
    plt.title(title)

    # contours
    m = [zz[i * nb:(i + 1) * nb] for i in range(nb)]
    plt.figure()
    # plt.contour(m)
    plt.contourf(m, 100)
    plt.title(title)


### plots the values of a function on a given interval
def curve(f, inter=[0, 1], title='', labels=[], nb=100):
    if np.ndim(f) == 1:
        funcs = f
    else:
        funcs = [f]
    if type(labels) == str:
        labels = [labels]
    x = [i * (inter[1] - inter[0]) / nb + inter[0] for i in range(nb + 1)]
    plt.figure()
    plt.title(title)
    if len(funcs) == len(labels):
        for ff in zip(*[funcs, labels]):
            y = [ff[0](i) for i in x]
            plt.plot(x, y, label=ff[1])
        plt.legend()
    else:
        for ff in funcs:
            y = [ff(i) for i in x]
            plt.plot(x, y)


### This function prints a table in latex language
def table_latex(l, xlabels=None, ylabels=None, title=None):
    if type(l[0]) == dict:
        keys = list(l[0].keys())
        l = [[line[k] for k in keys] for line in l]
        temp = keys
        if xlabels is not None:
            if type(xlabels) == dict:
                temp = [xlabels[k] for k in keys]
        xlabels = temp

    hei = len(l)
    wid = len(l[0])

    if ylabels is not None:
        if len(ylabels) != hei:
            print('wrong length for \'ylabels\' ')
            ylabels = None

    if xlabels is not None:
        if len(xlabels) != wid:
            print('wrong length for \'xlabels\' ')
            xlabels = None

    text = '\n\n\\begin{center}\n'
    if title is not None:
        text = '%s\\textbf{%s}\\\\ \n\\ \\\\ \n' % (text, title)
    text = '%s\\begin{tabular}{||' % text

    if ylabels is not None:
        text = '%s c |' % text

    for i in range(wid):
        text = '%s c ' % text
    text = '%s||}\n' % text
    text += '\\hline\n\\hline\n'

    if xlabels is not None:
        if ylabels is not None:
            text = '%s&' % text
        text += ' & '.join(['\\textbf{%s}' % el for el in xlabels])
        text = '%s \\\\ \n' % text
        text += '\\hline\n'

    for line in range(hei):
        if ylabels is not None:
            text = '%s\\textbf{%s} & ' % (text, ylabels[line])
        for i, el in enumerate(l[line]):
            if np.ndim(el) == 1 | (type(el) == set):
                text += '('
                for el_child in el:
                    text = '%s %s,' % (text, el_child)
                text = '%s) & ' % text[:-1]
            else:
                text = '%s%s ' % (text, el)
                if i != len(l[line]) - 1:
                    text += '& '
        text = '%s \\\\ \n' % text[:-1]
    text = '%s\\hline\n\\hline\n\\end{tabular} \n\\end{center}\n\n' % text
    print(text)
    return text


### This method plots a bar-plot in wich the values are in increasing order (it sorts xlabels as well)
def sorted_barplot(val, title=None, xlabels=None, min_zero=False):
    if xlabels is not None:
        wid = max([len(i) for i in xlabels]) / 10 * len(val)
    else:
        wid = 9

    if min_zero:
        m = min(val)
        val = [i - m for i in val]

    if type(val) == dict:
        keys = list(val.keys())
        val = [val[k] for k in keys]
        temp = keys
        if xlabels is not None:
            if type(xlabels) == dict:
                temp = [xlabels[k] for k in keys]
        xlabels = temp

    xvalues = list(range(len(val)))
    indices = sorted(xvalues, key=val.__getitem__)
    val_bis = list(map(val.__getitem__, indices))
    xlabels_bis = list(map(xlabels.__getitem__, indices))

    plt.figure(figsize=(wid, 6))
    if title is not None:
        plt.title(title)
    plt.bar(xvalues, val_bis)
    plt.xticks(np.array(xvalues) + 0.4, xlabels_bis)
    rcParams.update({'font.size': 9})


def hist3D(x, y, bins):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    hist, xedges, yedges = np.histogram2d(x, y, bins=bins)

    elements = (len(xedges) - 1) * (len(yedges) - 1)
    xpos, ypos = np.meshgrid(xedges[:-1], yedges[:-1])

    xpos = xpos.flatten()
    ypos = ypos.flatten()
    zpos = np.zeros(elements)
    dx = 1 / bins * np.ones_like(zpos)
    dy = dx.copy()
    dz = hist.flatten()

    return ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='b', zsort='average')
# ------------------------------ for distributions estimates -----------------------------------------------------------

### returns the pdf function of a gaussian fitted to vects
def create_gaussian_density(vects):
    dim = len(vects)
    covariance = np.cov(vects)
    means = np.mean(vects, axis=1)
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


### returns the empirical cumulative distribution function
def empirical_CDF(unifs):
    # checking the format of unifs: list of lists of lists
    msg = 'unifs must be a list of lists: unifs is a list of coordinates representing copula points'

    if not is_list_of(unifs, child=list)[0]:
        print(msg)
        return None

    dim = len(unifs)
    length = len(unifs[0])
    for j in range(1, dim):
        if len(unifs[j]) != length:
            print(msg)
            return None

    def CDF(pt):
        if len(pt) != dim:
            return None
        else:
            res = 0
            for i in zip(*unifs):
                b = True
                for j in range(dim):
                    if pt[j] < i[j]:
                        b = False
                if b:
                    res += 1
        return res / length

    return CDF


# Returns a function [(x,y)->int]
#   () f: the spline 2D-repartition function
# arguments
#   () (unif1,unif2): two vectors of same length (-> 2d points)
def compute_distribution(unif1, unif2):
    count = len(unif1)
    if count == len(unif1):
        dim = max(1, int(np.floor(np.sqrt(count / 16))))
        distr = [[0 for i in range(dim)] for j in range(dim)]
        points = [[[] for i in range(dim)] for j in range(dim)]
        temp1 = unif1.copy()
        temp2 = unif2.copy()
        for i in range(count):
            a = temp1.pop()
            b = temp2.pop()
            x = min(max(0, int(np.floor(a * dim))), dim - 1)
            y = min(max(0, int(np.floor(b * dim))), dim - 1)
            distr[x][y] = distr[x][y] + 1
            points[x][y].append((a, b))
        del (temp1, temp2)
        for i in range(dim - 1):
            distr[1][i + 1] += distr[1][i]
            distr[i + 1][1] += distr[i][1]
        for i in range(dim - 1):
            for j in range(dim - 1):
                distr[i + 1][j + 1] += distr[i][j + 1] + distr[i + 1][j] - distr[i][j]

        def f(a, b):
            x = min(max(0, int(np.floor(a * dim))), dim - 1)
            y = min(max(0, int(np.floor(b * dim))), dim - 1)
            res = distr[x][y]
            for t in points[x][y]:
                if (a >= t[0]) & (b >= t[1]):
                    res += 1
            return res / count

        return f
    else:
        print('unif1 and unif2 must be same length')
        return -1


### return the Kernel estimate of the density function: (a function)
# arguments:
#   () vects: (list of list) points of the distribution
#   () f: a function to use as kernel
#   () K: number of points to define the kernel's bandwidth
def kernel_estimate(vects, func=None, K=10):
    # a few characteristics of the distribution
    dim = len(vects)
    length = len(vects[0])
    vec = list(zip(*vects))
    cov = np.cov(vects)
    if dim > 1:
        cov /= np.linalg.det(cov)
    else:
        cov = np.matrix([[1]])

    def distance(a, b, m=cov):
        n = np.matrix(a) - np.matrix(b)
        return (n * m * np.transpose(n))[0, 0]

    # computing the maximum distance to the K nearest neighbours

    # clustering the data
    nb_cluster = int(np.sqrt(length))
    cluster = random.sample(vec, nb_cluster)
    cluster_pts = [[] for i in range(nb_cluster)]
    pts_cluster = []

    def compare(a, d):
        if d is None:
            return True
        else:
            return a < d

    def next_cluster(a, levelm=0, levelM=None):
        nearest = -1
        incr = 0
        d = None
        for j in cluster:
            d_temp = distance(i, j)
            if (levelm < d_temp) & compare(d_temp, d) & compare(d_temp, levelM):
                nearest = incr
                d = d_temp
            incr += 1
        return (nearest, d)

    # creating the clusters around selected points
    print('creating clusters')
    for i in vec:
        n = next_cluster(i)
        cluster_pts[n[0]].append(i)
        pts_cluster.append(n)

    # finding the nearest neighbours
    print('finding the nearest neighbours')
    neigh_dist = []
    incr_i = 0
    for i in vec:
        dist_temp = []
        incr = 0
        max_temp = 0
        ind_temp = 0
        while (len(dist_temp) < K) & (incr <= (K + 1)):
            if incr != incr_i:
                dist_temp.append(distance(vec[incr], i))
            incr += 1
        for k in range(K):
            if dist_temp[k] > max_temp:
                max_temp = dist_temp[k]
                ind_temp = k

        cl = pts_cluster[incr_i]
        while cl[0] != -1:
            for j in cluster_pts[cl[0]]:
                d = distance(i, j)
                if d < max_temp:
                    dist_temp[ind_temp] = d
                    max_temp = 0
                    for k in range(K):
                        if dist_temp[k] > max_temp:
                            max_temp = dist_temp[k]
                            ind_temp = k
            if cl[1] < max_temp * 2:
                cl = next_cluster(i, levelm=cl[1], levelM=max_temp * 2)
            else:
                cl = (-1, 0)
        neigh_dist.append(max_temp)
        incr_i += 1

    print('definition of the kernel estimate\n')
    # definition of kernel estimate
    final = list(zip(*[vec, neigh_dist]))

    def d(x):
        if func is None:
            res = 0
            for i in final:
                res += math.exp(-distance(x, i[0]) / i[1]) / math.sqrt(math.pi * i[1])
            return res / length
        else:
            if type(x) in {float, int}:
                x = [x]
            res = 0
            for i in final:
                res += func(x, i[0], i[1])
            return res / length

    return d


### return the expectancy and variance of a non stationary time-serie at a given moment
def weighted_moments(obs, band_width=20):
    band_width = int(band_width)

    length = len(obs)
    weights = [math.exp(-x * 2 / band_width ** 2) for x in range(-3 * band_width + 1, 3 * band_width)]
    E = []  # expectancy
    V = []  # variance

    for i in range(length):
        weight_sum = 0
        weight_square_sum = 0
        E_temp = 0
        V_temp = 0
        for j in range(-band_width + 1, band_width):
            if not 0 <= i + j < length:
                continue
            w = weights[j]
            obs_el = obs[i + j]
            weight_sum += w
            weight_square_sum += w ** 2
            E_temp += obs_el * w
            V_temp += obs_el ** 2 * w

        E_temp /= weight_sum
        V_temp -= weight_sum * E_temp ** 2
        V_temp *= weight_sum / (weight_sum ** 2 - weight_square_sum)

        E.append(E_temp)
        V.append(V_temp)

    return E, V


# -------------------------------- for distances btw distributions -----------------------------------------------------

"""
# Returns a float:
#   () res: the L1-distance between the two repartition function, as estimated from the points
# arguments
#   () (unif_a1,unif_a2): two vectors of same length (-> 2d points)
#   () (unif_b1,unif_b2): two other vectors of same length (-> 2d points)
def compute_distance_emd(a1, a2, b1, b2, precision=10):
    # with the earth mover distance
    counta = len(a1)
    if (len(a2) != counta):
        print('given vectors must be same length, here: %d and %d' % (len(a1), len(a2)))
        return -1
    countb = len(b1)
    if (len(b2) != countb):
        print('given vectors must be same length, here: %d and %d' % (len(b1), len(b2)))
        return -1
    prec2 = precision ** 2
    distance = np.array([[float(
        np.sqrt((i % precision - j % precision) ** 2 + (i // precision - j // precision) ** 2)) / precision for j in
                          range(precision ** 2)] for i in range(precision ** 2)])
    sa = np.array([0. for i in range(prec2)])
    sb = np.array([0. for i in range(prec2)])
    inva = float(1 / counta)
    for i in range(counta):
        index = min(int(np.floor(a1[i] * precision)), precision - 1) + precision * min(int(np.floor(a2[i] * precision)),
                                                                                       precision - 1)
        sa[index] = sa[index] + inva
    #
    invb = float(1 / countb)
    for i in range(countb):
        index = min(int(np.floor(b1[i] * precision)), precision - 1) + precision * min(int(np.floor(b2[i] * precision)),
                                                                                       precision - 1)
        sb[index] = sb[index] + invb
    return emd(sa, sb, distance)
"""


# uses a plugin to compute the emd (earth mover distance) between two lists of coordinate lists
def compute_emd(v1, v2, norm=2):
    length = len(v1)
    if len(v2) != length:
        print('vectors must be same length')
        return None
    width = len(v1[0])

    # initializing distance matrix
    m = np.identity(2 * length)
    for i in range(2 * length):
        for j in range(i, 2 * length):
            m[i, j] = 0
            if i != j:
                if i < length:
                    w1 = v1[i]
                else:
                    w1 = v2[i - length]
                if j < length:
                    w2 = v1[j]
                else:
                    w2 = v2[j - length]

                temp = 0

                for k in range(width):
                    temp += (w1[k] - w2[k]) ** norm
                m[i, j] = temp ** (1 / norm)
                m[j, i] = m[i, j]

    array1 = np.array([float(i // length) for i in range(2 * length)])
    array2 = np.array([1 - i for i in array1])
    return emd(array1, array2, m)


# returns a L-distance between two copulae:
# arguments:
#   () unif1,unif2: list of lists -> points of a copula
#   () precision
#   () rand: if true, using monte carlo method else a rectangular grid
#   () norm: computing the L_norm distance
def compute_distance_l(unif1, unif2, precision=30, rand=False, norm=2):
    dim = len(unif2)
    if len(unif1) != dim:
        print('the two distributions must have the same dimension')
        return None

    maximum = 1000000
    if precision ** dim > maximum:
        precision = int(maximum ** (1 / dim))

    F1 = empirical_CDF(unif1)
    F2 = empirical_CDF(unif2)
    res = 0

    if rand:
        test = [np.random.uniform(0, 1, precision ** dim).tolist() for i in range(dim)]
        for i in zip(*test):
            res += (F1(i) - F2(i)) ** norm
    else:
        inv = 1 / precision
        p = [precision ** (i + 1) for i in range(dim)]
        incr = [0 for i in range(dim)]
        test = [inv / 2 for i in range(dim)]
        for i in range(precision ** dim):
            res += (F1(test) - F2(test)) ** norm
            incr[0] += 1
            test[0] += inv
            for j in range(dim - 1):
                if incr[j] == p[j]:
                    incr[j] = 0
                    incr[j + 1] += 1
                    test[j] = inv / 2
                    test[j + 1] += inv

    return (res / precision ** dim) ** (1 / norm)

    # bin_edge=max(int((max(length/5,0))**(1/dim)),2)
    # dim_eff=bin_edge**dim
    # bins=[0 for j in range(dim_eff)]
    # bin_points=[[] for j in range(dim_eff)]
    # p=[bin_edge**i for i in range(dim)]
    # inv_edge=1/bin_edge
    #
    # for point in zip(*unifs):
    #     ind=0
    #     for i in range(dim):
    #         ind+=min(bin_edge-1,int(point[i]*bin_edge))*p[i]
    #     bins[ind]+=1
    #     bin_points[ind].append(point)
    # def f(x):
    #     return x/length
    # bins=list(map(f,bins))
    #
    # for i in range(dim_eff):
    #     r=[int(i//(j*bin_edge)) for j in p]
    #     temp=bins[i]
    #     for j in range(dim):
    #         if r[j]>0:
    #             temp+=bins[i-p[j]]
    #     bins[i]=temp
    #
    # def CDF(pt):
    #     if len(pt)!=dim:
    #         return None
    #     ind=0
    #     ind_inf=0
    #     for i in range(dim):
    #         ind+=min(bin_edge-1,int(point[i]*bin_edge))*p[i]
    #     for i in range(dim):
    #         coor=min(bin_edge-1,int(point[i]*bin_edge))*(p[i]-1)
    #         if coor<=0:
    #             ind_inf=0
    #             break
    #         else:
    #             ind_inf+=coor
    #
    #     res=bins[]


# returns the log likelihood of vect1,
def emp_log_likelihood(vect1, vect2, density1=None, density2=None, visualize=False):
    dim, dim2 = len(vect1), len(vect2)
    if (dim != dim2) | (dim < 1):
        raise (RuntimeError('the dimensions of the distributions must be equal and superior than 0'))

    def aux(vec1, vec2, den):
        len1, len2 = len(vec1[0]), len(vec2[0])
        if den is None:
            kde = stats.gaussian_kde(vec1)
            res = 0
            for pt in list(zip(*vec2)):
                res += math.log(kde(pt))
        else:
            print(vec2)
            den_list = den(vec2)
            print(den_list)
            res = np.sum(list(map(math.log, den_list)))
        res /= len2
        return res

    a, b = aux(vect1, vect2, density1), aux(vect2, vect1, density2)
    return ((a + b) / 2, a, b)


# returns the projection of vectors along axis [1,1..1],[-1,1,..1],[1,-1,..,1] etc as a list of lists
# arguments:
#   () copula: if True, it supposes that vects are uniforms and makes it so that the projections are in (0,1)
def tail_projection(vects, copula=True):
    dim = len(vects)

    nb_axes = 2 ** (dim - 1)

    res = []
    axis = [-1] * dim
    axis[-1] = 1
    denom = np.sqrt(2) ** (dim)

    for i in range(nb_axes):
        temp = []
        starting_corner = [(1 - j) / 2 for j in axis]
        if copula:
            offset = np.dot(starting_corner, axis) / denom
        else:
            offset = 0

        for pt in zip(*vects):
            value = np.dot(list(pt), axis) / denom - offset
            temp.append(value)
        res.append(temp)

        axis = table_increment(axis, mi=-1, ma=3, incr=2, start=0)

    return res


# returns the EMD (earth-mover-distance) between the first and last quantile of two list of points.
# vec1 and vec2 need not have the same length
def univariate_EMD_in_tails(vec1, vec2, quantile=0.1):
    a = sorted(vec1)
    b = sorted(vec2)

    nb_a = max(int(len(a) * quantile), 1)
    nb_b = max(int(len(b) * quantile), 1)

    factor = nb_a / nb_b

    remains = [1, factor]
    it_a = iter(a[:nb_a])
    it_b = iter(b[:nb_b])
    current = [next(it_a), next(it_b)]

    res = [0, 0]
    try:
        while True:

            if remains[0] < remains[1]:
                res[0] += abs(current[0] - current[1]) * remains[0]
                current[0] = next(it_a)
                remains[1] -= remains[0]
                remains[0] = 1
            else:
                res[0] += abs(current[0] - current[1]) * remains[1]
                current[1] = next(it_b)
                remains[0] -= remains[1]
                remains[1] = factor
    except StopIteration:
        pass

    remains = [1, factor]
    it_a = iter(a[-nb_a:])
    it_b = iter(b[-nb_b:])
    current = [next(it_a), next(it_b)]

    try:
        while True:

            if remains[0] < remains[1]:
                res[1] += abs(current[0] - current[1]) * remains[0]
                current[0] = next(it_a)
                remains[1] -= remains[0]
                remains[0] = 1
            else:
                res[1] += abs(current[0] - current[1]) * remains[1]
                current[1] = next(it_b)
                remains[0] -= remains[1]
                remains[1] = factor
    except StopIteration:
        pass

    res = [i / nb_a for i in res]

    return res


# returns the EMD distance between the projections of simulations items and past_obs
# also returns the quantile of the projections of obs as compared to the projections of simulations items
def compare_tails(simulations, past_obs, obs, quantile=0.1, visualize=False):
    proj_past_obs = tail_projection(past_obs)
    proj_obs = tail_projection([[i] for i in obs])

    fit = []
    rank = []

    for sim in simulations:

        proj_tp = tail_projection(sim)
        fit_tp = []
        for i in zip(*[proj_past_obs, proj_tp]):
            fit_tp.append(univariate_EMD_in_tails(i[0], i[1], quantile=quantile))
            if visualize:
                curve(stats.gaussian_kde(i[1]), inter=[-1000, 1000])

        fit.append(fit_tp)

        rank_tp = []
        for i in range(len(proj_obs)):
            rank_tp.append(empirical_CDF_scalar(proj_tp[i])(proj_obs[i][0]))

        rank.append(rank_tp)

    return fit, rank


# ------------------------------ for distributions: others -----------------------------------------------------------

# returns:
#   () points distributed according to a normal law, fitted to the input data
#   () the length of the output vectors
# arguments:
#   () vect1,vect2: specifies the 2 data range of same length to fit the normal law
def redistribute_gaussian(vect1, vect2, length=0, fit=True):
    if (length > 50000) or (length < 10):
        count = len(vect1)
    else:
        count = length
    if fit:
        cov = np.cov(np.matrix([vect1, vect2]))
        print(cov)
    else:
        cov = np.array([[1, 0], [0, 1]])
        print(cov)
    print('in gaussian, count=%d' % count)
    redistributed = np.transpose(np.random.multivariate_normal([0, 0], cov, count))

    return (redistributed[0], redistributed[1], count)


# returns pearson's rho (correlation of the ranks = correlation of the empirical copula)
# arguments:
#   () vects: list of coordinates lists
#   () holes: value of 'vects' that should be ignored (for missing data...)
def pearson_with_holes(vects, hole=None):
    dim = len(vects)

    res = np.identity(dim)
    for i in range(dim):
        for j in range(i + 1, dim):
            vec = list(zip(*[vects[i], vects[j]]))
            vec_bis = []
            incr = 0
            for k in vec:
                if (k[0] != hole) and (k[1] != hole):
                    vec_bis.append(k)
                    incr += 1
            if len(vec_bis) > 1:
                vec_bis = list(zip(*vec_bis))
                for a in [0, 1]:

                    ind = sorted(range(incr), key=vec_bis[a].__getitem__)
                    vec_bis[a] = [0 for i in range(incr)]
                    for k in range(incr):
                        vec_bis[a][ind[k]] = k
                cov = np.cov(vec_bis)
                res[i, j] = cov[0, 1] / math.sqrt(cov[0, 0] * cov[1, 1])
            else:
                res[i, j] = 0
            res[j, i] = res[i, j]

    return res


### returns a rough 2d integral of a pdf to check if this is equal to 1
# arguments:
#   () pdf: function: [a,b] -> [pdf(a,b)]
#   () mi,ma: min and max bounds
def test_pdf(pdf, mi=-1000, ma=1000, prec=30, scalar=None):
    s = 0

    if scalar is None:
        scalar = False
        if type(pdf([(ma - mi) / 2, (ma - mi) / 2])) in {float, np.float}:
            scalar = True

    fact = (ma - mi) / prec
    for i in range(prec):
        for j in range(prec):
            if scalar:
                s += pdf([(i + 0.5) * fact + mi, (j + 0.5) * fact + mi])
            else:
                s += pdf([(i + 0.5) * fact + mi, (j + 0.5) * fact + mi])[0]
    s *= fact ** 2
    return s


# ------------------------------ for marginals ------------------------------------------------------------------------



# finds quantiles, given a CDF
# l is the list specifying the desired quantiles
def find_quantiles(cdf, l, start=0, step=0.1, precision=0.00001):
    max_iter = 1000
    res = []
    cur = start
    length = len(l)
    indices = sorted(range(length), key=l.__getitem__)
    ordered_quantiles = [l[i] for i in indices]

    for i in ordered_quantiles:
        cur_step = step
        while cdf(cur) < i:
            cur += cur_step
        while cdf(cur) > i:
            cur -= cur_step
        diff = 1
        inter = [cur, cur + step]
        incr = 0
        while (incr < max_iter) & (abs(diff) > precision):
            incr += 1
            cur = (inter[0] + inter[1]) / 2
            diff = i - cdf(cur)
            if diff > 0:
                inter = [cur, inter[1]]
            else:
                inter = [inter[0], cur]
        res.append(cur)
        if incr >= max_iter:
            print('max number of iteration reached')

    ind_inverse = sorted(range(length), key=indices.__getitem__)
    result = [res[i] for i in ind_inverse]

    return result


# returns an empirical cdf function, given a list of observations
# if exp_tails, fits an exponential function at the extremities
def empirical_CDF_scalar(obs, exp_tails=True, mi=None, ma=None):
    # an auxiliar function: doesn't take min  or max bounds
    def aux(obs, exp_tails=False):
        length = len(obs)
        wei = []  # cdf difference between current and previous point
        pts = []  # points without doubles
        store = 0
        first = True
        for i in sorted(obs):
            if first:
                pts.append(i)
                wei.append(1)
                first = False
            else:
                if i == pts[-1]:
                    wei[-1] += 0.5
                    store += 0.5
                else:
                    pts.append(i)
                    wei.append(1 + store)
                    store = 0

        wei = [i / (length + 1) for i in wei]
        cum_wei = np.cumsum(wei)  # cdf value of the point

        # print(wei)
        # print(pts)
        # print(cum_wei)

        def f(x):
            ind = int(np.searchsorted(pts, x))

            if ind == 0:
                if exp_tails:
                    # exponential tail with matching derivative
                    return wei[0] * np.exp(-(pts[0] - x) * wei[1] / (wei[0] * (pts[1] - pts[0])))
                else:
                    # linear tail with matching derivative
                    return max(0, wei[0] * (1 - (pts[0] - x) * wei[1] / (wei[0] * (pts[1] - pts[0]))))
            elif ind == len(pts):
                if exp_tails:
                    return 1 - (1 - cum_wei[-1]) * np.exp(
                        - (x - pts[-1]) * wei[-1] / ((1 - cum_wei[-1]) * (pts[-1] - pts[-2])))
                return min(1, 1 - (1 - cum_wei[-1]) * (
                1 - (x - pts[-1]) * wei[-1] / ((1 - cum_wei[-1]) * (pts[-1] - pts[-2]))))
            else:
                a = (x - pts[ind - 1]) / (pts[ind] - pts[ind - 1])
                return cum_wei[ind - 1] * (1 - a) + cum_wei[ind] * a

        return f

    bi, ba = (mi is not None), (ma is not None)

    if bi:
        if mi > min(obs):
            raise (RuntimeError('minimum given (%f) is higher than some observations (%f)' % (mi, min(obs))))
    if ba:
        if ma < max(obs):
            raise (RuntimeError('maximum given (%f) is lower than some observations (%f)' % (ma, max(obs))))
    if len(set(obs)) < 2:
        raise (RuntimeError('impossible to compute a CDF with less than 2 points'))

    # taking min and max into account
    if bi or ba:
        obs_bis = []
        if bi:
            obs_bis.append(mi)
        obs_bis.extend(obs)
        if ba:
            obs_bis.append(ma)
        g = aux(obs_bis)

        def f(x):
            if bi and x < mi:
                return 0
            if ba and x > ma:
                return 1

            # creating linear CDF with two added observations: mi and ma
            res = g(x)
            # normalizing the interpolation to get 0 at mi and 1 at ma (unless doubles)
            if bi:
                res -= 1 / (len(obs_bis) + 1)
            res *= (len(obs_bis) + 1) / (len(obs) + 1)
            return res

        return f

    # no min and max given
    else:
        return aux(obs, exp_tails=exp_tails)


# returns the inverse of a cdf function, given a list of observations
# if log_tails, fits a log function at the extremities
# (this is the almost exact inverse of the previous)
def empirical_CDF_inv(obs, log_tails=True, mi=None, ma=None):
    # an auxiliar function: doesn't take min  or max bounds
    def aux(obs, log_tails=False):
        length = len(obs)
        wei = []  # cdf difference between current and previous point
        pts = []  # points without doubles
        store = 0
        first = True
        for i in sorted(obs):
            if first:
                pts.append(i)
                wei.append(1)
                first = False
            else:
                if i == pts[-1]:
                    wei[-1] += 0.5
                    store += 0.5
                else:
                    pts.append(i)
                    wei.append(1 + store)
                    store = 0

        wei = [i / (length + 1) for i in wei]
        cum_wei = np.cumsum(wei)  # cdf value of the point

        # print(wei)
        # print(pts)
        # print(cum_wei)

        def f(x):
            if (x > 1) or (x < 0):
                raise (RuntimeError('x should be comprised between 0 and 1'))
            if log_tails and ((x == 1) or (x == 0)):
                raise (RuntimeError('x should be comprised between 0 and 1'))

            ind = int(np.searchsorted(cum_wei, x))

            if ind == 0:
                if log_tails:
                    # log tail with matching derivative
                    return np.log(x / wei[0]) * wei[0] * (pts[1] - pts[0]) / wei[1] + pts[0]
                else:
                    # linear tail with matching derivative
                    return pts[0] + (wei[0] - x) / wei[1] * (pts[0] - pts[1])
            elif ind == len(pts):
                if log_tails:
                    return pts[-1] - np.log((1 - x) / (1 - cum_wei[-1])) / wei[-1] * (
                    (1 - cum_wei[-1]) * (pts[-1] - pts[-2]))
                else:
                    return pts[-1] + (x - cum_wei[-1]) * (pts[-1] - pts[-2]) / wei[-1]
            else:
                a = (x - cum_wei[ind - 1]) / (cum_wei[ind] - cum_wei[ind - 1])
                return pts[ind - 1] * (1 - a) + pts[ind] * a

        return f

    bi, ba = (mi is not None), (ma is not None)

    if len(set(obs)) < 2:
        raise (RuntimeError('impossible to compute a CDF with less than 2 points'))

    # taking min and max into account
    if bi or ba:
        obs_bis = []
        if bi:
            obs_bis.append(mi)
        obs_bis.extend(obs)
        if ba:
            obs_bis.append(ma)
        g = aux(obs_bis)
        length_bis = (len(obs_bis) + 1)
        ratio = (len(obs) + 1) / length_bis

        def f(x):
            x_bis = x * ratio
            if bi:
                x_bis += 1 / length_bis
            return g(x_bis)

        return f

    # no min and max given
    else:
        return aux(obs, log_tails=log_tails)


### returns a function that computes the approximate CDFs (cumulative density functions) of the marginals of a distribution.
### the approximate is using epi-splines
# arguments:
#   () vects: list of coordinate lists of the observations
#   () seg_N : argument for the epi-spline approximation
# returns:
#   () res: a list of CDF functions
def marginals_cdf(vects, seg_N=10, interpolate='linear'):
    res = []
    dim = len(vects)
    if dim == 0:
        raise (RuntimeError('dimension of array must be >0'))
    if interpolate == 'epi-spline':
        bins = 32
        for i in range(dim):
            h = np.histogram(vects[i], bins=bins)
            xx = []
            for j in range(bins):
                xx.append(float(h[1][j] + h[1][j + 1]) / 2)
            yy = [h[0][0]]
            for j in range(bins - 1):
                yy.append(float(yy[j] + h[0][j + 1]))

            dico = find_spline(xx, yy, seg_N=seg_N, visualize=False, positiveness_constraint=True,
                               increasingness_constraint=True)
            f = create_spline(dico)
            df = create_spline_d(dico)

            mi, ma = f(h[1][0]), f(h[1][-1])
            div = 1 / (ma + mi)
            ma *= div
            mi *= div
            dmi, dma = df(h[1][0]) * div, df(h[1][-1]) * div

            if (mi > 0) and (dmi > 0) and (dma > 0):
                mi_a = dmi / mi
                mi_b = mi / (math.exp(mi_a * h[1][0]))
                ma_a = dma / (1 - ma)
                ma_b = (1 - ma) / math.exp(-ma_a * h[1][-1])

                def temp(x):
                    if x < h[1][0]:
                        return mi_b * math.exp(mi_a * x)
                    elif x > h[1][-1]:
                        return 1 - ma_b * math.exp(-ma_a * x)
                    else:
                        return f(x) * div

                res.append(temp)
            else:
                def temp(x):
                    if x < h[1][0]:
                        return 0
                    elif x > h[1][-1]:
                        return 1
                    else:
                        return (f(x) - mi / div) / (ma - mi) * div

                res.append(temp)
    else:
        for i in vects:
            res.append(empirical_CDF_scalar(i, mi=-3000, ma=3000))

    return res


### similarily, computes the inverses of the marginals'CDFs
def marginals_cdf_inv(vects, seg_N=10, limits=None, interpolate='linear'):  # 'epi-spline'):
    bins = 100
    res = []
    dim = len(vects)
    limited = isinstance(limits, tuple) and (len(limits) == 2)

    if dim == 0:
        print('dimension of array must be >0')
        return None

    for i in range(dim):
        if interpolate == 'epi-spline':
            # creating a histogram of the distribution
            hist, bin_edges = np.histogram(vects[i], bins=bins)
            xx = []
            for j in range(bins):
                xx.append(float(bin_edges[j] + bin_edges[j + 1]) / 2)
            yy = [hist[0]]

            div = 0
            for j in range(bins - 1):
                div = float(yy[j] + hist[j + 1])
                yy.append(div)

            div -= yy[0]
            print(div)

            # creating the spline estimator
            dico = find_spline(yy, xx, seg_N=seg_N,
                               visualize=False)  # ,positiveness_constraint=True,increasingness_constraint=True)
            f = create_spline(dico)

            # if the support of the distribution is known
            if limited:
                def temp(a):
                    return (min(limits[1], max(limits[0], f(a * div))))

                res.append(temp)
            else:
                def temp(a):
                    return (f(a * div))

                res.append(temp)

        else:
            if limited:
                res.append(empirical_CDF_inv(vects[i], mi=limits[0], ma=limits[1]))
            else:
                res.append(empirical_CDF_inv(vects[i]))

    return res


### similarily, computes the marginals' PDFs
def marginals_pdf(vects, seg_N=8, interpolate='kernel'):
    length = len(vects[0])
    bins = int(math.ceil(math.sqrt(length)))
    res = []
    dim = len(vects)
    if dim == 0:
        raise (RuntimeError('dimension of array must be >0'))
    for i in range(dim):

        if interpolate == 'kernel':
            res.append(stats.gaussian_kde(vects[i]))
            continue

        h = np.histogram(vects[i], bins=bins)
        xx = []
        for j in range(bins):
            xx.append(float(h[1][j] + h[1][j + 1]) / 2)
        norm = length * (xx[-1] - xx[0]) / bins
        yy = [i / norm for i in h[0]]

        if interpolate == 'epi-spline':
            dico = find_spline(xx, yy, seg_N=seg_N, visualize=False, positiveness_constraint=True,
                               increasingness_constraint=True)
            f = create_spline(dico)
        elif interpolate == 'linear':

            # linear interpolation with exponential tails
            decr = bins / (xx[-1] - xx[0])

            def f(x, xx=xx, yy=yy, decr=decr, bins=bins):
                ind = np.searchsorted(xx, x)
                if ind == 0:
                    return yy[0] * math.exp(-(xx[0] - x) * decr)
                elif ind == bins:
                    return yy[-1] * math.exp(-(x - xx[-1]) * decr)
                else:
                    return (yy[ind - 1] * (xx[ind] - x) + yy[ind] * (x - xx[ind - 1])) / (xx[ind] - xx[ind - 1])

        else:
            f = interpol_simple(xx, yy)
            raise (RuntimeWarning('\'interpolate\' argument is not valid'))

        # normalizing the pdf
        mi, ma = min(vects[i]), max(vects[i])
        inter_length = ma - mi
        mi -= 0.1 * inter_length
        ma += 0.1 * inter_length
        fact = (ma - mi) / (10 * bins)
        div = 0
        for i in range(10 * bins):
            div += f(mi + (i + 0.5) * fact)
        div *= fact

        def f_normalized(x, f=f):
            return f(x) / div

        res.append(f_normalized)

    return res


### using previous functions, computes marginals to create distribution points from copula points
def copula_to_distribution(vects, visualize=False):
    dim = len(vects)
    f_inv = marginals_cdf_inv(vects)

    def func(unifs, visu=visualize):
        if len(unifs) != dim:
            raise RuntimeError('vects and unifs must have same dimensions')
        res = [[f_inv[i](x) for x in unifs[i]] for i in range(dim)]

        if visu:
            plt.figure()
            plt.plot(vects[0], vects[1], '.')
            plt.title('empirical')
            plt.figure()
            plt.plot(res[0], res[1], '.')
            plt.title('approximate')
        return res

    return func


### given a distribution and the densities of various copulae, returns the densities of the corresponding distributions
def copula_to_densities(vects, densities, visualize=False, log_return=False):
    return_list = True
    if callable(densities):
        return_list = False
        densities = [densities]
    dim = len(vects)
    mar_pdf = marginals_pdf(vects, interpolate='kernel')
    mar_cdf = marginals_cdf(vects, interpolate='linear')

    if visualize:
        mi = min(min(vects))
        ma = max(max(vects))
        for i in mar_pdf:
            curve(i, inter=[mi, ma])
        for i in mar_cdf:
            curve(i, inter=[mi, ma])

    res = []
    for i in densities:
        if i is not None:
            def f(x, den=i, log_return=log_return):
                if type(x[0]) in {float, int, np.float, np.int}:
                    x = [[j] for j in x]
                pdf = [list(map(mar_pdf[j], x[j])) for j in range(dim)]
                cdf = [list(map(mar_cdf[j], x[j])) for j in range(dim)]
                if log_return:
                    prod = [np.sum(list(map(math.log, j))) for j in zip(*pdf)]
                    return [j[0] + math.log(j[1]) for j in zip(*[prod, den(cdf)])]
                else:
                    prod = [np.prod(j) for j in zip(*pdf)]
                    return [j[0] * j[1] for j in zip(*[prod, den(cdf)])]

            res.append(f)
        else:
            res.append(None)

    if return_list:
        return res
    else:
        return res[0]


### given a distribution's observations, and a list of pdf approximates for the distribution,
###  returns the densities of the corresponding copula
def distribution_to_copula_densities(vects, densities, log_return=False):
    return_list = True
    if callable(densities):
        return_list = False
        densities = [densities]
    dim = len(vects)
    mar_pdf = marginals_pdf(vects, interpolate='kernel')
    mar_cdf_inv = marginals_cdf_inv(vects, interpolate='linear')

    res = []
    for i in densities:
        if i is not None:
            def f(x, den=i, log_return=log_return):
                if type(x[0]) in {float, int, np.float, np.int}:
                    x = [[j] for j in x]
                cdf = [list(map(mar_cdf_inv[j], x[j])) for j in range(dim)]
                pdf = [list(map(mar_pdf[j], cdf[j])) for j in range(dim)]
                if log_return:
                    prod = [np.sum(list(map(math.log, j))) for j in zip(*pdf)]
                    return [math.log(j[1]) - j[0] for j in zip(*[prod, den(cdf)])]
                else:
                    prod = [np.prod(j) for j in zip(*pdf)]
                    return [j[1] / j[0] for j in zip(*[prod, den(cdf)])]

            res.append(f)
        else:
            res.append(None)

    if return_list:
        return res
    else:
        return res[0]


# returns the rank of vects normalized by the length.
# If rand, it will simulate uniform points by replacing the rank by the value of the n^th point
def uniforms(vects, rand=True):
    length = len(vects[0])
    res = []
    for vect in vects:
        ranks = stats.rankdata(vect)
        if (rand):
            unifs = sorted(np.random.rand(length))
            ## Sort uniformly distributed points by ranks
            res.append([unif for (rank, unif) in sorted(zip(ranks, unifs), key=lambda pair: pair[0])])
        else:
            res.append((ranks / (length + 1)).tolist())

    return res


# ------------------------------ utilities for matrices ------------------------------------------------

def conditional_matrix_gaussian(M):
    res = []
    length = len(M)
    for i in range(length):
        m11 = np.matrix(M[i, i])
        m12 = np.zeros([1, length - 1])
        m22 = np.zeros([length - 1, length - 1])
        inc1 = 0
        for j in range(length - 1):
            if inc1 == i:
                inc1 += 1
            m12[0, j] = M[i, inc1]
            inc2 = 0
            for k in range(length - 1):
                if inc2 == i:
                    inc2 += 1
                m22[j, k] = M[inc1, inc2]
                inc2 += 1
            inc1 += 1
        m12 = np.matrix(m12)
        m22 = np.matrix(m22)

        print('m11: %r, \n m12: %r \n m22: %r' % (m11, m12, m22))
        res.append(float(m11 - m12 * m22 ** (-1) * np.transpose(m12)))

    return res


# ------------------------------ utilities for solar hour ------------------------------------------------

sun = ephem.Sun()
california = ephem.Observer()
california.lat = '37'
california.lon = '-122'


### returns the sunrise and sunset time (in Californian time)
def sun_rise_set(date):
    california.date = date
    day = int(california.date % 365.25)
    day = day - day % 7

    offset = -8
    if 67 <= day < 305:
        offset += 1
    offset /= 24
    california.date -= offset
    rise = (california.next_rising(sun) + offset + 0.5) % 1 * 24
    set = (california.next_setting(sun) + offset + 0.5) % 1 * 24
    # print('%d:%d %d:%d'% (rise//1,rise%1*60//1,set//1,set%1*60//1))
    return (rise, set)


### This function prepares the solar data to take into account the changing daily pattern of this resource (-> solar hours)
# returns a tuple containing:
#   () date: ordered dates
#   () hour_sol: the solar hours
#   () the corresponding observation
# arguments:
# options (for the epi spline interpolation)
# l: a dictionary containing the data at keys: 'act' (actuals), 'for' (forecasts) and 'date'
def prepare_solar(l, visualize=False):
    length = len(l['date'])
    error = [actual - forecast for actual, forecast in zip(l['act'], l['for'])]
    hour_sol = []
    date = []

    # renormalizing the errors by square root of the maximum forecast (proportional to the number of generators...)
    error_temp = []
    actuals = []
    forecasts = []
    f_max = 1
    forecast_max = []

    nb_slices = 32

    # Creates a list of maximums of f[0:i] for i in len(f)
    for forecast in l['for']:
        f_max = max(f_max, forecast)
        forecast_max.append(f_max)

    for i in range(50, length):
        time = l['date'][i]
        err = error[i] / math.sqrt(forecast_max[i])

        # estimating the solar hour
        rise, sete = sun_rise_set('%d-%d-%d' % (time.year, time.month, time.day))
        sol_hour = ((time.hour + time.minute / 60 + time.second / 3600) - rise) / (sete - rise) * 12

        if 0.1 < sol_hour < 11.9:
            date.append(str(time))
            hour_sol.append(sol_hour)
            error_temp.append(err)
            forecasts.append(l['for'][i])
            actuals.append(l['act'][i])

    obs = error_temp.copy()

    # estimating the mean square error as a function of time
    length1 = len(hour_sol)
    variance = []
    incr = 0
    nb_obs = length1 // nb_slices
    x = [0] * (nb_slices + 1)
    slices = [[] for _ in range(nb_slices + 1)]
    indices = sorted(range(length1), key=hour_sol.__getitem__)
    obs = [obs[i] for i in indices]
    hour_sol_sorted = [hour_sol[i] for i in indices]

    for error, solar_hour in zip(obs, hour_sol_sorted):
        slices[incr // nb_obs].append(error)
        x[incr // nb_obs] += solar_hour
        incr += 1

    for i, slice in enumerate(slices[:nb_slices]):
        len_temp = len(slice)
        sd = math.sqrt(np.var(slice))
        variance.append(sd)
        if (i == 0) or (i == nb_slices - 1):
            variance.append(sd)
        x[i] /= len_temp
    x = [0] + x[:-1] + [12]

    # creating the spline estimator of the mean square error
    inter = create_spline(find_spline(x, variance, visualize=False, seg_N=nb_slices // 2))

    # normalizing the error
    normed_error = []
    for solar_hour, error in zip(hour_sol, error_temp):
        normed_error.append(error / inter(solar_hour))

    if visualize:
        plt.figure()
        plt.plot(hour_sol, error_temp, '.')
        abscisses = [i / 10 for i in range(120)]
        y = [inter(i) for i in abscisses]
        plt.plot(abscisses, y, color='red')

        plt.figure()
        plt.plot(hour_sol, normed_error, '.')
        plt.draw()

    return {'date': date, 'hour_sol': hour_sol, 'error': normed_error, 'for': forecasts, 'act': actuals,
            'var_function': inter}


def intersect_dates(intervals, current, width, visualize=False):
    if isinstance(current, str):
        current = dt.parse(current)
    if np.ndim(intervals[0]) == 0:
        intervals = [intervals]
    if isinstance(intervals[0][0], str):
        intervals = [(dt.parse(a), dt.parse(b)) for (a, b) in intervals]

    incr = 2
    if visualize:
        plt.figure()
        plt.plot([incr, incr + 0.03, incr + .03, incr], [0, 0, 0.03, 0.03])

    res = []
    for start, end in intervals:
        mi = (start - current).days
        ma = (end - current).days
        print(mi, ma)

        if visualize:
            plt.plot([mi, ma], [incr - 1, incr - 1], c='blue')

        if 0 <= (mi + width) % 365.25 <= 2 * width:
            temp = int(((mi + width) // 365.25) * 365.25 + width)
            res.append((str(start), str(current + relativedelta(days=min(temp, ma)))))

            if visualize:
                print('1: %r' % [mi, min(temp, ma)])
                plt.plot([mi, min(temp, ma)], [incr, incr], c='red')
                incr += 0.1

        if 0 <= (ma + width) % 365.25 <= 2 * width:
            temp = int(((ma + width) // 365.25) * 365.25 - width)
            if temp > mi:
                res.append((str(current + relativedelta(days=temp)), str(end)))
                if visualize:
                    print('2: %r' % [temp, ma])
                    plt.plot([temp, ma], [incr, incr], c='red')
                    incr += 0.1

        for j in range(int(math.ceil((mi + width) / 365.25)), 1 + int(math.floor((ma - width) / 365.25))):
            res.append((str(relativedelta(days=int(j * 365.25) - width) + current),
                        str(relativedelta(days=int(j * 365.25) + width) + current)))
            if visualize:
                print('3: %r' % [int(j * 365.25) - width, int(j * 365.25) + width])
                plt.plot([int(j * 365.25) - width, int(j * 365.25) + width], [incr, incr], c='red')
                incr += 0.1

    return res


# dr_par=[]
#         cur=dt.parse(i[2]).day
#         for j in date_range:
#             mi,ma=(j[0]-i[2]).days,(j[1]-i[2]).days # difference between the date and the date range bounds (negative)
#             # mi,ma=j[0].day,j[1].day
#             do_loop=True
#             if not 0<=(-win_days-mi)%365.25<=2*win_days:
#                 temp=int(mi+(win_days-mi)%365.25 +2*win_days)
#                 if temp<ma:
#                     dr_par.append((mi,ma))
#                     do_loop=False
#                 else:
#                     dr_par.append((mi,temp))
#
#
#
#
#
#             if do_loop:
#                 start=(mi+win_days)//365.25
#                 end=(ma-win_days)//365.25
#                 for k in range(start,end):
#                     dr_par.append(int(start*365.25))


# --------------------------------- other -------------------------------------------------------------------------------

### this function returns a list used to go through a multidimensional grid on [mi,ma]^n, with increment 'incr
### When it has went through all the grid, it returns None
def table_increment(t, mi=0, ma=2, incr=1, length=None, start=0):
    if length is None:
        length = len(t)
    if t[start] + incr < ma:
        t[start] += incr
        return t
    elif length > start + 1:
        temp = table_increment(t, mi=mi, ma=ma, incr=incr, length=length, start=start + 1)
        if temp is not None:
            temp[start] = mi
        return temp
    else:
        return None
