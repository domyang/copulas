import matplotlib.pyplot as plt
import math
import random
import numpy as np
from scipy import stats
import scipy.special as spe
from scipy.optimize import root
import scipy.optimize as opt
import subsetModel as sub
import utilities as ut


### upper class for copula models
#   () val: Points of this copula
#   () f:   A function (of 'count') returning 'count' points distributed according to the copula
#   () par: The parameters of the copula
#   () pdf: The pdf of the copula
class cop_mod(object):
    dim = 0
    length = 0
    name = ''
    ACR = ''
    val = []
    par = {}

    def simulate(self, count):
        return [[]]

    def pdf(self, x):
        return 1

    def pprint(self):
        s = '\n### Copula Model: ' + self.name + ' ###\n\n'
        s += 'val: (%d,%d)\n' % (self.dim, self.length)
        for i in range(self.dim):
            s += '    unifs[%d]: [' % i
            for j in range(6):
                s += ' %.3f,' % self.val[i][j]
            s = '%s...]\n' % s[:-1]
        s += 'parameters: %r' % self.par
        print(s)

    def plot2d(self, axe=[0, 1], method={'points'}, nbPoints=None):
        if (self.dim > 1) & (0 <= min(axe) < max(axe) < self.dim):

            if type(nbPoints) == int:
                if 0 < nbPoints < self.length:
                    vals = []
                    for i in self.val:
                        temp = []
                        for j in range(nbPoints):
                            temp.append(i[j])
                        vals.append(temp)
            else:
                vals = self.val

            if len(method.intersection({'contour', 'distribution3D'})) > 0:
                def f(x, y, c):
                    return math.exp((-(x[0] - y[0]) ** 2 - (x[1] - y[1]) ** 2) / 0.01) / (math.pi * 0.01)

                k = ut.kernel_estimate([vals[axe[0]], vals[axe[1]]], func=f)
                if 'contour' in method:
                    fig = plt.figure()
                    plt.title(self.name)
                    z = []
                    for i in range(21):
                        temp = []
                        for j in range(21):
                            temp.append(k([i / 20, j / 20]))
                        z.append(temp)
                    plt.contour(z)
                if 'distribution3D' in method:
                    fig = plt.figure()
                    ax = fig.add_subplot(111, projection='3d')
                    ax.set_title(self.name)
                    x = [(i // 21) / 20 for i in range(21 ** 2)]
                    y = [(i % 21) / 20 for i in range(21 ** 2)]
                    z = []
                    for i in range(21 ** 2):
                        z.append(k([x[i], y[i]]))
                    ax.plot_trisurf(x, y, z)

            if 'points' in method:
                plt.figure()
                plt.title(self.name)
                plt.plot(vals[axe[0]], vals[axe[1]], '.')

    def plot3d(self, axe=[0, 1, 2], nbPoints=None):
        if (self.dim > 2) & (0 <= min(axe) < max(axe) < self.dim):

            if type(nbPoints) == int:
                if 0 < nbPoints < self.length:
                    vals = []
                    for i in self.val:
                        temp = []
                        for j in range(nbPoints):
                            temp.append(i[j])
                        vals.append(temp)
            else:
                vals = self.val

            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.set_title(self.name)
            ax.scatter(vals[axe[0]], vals[axe[1]], vals[axe[2]])


### instead of replacing the values of the observation by their ranks to get an 'empirical copula',
### replacing them by the ordered values of uniforms
class cop_emp_redistributed(cop_mod):
    ACR = 'ER'

    def __init__(self, vects, model=None):
        length = len(vects[0])
        dim = len(vects)
        self.dim = dim
        self.length = length

        def distribute(count):
            nb = int(math.ceil(count / length))
            temp1 = list(zip(*ut.uniforms(vects, rand=True)))
            fin_nb = count % length
            if fin_nb == 0:
                fin_nb = length
            temp2 = list(zip(*random.sample(temp1, fin_nb)))
            val_fin = []
            for i in temp2:
                val_fin.append(list(i))

            if nb > 1:
                val = ut.uniforms(vects, rand=True)
                for i in range(nb - 1):
                    if i < nb - 2:
                        val_temp = ut.uniforms(vects, rand=True)
                    else:
                        val_temp = val_fin
                    for j in range(dim):
                        val[j].extend(val_temp[j])
                return val
            return val_fin

        if model is None:
            self.name = 'empirical redistributed'
        else:
            self.name = model.name + ' redistributed'
        self.simulate = distribute
        self.val = distribute(length)

        # functions to estimate density
        if model is None:
            temp_density = stats.gaussian_kde(vects)

            def prev_density(x):
                return [j / length for j in temp_density(x)]
        else:
            prev_density = model.pdf

        mar_density = ut.marginals_pdf(vects, interpolate='epi-spline')
        mar_inv = ut.marginals_cdf_inv(vects)

        self.prev_density = prev_density
        self.mar_density = mar_density
        self.mar_inv = mar_inv

        # estimating density
        def density(x, vects=vects):
            if np.ndim(x[0]) == 0:
                x = [x]
            else:
                x = list(zip(*x))

            res = []
            for i in x:
                v = [self.mar_inv[j](i[j]) for j in range(dim)]
                d = [self.mar_density[j](v[j]) for j in range(dim)]

                # print('v: %r'%v);print('d: %r'%d)

                temp = 1 / np.prod(d) * self.prev_density(v)[0]
                res.append(temp)
            return res

        # quick fix to make sure the density integrates to 1
        if dim < 4:
            nb = 20
        else:
            nb = max(5, int(10000 ** (1 / dim)))

        nb_inv = 1 / nb
        div = 0
        cur = [nb_inv / 2 for i in range(dim)]
        incr = 0
        while cur is not None:
            incr += 1
            div = div + density(cur)[0]
            cur = ut.table_increment(cur, mi=nb_inv / 2, ma=1, incr=nb_inv)
        div = div / incr

        def density_bis(x):
            return [i / div for i in density(x)]

        self.pdf = density_bis


### returns a 'gaussian copula' fitted to the given points
class cop_gaussian(cop_mod):
    ACR = 'GA'

    def __init__(self, uniforms, length=0):
        if not (ut.good_format(uniforms)[1]):
            raise (RuntimeError('uniforms should be a list of same length list'))

        dim = len(uniforms)
        if (length > 50000) | (length < 10):
            count = len(uniforms[0])
        else:
            count = length

        a = stats.norm([0], [1])
        vects = [a.ppf(coor) for coor in uniforms]

        cov = norm_matrix(np.cov(np.matrix(vects)))
        cov = np.cov(np.matrix(uniforms))  # spearman's rho
        for i in range(dim):
            for j in range(dim):
                cov[i, j] = 2 * math.sin(math.pi / 6 * cov[i, j])  # linear correlation (formula for elliptical copulae)
        cov = norm_matrix(cov)

        print('in gaussian, count=%d' % count)

        def simulate(x):
            redistributed = np.transpose(np.random.multivariate_normal([0 for i in range(dim)], cov, x))
            return ut.uniforms(redistributed.tolist())

        self.dim = dim
        self.length = len(uniforms[0])
        self.name = 'gaussian'
        self.val = simulate(count)
        self.par = {'cov': cov}
        self.simulate = simulate

        # density function definition

        def density(unifs, vects=None):
            if np.ndim(unifs[0]) == 0:
                unifs = [[i] for i in unifs]
            vecs = []
            marginal_density = []
            for i in range(dim):
                a = stats.norm([0], np.sqrt(cov[i, i]))
                if vects is None:
                    temp = a.ppf(unifs[i])
                else:
                    temp = vects[i]
                marginal_density.append(a.pdf(temp))
                vecs.append(temp)
            res = []
            vecs = list(zip(*vecs))
            marginal_density = list(zip(*marginal_density))
            factor = np.sqrt((2 * math.pi) ** dim * np.linalg.det(cov))
            cov_inv = cov ** (-1)
            for i in list(zip(*[vecs, marginal_density])):
                v = np.matrix(i[0])
                m = np.prod(i[1])
                if m > 0:
                    temp = 1 / (m * factor) * math.exp(-0.5 * v * cov_inv * np.transpose(v))
                else:
                    temp = 0
                res.append(temp)
            return res

        self.pdf = density


### returns a 'student copula' fitted to the given points
class cop_student(cop_mod):
    ACR = 'ST'

    def __init__(self, uniforms, vv=None):
        length, format = ut.good_format(uniforms)
        if not (format):
            raise (RuntimeError('uniforms do not have the good format'))

        dim = len(uniforms)

        # computing correlation
        a = stats.norm([0], [1])
        vects = [a.ppf(coor) for coor in uniforms]
        dim = len(uniforms)
        cor = norm_matrix(np.cov(np.matrix(vects)))

        cor_inv = np.linalg.inv(cor)

        # definition of the log-likelihood
        def log_likelihood(unifs, vv=None):
            if unifs is None:
                unifs = uniforms
            not_list = np.ndim(unifs[0]) == 0
            if not_list:
                unifs = [[i] for i in unifs]

            if vv is None:
                v = self.par['v']
            else:
                v = vv

            sigma = norm_matrix(cor)
            sigma_inv = cor_inv

            # back in the real distribution space...
            vecs = []
            marginal_density = []
            t_obj = stats.t(v)

            for i in range(dim):
                temp = t_obj.ppf(unifs[i])
                marginal_density.append(t_obj.pdf(temp))
                vecs.append(temp)
            vecs = list(zip(*vecs))
            marginal_density = list(zip(*marginal_density))

            # factor of the density in the real space
            factor = math.log(spe.gamma((v + dim) / 2)) \
                     - math.log(spe.gamma(v / 2)) \
                     - 0.5 * dim * math.log((math.pi * v)) \
                     - 0.5 * math.log(np.linalg.det(sigma))

            # computing the log likelihood
            res = []
            for pt in list(zip(*[vecs, marginal_density])):
                x = np.matrix(pt[0])
                try:
                    m = np.sum([math.log(density) for density in pt[1]])
                    temp = -m + factor + math.log(1 + x * sigma_inv * np.transpose(x) / v) * (-(v + dim) / 2)
                except:
                    print('error computing the density of student copula: pdf value is too low')
                    temp = 1
                res.append(temp)

            return sum(res)

        # maximizing log likelihood to find the number of degree of freedom
        print('maximizing log likelihood to find the number of degree of freedom')
        if vv is None:
            opt_result = opt.minimize(lambda x: -(log_likelihood(uniforms, vv=x[0])), [2.1],
                                      bounds=[(2.001, 10)], method='L-BFGS-B')
            if opt_result['success']:
                v = opt_result['x'][0]
            else:
                print('Warning: optimization did not converge')
                v = 3
        else:

            v = vv
        v = max(v, 2.00001)
        s = cor * (v - 2) / v
        par = {'v': v, 'sigma': s}

        # computing t distribution (see " https://en.wikipedia.org/wiki/Multivariate_t-distribution")
        def distribute(nb, vv=v):
            sigma = cor * (vv - 2) / vv
            Y = np.transpose(np.random.multivariate_normal([0 for i in range(dim)], sigma, nb))
            u = np.sqrt(vv / np.random.chisquare(vv, nb))
            temp = []
            for yy in Y:
                temp.append(list(np.array(yy) * u))
            return (ut.uniforms(temp))

        self.dim = dim
        self.length = length
        self.name = 'student'
        self.val = distribute(length)
        self.par = par
        self.simulate = distribute
        self.log_likelihood = log_likelihood

        # density function definition

        def density(unifs):
            if np.ndim(unifs[0]) == 0:
                unifs = [[i] for i in unifs]

            v = self.par['v']
            sigma = norm_matrix(self.par['sigma'])
            std = math.sqrt(v / (v - 2))

            vecs = []
            marginal_density = []
            t_obj = stats.t(v)
            for i in range(dim):
                temp = t_obj.ppf(unifs[i])
                marginal_density.append(t_obj.pdf(temp))
                vecs.append(temp)
            res = []

            vecs = list(zip(*vecs))
            marginal_density = list(zip(*marginal_density))
            factor = spe.gamma((v + self.dim) / 2) / (
            spe.gamma(v / 2) * math.sqrt((math.pi * v) ** self.dim * np.linalg.det(sigma)))
            cov_inv = sigma ** (-1)

            for i in list(zip(*[vecs, marginal_density])):
                x = np.matrix(i[0])
                m = np.prod(i[1])
                if m > 0:
                    temp = 1 / m * factor * math.exp(
                        math.log(1 + x * cov_inv * np.transpose(x) / v) * (-(v + self.dim) / 2))
                else:
                    temp = 0
                res.append(temp)
            return res

        self.pdf = density


### returns a customized copula:
### it will 'squeeze' a given model to fit a few functions of interest, like tail functions,
### diagonal and anti diagonal distributions
class cop_customized(cop_mod):
    ACR = 'CU'

    nb_points_returned = []
    pdf_list = []
    uniform_marginals = True
    K = []  # measuring how much of the points are rejected in customized

    # arguments:
    #   () model: a copula model (subclass of cop_mod) (ex: cop_gaussian -- not cop_gaussian(unifs))
    #   () unifs: uniforms -- empirical copula
    #   () subsetFunctions: a list of subclasses of SubsetModel
    #   () fact: for each computed point, we will compute 'fact' points from model
    #               'fact' must be >1. The greater it is, the more accurate, but the longer computation time too.
    def __init__(self, model, unifs, subsetFunctions=None, fact=5, visualize=False, redistribute=False, testMode=False):

        if subsetFunctions is None:
            subsetFunctions = sub.all_subsets(len(unifs))

        model_created = model(unifs)
        subset_created = []

        def distribute(count, visu=False):
            count_init = fact * count + 2 * fact * int(math.sqrt(count)) + 250
            density_multiplier = len(unifs[0]) * fact / count_init

            vec = model_created.simulate(count_init)
            dim = len(vec)

            # creating subset models and corresponding densities for empirical and model
            print('creating subset models and corresponding densities for empirical and model')
            f_emp = []
            f = []
            inverse = []
            intervals_updated = []
            densities_added_to_model = []
            nb_points_added = 0

            for subset in subsetFunctions:
                sset_emp = subset(unifs)
                sset = subset(vec, model=model)
                if len(subset_created) < len(subsetFunctions):
                    subset_created.append(sset)
                print('\n  # %s #\n' % sset_emp.name)
                f_emp_tp = sset_emp.density()
                f_tp_tp = sset.density()
                f_tp = []

                # scaling density function
                print('scaling density function')

                for func in f_tp_tp:
                    def func2(d, arg_func=func):
                        return fact * arg_func(d)

                    f_tp.append(func2)

                if visu:
                    for func3 in zip(*[f_tp, f_emp_tp]):
                        x = np.arange(0, 1, 0.05)
                        y1 = [func3[0](i) for i in x]
                        y11 = [func3[0](i) / fact for i in x]
                        y2 = [func3[1](i) for i in x]
                        plt.figure()
                        plt.title('densities comparison for %s' % sset_emp.name)
                        plt.plot(x, y1, c='blue', label='model')
                        plt.plot(x, y11, c='green', label='model-normalized')
                        plt.plot(x, y2, c='red', label='empirical')

                # adding points to 'model' points when density is inferior to empirical density
                print('adding points to \'model\' points when density is inferior to empirical density')
                f_tp_updated = []

                for ii in range(len(f_emp_tp)):
                    inters = []
                    added_real_density = []
                    added_subset_density = []
                    nb_points_in_quadrant = sset.quad_points[ii]
                    nb_points_added_in_quadrant = 0

                    for interval in find_under(f_tp[ii], f_emp_tp[ii], mi=0, ma=1):

                        if interval[0] >= interval[1]:
                            continue
                        # M is the number of points to add so that the subset 'densities' of the model are higher
                        # than those of the empirical.
                        area = sset.area(interval)
                        M = int(math.ceil(
                            max(0, -find_min(f_tp[ii], f_emp_tp[ii], mi=interval[0], ma=interval[1])[0]) * area *
                            sset.quad_points[ii]))
                        if M == 0:
                            continue
                        nb_points_added_in_quadrant += M

                        # creating M uniformly distributed points in these area and adding these to the model's points
                        vec_tp = list(zip(*sset.uniform(ii, interval, M)))
                        if vec_tp != []:
                            for j in range(dim):
                                vec[j].extend(list(vec_tp[j]))

                        inters.append(interval)
                        added_subset_density.append(M / ((interval[1] - interval[0]) * nb_points_in_quadrant))
                        added_real_density.append(M / area)

                    nb_points_added += nb_points_added_in_quadrant
                    densities_added_to_model.append(added_real_density)
                    intervals_updated.append(inters)

                    # updating 'model' density to take the previous modification into account
                    print('updating \'model\' density to take the previous modification into account')
                    if len(inters) > 0:
                        def f_tp_modified(d, arg_ii=ii, arg_f_tp=f_tp, arg_inters=inters,
                                          arg_added_density=added_subset_density,
                                          arg_added_points=nb_points_added_in_quadrant,
                                          arg_nb_pt_quadrant=nb_points_in_quadrant):

                            ratio = arg_added_points / (arg_nb_pt_quadrant + arg_added_points)
                            res = arg_f_tp[arg_ii](d) * (1 - ratio)

                            for j in range(len(arg_inters)):
                                if arg_inters[j][0] < d < arg_inters[j][1]:
                                    res += arg_added_density[j] * ratio
                            return res

                        f_tp_updated.append(f_tp_modified)
                    else:
                        f_tp_updated.append(f_tp[ii])

                inverse.append((sset.inverse, len(f)))
                f_emp.extend(f_emp_tp)
                f.extend(f_tp_updated)

            print('total number of points added: %d' % nb_points_added)

            # looking at how the various functions are correlated (kind of dimensionality?)
            # computing correlation ('cosinuses') for the inverse of the points
            print('computing correlation between the arguments of each functions of interest')

            # computing the arguments of each 1-dimensional function of interest
            nb_func = len(f)
            inv_points = []

            def compute_inverse(pt):
                inv_i = [None for i in range(nb_func)]
                for j in inverse:
                    temp = j[0](pt)
                    inv_i[j[1] + temp[0]] = temp[1]
                return inv_i

            model_updated_points = list(zip(*vec))
            for i in model_updated_points:
                inv_i = compute_inverse(i)
                inv_points.append(inv_i)

            # computing cosinus and sinus
            cosinuses = abs(ut.pearson_with_holes(list(zip(*inv_points))))
            sinuses = [[0 for i in range(nb_func)] for j in range(nb_func)]

            for i in range(nb_func):
                for j in range(i + 1, nb_func):
                    temp = math.sqrt(1 - cosinuses[i, j] ** 2)
                    sinuses[i][j] = temp
                    sinuses[j][i] = temp
            sinuses = np.matrix(sinuses)

            # computing density estimates for each points (for both empirical copula and modelled copula)
            print('computing density estimates')

            def compute_density(pt, func):
                indices = []
                nb_temp = 0
                for i in range(nb_func):
                    if pt[i] is not None:
                        nb_temp += 1
                        indices.append(i)

                if nb_temp > 0:
                    res = 0
                    for i in indices:
                        prod = func[i](pt[i])
                        for j in indices:
                            if j != i:
                                prod *= func[j](pt[j]) * sinuses[i, j] + cosinuses[i, j]
                        res += prod
                    res /= nb_temp
                else:
                    res = 1
                return res

            # trying not to exclude too many points: K is the max ratio between empirical 'density' and calculated 'density'
            # the previous addition of points ensures that 0<K<1 but it is most of the time well below 1
            # (especially for high dimensions of subsetfunction)
            def div_safe(a, b):
                if b > 0:
                    return a / b
                else:
                    return 1

            K = 0
            density_ratio = []
            for pts in inv_points:
                ratio = div_safe(compute_density(pts, f_emp), compute_density(pts, f))
                density_ratio.append(ratio)
                K = max(ratio, K)
            self.K.append(K)
            print('K: %.5f' % K)

            # creating the density function if the points were kept with a probability proportional to density_ratio
            def density(x):
                if np.ndim(x[0]) == 0:
                    x = [[i] for i in x]
                res = []
                cste = np.mean(density_ratio) / K  # proportion of kept points
                added_pt_ratio = nb_points_added / (len(model_updated_points))

                for pt in list(zip(*x)):
                    inv = compute_inverse(pt)
                    # computing the pdf value of the 'updated' model (taking point additions into account)
                    res_tp = model_created.pdf(pt)[0] * (1 - added_pt_ratio)
                    incr = 0
                    for j in inv:
                        if j is None:
                            continue
                        for interval in zip(*[intervals_updated[incr], densities_added_to_model[incr]]):
                            if interval[0][0] < j < interval[0][1]:
                                res_tp += interval[1] / nb_points_added * added_pt_ratio
                        incr += 1

                    temp = (compute_density(inv, f_emp), compute_density(inv, f))
                    # probability of being kept in the first rejection algoritm
                    res_tp *= div_safe(temp[0], temp[1]) / (K * cste)
                    res.append(res_tp)

                return res

            def rejection_algo(points, probability, max_probability):
                inv_max_probability = 1 / max_probability
                incr = 0
                sampled_points = []
                points_plus_probability = list(zip(*[points, probability]))
                random.shuffle(points_plus_probability)

                for pt in points_plus_probability:
                    if random.random() < pt[1] * inv_max_probability:
                        sampled_points.append(pt[0])
                        incr += 1
                        if incr >= count:
                            break
                return sampled_points

            ### making sure the density has uniform marginals
            if not redistribute:
                normed = norm_distribution(density, dim)
            else:
                normed = None

            if normed is None:
                self.uniform_marginals = False

                # sampling points with the rejection algorithm
                sampled_points = rejection_algo(model_updated_points, density_ratio, K)
                new_density = density

            else:
                new_density, marginal_multipliers = normed[0], normed[1]
                new_density_ratio = []
                for pts in zip(*[model_updated_points, density_ratio]):
                    ratio = np.prod([marginal_multipliers[i](pts[0][i]) for i in range(dim)]) * pts[1]
                    new_density_ratio.append(ratio)
                new_K = max(new_density_ratio)
                print('new_K: %.5f' % new_K)

                # sampling points with the rejection algorithm
                sampled_points = rejection_algo(model_updated_points, new_density_ratio, new_K)

            # for visualization purpose
            if visu & (dim == 2):
                nb_pt = 20
                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')
                ax.set_title('probability of the points being kept')
                xx = [(i // (nb_pt + 1)) / nb_pt for i in range((nb_pt + 1) ** 2)]
                yy = [(i % (nb_pt + 1)) / nb_pt for i in range((nb_pt + 1) ** 2)]
                xxx = [((i // (nb_pt + 1)) / nb_pt, (i % (nb_pt + 1)) / nb_pt) for i in range((nb_pt + 1) ** 2)]
                zz = []
                for j in xxx:
                    zz.append(div_safe(compute_density(compute_inverse(j), f_emp),
                                       (compute_density(compute_inverse(j), f) * K)))
                ax.plot_trisurf(xx, yy, zz)

                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')
                ax.set_title('density emp')
                zz = []
                for j in xxx:
                    zz.append(compute_density(compute_inverse(j), f_emp))
                ax.plot_trisurf(xx, yy, zz)

                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')
                ax.set_title('density model')
                zz = []
                for i in xxx:
                    zz.append(compute_density(compute_inverse(j), f))
                ax.plot_trisurf(xx, yy, zz)

            if visu:
                print('cosinuses: \n%r' % cosinuses)
                print('sinuses: \n%r' % sinuses)

            return sampled_points, new_density

        # looping until enough points are generated
        def distribute2(count, visualize=False):

            print('\n### distribute ###\n\n')

            total_length = 0
            length_per_call = []
            it = 1
            max_iter = 20
            result = []
            self.uniform_marginals = True

            while (total_length < count) & (it < max_iter):

                # trying to create the right number of points:
                if it == 1:
                    count_bis = max(300, count)
                else:
                    count_bis = max(int(1.2 * (count - total_length) * max(300, count) / length_per_call[0]), 300)

                # creating the points using distribute
                points, density = distribute(count_bis, visu=visualize)
                length = len(points)

                self.pdf_list.append(density)
                self.nb_points_returned.append(min(length, count - total_length))

                length_per_call.append(length)

                if total_length + length < count:
                    result.extend(points)
                    total_length += length
                else:
                    result.extend(points[:(count - total_length)])
                    total_length = count
                it += 1

            it -= 1
            print('number of calls: %d  (nb pts: %r)\n\n' % (it, length_per_call))

            # creating the final pdf
            if it == 1:
                self.pdf = self.pdf_list[0]
            else:
                coeffs = [i / count for i in self.nb_points_returned]

                def final_pdf(x):
                    pdf_returned = [l(x) for l in self.pdf_list]
                    res_tp = []
                    for den in zip(*pdf_returned):
                        res_tp.append(np.sum([den[i] * coeffs[i] for i in range(it)]))
                    return res_tp

                self.pdf = final_pdf

            # transforming from a list of points to a list of coordinate lists
            result = list(zip(*result))
            for i in range(len(result)):
                result[i] = list(result[i])

            name = model_created.name + '_customized'
            for sset in subset_created:
                name += '_' + sset.acr
            # if the 'norm_distribution' function returned 'None', redistributing the points according to their ranks
            # to get a copula (with uniform marginals)
            if not self.uniform_marginals:
                cop_red = cop_emp_redistributed(result, self)
                result = cop_red.val
                self.pdf = cop_red.pdf
                name += '_not_normed_redistributed'
                self.name = name
            else:
                self.name = name

            return (result)

        self.dim = len(unifs)
        self.length = len(unifs[0])

        self.simulate = distribute2
        self.val = distribute2(self.length, visualize=visualize)


### return a uniform copula for comparison purposes
class cop_uniform(cop_mod):
    ACR = 'UN'

    def __init__(self, unifs):
        if unifs is None:
            unifs = [[0.5], [0.5]]
        self.dim = len(unifs)
        self.length = len(unifs[0])

    def simulate(self, count):
        return np.random.uniform(0, 1, (self.dim, count))

    def pdf(self, x):
        if len(np.shape(x)) == 1:
            x = [x]
        return [1 for i in range(len(x[0]))]


### returns the independent combination of various copulae
# arguments:
#   () models: a list of copula models to be used on subsets of the data (ex: [cop_gaussian,cop_student])
#   () separators: a list of indexes specifying the different subsets (ex: [2] for the precedent example)
class cop_independent_combination(cop_mod):
    ACR = 'IC'

    def __init__(self, uniforms, models, separators):

        self.dim = dim = len(uniforms)
        self.length = length = len(uniforms[0])
        self.nb_copulae = nb_copulae = len(separators) + 1

        sep = [0]
        sep.extend(separators)
        sep.append(dim)
        self.sep = sep

        if min([sep[i + 1] - sep[i] for i in range(nb_copulae)]) < 1:
            raise (RuntimeError('A copula doesn\'t get data'))

        if len(models) != nb_copulae:
            raise (RuntimeError('length of \'models\' and \'separators\' don\'t correspond'))

        copulae = []

        for i in range(nb_copulae):
            if sep[i + 1] - sep[i] == 1:
                copulae.append(cop_uniform(uniforms[sep[i]:sep[i + 1]]))
            else:
                copulae.append(models[i](uniforms[sep[i]:sep[i + 1]]))

        name = 'combination'
        for cop in copulae:
            name = '%s_%s' % (name, cop.ACR)

        self.name = name
        self.copulae = copulae

    def pdf(self, unifs):
        sep = self.sep
        densities = [self.copulae[i].pdf(unifs[sep[i]:sep[i + 1]]) for i in range(self.nb_copulae)]
        return [float(np.prod(i)) for i in zip(*densities)]

    def simulate(self, count):
        sims = []
        for cop in self.copulae:
            sims.extend(cop.simulate(count))
        return sims


# -----------------------------------------------------------------------------------------------------------------------





# returns the kendall matrix associated with these vectors
def kendall(vects):
    length, format = ut.good_format(vects)
    if not format:
        return None

    def sign(x):
        if (x > 0):
            return 1
        elif x < 0:
            return -1
        else:
            return 0

    dim = len(vects)
    ken = np.identity(dim)
    for v1 in zip(*vects):
        for v2 in zip(*vects):
            for j in range(dim):
                for k in range(j + 1, dim):
                    ken[j, k] += sign((v2[j] - v1[j]) * (v2[k] - v1[k]))
    for j in range(dim):
        for k in range(j + 1, dim):
            ken[j, k] = (ken[j, k]) / (length * (length - 1))
            ken[k, j] = ken[j, k]
    return ken


# follows the gradient function 'df' to find the coordinate ('par') corresponding to the highest value of a corresponding "f"
def grad_ascent(df, par, step=1., prec=0.02):
    print('\ngrad ascent to find the parameters')
    if np.ndim(par) == 0:
        par = np.array([par])
    elif not type(par) == np.array:
        par = np.array(par)
    max_iter = 100
    incr = 0
    last_val = None

    while (step > prec) & (incr < max_iter):
        df_val = df(par)
        norm = np.sqrt(sum(df_val ** 2))
        if (norm > 0):
            df_val = df_val * (step / norm)
        if last_val is not None:
            if sum(df_val * last_val) <= 0:
                step /= 1.414
        par = par + df_val
        incr += 1
        last_val = df_val
        # print('df_val: %r'%df_val)
        # print('par: %r'%par.tolist())

    if (incr == max_iter):
        print('Warning: reached maximum number of iterations in grad_descent')
    return par


# assuming f1(0)>f2(0) both increasing function on [0,0.5], finds a point just after both function cross for the last time
def find_crossing(f1, f2, precision=0.0002):
    curr = 0.5
    step = 0.02
    incr = 0
    while (step > precision) & (incr < 50):
        incr += 1
        if (curr - step) <= 0:
            return 0
        if f1(curr - step) > f2(curr - step):
            curr -= step
        else:
            step /= 2
    return (curr)


# find crossing of "f" (-f2 when "f2" is given) and "level" in the range("mi","ma"), with precision = "precision"
# It migth not find changes of sign within distance smaller than "step_wid" (>0)
def find_all_crossing(f, f2=None, level=0, mi=0, ma=1, precision=0.0001, step_wid=0.01):
    curr = mi
    if f2 is not None:
        def temp(x):
            return f(x) - f2(x)
    else:
        temp = f
    res = []
    sign = 1
    if temp(mi) < level:
        sign = -1
    if temp(mi) == level:
        step = step_wid
        for k in range(10):
            val = temp(mi + step) - level
            if val != 0:
                if val < 0:
                    sign = -1
                break
            step /= 2

    while curr < ma - precision:
        incr = 0
        step = 0.01
        while (step > precision) & (incr < 50):
            incr += 1
            if (curr + step) >= ma:
                break
            if sign * temp(curr + step) > level:
                curr += step
            else:
                step /= 2
        if (curr + step) < ma:
            res.append(curr)
        curr += step
        sign *= -1
    return (res)


# find intervalls where "f" falls under "level"
def find_under(f, f2=None, level=0, mi=0, ma=1, precision=0.0001, step_wid=0.01):
    if f2 is not None:
        def f_tp(x):
            return f(x) - f2(x)
    else:
        f_tp = f
    crossing = find_all_crossing(f, f2=f2, level=level, mi=mi, ma=ma, precision=precision, step_wid=step_wid)
    if crossing[-1] != ma:
        crossing.append(ma)
    if crossing[0] != mi:
        temp = [mi]
        temp.extend(crossing)
        crossing = temp
    res = []
    for i in range(len(crossing) - 1):
        if f_tp((crossing[i] + crossing[i + 1]) / 2) < level:
            res.append((crossing[i], crossing[i + 1]))
    return res


# finds the minimum value of f on a grid with precision 0.01
def find_min(f, f2=None, precision=100, mi=0, ma=1):
    if f2 is not None:
        def temp(x):
            return f(x) - f2(x)
    else:
        temp = f
    incr = mi
    step = (ma - mi) / precision
    res = [temp(mi) + 0.01, None]
    for i in range(precision):
        val = temp(incr)
        if val < res[0]:
            res = [val, incr]
        incr += step
    return res


# normalizes a correlation matrix
def norm_matrix(m):
    format = False
    dim = 0
    if m is not None:
        dim = len(m)
        if dim > 0:
            if m.shape == (dim, dim):
                format = True
    if not format:
        print('norm_matrix requires a square matrix as argument')
        return None

    allPos = True
    for i in range(dim):
        for j in range(dim):
            allPos &= float(m[i, j]) >= 0
            if not allPos:
                break

    mbis = None
    if not allPos:
        eig, rmat = np.linalg.eigh(m)
        rmat = np.matrix(rmat)
        temp = np.matrix(np.identity(dim))
        for i in range(dim):
            if eig[i] > 0:
                temp[i, i] = eig[i]
            else:
                temp[i, i] = 0.001
        mbis = (rmat ** (-1)) * temp * rmat
    else:
        mbis = m.copy()

    d = np.diag(mbis).tolist()
    for i in range(dim):
        for j in range(dim):
            mbis[i, j] /= np.sqrt(abs(d[i] * d[j]))

    return np.matrix(mbis)


### This function is primarily written to be used in norm_distribution
### arguments:
#   () dim is the dimension of the distribution space
#   () nb is the number of subdivisions per axis
#   () val are cdf values from a multivariate distribution at regular intervals (list of length nb^dim)
#   () x is a set of marginals multiplier, (list of length nb*dim), starting with first axis multipliers, then second...
### returns the discrete marginals values of the multivariate distribution obtained by multiplying the marginals
### by the marginals' multipliers 'x'
def discrete_marginals_multiplied(val, x, dim, nb):
    # compute the updated discrete marginal of the first axis
    def sum_except_first_axis(val_tp, x_tp, dim_tp):
        print(dim_tp)
        if dim_tp == 1:
            return val_tp
        else:
            incr = [0, 0]
            res_tp = []
            temp = [0 for i in range(nb)]
            for i in val_tp:
                temp[incr[0]] = temp[incr[0]] + i * x_tp[incr[0]]
                incr[0] = incr[0] + 1
                if incr[0] == nb:
                    incr[1] += 1
                    incr[0] = 0
                if incr[1] == nb:
                    incr[1] = 0
                    # print('temp: %r'%temp)
                    res_tp.extend(temp.copy())
                    temp = [0 for i in range(nb)]
            return sum_except_first_axis(res_tp, x_tp[nb:], dim_tp - 1)

    # this recursive function will compute the sum of the numbers on the hyperplan x1=1,
    # multiplied by the discrete marginals specified in x_tp (...except that corresponding to x1, the followed axis)
    # and then x1=2,x1=3...x2=1,.. xdim=1,..xdim=nb
    def sum_along_tree(val_tp, x_tp, dim_tp):
        print('length val: %d, length x : %d, dim: %d' % (len(val_tp), len(x_tp), dim_tp))
        if dim_tp == 1:
            return val_tp
        else:
            incr = 0
            res_tp = []
            temp = 0
            for i in val_tp:
                temp = temp + i * x_tp[incr]
                incr = incr + 1
                if incr == nb:
                    incr = 0
                    res_tp.append(temp)
                    temp = 0

            res = sum_except_first_axis(val_tp, x_tp[nb:], dim_tp)
            res.extend(sum_along_tree(res_tp, x_tp[nb:], dim_tp - 1))
            return res

    final_res = sum_along_tree(val, x.copy(), dim)
    final_res = [i[0] * i[1] for i in zip(*[final_res, x])]
    return final_res


### given an approximate copula density with non-uniform marginals,
### (ie a distribution pdf on [0,1]^n with marginals that vary less than an order of magnitude)
### This function computes a new distribution with uniform marginals by multiplying the marginals with 1d-functions
#
### Arguments:
#   () pdf: the approximate copula density
#   () dim: the dimension of the space
### Returns:
#   () None if unable to compute
#   () (new_pdf,marginals_multiplier) if everything was fine, marginals_multiplier being a list of functions [0,1]->R+
def norm_distribution(pdf, dim):
    if dim >= 10:
        print('dimensions above 10 requires too much computation power')
        return None

    if dim < 4:
        nb = 20
    else:
        nb = max(5, int(10000 ** (1 / dim)))

    xx = [(i + 0.5) / nb for i in range(nb)]
    increment = 1 / nb

    # pdf_val contains values of the pdf on a regular grid
    pdf_val = []
    current = [0.5 / nb for i in range(dim)]
    while current is not None:
        pdf_val.append(pdf(current)[0])
        current = ut.table_increment(current, mi=increment / 2, ma=1, incr=increment, length=dim)

    # f returns:
    # 'objective': the sum of pdf values multiplied by the marginals multiplier 'x'
    #              across the hyperplans X1=xx[0], X1=xx[1]...Xdim=xx[dim*nb],.. Xdim=xx[dim*nb + nb-1]
    #              minus 'sum_to' (what uniform marginals should give)
    # 'jacobian':  the jacobian of objective
    #
    sum_to = nb ** (dim - 1)

    def f(x):
        # sum_returned=discrete_marginals_multiplied(pdf_val,x,dim,nb)
        resulting_marginals = np.zeros(nb * dim)
        jacobian = np.matrix(np.zeros((nb * dim, nb * dim)))

        current = [0 for i in range(dim)]
        multipliers = np.ones(dim)
        for pt in pdf_val:

            for i in range(dim):
                multipliers[i] = x[i * nb + current[i]]
            pt = pt * np.prod(multipliers)

            for d in range(dim):
                for e in range(dim):
                    jacobian[e * nb + current[e], d * nb + current[d]] = jacobian[e * nb + current[e], d * nb + current[
                        d]] + pt

            current = ut.table_increment(current, mi=0, ma=nb, incr=1, length=dim)

        for i in range(nb * dim):
            resulting_marginals[i] = jacobian[i, i]
            for j in range(nb * dim):
                jacobian[i, j] = jacobian[i, j] / x[i]

        objective = [i - sum_to for i in resulting_marginals]

        return objective, jacobian

    # using scipy optimizer to find the root of f
    start = np.ones(nb * dim)
    # print('start=%r'%(start))
    marginals_multiplier_values = root(f, start, jac=True, tol=0.0001)

    # returns 'None' if solution didn't converge
    if not marginals_multiplier_values['success']:
        print('solution did not converge in \'norm_distribution\'')
        print(marginals_multiplier_values['message'])
        return None
    else:
        marginals_multiplier_values = marginals_multiplier_values['x']

    marginals_multiplier = []
    for i in range(dim):
        yy = list(marginals_multiplier_values[nb * i:nb * (i + 1)])
        marginals_multiplier.append(ut.create_spline(ut.find_spline(xx, yy, positiveness_constraint=True)))

    # multiplies the previous pdf with the computed marginals
    def new_pdf(pts):
        if np.ndim(pts[0]) == 0:
            pts = [[i] for i in pts]

        prod = [np.prod([marginals_multiplier[i](point[i]) for i in range(dim)]) for point in zip(*pts)]
        res = [i[0] * i[1] for i in zip(*[prod, pdf(pts)])]

        return res

    return new_pdf, marginals_multiplier
