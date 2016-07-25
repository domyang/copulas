import matplotlib.pyplot as plt
import math
import time
import random
import numpy as np
from scipy import stats
import scipy.special as spe
import scipy.optimize as optimize
from scipy.integrate import quad
import utilities as ut

# -------------------------------------- 2 types of vines ------------------------------------------------------------
epsilon = 0.000001


class conditional(object):
    dim = 1
    length = 0
    name = ''
    val = []
    cop_mod = None
    copulae = [[]]

    mode = 'conditional'

    def update_parameter(self, par):
        pass

    def pdf(self, unifs):
        return [1]

    def log_likelihood(self, par, arg_unifs=[[]]):
        return 0

    def simulate(self, count):
        return []

    def print_par(self):
        for i in self.copulae:
            for j in i:
                print(j.par)

    def print_names(self):
        for i in self.copulae:
            for j in i:
                print(j.name)

    def pprint(self):
        pass

    def conditionnalPDF(self, obs):
        def pdf(x):
            if 0 < x < 1:
                return 1
            else:
                return 0

        return pdf


class C_vine(conditional):
    name = 'canonical vine'
    ACR = 'CV'

    def __init__(self, uniforms, list_models=None, add_name=''):

        # arranging the coordinates for the canonical-vine decomposition: n-1,n-2,..,1,n so that the unknown Xn
        # will be estimated using Xn-1 first, then Xn-2 ...etc
        unifs = uniforms[:-1]
        unifs.reverse()
        unifs.append(uniforms[-1])

        dim = len(unifs)
        if list_models is None:
            list_models = [[cop2d_gaussian, cop2d_student, cop2d_clayton, cop2d_frank, cop2d_gumbel] for i in
                           range(dim - 1)]
        elif not type(list_models[0]) in {tuple, list, np.array}:
            list_models = [list_models for i in range(dim - 1)]

        copulae = []
        obs = unifs

        # creating the vine copula and selecting the models for each pair
        print('creating the vine copula and selecting the models for each pair')
        for step in range(dim - 1):

            print('step: %d' % step)

            copulae_for_step = []
            sum_lld = 0

            # selecting and fitting the best copula model for each pair of observations
            for pair in range(1, dim - step):

                best_copula = None
                best_lld = 0

                for cop in list_models[step]:
                    try:
                        copula_tp = cop([obs[0], obs[pair]])
                        lld_tp = copula_tp.log_likelihood([obs[0], obs[pair]])

                        if best_copula is None:
                            best_lld = lld_tp
                            best_copula = copula_tp
                        else:

                            if lld_tp > best_lld:
                                best_copula = copula_tp
                                best_lld = lld_tp

                    except:
                        print("could not estimate log likelihood of %s" % cop.name)

                sum_lld += best_lld
                copulae_for_step.append(best_copula)

            # print('sum lld: %f' % sum_lld)
            copulae.append(copulae_for_step)

            # computing the observations for the next step

            obs_for_next_step = []
            for pair in range(1, dim - step):
                obs_for_next_step.append(copulae_for_step[pair - 1].conditionalCDF([obs[0], obs[pair]], 0))

            obs = obs_for_next_step

        self.dim = dim
        self.copulae = copulae
        self.uniforms = uniforms
        self.name += add_name

        # defining the log likelihood function for our problem

        par_start = []
        for step in range(dim - 1):
            for pair in range(1, dim - step):
                par_tp = copulae[step][pair - 1].parameter_list()
                par_start.extend(par_tp)
        self.par_start = par_start

        # function to update the parameters of the copula

    def update_parameter(self, par):
        cur = 0
        for step in range(self.dim - 1):
            for pair in range(1, self.dim - step):
                length_par = len(self.copulae[step][pair - 1].parameter_list())
                self.copulae[step][pair - 1].assign_parameter(par[cur:cur + length_par])
                cur = cur + length_par


                # returns the log-likelihood

    def log_likelihood(self, par, arg_unif=None, restricted=False):

        if arg_unif is None:
            arg_unif = self.uniforms

        not_list = np.ndim(arg_unif[0]) == 0
        if not_list:
            arg_unif = [[i] for i in arg_unif]

        unifs = arg_unif[:-1]
        unifs.reverse()
        unifs.append(arg_unif[-1])

        if par is not None:
            self.update_parameter(par)

        obs = unifs
        res = 0

        # computing the log-likelihood
        for step in range(self.dim - 1):
            # print('step: %d, \no1: %r , \no2:%r'%(step,obs[0][:10],obs[-1][:10]))
            if restricted:
                res = res + self.copulae[step][-1].log_likelihood([obs[0], obs[-1]])
            else:
                for pair in range(1, self.dim - step):
                    # print(' we add: %r'%(self.copulae[step][pair-1].log_likelihood([obs[0],obs[pair]])))
                    res = res + self.copulae[step][pair - 1].log_likelihood([obs[0], obs[pair]])

            obs_for_next_step = []
            for pair in range(1, self.dim - step):
                obs_for_next_step.append(self.copulae[step][pair - 1].conditionalCDF([obs[0], obs[pair]], 0))
            obs = obs_for_next_step

        return res

    # returns the density
    def pdf(self, unifs):
        not_list = np.ndim(unifs[0]) == 0
        if not_list:
            unifs = [[i] for i in unifs]

        res = list(map(math.exp, [self.log_likelihood(None, arg_unif=[[i] for i in x]) for x in zip(*unifs)]))

        return res

    # print('maximizing the log-likelihood')
    # par=optimize.minimize(lambda x: -log_likelihood(x),par_start.copy())

    def conditionalPDF(self, observed):
        if (min(observed) <= 0) | (max(observed) >= 1):
            raise (RuntimeError('\'observed\' must be in the unit cube'))

        arrange_observed = observed
        arrange_observed.reverse()

        def pdf(x):

            not_list = np.ndim(x) == 0
            if not_list:
                x = [x]

            if (max(x) >= 1) | (min(x) <= 0):
                raise (RuntimeError('\'observed\' must be in the unit cube'))

            arguments = []
            for xx in x:
                temp = arrange_observed.copy()
                temp.append(xx)
                arguments.append(temp)

            return [math.exp(-self.log_likelihood(None, unifs=arg)) for arg in arguments]

        return pdf

    def partial(self, unifs, a, b):
        if (b >= a) | (not (type(b) == type(a) == int)) | b < -1:
            raise (RuntimeError('Wrong argument in partial'))
        if b == -1:
            return unifs[a]
        else:
            u1 = self.partial(unifs, a, b - 1)
            u2 = self.partial(unifs, b, b - 1)
            return self.copulae[b][a - b - 1].conditionalCDF([u1, u2], coor=1)

    def simulate(self, nb):
        unifs = np.random.uniform(0, 1, (self.dim, nb))

        res = [unifs[0]]
        for i in range(1, self.dim):
            temp = [self.partial(res, j, j - 1) for j in range(i)]

            def cdf(x, pt, temp=temp):
                res_tp = [x]
                for j in range(i):
                    res_tp = self.copulae[j][i - j - 1].conditionalCDF([res_tp, [temp[j][pt]]], coor=1)
                return res_tp[0]

            res_tp = [optimize.brentq(lambda x: cdf(x, pt) - unifs[i][pt], epsilon, 1 - epsilon, xtol=epsilon / 2) for
                      pt in
                      range(nb)]
            res.append(res_tp)

        return res


class D_vine(conditional):
    name = 'D_vine'
    ACR = 'DV'

    def __init__(self, uniforms, list_models=None, rearrange=False, add_name=''):

        # arranging the coordinates for the canonical-vine decomposition: n-1,n-2,..,1,n so that the unknown Xn
        # will be estimated using Xn-1 first, then Xn-2 ...etc
        unifs = uniforms

        dim = len(unifs)

        if list_models is None:
            list_models = [[cop2d_gaussian, cop2d_student, cop2d_clayton, cop2d_frank, cop2d_gumbel] for i in
                           range(dim - 1)]
        elif np.ndim(list_models[0]) == 0:
            list_models = [list_models for i in range(dim - 1)]

        copulae = []
        obs1 = unifs[:-1]
        obs2 = unifs[1:]

        self.dump = []

        # creating the vine copula and selecting the models for each pair
        print('creating the vine copula and selecting the models for each pair')
        for step in range(dim - 1):

            print('step: %r' % step)

            copulae_for_step = []

            # selecting and fitting the best copula model for each pair of observations
            for pair in range(dim - step - 1):

                pair_obs = [obs1[pair], obs2[pair]]
                best_copula = max(list_models[step], key=lambda cop: cop(pair_obs).log_likelihood(pair_obs))
                best_lld = 0

                # print(list_models)
                # for cop in list_models[step]:
                #     # try:
                #     copula_tp = cop([obs1[pair], obs2[pair]])
                #     lld_tp = copula_tp.log_likelihood([obs1[pair], obs2[pair]]) - len(copula_tp.parameter_list())
                #     # print('\n### step: %d copula: %s, log: %r, par: %r\n'%(step,copula_tp.name,lld_tp,copula_tp.par))
                #
                #     if best_copula is None:
                #         best_lld = lld_tp
                #         best_copula = copula_tp
                #     else:
                #         if lld_tp > best_lld:
                #             best_copula = copula_tp
                #             best_lld = lld_tp
                #             # except:
                #             #     print([obs1[pair], obs2[pair]])
                #             #     print("could not estimate log likelihood of %dth copula" % pair)

                copulae_for_step.append(best_copula)

            copulae.append(copulae_for_step)

            # computing the observations for the next step

            obs_for_next_step_1 = []
            obs_for_next_step_2 = []

            for pair in range(0, dim - step - 2):
                obs_for_next_step_1.append(copulae[step][pair].conditionalCDF([obs1[pair], obs2[pair]], 1))

            for pair in range(1, dim - step - 1):
                obs_for_next_step_2.append(copulae[step][pair].conditionalCDF([obs1[pair], obs2[pair]], 0))

            if rearrange:
                obs1 = []
                obs2 = []
                for obs_for_copula in zip(*[obs_for_next_step_1, obs_for_next_step_2]):
                    unif_tp = ut.uniforms(obs_for_copula, rand=False)
                    obs1.append(unif_tp[0])
                    obs2.append(unif_tp[1])
            else:
                obs1 = obs_for_next_step_1
                obs2 = obs_for_next_step_2

            # print('obs: %r %r'%(obs1,obs2))

            self.dump.append((obs1, obs2))

        name = 'D_vine'
        for i in copulae:
            for j in i:
                name = '%s_%s' % (name, j.ACR)
            name += '|'
        name += add_name

        self.dim = dim
        self.copulae = copulae
        self.rearrange = rearrange
        self.uniforms = uniforms
        self.name = name

        # defining the log likelihood function for our problem

        par_start = []
        for step in range(dim - 1):
            for pair in range(1, dim - step):
                par_tp = copulae[step][pair - 1].parameter_list()
                par_start.extend(par_tp)
        self.par_start = par_start

    # function to update the parameters of the copula
    def update_parameter(self, par):
        cur = 0
        for step in range(self.dim - 1):
            for pair in range(1, self.dim - step):
                length_par = len(self.copulae[step][pair - 1].parameter_list())
                self.copulae[step][pair - 1].assign_parameter(par[cur:cur + length_par])
                cur = cur + length_par

    # returns the log-likelihood
    def log_likelihood(self, par, arg_unif=None):

        if arg_unif is None:
            arg_unif = self.uniforms

        not_list = np.ndim(arg_unif[0]) == 0
        if not_list:
            arg_unif = [[i] for i in arg_unif]

        if par is not None:
            self.update_parameter(par)

        obs1 = arg_unif[:-1]
        obs2 = arg_unif[1:]
        res = 0

        # computing the log-likelihood
        for step in range(self.dim - 1):
            for pair in range(self.dim - step - 1):
                res = res + self.copulae[step][pair].log_likelihood([obs1[pair], obs2[pair]])
            # creating the conditional values of the CDF
            obs_for_next_step_1 = []
            obs_for_next_step_2 = []

            for pair in range(0, self.dim - step - 2):
                obs_for_next_step_1.append(self.copulae[step][pair].conditionalCDF([obs1[pair], obs2[pair]], 1))

            for pair in range(1, self.dim - step - 1):
                obs_for_next_step_2.append(self.copulae[step][pair].conditionalCDF([obs1[pair], obs2[pair]], 0))

            if self.rearrange:
                obs1 = []
                obs2 = []
                for obs_for_copula in zip(*[obs_for_next_step_1, obs_for_next_step_2]):
                    unif_tp = ut.uniforms(obs_for_copula, rand=False)
                    obs1.append(unif_tp[0])
                    obs2.append(unif_tp[1])
            else:
                obs1 = obs_for_next_step_1
                obs2 = obs_for_next_step_2

        return res

    # returns the density
    def pdf(self, unifs):
        not_list = np.ndim(unifs[0]) == 0
        if not_list:
            unifs = [[i] for i in unifs]

        res = list(map(math.exp, [self.log_likelihood(None, arg_unif=[[i] for i in x]) for x in zip(*unifs)]))

        return res

    def partial(self, unifs, first, last, which=0):
        height = last - first
        if not ((type(first) == type(last) == int) & (first <= last)):
            raise (RuntimeError('wrong argument in partial'))
        if first == last:
            return unifs[first]
        else:
            a = self.partial(unifs, first, last - 1, which=1)
            b = self.partial(unifs, first + 1, last, which=0)
            return self.copulae[height - 1][first].conditionalCDF([a, b], coor=which)

    def simulate(self, nb):
        unifs = [list(i) for i in np.random.uniform(0, 1, (self.dim, nb))]

        res = [unifs[0]]
        for i in range(1, self.dim):
            temp = [self.partial(res, j, i - 1, which=1) for j in range(i)]

            def cdf(x, temp=temp):
                res_tp = x
                for j in range(i):
                    res_tp = self.copulae[j][i - j - 1].conditionalCDF([temp[i - j - 1], res_tp], coor=0)
                return res_tp

            res_tp = ut.simultaneous_quantile(cdf, unifs[i])
            res.append(res_tp)

        return res


# ------------------------------------- defining 2d-copulae -----------------------------------------------------------


### upper class for copula models
#   () val: Points of this copula
#   () simulate:   A function (of 'count') returning 'count' points distributed according to the copula
#   () par: The parameters of the copula
#   () pdf: The pdf of the copula
class cop2d(object):
    length = 0
    name = ''
    ACR = ''
    par = {}
    lld = 0
    nbPar = 0

    # simulates points from ths copula
    def simulate(self, count):
        return [[]]

    # pdf of th ecopula
    def pdf(self, x):
        return 1

    ### returns the conditional CDF with respect to the 'coor' coordinate:
    # arguments:
    #   () unifs: a list of two coordinate lists
    #   () coor: 0 or 1: we will compute the CDF knowing unifs[coor] of unif[other]
    #   P(X[coor] <= c[coor] | X[other] = c[other])
    def conditionalCDF(self, unifs, coor):
        return [1]

    ### returns the sum of the log-likelihood of unifs with the parameters theta
    def log_likelihood(self, unifs, theta=None):
        return 0

    ## returns a list of the parameters
    def parameter_list(self):
        return []

    ## modifies the parameters of the copula
    def assign_parameter(self, theta):
        pass

    ## print method for 2d copulae
    def pprint(self):
        s = '\n### Copula Model: ' + self.name + ' ###\n\n'
        s += 'parameters: %r' % self.par
        print(s)

    def plot_points(self, nb=1800):
        v = self.simulate(nb)
        plt.figure()
        plt.plot(v[0], v[1], '.')
        plt.title('%s copula' % self.name)

    def plot_pdf(self):
        ut.curve3d(lambda x: self.pdf([[x[0]], [x[1]]])[0], title='%s copula' % self.name)


# -------------------------------------- 2d copulae models ------------------------------------------------------------

## returns a 'gaussian copula' fitted to the given points
class cop2d_gaussian2(cop2d):
    name = 'gaussian'
    ACR = 'GA'

    def __init__(self, uniforms):

        # Inference
        if uniforms is None:
            uniforms = [[0.5, 0.6], [0.5, 0.4]]

        if not (ut.good_format(uniforms)[1]):
            print(uniforms)
            raise (RuntimeError('uniforms should be a list of same length list'))

        a = stats.norm([0], [1])
        vects = [[a.ppf(j)[0] for j in coor] for coor in uniforms]

        cov = norm_matrix(np.cov(np.matrix(vects)))

        self.par['cov'] = cov

        self.uniforms = uniforms
        self.vects = vects
        self.length = len(uniforms[0])
        self.name = 'gaussian'
        self.nbPar = 1
        self.bounds = bounds = [(epsilon - 1, 1 - epsilon)]

        opt_result = optimize.minimize_scalar(lambda x: -self.log_likelihood(uniforms, x),
                                              bounds=bounds[0], method='bounded')
        if opt_result['success']:
            # print('old parameter: %f, new: %f ' % (self.par['cov'][0, 1], opt_result['x']))
            self.par['corr'] = opt_result['x']
        else:
            print('Maximization of log likelihood was unsuccessful in cop2d_gaussian')

        # defining the conditional cdf
        self.conditionalCDF = self.createConditionalCDF()

    def pdf(self, unifs):
        """Computes pdf for a list of lists"""
        if np.ndim(unifs) < 2:
            unifs = [[i] for i in unifs]
        points = zip(*unifs)
        return [self.single_pdf(u, v, self.par['corr']) for (u, v) in points]

    # Simulation
    def simulate(self, x):
        redistributed = np.transpose(np.random.multivariate_normal([0] * 2, self.par['cov'], x))
        return ut.uniforms(redistributed.tolist())

    # density function
    def single_pdf(self, u, v, corr=None):

        if corr is None:
            corr = self.par['corr']

        first = 1 / np.sqrt(1 - corr ** 2)
        x1, x2, = stats.norm.ppf(u), stats.norm.ppf(v)
        arg = -(corr ** 2 * (x1 ** 2 + x2 ** 2) - 2 * corr * x1 * x2) / (2 * (1 - corr ** 2))
        return first * np.exp(arg)

    # partial CDF: dC/dx
    # coor specifies along which coordinate (x0 or x1) it is differentiated
    def createConditionalCDF(self):
        cov = self.par['cov']

        new_var = (1 - cov[0, 1] ** 2)
        self.norm_obj0 = stats.norm(0, 1)
        self.norm_obj1 = stats.norm(0, np.sqrt(new_var))
        self.last_par_1 = cov[0, 1]

        def conditionalCDF(unifs, coor):

            if not self.last_par_1 == self.par['cov'][0, 1]:
                new_var = (1 - cov[0, 1] ** 2)
                self.norm_obj1 = ut.fast_stats(stats.norm([0], np.sqrt(new_var)))
                self.last_par_1 = cov[0, 1]

            not_list = np.ndim(unifs) < 2
            if not_list:
                unifs = [[i] for i in unifs]

            if (coor != 0) & (coor != 1):
                raise (RuntimeError('coor specifies the coordinate: it must be either 0 or 1'))
            if coor == 1:
                other = 0
            else:
                other = 1

            # the corresponding points in the 'real space' (ie. where the are distributed according to a normal law))
            vecs = [self.norm_obj0.ppf(unifs[i]) for i in range(2)]
            temp = [i[other] - cov[0, 1] * i[coor] for i in zip(*vecs)]
            res = self.norm_obj1.cdf(temp)
            return res


        return conditionalCDF

    # log-likelihood
    def log_likelihood(self, unifs, corr=None):
        if corr is None:
            corr = self.par['corr']
        total = 0
        for (u, v) in zip(*unifs):
            first = (-1/2) * np.log(1 - corr ** 2)
            x1 = stats.norm.ppf(u)
            x2 = stats.norm.ppf(v)
            second = - (corr ** 2 * (x1 ** 2 + x2 ** 2) - 2 * corr * x1 * x2) / (2 * (1 - corr ** 2))
            total += first + second
        return total

    def parameter_list(self):
        return [self.par['cov'][0, 1]]

    def assign_parameter(self, theta):
        self.par['cov'] = np.matrix([[1, theta[0]], [theta[0], 1]])

### returns a 'gaussian copula' fitted to the given points
class cop2d_gaussian(cop2d):
    name = 'gaussian'
    ACR = 'GA'

    def __init__(self, uniforms):

        # Inference
        if uniforms is None:
            uniforms = [[0.5, 0.6], [0.5, 0.4]]

        if not (ut.good_format(uniforms)[1]):
            print(uniforms)
            raise (RuntimeError('uniforms should be a list of same length list'))

        a = stats.norm([0], [1])
        vects = [[a.ppf(j)[0] for j in coor] for coor in uniforms]

        cov = norm_matrix(np.cov(np.matrix(vects)))

        self.uniforms = uniforms
        self.vects = vects
        self.length = len(uniforms[0])
        self.name = 'gaussian'
        self.nbPar = 1
        self.par = {'cov': cov}
        self.bounds = bounds = [(epsilon - 1, 1 - epsilon)]

        opt_result = optimize.minimize_scalar(self.obs_likelihood, bounds=bounds[0], method='bounded')
        if opt_result['success']:
            # print('old parameter: %f, new: %f ' % (self.par['cov'][0, 1], opt_result['x']))
            self.par['cov'][0, 1] = opt_result['x']
            self.par['cov'][1, 0] = opt_result['x']
        else:
            print('Maximization of log likelihood was unsuccessful in cop2d_gaussian')

        # defining the conditional cdf
        self.createConditionalCDF()

    # Simulation
    def simulate(self, x):
        redistributed = np.transpose(np.random.multivariate_normal([0 for i in range(2)], self.par['cov'], x))
        return ut.uniforms(redistributed.tolist())

    # density function
    def pdf(self, unifs):

        cov = self.par['cov']
        not_list = np.ndim(unifs[0]) == 0
        if not_list:
            unifs = [[i] for i in unifs]

        # the corresponding points in the 'real space' (ie. where the are distributed according to a normal law))
        vecs = []

        marginal_density = []
        for i in range(2):
            a = stats.norm([0], np.sqrt(cov[i, i]))
            temp = [float(a.ppf(j)) for j in unifs[i]]
            marginal_density.append([float(a.pdf(j)) for j in temp])
            vecs.append(temp)
        res = []

        vecs = list(zip(*vecs))
        marginal_density = list(zip(*marginal_density))
        factor = np.sqrt((2 * math.pi) ** 2 * np.linalg.det(cov))
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

    # partial CDF: dC/dx
    # coor specifies along which coordinate (x0 or x1) it is differentiated
    def createConditionalCDF(self):
        cov = self.par['cov']

        new_var = (1 - cov[0, 1] ** 2)
        self.norm_obj0 = stats.norm(0, 1)
        self.norm_obj1 = stats.norm(0, np.sqrt(new_var))
        self.last_par_1 = cov[0, 1]

        def conditionalCDF(unifs, coor):

            if not self.last_par_1 == self.par['cov'][0, 1]:
                new_var = (1 - cov[0, 1] ** 2)
                self.norm_obj1 = ut.fast_stats(stats.norm([0], np.sqrt(new_var)))
                self.last_par_1 = cov[0, 1]

            not_list = np.ndim(unifs) < 2
            if not_list:
                unifs = [[i] for i in unifs]

            if (coor != 0) & (coor != 1):
                raise (RuntimeError('coor specifies the coordinate: it must be either 0 or 1'))
            if coor == 1:
                other = 0
            else:
                other = 1

            # the corresponding points in the 'real space' (ie. where the are distributed according to a normal law))
            vecs = [self.norm_obj0.ppf(unifs[i]) for i in range(2)]
            temp = [i[other] - cov[0, 1] * i[coor] for i in zip(*vecs)]
            res = self.norm_obj1.cdf(temp)
            return res

        self.conditionalCDF = conditionalCDF
        return conditionalCDF

    # log-likelihood
    def log_likelihood(self, unifs, theta=None, vects=None):
        if unifs is None:
            unifs = self.uniforms
        not_list = np.ndim(unifs[0]) == 0
        if not_list:
            unifs = [[i] for i in unifs]

        use_theta = False
        if theta is not None:
            if -1 < theta < 1:
                use_theta = True
                # else:
                #     raise(RuntimeError('theta must be in the interval (-1,1)'))

        if not use_theta:
            theta = self.par['cov'][0, 1]

        if not -1 < theta < 1:
            return self.log_likelihood(unifs, theta=0.99, vects=vects) * 10 * theta ** 2
        else:

            a = stats.norm([0], [1])

            log0, log1 = 0, 0
            if vects is None:
                for i in zip(*unifs):
                    i = [a.ppf(i[0])[0], a.ppf(i[1])[0]]
                    log0 = log0 + i[0] ** 2 + i[1] ** 2
                    log1 = log1 + i[0] * i[1]
            else:
                for i in zip(*vects):
                    log0 = log0 + i[0] ** 2 + i[1] ** 2
                    log1 = log1 + i[0] * i[1]

            log0 = -log0 * theta ** 2 / (1 - theta ** 2) / 2
            log1 = log1 * theta / (1 - theta ** 2)

            return log0 + log1 - len(unifs[0]) * math.log(1 - theta ** 2) / 2

    def parameter_list(self):
        return [self.par['cov'][0, 1]]

    def assign_parameter(self, theta):
        self.par['cov'] = np.matrix([[1, theta[0]], [theta[0], 1]])

    # maximize log-likelihood
    def obs_likelihood(self, theta):
        if np.ndim(theta) > 0:
            theta = theta[0]

        if -1 < theta < 1:
            return -self.log_likelihood(self.uniforms, theta=theta, vects=self.vects)
        else:
            return 10 * self.obs_likelihood(0.99) * theta ** 2


### returns a 'student copula' fitted to the given points
class cop2d_student(cop2d):
    name = 'student'
    ACR = 'ST'

    def __init__(self, uniforms, precise=False):
        length, format = ut.good_format(uniforms)
        if not (format):
            print(uniforms)
            raise (RuntimeError('uniforms do not have the good format'))

        if len(uniforms) != 2:
            raise (RuntimeError('uniforms should be of dimension 2'))

        # computing correlation
        a = stats.norm([0], [1])
        vects = [[a.ppf(j)[0] for j in coor] for coor in uniforms]
        dim = len(uniforms)
        cor = norm_matrix(np.cov(np.matrix(vects)))
        cor_inv = np.linalg.inv(cor)

        self.cor = cor
        self.precise = precise
        self.uniforms = uniforms
        self.dim = 2
        self.name = 'student'
        self.length = length
        self.bounds = [((-1 + 9 * cor[0, 1]) / 10, (1 + 9 * cor[0, 1]) / 10), (2 + epsilon, 10)]

        theta = self.find_parameter()

        # print(theta)
        v = max(theta[1], 2 + epsilon)
        s = np.matrix([[1, theta[0]], [theta[0], 1]])  # no need to multiply by factor (v-2)/v for the uniform
        par = {'v': v, 'sigma': s}

        self.par = par
        self.createConditionalCDF()

    # definition of the log_likelihood
    def log_likelihood(self, unifs, theta=None):

        if unifs is None:
            unifs = self.uniforms
        not_list = np.ndim(unifs[0]) == 0
        if not_list:
            unifs = [[i] for i in unifs]

        if theta is None:
            v = self.par['v']
            sigma = norm_matrix(self.par['sigma'])
        else:
            # print(theta)
            v = theta[1]
            sigma = np.matrix([[1, theta[0]], [theta[0], 1]])

        sigma_inv = sigma ** (-1)

        # back in the real distribution space...
        vecs = []
        marginal_density = []
        t_obj = stats.t(v)

        for i in range(self.dim):
            temp = [float(j) for j in t_obj.ppf(unifs[i])]
            marginal_density.append([float(j) for j in t_obj.pdf(temp)])
            vecs.append(temp)
        vecs = list(zip(*vecs))
        marginal_density = list(zip(*marginal_density))

        # factor of the density in the real space
        factor = math.log(spe.gamma((v + self.dim) / 2)) \
                 - math.log(spe.gamma(v / 2)) \
                 - 0.5 * self.dim * math.log((math.pi * v)) \
                 - 0.5 * math.log(np.linalg.det(sigma))

        # computing the log likelihood
        res = []
        for pt in list(zip(*[vecs, marginal_density])):
            x = np.matrix(pt[0])
            try:
                m = np.sum([math.log(density) for density in pt[1]])
                temp = -m + factor + math.log(1 + x * sigma_inv * np.transpose(x) / v) * (-(v + self.dim) / 2)
            except:
                print(self.par)
                print('error computing the density of student copula: pdf value is too low')
                temp = 1
            res.append(temp)

        return sum(res)

    # parameters of the student copula

    # maximizing log-likelihood to find the parameters
    def find_parameter(self):
        cor = self.cor
        if self.precise:
            opt_result = optimize.minimize(lambda x: -(self.log_likelihood(self.uniforms, theta=x)), [cor[0, 1], 2.1],
                                           bounds=self.bounds, method='L-BFGS-B')
        else:
            opt_result = optimize.minimize(lambda x: -(self.log_likelihood(self.uniforms, theta=[cor[0, 1], x[0]])),
                                           [2.1],
                                           bounds=[(2 + epsilon, 10)], method='L-BFGS-B')

        if opt_result['success']:
            if self.precise:
                theta = opt_result['x']
            else:
                theta = [cor[0, 1]]
                theta.append(float(opt_result['x']))
        else:
            print('Warning: optimization did not converge')
            theta = [cor[0, 1], 5]
        return theta

    # computing t distribution (see " https://en.wikipedia.org/wiki/Multivariate_t-distribution")
    def simulate(self, nb, theta=None):
        if theta is None:
            theta = [self.par['sigma'][0, 1], self.par['v']]
        s = np.matrix([[1, theta[0]], [theta[0], 1]])
        Y = np.transpose(np.random.multivariate_normal([0 for i in range(self.dim)], s, nb))
        u = np.sqrt(theta[1] / np.random.chisquare(theta[1], nb))
        temp = []
        for yy in Y:
            temp.append(list(np.array(yy) * u))
        return (ut.uniforms(temp))

    # density function definition
    def pdf(self, unifs):
        not_list = np.ndim(unifs[0]) == 0
        if not_list:
            unifs = [[i] for i in unifs]

        v = self.par['v']
        sigma = norm_matrix(self.par['sigma'])
        sigma_inv = sigma ** (-1)

        # back in the real distribution space...
        vecs = []
        marginal_density = []
        t_obj = stats.t(v)
        for i in range(self.dim):
            temp = [float(j) for j in t_obj.ppf(unifs[i])]
            marginal_density.append([float(j) for j in t_obj.pdf(temp)])
            vecs.append(temp)
        res = []
        vecs = list(zip(*vecs))
        marginal_density = list(zip(*marginal_density))

        # factor of the density in the real space
        factor = spe.gamma((v + self.dim) / 2) / (
            spe.gamma(v / 2) * math.sqrt((math.pi * v) ** self.dim * np.linalg.det(sigma)))

        for i in list(zip(*[vecs, marginal_density])):
            x = np.matrix(i[0])
            m = np.prod(i[1])
            if m > 0:
                temp = 1 / m * factor * math.exp(
                    math.log(1 + x * sigma_inv * np.transpose(x) / v) * (-(v + self.dim) / 2))
            else:
                # print(self.par)
                # print('error computing the density of student copula: pdf value is too low')
                temp = 0
            res.append(temp)

        return res

    # estimating the conditional distributions
    def createConditionalCDF(self):
        v = self.par['v']
        sigma = self.par['sigma']
        self.condCDFpar = [v, sigma[0, 1]]

        factor1 = sigma[0, 1] / 1
        factor2 = 1 - sigma[0, 1] / 1 * sigma[0, 1]
        variance = math.sqrt(v / (v - 2))

        t_obj = stats.t(v)
        new_t_obj = stats.t(v + 1)

        def conditionalCDF(unifs, coor):
            if not ((self.par['v'] == self.condCDFpar[0]) & (self.par['sigma'][0, 1] == self.condCDFpar[1])):
                return self.createConditionalCDF()(unifs, coor)

            not_list = np.ndim(unifs) < 2
            if not_list:
                unifs = [[i] for i in unifs]

            if (coor != 0) & (coor != 1):
                raise (RuntimeError('coor specifies the coordinate: it must be either 0 or 1'))
            if coor == 0:
                other = 1
            else:
                other = 0

            vecs = np.array([t_obj.ppf(unifs[i]) for i in range(2)]) * variance
            new_sigmas = (v + vecs[coor] ** 2) / (v + 1) * factor2
            new_obs = (vecs[other] - factor1 * vecs[coor]) / np.sqrt(new_sigmas)
            probas = new_t_obj.cdf(new_obs)

            return probas

        self.conditionalCDF = conditionalCDF
        return conditionalCDF

    def parameter_list(self):
        return [self.par['sigma'][0, 1], self.par['v']]

    def assign_parameter(self, theta):
        self.par['sigma'] = np.matrix([[1, theta[0]], [theta[0], 1]])
        self.par['v'] = theta[1]


### return a uniform copula for comparison purpose
class cop2d_uniform(cop2d):
    name = 'uniform'
    ACR = 'UN'

    def __init__(self, unifs):
        if unifs is None:
            unifs = [[0.5], [0.5]]
        self.dim = len(unifs)
        self.length = len(unifs[0])
        self.bounds = []

    def simulate(self, count):
        return np.random.uniform(0, 1, (self.dim, count))

    def pdf(self, x):
        if len(np.shape(x)) == 1:
            x = [x]
        return [1 for i in range(len(x[0]))]

    def conditionalCDF(self, unifs, coor):
        if coor == 0:
            other = 1
        else:
            other = 0
        return unifs[other]

    def log_likelihood(self, unifs, theta=None):
        return 0

    def parameter_list(self):
        return []

    def assign_parameter(self, theta):
        pass


######################
### From dominique ###
######################




# A class which returns a copula which maximizes the log likelihood of the weighted sums of copulas
# Max(x1 * c1 + x2 * c2 + x3 * c3 + ... + xn * cn)
# where x1 + x2 + x3 + ... + xn = 1,
#     x1, x2, ... , xn >= 0
#
# To initialize pass in a list of points, and a list of copula models to test
class WeightedCopula(cop2d):
    ACR = 'WE'

    def __init__(self, vects, models, precise=True, max_models=10):
        dim = len(vects)

        if dim < 2:
            raise ValueError('Data must be of at least dimension 2')
        elif dim > 2:
            vects = vects[:2]  # Using only first two rows

        count = len(models)
        copulas = [copula(vects) for copula in models]
        name = 'weighted'
        if not precise:
            name += '_simple'
        for i in copulas:
            name = '%s_%s' % (name, i.ACR)
            self.ACR += '-' + i.ACR

        if count > max_models:
            lld = [len(cop.parameter_list()) - cop.log_likelihood(vects) for cop in copulas]
            # print(lld)
            indexes = sorted(range(count), key=lld.__getitem__)[:max_models]
            copulas = list(map(copulas.__getitem__, indexes))
            models = list(map(models.__getitem__, indexes))
            count = max_models

        self.models = models
        self.copulas = copulas
        self.count = count
        self.uniforms = vects

        self.name = name
        '''
        if not precise:
            cons = {'type': 'eq', 'fun': lambda x: sum(x) - 1}
            bnds = ((0, 1),) * count

            guess = [1 / count] * count

            def func(xs):
                trans = [x / sum(xs) for x in xs]
                return -self.log_likelihood(None, theta=trans) + (1 - sum(xs)) ** 2

            opt_results = optimize.minimize(lambda x: -self.log_likelihood(None, theta=x), guess, method='SLSQP',
                                            bounds=bnds, constraints=cons)
            if opt_results['success']:
                self.par = opt_results['x']
            else:
                raise RuntimeError('No optimal weights found')
        else:
        '''
        par_indexes = [count]
        guess = [1 / count] * count
        bounds = [(0, 1), ] * count
        optimize_too = []  # the copulas index whose parameter will be optimized along the weights

        if precise:
            for i in range(count):
                cop = self.copulas[i]
                if not cop.ACR in {'GA', 'ST'}:
                    optimize_too.append(i)
                    par_tp = cop.parameter_list()
                    par_indexes.append(par_indexes[-1] + len(par_tp))
                    guess.extend(par_tp)
                    bounds.extend(cop.bounds)

        self.optimize_too = optimize_too
        self.last_pdfs = [None] * self.count
        self.last_par = [np.inf] * self.count
        self.par_indexes = par_indexes

        def log_likelihood_opt(par):
            a = 10
            for i in range(len(optimize_too)):
                self.copulas[optimize_too[i]].assign_parameter(par[par_indexes[i]:par_indexes[i + 1]])

            sum_coeffs = sum(par[:count])
            res_tp = -self.log_likelihood(None, theta=par)
            res_tp += a * (1 - sum_coeffs) ** 2

            return res_tp

        self.log_likelihood_opt = log_likelihood_opt

        opt_results = optimize.minimize(log_likelihood_opt, guess, bounds=bounds)
        if opt_results['success']:
            sum_coeffs = sum(opt_results['x'][:count])
            for i in range(count):
                opt_results['x'][i] /= sum_coeffs
            self.par = opt_results['x']
        else:
            raise RuntimeError('No optimal weights found')

    def single_pdf(self, x, weights=None):
        """Returns a single pdf given a list x1, x2, ..., xn of weights and a point x"""
        if weights is None:
            weights = self.par[:self.count]
        assert len(weights) == len(self.models)

        return sum(weight * copula.pdf(x)[0] for (weight, copula) in zip(weights, self.copulas))

    def log_likelihood(self, unifs, theta=None):
        unifs_none = unifs is None
        if unifs_none:
            unifs = self.uniforms
        if theta is None:
            theta = self.par

        # checking if the pdf values have already been computed
        already_computed = [(self.last_pdfs[i] is not None) & unifs_none for i in range(self.count)]

        for i in range(len(self.optimize_too)):
            index = self.optimize_too[i]
            start, end = self.par_indexes[i], self.par_indexes[i + 1]
            last_parameter = [j for j in self.last_par[start:end]]
            new_parameter = [j for j in theta[start:end]]
            already_computed[index] &= all([j[0] == j[1] for j in zip(*[last_parameter, new_parameter])])

        pdfs_tp = []
        for i in range(self.count):
            if already_computed[i]:
                pdfs_tp.append(self.last_pdfs[i])
            else:
                pdfs_tp.append(self.copulas[i].pdf(unifs))

        if unifs_none:
            self.last_par = theta.copy()
            self.last_pdfs = pdfs_tp

        # computing the log-likelihood
        pdfs = zip(*pdfs_tp)
        sum_coeffs = sum(theta[:self.count])
        coeffs = [i / sum_coeffs for i in theta[:self.count]]
        res = sum(np.log(sum(weight * pdf for (weight, pdf) in zip(coeffs, point))) for point in pdfs)
        return res

    def pdf(self, unifs):
        if np.ndim(unifs) < 2:
            unifs = [[i] for i in unifs]
        points = zip(*unifs)
        return [self.single_pdf(point) for point in points]

    def conditionalCDF(self, unifs, coor):
        if np.ndim(unifs) < 2:
            unifs = [[i] for i in unifs]
        return list(
            np.sum([np.array(self.copulas[i].conditionalCDF(unifs, coor)) * self.par[i] for i in range(self.count)], 0))

    '''
    def simulate(self, n):
        """Simulates n random numbers from weighted distribution, only works in 2D now"""
        curr_U, curr_V = np.zeros(n), np.zeros(n)
        for i, copula in enumerate(self.copulas):
            U, V = np.array(copula.simulate(n))
            curr_U += self.par[i] * U
            curr_V += self.par[i] * V
        return list(curr_U), list(curr_V)
    '''

    def simulate(self, count):
        nb_models = len(self.copulas)
        which = []
        bounds = np.cumsum(self.par[:count])
        nb_simulations = [0 for i in range(nb_models)]
        for i in np.random.uniform(0, 1, count):
            for j in range(nb_models):
                if i <= bounds[j]:
                    which.append(j)
                    nb_simulations[j] = nb_simulations[j] + 1
                    break

        # print(which)
        # print(nb_simulations)
        # print([len(self.copulas[j].simulate(nb_simulations[j])[0]) for j in range(nb_models)])
        simulated = [iter(list(zip(*self.copulas[j].simulate(nb_simulations[j])))) for j in range(nb_models)]

        res = []
        for i in which:
            res.append(next(simulated[i]))
        return [list(i) for i in zip(*res)]

    def get_names(self):
        return [self.copulas[i].ACR for i in range(self.count) if self.par[i] != 0]

    def parameter_list(self):
        return self.par

    def assign_parameter(self, theta):
        self.par = theta

    __call__ = pdf


### returns a 'clayton copula' fitted to the given points
class cop2d_clayton(cop2d):
    name = 'clayton'
    ACR = 'CL'

    def __init__(self, vects, length=0):
        dim = len(vects)

        if dim < 2:
            raise ValueError('Data must be of at least dimension 2')
        elif dim > 2:
            vects = vects[:2]  # Using only first two rows
        self.uniforms = vects
        self.name = 'clayton'
        self.bounds = bounds = [(epsilon, 10)]

        tau = stats.kendalltau(vects[0], vects[1])[0]
        # Calculate Theta
        theta = 2 * tau / (1 - tau)
        opt_result = optimize.minimize_scalar(lambda x: -self.log_likelihood(self.uniforms, theta=x),
                                              bounds=bounds[0], method='bounded')
        if opt_result['success']:
            theta = opt_result['x']
        else:
            theta = max(epsilon, theta)

        self.dim = dim
        self.length = len(vects[0])
        self.par = {'theta': theta}

    def cdf(self, u, v):
        """Computes theta at a single point"""
        theta = self.par['theta']
        return (u ** (-theta) + v ** (-theta) - 1) ** (-1 / theta)

    def single_pdf(self, u, v, theta=None):
        """computes pdf for a single point given a parameter theta"""
        if theta is None:
            theta = self.par['theta']
        return (theta + 1) * (u ** (-theta) + v ** (-theta) - 1) ** ((-1 / theta) - 2) * u ** (-theta - 1) * v ** (
            -theta - 1)

    def pdf(self, unifs):
        """Computes pdf for a list of lists"""
        theta = self.par['theta']
        if np.ndim(unifs) < 2:
            unifs = [[i] for i in unifs]
        points = zip(*unifs)
        return [self.single_pdf(u, v, theta) for (u, v) in points]

    def conditional(self, u, v):
        """Calculated based on first order derivative of cdf, P(v <= v| u=u)"""
        theta = self.par['theta']
        return (u ** (-theta) + v ** (-theta) - 1) ** (-1 / theta - 1) * u ** (-theta - 1)

    def conditionalCDF(self, unifs, coor):
        not_list = np.ndim(unifs[0]) == 0
        if not_list:
            unifs = [[i] for i in unifs]

        if coor == 0:
            other = 1
        else:
            other = 0

        return [self.conditional(i[coor], i[other]) for i in zip(*unifs)]

    def _ll(self, theta):
        """Alternate formulation of log-likelihood, unsure if they are the same"""
        n = self.length
        m = self.dim
        first = n * (m * np.log(theta) + np.log(spe.gamma(1 / theta + m)) - np.log(spe.gamma(1 / theta)))
        second = -(theta + 1) * sum(sum(x for x in vect) for vect in self.uniforms)
        third = -(1 / theta + m) * sum(np.log(sum(x ** (-theta) for x in vect)) - m + 1 for vect in self.uniforms)
        return first + second + third

    def log_likelihood(self, unifs, theta=None):
        if theta is None:
            theta = self.par['theta']
        if unifs is None:
            unifs = self.uniforms
        return sum(np.log(self.single_pdf(u, v, theta)) for (u, v) in zip(*unifs))

    def simulate(self, n):
        """Generates based on inverse conditional distribution"""
        theta = self.par['theta']
        U = np.random.uniform(0, 1, n)
        W = np.random.uniform(0, 1, n)
        V = U * (W ** (theta / -(1 + theta)) - 1 + U ** theta) ** (-1 / theta)
        return U, V

    def parameter_list(self):
        return [self.par['theta']]

    def assign_parameter(self, theta):
        self.par['theta'] = theta[0]


### returns a 'Gumbel Hougard copula' fitted to the given points
# WARNING: this should not be mistaken for an EFGM ('Eyraud Farlie Gumbel Morgenstern') copula
class cop2d_gumbel(cop2d):
    name = 'gumbel'
    ACR = 'GU'

    def __init__(self, vects, length=0):
        dim = len(vects)

        if dim < 2:
            raise ValueError('Data must be of at least dimension 2')
        elif dim > 2:
            vects = vects[:2]  # Using only first two rows
        self.uniforms = vects
        self.name = 'gumbel'
        self.bounds = bounds = [(1 + epsilon, 10)]

        tau = stats.kendalltau(vects[0], vects[1])[0]
        # Calculate Theta
        theta = 1 / (1 - tau)
        opt_result = optimize.minimize_scalar(lambda x: -self.log_likelihood(self.uniforms, theta=x),
                                              bounds=bounds[0], method='bounded')
        if opt_result['success']:
            theta = opt_result['x']
        else:
            theta = max(1 + epsilon, theta)

        self.dim = dim
        self.length = len(vects[0])
        self.par = {'theta': theta}

    def cdf(self, u, v, theta):
        """Computes cdf for a single point"""
        interior = -((-np.log(u)) ** theta + (-np.log(v)) ** theta) ** (1 / theta)
        return np.exp(interior)

    def single_pdf(self, u, v, theta=None):
        """Computes pdf for a single point given parameter theta"""
        if theta is None:
            theta = self.par['theta']
        alpha = 1 / theta
        x = self.gen(u, theta) + self.gen(v, theta)
        f1 = alpha
        f2 = np.exp(-x ** alpha)
        f3 = x ** (-2 + alpha)
        f4 = (alpha * (x ** alpha - 1) + 1)
        f5 = theta * (-np.log(u)) ** (theta - 1) * theta * (-np.log(v)) ** (theta - 1) * (1 / (u * v))
        return f1 * f2 * f3 * f4 * f5

    def pdf(self, unifs):
        """computes pdf for a list of lists [[x1, x2, ..., xn], [y1, y2, ..., yn"""
        if np.ndim(unifs) < 2:
            unifs = [[i] for i in unifs]
        points = zip(*unifs)
        return [self.single_pdf(u, v, self.par['theta']) for (u, v) in points]

    def conditional(self, u, v):
        theta = self.par['theta']
        first = np.exp(-((-np.log(u)) ** theta + (-np.log(v)) ** theta) ** (1 / theta))
        second = ((-np.log(u)) ** theta + (-np.log(v)) ** theta) ** (1 / theta - 1)
        third = (-np.log(u)) ** (theta - 1) * (1 / u)
        return first * second * third

    def conditionalCDF(self, unifs, coor):
        not_list = np.ndim(unifs[0]) == 0
        if not_list:
            unifs = [[i] for i in unifs]

        if coor == 0:
            other = 1
        else:
            other = 0

        return [self.conditional(i[coor], i[other]) for i in zip(*unifs)]

    def gen(self, x, theta):
        return (-np.log(x)) ** theta

    def _ll(self, u, v, theta):
        """Alternate formulation of log_likelihood, should be equivalent to log(single_pdf)"""
        alpha = 1 / theta
        x = self.gen(u, theta) + self.gen(v, theta)
        f1 = alpha
        f2 = np.exp(-x ** alpha)
        f3 = x ** (-2 + alpha)
        f4 = (alpha * (x ** alpha - 1) + 1)
        f5 = theta * (-np.log(u)) ** (theta - 1) * theta * (-np.log(v)) ** (theta - 1) * (1 / (u * v))
        return np.log(f1) + np.log(f2) + np.log(f3) + np.log(f4) + np.log(f5)

    def log_likelihood(self, unifs, theta=None):
        if theta is None:
            theta = self.par['theta']
        if unifs is None:
            unifs = self.uniforms
        return sum(np.log(self.single_pdf(u, v, theta)) for (u, v) in zip(*unifs))

    def simulate(self, n):
        """Generates first a sample from a stable distribution st(1/theta, 1, gamma, 0) according to Wikipedia
            Then performs transform according to """
        theta = self.par['theta']
        gamma = np.cos(np.pi / (2 * theta)) ** theta
        alpha = 1 / theta
        zeta = -np.tan(np.pi * alpha / 2)
        xi = (1 / alpha) * np.arctan(-zeta)
        u1 = np.random.uniform(-np.pi / 2, np.pi / 2, n)
        u2 = np.random.uniform(0, 1, n)
        W = -np.log(u2)
        first = (1 + zeta ** 2) ** (1 / (2 * alpha))
        second = np.sin(alpha * (u1 + xi)) / np.cos(u1) ** (1 / alpha)
        third = (np.cos(u1 - alpha * (u1 + xi)) / W) ** ((1 - alpha) / alpha)
        X = first * second * third  # X is distributed according to st(1/theta, 1, 0, 0)
        scaled = gamma * X

        x1 = np.random.uniform(0, 1, n)
        x2 = np.random.uniform(0, 1, n)
        t1 = -np.log(x1) / scaled
        t2 = -np.log(x2) / scaled
        U = np.exp(-t1 ** (1 / theta))
        V = np.exp(-t2 ** (1 / theta))
        return U, V

    def parameter_list(self):
        return [self.par['theta']]

    def assign_parameter(self, theta):
        self.par['theta'] = theta[0]


### returns a 'franck copula' fitted to the given points
class cop2d_frank(cop2d):
    name = 'frank'
    ACR = 'FR'

    def __init__(self, vects, length=0):
        dim = len(vects)

        if dim < 2:
            raise ValueError('Data must be of at least dimension 2')
        elif dim > 2:
            vects = vects[:2]  # Using only first two rows

        self.uniforms = vects
        self.name = 'frank'

        self.tau = tau = stats.kendalltau(vects[0], vects[1])[0]
        # Calculate Theta
        theta = self.find_root()
        # print(theta)
        if theta > 0:
            self.bounds = bounds = [(epsilon, 40)]
            opt_result = optimize.minimize_scalar(lambda x: -self.log_likelihood(vects, theta=x), bounds=bounds[0],
                                                  method='bounded')
            if opt_result['success']:
                theta = opt_result['x']
            else:
                print('optimiztion did not converge')
                theta = max(epsilon, theta)
        else:
            self.bounds = bounds = [(-40, -epsilon)]
            opt_result = optimize.minimize_scalar(lambda x: -self.log_likelihood(vects, theta=x),
                                                  bounds=bounds[0],
                                                  method='bounded')
            if opt_result['success']:
                theta = opt_result['x']
            else:
                print('optimiztion did not converge')
                theta = min(-epsilon, theta)

        self.dim = dim
        self.length = len(vects[0])
        self.par = {'theta': theta}

    def cdf(self, u, v):
        """Computes cdf at a single point"""
        theta = self.par['theta']
        return -(1 / theta) * np.log(1 + ((np.exp(-theta * u) - 1) * (np.exp(-theta * v) - 1)) / (np.exp(-theta) - 1))

    def pdf(self, unifs):
        """Computes pdf for a list of lists"""
        if np.ndim(unifs) < 2:
            unifs = [[i] for i in unifs]
        points = zip(*unifs)
        return [self.single_pdf(u, v, self.par['theta']) for (u, v) in points]

    def single_pdf(self, u, v, theta=None):
        """Computes pdf for a single point"""
        if theta is None:
            theta = self.par['theta']
        if theta < 0:
            numer = theta * (1 - np.exp(-theta)) * np.exp(-theta * (u + v))
            denom = (1 - np.exp(-theta) - (1 - np.exp(-theta * u)) * (1 - np.exp(-theta * v))) ** 2
            return numer / denom
        else:
            """
            x = self.gen(u, theta) + self.gen(v, theta)
            phi = ((theta * np.exp(-theta * u))/(np.exp(-theta * u) - 1) *
            (theta * np.exp(-theta * v)) / (np.exp(-theta * v) - 1))
            y = (1 / (1 + np.exp(-x) * (np.exp(-theta) - 1)))
            return (1/theta) * y * (y - 1) * phi
            """
            numer = (np.exp(theta) - 1) * theta * np.exp(theta * (u + v + 1))
            denom = (np.exp(theta) - np.exp(theta + theta * u) - np.exp(theta + theta * v) + np.exp(
                theta * (u + v))) ** 2
            return numer / denom

    def conditional(self, u, v):
        theta = self.par['theta']
        numer = (np.exp(-theta * u) * (np.exp(-v * theta) - 1))
        denom = (np.exp(-theta) - 1 + (np.exp(-theta * u) - 1) * (np.exp(-theta * v) - 1))
        return numer / denom

    def conditionalCDF(self, unifs, coor):
        not_list = np.ndim(unifs[0]) == 0
        if not_list:
            unifs = [[i] for i in unifs]

        if coor == 0:
            other = 1
        else:
            other = 0

        return [self.conditional(i[coor], i[other]) for i in zip(*unifs)]

    def gen(self, x, theta):
        return -np.log((np.exp(-theta * x) - 1) / (np.exp(-theta) - 1))

    def log_likelihood(self, unifs, theta=None):
        if theta is None:
            theta = self.par['theta']
        if unifs is None:
            unifs = self.uniforms
        return sum(np.log(self.single_pdf(u, v, theta)) for (u, v) in zip(*unifs))

    def debye(self, x):
        if x == 0: return 1
        return (1 / x) * quad(lambda t: t / (np.exp(t) - 1), 0, x)[0]

    def frank_fun(self, x):
        """Function for determining theta based on relation to tau"""
        return 1 - 4 * (1 / x) * (1 - self.debye(x)) - self.tau

    def find_root(self):
        if self.tau == 0: return 0
        start_a, start_b = -10, 10
        while self.frank_fun(start_a) > 0:
            start_a *= 10
        while self.frank_fun(start_b) < 0:
            start_b *= 10
        return optimize.brentq(self.frank_fun, start_a, start_b)

    def simulate(self, n):
        """Generates n random variables based on inverse transform"""
        theta = self.par['theta']
        U = np.random.uniform(0, 1, n)
        W = np.random.uniform(0, 1, n)
        V = -np.log(
            ((1 - W) / W * np.exp(-U * theta) + np.exp(-theta)) / (1 + ((1 - W) / W) * np.exp(-U * theta))) / theta
        return U, V

    def parameter_list(self):
        return [self.par['theta']]

    def assign_parameter(self, theta):
        self.par['theta'] = theta[0]


# ---------------------------------------------- utils ------------------------------------------------------------------

# Makes sure that the correlation matrix is positive:
# if not, computes the eigenvalue form and replace the negative eigenvalues by a small epsilon
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
