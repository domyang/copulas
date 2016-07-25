import numpy as np
import vines
import copula_analysis as ca

def pos_semidef(dim):
    matrix = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            matrix[i, j] = matrix[j, i] = np.random.uniform(0, 10)
    return np.dot(matrix, matrix.T)


def uniform_sample(size=1000, bounds=None, dim=2):
    if bounds is None:  # if no bounds passed in, defaults to unit cube
        bounds = [[0, 1]] * dim

    sample = [np.random.uniform(bound[0], bound[1], size) for bound in bounds]
    return np.array(sample).T


def student_sample(size=1000, df=None, mean=None, cov=None, dim=2):
    if df is None:
        df = np.random.randint(3, 50)
    if mean is None:
        mean = np.random.uniform(-10, 10, dim)
    if cov is None:
        cov = pos_semidef(dim)

    y = np.random.multivariate_normal(np.zeros(dim), cov, size)
    u = np.random.chisquare(df, size)

    x = mean + (y.T * np.sqrt(df / u)).T
    return x


def gaussian_sample(size=1000, mean=None, cov=None, dim=2):
    if cov is None:
        cov = pos_semidef(dim)
    if mean is None:
        mean = np.random.uniform(-10, 10, dim)
    return np.random.multivariate_normal(mean, cov, size)

if __name__ == '__main__':
    gu = vines.cop2d_gumbel
    ga = vines.cop2d_gaussian2
    st = vines.cop2d_student
    fr = vines.cop2d_frank
    un = vines.cop2d_uniform
    cl = vines.cop2d_clayton

    list1 = [gu, ga, st, fr, un, cl]

    s_samples = []
    g_samples = []
    u_samples = []
    for _ in range(10):
        s_samples.append(student_sample(20000).T.tolist())
        u_samples.append(uniform_sample(20000).T.tolist())
        g_samples.append(gaussian_sample(20000).T.tolist())

    g_results = []
    s_results = []
    u_results = []

    for sample in s_samples:
        copula = ca.fake_copulaManager(sample)
        s_results.append(ca.test_models(copula, list_models=list1))

    for sample in u_samples:
        copula = ca.fake_copulaManager(sample)
        u_results.append(ca.test_models(copula, list_models=list1))

    for sample in g_samples:
        copula = ca.fake_copulaManager(sample)
        g_results.append(ca.test_models(copula, list_models=list1))