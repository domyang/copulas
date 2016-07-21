import numpy as np


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
    U = uniform_sample()
    S = student_sample()
    G = gaussian_sample()
