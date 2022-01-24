import numpy as np
from scipy.interpolate import RectBivariateSpline
import matplotlib.pyplot as plt;plt.ion()


# input data
data = np.array([[0, 1, 2, 3, 4],
                 [1, 2, 3, 4, 5],
                 [2, 3, 4, 5, 6],
                 [3, 4, 5, 6, 7],
                 [4, 5, 6, 7, 8]])

upper_left = 2, 0
upper_right = 0, 2
lower_right = 2, 4   # note that I swapped this
lower_left = 4, 2    # and this
n_steps = 3, 3


# build interpolator
m, n = data.shape
x, y = np.arange(m), np.arange(n)

interpolator = RectBivariateSpline(x, y, data)

# build grid
ul,ur,ll,lr = map(np.array, (upper_left,upper_right,lower_left,lower_right))
assert np.allclose(ul + lr, ur + ll)    # make sure edges are parallel

x, y = ul[:, None, None] \
       + np.outer(ll-ul, np.linspace(0.0, 1.0, n_steps[0]))[:, :, None] \
       + np.outer(ur-ul, np.linspace(0.0, 1.0, n_steps[1]))[:, None, :]

breakpoint()
