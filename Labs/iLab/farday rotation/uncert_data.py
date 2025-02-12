import numpy as np

import matplotlib as plt

def uncert_mean(vals):
    mean = np.mean(vals)
    uncert = np.sqrt(np.sum((1/vals)**2))*mean
    return np.array([mean,uncert])



meas_1 = np.array([315,
                 308,
                 297,
                 366,
                 370,
                 383])

meas_2 = np.array([322,
                   315,
                   303,
                   366,
                   373,
                   387])

meas_3 = np.array([308,
                   301,
                   289,
                   360,
                   366,
                   377])

print(uncert_mean(meas_1))

print(uncert_mean(meas_2))

print(uncert_mean(meas_3))


