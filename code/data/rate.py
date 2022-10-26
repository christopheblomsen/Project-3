import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt


"""We define two functions since we need all simulation done before we can plot the error
so max_error is here to find the maximum deviation between the analytical and numerical for both methods
and then rel_error will just compute the error"""

def max_error(r_analytical, r_numerical):
   """r_analytical and r_numerical are (3, N) so we compute
   the norm for each column returning a N size array"""
   return np.max(np.linalg.norm(r_analytical - r_numerical, axis=0))


def rate(del_max, h_k):

    rate = 0

    """two arrays with all the values of interest for the different step size"""
    for k in range(1, len(del_max)):
        log_err = np.log(del_max[k] / del_max[k-1])
        log_step = np.log(h_k[k] / h_k[k-1])

        rate += log_err / log_step

    return 1/3 * rate


def rel_err(r_true, r_num):
   N = np.shape(r_true)[1]
   rel_err_arr = np.zeros(N)
   abs_err = np.linalg.norm(r_true - r_num[:, :, 0], axis = 0)
   abs_true = np.linalg.norm(r_true, axis = 0)
   rel_err = abs_err/abs_true
   return rel_err

##########calculate analytical solution


x0 = 20
z0 = 20
vy0 = 25

B0 = 9.65e1
V0 = 2.41e6
m = 40.08
d = 500
q = 1

w_z = np.sqrt(2 * q * V0 / (m * d**2))
w_0 = B0*q / m
w_plus = (w_0 + np.sqrt(w_0**2 - 2*w_z**2)) / 2
w_minus = (w_0 - np.sqrt(w_0**2 - 2*w_z**2)) / 2

A_plus = (vy0 + w_minus * x0) / (w_minus - w_plus)
A_minus = - (vy0 + w_plus * x0) / (w_minus - w_plus)

R_plus = np.abs(A_plus + A_minus)
R_minus = np.abs(A_plus - A_minus)

r_analytical = []
time_4 = np.linspace(0, 50, 4000)
time_8 = np.linspace(0, 50, 8000)
time_16 = np.linspace(0, 50, 16000)
time_32 = np.linspace(0, 50, 32000)
time_list = [time_4, time_8, time_16, time_32]
for time in  time_list:
   x = A_plus * np.cos(-w_plus * time) + A_minus * np.cos(-w_minus * time)
   y = A_plus * np.sin(-w_plus * time) + A_minus * np.sin(-w_minus * time)
   z = z0 * np.cos(w_z * time)
   r_analytical.append(np.array([x, y, z]))




   ################### comparing to num solution



methods = ['RK4_', 'FE_']
iterations = ['4000', '8000', '16000', '32000']
it_int = [4000, 8000, 16000, 32000]
rel = []
time = []

for i in range(3):
    time.append(np.linspace(0, 50, it_int[i]))
    r = pa.cube()
    r.load('position_RK4_'+ iterations[i] + '.bin')
    r = np.array(r)
    x = r[:, 0]
    y = r[:, 1]
    z = r[:, 2]
    r = np.array([x, y, z])
    rel.append(rel_err(r_analytical[i], r))
    plt.plot(time[i], rel[i])

plt.show()
