from tkinter import S
import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

plt.rc('text', usetex=True)
plt.rc('font', family='Helvetica')
plt.style.use('seaborn')

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
   abs_err = np.linalg.norm(r_true[0] - r_num[0], axis = 0)
   abs_true = np.linalg.norm(r_true, axis = 0)
   rel_err = abs_err/abs_true
   return rel_err



###############Import sim data##########################

xrk4000 = pa.cube()
xrk4000.load('position_RK4_4000.bin')
xrk4000 = np.array(xrk4000)[:,:,0].T

xrk8000 = pa.cube()
xrk8000.load('position_RK4_8000.bin')
xrk8000 = np.array(xrk8000)[:,:,0].T

xrk16000 = pa.cube()
xrk16000.load('position_RK4_16000.bin')
xrk16000 = np.array(xrk16000)[:,:,0].T

xrk32000 = pa.cube()
xrk32000.load('position_RK4_32000.bin')
xrk32000 = np.array(xrk32000)[:,:,0].T

xFE4000 = pa.cube()
xFE4000.load('position_FE_4000.bin')
xFE4000 = np.array(xFE4000)[:,:,0].T

xFE8000 = pa.cube()
xFE8000.load('position_FE_8000.bin')
xFE8000 = np.array(xFE8000)[:,:,0].T

xFE16000 = pa.cube()
xFE16000.load('position_FE_16000.bin')
xFE16000 = np.array(xFE16000)[:,:,0].T

xFE32000 = pa.cube()
xFE32000.load('position_FE_32000.bin')
xFE32000 = np.array(xFE32000)[:,:,0].T


###############analytical solution#####################

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

# print(len(r_analytical))
# To check if our solution is correct 
plt.figure()
plt.scatter(xFE32000[0,:], xFE32000[1, :], s=.7)
plt.scatter(r_analytical[3][0], r_analytical[3][1], s=.7)
plt.axis('equal')


rel_err_RK4_4000 = rel_err(r_analytical[0], xrk4000)
rel_err_RK4_8000 = rel_err(r_analytical[1], xrk8000)
rel_err_RK4_16000 = rel_err(r_analytical[2], xrk16000)
rel_err_RK4_32000 = rel_err(r_analytical[3], xrk32000)

rel_err_FE_4000 = rel_err(r_analytical[0], xFE4000)
rel_err_FE_8000 = rel_err(r_analytical[1], xFE8000)
rel_err_FE_16000 = rel_err(r_analytical[2], xFE16000)
rel_err_FE_32000 = rel_err(r_analytical[3], xFE32000)


fig, [ax1, ax2] = plt.subplots(2, 1, figsize=(12, 10))
ax1.plot(time_list[0], rel_err_RK4_4000, label = 'n = 4000') # log2 on the y-axis
ax1.plot(time_list[1], rel_err_RK4_8000, label = 'n = 8000')
ax1.plot(time_list[2], rel_err_RK4_16000, label = 'n = 16000')
ax1.plot(time_list[3],rel_err_RK4_32000, label = 'n = 32000')
<<<<<<< HEAD
ax1.set_ylabel('Relative error [-]')
ax1.set_xlim(0, 55)
=======
ax1.set_ylabel('log Relative error [-]')
ax1.set_xlim(0, 55)
ax1.set_yscale('log')
>>>>>>> refs/remotes/origin/main
ax1.set_title('Runge Kutta 4')

ax2.plot(time_list[0], rel_err_FE_4000, label = 'n = 4000')
ax2.plot(time_list[1], rel_err_FE_8000, label = 'n = 8000')
ax2.plot(time_list[2], rel_err_FE_16000, label = 'n = 16000')
ax2.plot(time_list[3], rel_err_FE_32000, label = 'n = 32000')
ax2.set_xlabel('t [$\mu$s]')
ax2.set_xlim(0, 55)
<<<<<<< HEAD
ax2.set_ylabel('Relative error [-]')
ax2.set_title('Forward Euler')
=======
ax2.set_ylabel('log Relative error [-]')
ax2.set_title('Forward Euler')
ax2.set_yscale('log')
>>>>>>> refs/remotes/origin/main
plt.legend()
plt.savefig('convergence.pdf')
plt.show()
#######################calculating r_err########################
del_max_RK4 = []
del_max_RK4.append(max_error(r_analytical[0], xrk4000))
del_max_RK4.append(max_error(r_analytical[1], xrk8000))
del_max_RK4.append(max_error(r_analytical[2], xrk16000))
del_max_RK4.append(max_error(r_analytical[3], xrk32000))


del_max_FE = []
del_max_FE.append(max_error(r_analytical[0], xFE4000))
del_max_FE.append(max_error(r_analytical[1], xFE8000))
del_max_FE.append(max_error(r_analytical[2], xFE16000))
del_max_FE.append(max_error(r_analytical[3], xFE32000))

n_k = np.array([4000, 8000, 16000, 32000])
h_k = 50/n_k

rate_RK4 = rate(del_max_RK4, h_k)
rate_FE = rate(del_max_FE, h_k)
print(f"Convergence rate_RK4: {rate_RK4}")
print(f"Convergence rate_FE: {rate_FE}")
plt.plot(del_max_FE, "^-", label='FE')
plt.plot(del_max_RK4, "*-", label='RK4')
plt.savefig('del_max.pdf')
plt.legend()
plt.show()
