import pyarma as pa
import matplotlib.pyplot as plt
import numpy as np

# Graph parameters
plt.style.use('seaborn')
plt.rc('font', family='Helvetica')
plt.rc('text', usetex=True)


def plot_load(f, filename_O, filename_N, linestyle = None, label = None):
    N = pa.mat()
    omega = pa.mat()
    N.load(filename_N)
    omega.load(filename_O)
    
    omega = np.array(omega)
    N = np.array(N)
    
    plt.plot(omega, N, label = label, linestyle = linestyle)
    plt.xlabel(r"$\omega_{\rm V} [\rm MHz]$", fontsize=14)
    plt.xlim(0, 2.8)
    plt.ylabel(r"$N_{\rm trapped} / N_{\rm escaped}$", fontsize=14)
        
# Plotting for part one of the simualtion in task 9:
l_ = ["-", ":", "--"]
f_ = ["0.100000", "0.400000", "0.700000"]

# Plots just with N = 32000
N_ = "32000"
for f, l in zip(f_, l_):
    filename_N = f"fraction_left_f_{f}_N_{N_}_.bin"
    filename_O = f"omega_f_{f}_N_32000_.bin"
    plot_load(f, filename_O, filename_N, linestyle = l, label = f"f={float(f):0.1f}")


#add analytical values
omega = 2.3033
omega_ = 1.044
plt.plot(omega, 0, 'o', color = 'k')
plt.plot(omega_*2, 0, 'o', color = 'k')
plt.legend()
plt.savefig("fraction_left_omega.pdf")
plt.show()


# Plotting for part two of the simualtion in task 9:
inter = ["interaction", "no_interaction_"]
linestyles = ["-", ":"]
f = "0.100000"


'''When we get the data '''
# for i, l in zip(inter, linestyles):
for i in range(2):
    O = inter[i]
    l = linestyles[i]
    filename_O = f"fraction_left_fine_grain_f_0.100000_{O}.bin"
    filename_N = f"omega_fine_grain_f_0.100000_{O}.bin"
    plot_load(f, filename_N, filename_O, linestyle = l, label = O)

plt.legend()
plt.show()


"""
Found out that the dip is only visible for N = 32000
"""
# Plots for all N: 4000, 8000, 16000, 32000
# N_ = ["4000", "8000", "16000", "32000"]
# for _ in N_:
#     for f, l in zip(f_, l_):
#         plot_load(f, linestyle = l, N_ = _)

#     plt.title(f"N = {_}")
#     plt.legend()
#     plt.show()
