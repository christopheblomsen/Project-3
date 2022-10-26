import pyarma as pa
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys

# General parameters
# plt.style.use('science')
plt.rc('font', family='Helvetica')
plt.rc('text', usetex=True)

#Analytical solution
B0 = 9.64852558e1 
V0 = 2.41e6
d = 500
masse = 40.8
q = 1

x0 = 20
z0 = 20
vy0 = 25

omega_0 = q*B0/masse
omega_z = np.sqrt(2*q*V0 / (masse*d**2))

omega_minus = (omega_0 - np.sqrt(omega_0**2 - 2 * omega_z**2))/2
omega_plus = (omega_0 + np.sqrt(omega_0**2 - 2 * omega_z**2))/2

A_plus = (vy0 + omega_minus*x0) / (omega_minus - omega_plus)
A_minus = -(vy0 + omega_plus*x0) / (omega_minus - omega_plus)

t = np.linspace(0, 50, 4000)

x = A_plus * np.cos(omega_plus*t) + A_minus * np.cos(omega_minus*t)
y = - A_plus * np.sin(omega_plus*t) - A_minus * np.sin(omega_minus*t)
z = z0 * np.cos(omega_z * t)

# Data
pos = pa.cube()
vel = pa.cube()
time = pa.mat()
# pos.load('pos_interaction_4000.bin')
pos.load(sys.argv[1])
vel.load(sys.argv[2])
time.load(sys.argv[3])
time = np.array(time)

pos = np.array(pos)
print(np.size(pos))
z1 = pos[:, 2, 0]
# for i in range(len(time)):
#     z1[i] = pos[2, i, 0]

# pos = np.array(pos)
# print(f'Size: {np.size(pos)}')
# print(pos[0, 0, 0])
# print(np.size(pos))

print(f'analytical: {z[0:10]}')
print(f'simulation: {z1[0:10]}')
cut = 1500

plt.figure()
plt.plot(time[:cut], z1[:cut], label='Simulation')
plt.plot(t, z, label='Analytical')
plt.hlines(d/2, 0, time[cut], ls = '--', label = 'trap perimeter')
plt.hlines(-d/2, 0, time[cut], ls = '--')
plt.xlabel('Time [s]')
plt.ylabel('position [m]')
plt.legend()
plt.show()

# plt.figure()
# plt.plot(time, pos[:, 0, 0], label='x')
# plt.plot(t, x, label='x_analytical')
# plt.legend()

# plt.figure()
# plt.plot(time, pos[:, 1, 0], label='y')
# plt.plot(t, y, label='y_analytical')
# plt.legend()

cut = 800
plt.figure()
plt.plot(x, y, label='analytical')
plt.plot(pos[:cut,0,0], pos[:cut,1,0], label = 'numerical')
plt.hlines(d/2, -d/2, d/2, label = 'trap perimeter')
plt.hlines(-d/2, -d/2, d/2)
plt.vlines(d/2, -d/2, d/2)
plt.vlines(-d/2, -d/2, d/2)
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')
plt.show()

# for i in range(3):
#     plt.plot(time, pos[:, i, 0], 'b'))
#     plt.xlabel(r'Time [$\mu s$]')
#     plt.ylabel(r'Position [$\mu m$]')
#     plt.show()

# plt.plot(pos[:, 0, 0], pos[:, 1, 0])
# plt.show()
