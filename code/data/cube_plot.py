import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt

r_RK = pa.cube()
v_RK = pa.cube()
t = pa.mat()

r_RK.load("position_RK4_16000.bin")
v_RK.load("velocity_RK4_16000.bin")
t.load("time_16000.bin")


r_RK = np.array(r_RK)
v_RK = np.array(v_RK)
time = np.array(t)

x_RK, y_RK, z_RK = r_RK[:, 0, 0], r_RK[:, 1, 0], r_RK[:, 2, 0]
vx_RK, vy_RK, vz_RK = v_RK[:, 0, 0], v_RK[:, 1, 0], v_RK[:, 2, 0]

r_FE = pa.cube()
v_FE = pa.cube()

r_FE.load("position_FE_16000.bin")
v_FE.load("velocity_FE_16000.bin")
# t.load("time.bin")


r_FE = np.array(r_FE)
v_FE = np.array(v_FE)
# time = np.array(t)

x_FE, y_FE, z_FE = r_FE[:, 0, 0], r_FE[:, 1, 0], r_FE[:, 2, 0]
vx_FE, vy_FE, vz_FE = v_FE[:, 0, 0], v_FE[:, 1, 0], v_FE[:, 2, 0]

print(z_RK[0:5])

# Analytical solution
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

x = A_plus * np.cos(-w_plus * time) + A_minus * np.cos(-w_minus * time)
y = A_plus * np.sin(-w_plus * time) + A_minus * np.sin(-w_minus * time)
z = z0 * np.cos(w_z * time)

circle1 = plt.Circle((0, 0), R_plus, color='k', fill=False, label=r"$|R_+ - R_-$|")
circle2 = plt.Circle((0, 0), R_minus, color='g', fill=False,  label=r"$R_+ + R_-$")

plt.plot(time, z_RK, label="RK4")
plt.plot(time, z_FE, label="FE")
plt.plot(time, z, label="Analytical")
plt.legend()
plt.show()

print(f'FE: {z_FE[0:5]}')
