import numpy as np
import matplotlib.pyplot as plt

# data to plot
x_FE = np.loadtxt('x_FE.txt')
y_FE = np.loadtxt('y_FE.txt')
z_FE = np.loadtxt('z_FE.txt')

x_rk4 = np.loadtxt('x_RK4.txt')
y_rk4 = np.loadtxt('y_RK4.txt')
z_rk4 = np.loadtxt('z_RK4.txt')

x2_rk4 = np.loadtxt('x2_RK4.txt')
y2_rk4 = np.loadtxt('y2_RK4.txt')
z2_rk4  = np.loadtxt('z2_RK4.txt')

vx_1 = np.loadtxt('v_x_RK4.txt')
vy_1 = np.loadtxt('v_y_RK4.txt')
vz_1 = np.loadtxt('v_z_RK4.txt')

vx_2 = np.loadtxt('v_x2_RK4.txt')
vy_2 = np.loadtxt('v_y2_RK4.txt')
vz_2 = np.loadtxt('v_z2_RK4.txt')

x_no = np.loadtxt('x_RK4_no.txt')
y_no = np.loadtxt('y_RK4_no.txt')
z_no = np.loadtxt('z_RK4_no.txt')

vx_no = np.loadtxt('v_x_RK4_no.txt')
vy_no = np.loadtxt('v_y_RK4_no.txt')
vz_no = np.loadtxt('v_z_RK4_no.txt')

x2_no = np.loadtxt('x2_RK4_no.txt')
y2_no = np.loadtxt('y2_RK4_no.txt')
z2_no = np.loadtxt('z2_RK4_no.txt')

vx2_no = np.loadtxt('v_x2_RK4_no.txt')
vy2_no = np.loadtxt('v_y2_RK4_no.txt')
vz2_no = np.loadtxt('v_z2_RK4_no.txt')

xf = np.loadtxt('xf.txt')
yf = np.loadtxt('yf.txt')
zf = np.loadtxt('zf.txt')


time = np.linspace(0, 50, 50000)

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

circle1 = plt.Circle((0, 0), R_plus, color='k', fill=False, label= "Lower bound")
circle2 = plt.Circle((0, 0), R_minus, color='g', fill=False,  label="Upped bound")

r_1 = np.sqrt(x_rk4**2 + y_rk4**2 + z_rk4**2)
r_2 = np.sqrt(x2_rk4**2 + y2_rk4**2 + z2_rk4**2)
r = np.abs(r_1 - r_2)


########################################
#used figures report in

################harm oscillator z-t
omega_z = 0.693569
plt.figure()
plt.xlabel('t [$\mu$s]')
plt.ylabel('z [$\mu$m]')
#plt.plot(time, z, color = 'g', label='Analytical')
plt.plot(time, z_no, label='RK4', color = 'k')
plt.vlines(0, -20, 20, ls = '--', color = 'r')
plt.vlines(2*np.pi/omega_z, -20, 20, ls = '--', color = 'r')
plt.vlines(4*np.pi/omega_z, -20, 20, ls = '--', color = 'r')
plt.vlines(6*np.pi/omega_z, -20, 20, ls = '--', color = 'r')
plt.vlines(8*np.pi/omega_z, -20, 20, ls = '--', color = 'r')
plt.vlines(10*np.pi/omega_z, -20, 20, ls = '--', color = 'r')
plt.axis('equal')
plt.legend()
#plt.savefig('plots/z_t_.pdf')



###########path in 3D
fig = plt.figure()
fig.suptitle('Particle trajectory including interaction')
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x_rk4, y_rk4, z_rk4, s=.2,  label='Particle 1')
ax.scatter(x2_rk4, y2_rk4, z2_rk4,  s=.2, label='Particle 2')
ax.set_xlabel('x [$\mu$m]')
ax.set_ylabel('y [$\mu$m]')
ax.set_zlabel('z [$\mu$m]')
ax.legend()
#plt.savefig('3d_int.pdf')



fig = plt.figure()
fig.suptitle('Particle trajectory without interaction')
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x_no, y_no, z_no, s=.2,  label='particle 1')
ax.scatter(x2_no, y2_no, z2_no, s=.2,  label='particle 2')
ax.set_xlabel('x [$\mu$m]')
ax.set_ylabel('y [$\mu$m]')
ax.set_zlabel('z [$\mu$m]')
ax.legend()
#plt.savefig('3d_noint.pdf')



####xy plane with and without interaction
fig = plt.figure()
ax = plt.axes()
plt.title('Including particle interactions')
ax.scatter(x_rk4, y_rk4, s=.2, label='Particle1')
ax.scatter(x2_rk4, y2_rk4, s=.2, label='Particle2')
ax.add_patch(circle1)
ax.add_patch(circle2)
ax.set_xlabel('x [$\mu$m]')
ax.set_ylabel('y [$\mu$m]')
plt.axis('equal')
plt.legend()
#plt.savefig('xy_plane_int.pdf')



fig = plt.figure()
ax = plt.axes()
plt.title('No particle interactions')
ax.scatter(x_no, y_no, s=.2, label='Particle1')
ax.scatter(x2_no, y2_no, s=.2, label='Particle2')
#ax.add_patch(circle1)  #the circles must be added to one figure at the time
#ax.add_patch(circle2)
ax.set_xlabel('x [$\mu$m]')
ax.set_ylabel('y [$\mu$m]')
plt.legend()
plt.axis('equal')
#plt.savefig('xy_plane_noint.pdf')


####xz plane
fig = plt.figure()
ax = plt.axes()
plt.title('Including particle interactions')
ax.scatter(x_rk4, z_rk4, s=.2, label='Particle1')
ax.scatter(x2_rk4, z2_rk4, s=.2, label='Particle2')
ax.set_xlabel('x [$\mu$m]')
ax.set_ylabel('z [$\mu$m]')
plt.axis('equal')
plt.legend()
#plt.savefig('xz_plane_int.pdf')


fig = plt.figure()
ax = plt.axes()
plt.title('No particle interactions')
ax.scatter(x_no, z_no, s = .2, label='Particle1')
ax.scatter(x2_no, z2_no,  s = .2, label='Particle2')
ax.set_xlabel('x [$\mu$m]')
ax.set_ylabel('z [$\mu$m]')
plt.legend()
plt.axis('equal')
#plt.savefig('xz_plane_noint.pdf')



###########phase space x
fig, [ax1, ax2] = plt.subplots(2, 1)
fig.suptitle('x - Phase space plot with interaction')
ax1.scatter(x_rk4, vx_1, c = r , s=.2, cmap = 'bwr_r')
ax2.scatter(x2_rk4, vx_2, c = r, s=.2, cmap = 'bwr_r')
ax1.set_title('Particle 1')
ax2.set_title('Particle 2')
ax1.set_xlabel('x [$\mu$m]')
ax1.set_ylabel('$\dot{x}$  [$\mu$m/s]')
ax2.set_xlabel('x [$\mu$m]')
ax2.set_ylabel('$\dot{x}$  [$\mu$m/s]')
ax1.axis('equal')
ax2.axis('equal')
plt.tight_layout()
#plt.savefig('x_phase_int.pdf')


plt.figure()
plt.title('x - Phase space plot without interaction')
plt.scatter(x_no, vx_no, s=.2,  label='Particle1')
plt.scatter(x2_no, vx2_no, s=.2,  label='Particle2')
plt.xlabel('x [$\mu$m]')
plt.ylabel('$\dot{x}$  [$\mu$m/s]')
plt.legend()
plt.axis('equal')
#plt.savefig('x_phase_noint.pdf')


####phase space z
fig, [ax1, ax2] = plt.subplots(1, 2)
fig.suptitle('z - Phase space plot with interaction')
phase_1 = ax1.scatter(z_rk4, vz_1,  s=.2, c = r, label='Particle 1', cmap = 'bwr_r')
ax1.text(0, 20, 'Particle 1')
phase_2 = ax2.scatter(z2_rk4, vz_2,  s=.2, c = r, label='Particle2', cmap = 'bwr_r')
ax2.text(0, 20, 'Particle 2')
ax1.set_xlabel('z [$\mu$m]')
ax1.set_ylabel('$\dot{z}$ [$\mu$m/s]')
ax2.set_xlabel('z [$\mu$m]')
ax2.set_ylabel('$\dot{z}$ [$\mu$m/s]')
ax1.axis('equal')
ax2.axis('equal')
fig.colorbar(phase_2)
plt.tight_layout()
#plt.savefig('z_phase_int.pdf')


plt.figure()
plt.title('z- Phase space plot without interaction')
plt.scatter(z_no, vz_no,  s=.2, label='Particle1')
plt.scatter(z2_no, vz2_no, s=.2,  label='Particle2')
plt.xlabel('z [$\mu$m]')
plt.ylabel('$\dot{z}$ [$\mu$m/s]')
plt.legend()
plt.axis('equal')
#plt.savefig('z_phase_noint.pdf')


###3d phase space x, y, v_x

fig = plt.figure()
fig.suptitle('Particle trajectory without interaction')
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x_no, y_no, vx_no, s=.2,  label='particle 1')
ax.set_xlabel('x [$\mu$m]')
ax.set_ylabel('y [$\mu$m]')
ax.set_zlabel('$\dot{x}$ [$\mu$m/s]')
ax.legend()
#plt.savefig('3d_phase.pdf')



plt.show()
########################other figures

# plt.figure()
# plt.plot(time, z_FE, label='FE')
# plt.plot(time, z, label='Analytical')
# plt.axis('equal')
# plt.legend()

# plt.figure()
# plt.plot(time, x_rk4, label='RK4')
# #plt.plot(time, x, label='Analytical')
# plt.axis('equal')
# plt.legend()
# plt.show()

# plt.figure()
# plt.plot(time, y_rk4, label='RK4')
# plt.plot(time, y, label='Analytical')
# plt.axis('equal')
# plt.legend()

# fig = plt.figure()
# ax = plt.axes()
# ax.plot(time, z_rk4)
# ax.plot(time, z2_rk4)
# ax.plot(time, r, ls = '--')
# ax.plot(time[31000], r[31000], 'x')
# plt.show()
