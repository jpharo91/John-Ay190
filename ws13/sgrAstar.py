#!/usr/bin/env python

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as mpl3d

# global constants
ggrav = 6.67e-8
msun  = 1.99e33
seconds_per_year = 24.*3600*365 # roughly
cm_per_pc = 3.1e18
distance_to_sgrAstar = 8e3 * cm_per_pc

# system parameters
initial_data_file = "sun_earth.asc"
other_data_file = "sgrAstar.asc"
distance_unit_to_cm = 1.197e17
time_unit_to_s = seconds_per_year
mass_unit_to_g = msun
Nsteps = 10**4
t0 = 0
t1 = 100 * seconds_per_year
dt = (t1-t0)/Nsteps

final_data_file = "final_positions_agrA.asc"

def NbodyRHS(u,mass,time):
    x = list(u[:,0])
    y = list(u[:,1])
    z = list(u[:,2])
    vx = list(u[:,3])
    vy = list(u[:,4])
    vz = list(u[:,5])
    ax = []
    ay = []
    az = []

    for i, particle in enumerate(mass):
        #x_i = x.pop(i)
        #y_i = y.pop(i)
        #z_i = z.pop(i)
        ax_i = 0
        ay_i = 0
        az_i = 0
        for j in range(len(x)-1):
            if j != i:
                r_ij2 = np.abs(x[i]-x[j])**2 + np.abs(y[i]-y[j])**2 + np.abs(z[i]-z[j])**2
                dr = np.array([x[i] - x[j], y[i] - y[j], z[i] - z[j]])
                r_hat = dr/np.linalg.norm(dr)
                ax_i += -ggrav*mass[j]*r_hat[0]/r_ij2
                ay_i += -ggrav*mass[j]*r_hat[1]/r_ij2
                az_i += -ggrav*mass[j]*r_hat[2]/r_ij2
                #print ax_i, ay_i, az_i
        ax.append(ax_i)
        ay.append(ay_i)
        az.append(az_i)
        #x.insert(i, x_i)
        #y.insert(i, y_i)
        #z.insert(i, z_i)

    return np.array((vx, vy, vz, ax, ay, az)).transpose()

def NbodyRK4(u,mass,time,dt):
    k1 = np.zeros(3)
    k2 = np.zeros(3)
    k3 = np.zeros(3)
    k4 = np.zeros(3)

    k1 = dt*NbodyRHS(u,mass,time)
    k2 = dt*NbodyRHS(u + 0.5 * k1, mass, time + dt / 2.0)
    k3 = dt*NbodyRHS(u + 0.5 * k2, mass, time + dt / 2.0)
    k4 = dt*NbodyRHS(u + k3, mass, time + dt)
    #print k1, k2, k3, k4

    return u + (1.0/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4)

def TotalEnergy(u,mass,time):
    x = list(u[:,0])
    y = list(u[:,1])
    z = list(u[:,2])
    vx = list(u[:,3])
    vy = list(u[:,4])
    vz = list(u[:,5])

    E_kin = sum(0.5*mass*np.abs(vx)*np.abs(vx)) + sum(0.5*mass*np.abs(vy)*np.abs(vy)) + sum(0.5*mass*np.abs(vz)*np.abs(vz))

    E_pot = 0

    for i, particle in enumerate(mass):

        for j in range(len(x)-1):

            if j != i:

                r_ij = np.abs(x[i]-x[j])**2 + np.abs(y[i]-y[j])**2 + np.abs(z[i]-z[j])
                E_pot += -ggrav*particle*mass[j]/r_ij

    E_tot = E_kin + E_pot
    return E_tot

# main program
plt.ion()

(x,y,z,vx,vy,vz,mass) = np.loadtxt(other_data_file, unpack = True)

# convert from unitis in initial data file to cgs
x *= distance_unit_to_cm
y *= distance_unit_to_cm
z *= distance_unit_to_cm
vx *= distance_unit_to_cm / time_unit_to_s
vy *= distance_unit_to_cm / time_unit_to_s
vz *= distance_unit_to_cm / time_unit_to_s
mass *= mass_unit_to_g

xmin = np.amin(x)
xmax = np.amax(x)
ymin = np.amin(y)
ymax = np.amax(y)
zmin = np.amin(z)
zmax = np.amax(z)
rmax = 2.5*max(abs(xmin),abs(xmax),abs(ymin),abs(ymax),abs(zmin),abs(zmax))

# use a single state vector to simplify the ODE code
# indices:
# u[:,0] = x
# u[:,1] = y
# u[:,2] = z
# u[:,3] = vx
# u[:,4] = vy
# u[:,5] = vz
u = np.array((x,y,z,vx,vy,vz)).transpose()
energies = []

for it in range(0, Nsteps):
    time = t0 + it * dt
    u = NbodyRK4(u,mass,time,dt)
    if it % max(1,Nsteps/100) == 0:
      print "it = %d, time = %g years, energy = %g" % \
            (it, time / seconds_per_year, TotalEnergy(u,mass,time))
      energies.append(TotalEnergy(u,mass,time))
      plt.clf()
      fig = plt.gcf()
      ax = mpl3d.Axes3D(fig)
      ax.scatter(u[:,0],u[:,1],u[:,2])
      ax.set_xlim((-rmax,rmax))
      ax.set_ylim((-rmax,rmax))
      ax.set_zlim((-rmax,rmax))
      plt.draw()

energy_dif = max(energies) - min(energies)
energies2 = np.abs(energies - energies[-1])/energies
years = range(len(energies2))
plt.plot(years, energies2, "b")
plt.savefig("Energy-sgrA.pdf")
print energy_dif
# output result
file_header = "1:x 2:y 3:z 4:vx 5:vy 6:vz 7:mass"
np.savetxt(final_data_file, u)#, header=file_header)
