#https://cdn.fbsbx.com/v/t59.2708-21/285168815_715364306401213_6117355139697281868_n.py/dm4bem.py?_nc_cat=102&ccb=1-7&_nc_sid=0cab14&_nc_ohc=9CWZYxLoeecAX-GhwXt&_nc_ht=cdn.fbsbx.com&oh=03_AVJBYYYMBiZJStqEmTTxSDK9VDM7Z5-oZy1xYfHvb26VFQ&oe=62A16AD9&dl=1# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.

@author:ambrebanzet
"""
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import dm4bem

Kp = 1e4    

S_wall_W=6*3
S_wall_S=2*2*3
S_wall_E=5*3
S_wall_N=1*3

S_door=1*3

S_window_1=2*3
S_window_2=2*3

S_toilet_door=1*3
S_toilet_room_wall=3*(2+4)

S_toilet_wall_E=3*3
S_toilet_wall_N=4*3

S_c=S_wall_W+S_wall_S+S_wall_E+S_wall_N+S_toilet_wall_E+S_toilet_wall_N
S_g=S_window_1+S_window_2

wall = {'lambda': [1.4, 3e-2, 1.2],  
        'rho*c': [2000000, 26000, 1000000],        
        'w': [0.15, 0.08, 0.02] }
        
wall = pd.DataFrame(wall, index=['Concrete', 'Insulation', 'Glass'])

hout=20
hin=4

rho_air=.2                      # kg/m³
c_air=1000               # J/kg.K

Va = (6*8)*3          # m³ volume of air
ACH = 1             # air changes per hour
Va_dot = ACH * Va / 3600    # m³/s air infiltration
#matriceG


ε_wLW = 0.9     # long wave wall emmisivity (concrete)
α_wSW = 0.2     # absortivity white surface
ε_gLW = 0.9     # long wave glass emmisivity (glass pyrex)
τ_gSW = 0.83    # short wave glass transmitance (glass)
α_gSW = 0.1     # short wave glass absortivity

σ = 5.67e-8     # W/m².K⁴ Stefan-Bolzmann constant

Fwg = 1 / 5     # view factor wall - glass

Tm = 20 + 273   # mean temp for radiative exchange

GLW1 = ε_wLW / (1 - ε_wLW) * (S_wall_W+S_wall_S+S_wall_E+S_wall_N+S_toilet_wall_E+S_toilet_wall_N) * 4 * σ * Tm**3
GLW2 = Fwg * (S_wall_W+S_wall_S+S_wall_E+S_wall_N+S_toilet_wall_E+S_toilet_wall_N) * 4 * σ * Tm**3
GLW3 = ε_gLW / (1 - ε_gLW) * (S_wall_W+S_wall_S+S_wall_E+S_wall_N+S_toilet_wall_E+S_toilet_wall_N)* 4 * σ * Tm**3
# long-wave exg. wall-glass
GLW = 1 / (1 / GLW1 + 1 / GLW2 + 1 / GLW3)

Gv = Va_dot * rho_air * c_air


g=np.zeros(49)
g[0]=hout*S_wall_W
g[1]=wall['lambda'][0]/wall['w'][0]*S_wall_W
g[2]=g[1]
g[3]=wall['lambda'][1]/wall['w'][1]*S_wall_W
g[4]=g[3]
g[5]=hin*S_wall_W

g[6]=hout*S_wall_S
g[7]=wall['lambda'][0]/wall['w'][0]*S_wall_S
g[8]=g[7]
g[9]=wall['lambda'][1]/wall['w'][1]*S_wall_S
g[10]=g[9]
g[11]=hin*S_wall_S

g[12]=hin*S_wall_W
g[13]=g[12]
g[14]=wall['lambda'][1]/wall['w'][1]*S_wall_N
g[15]=g[14]
g[16]=hin*S_wall_W

g[17]=hout*S_wall_E
g[18]=wall['lambda'][0]/wall['w'][0]*S_wall_E
g[19]=g[18]
g[20]=wall['lambda'][1]/wall['w'][1]*S_wall_E
g[21]=g[20]
g[22]=hin*S_wall_E

g[23]=S_door*hout

g[24]=2 #VALEUR A COMPLETER 

g[25]=hout*S_window_1
g[26]=hin*S_window_1

g[27]=hout*S_window_2
g[28]=hin*S_window_2

g[29]=S_toilet_door*hin

g[30]=hin*S_toilet_room_wall
g[31]=(wall['lambda'][0]/wall['w'][0])*S_toilet_room_wall
g[32]=g[31]
g[33]=(wall['lambda'][0]/wall['w'][0])*S_toilet_room_wall
g[34]=g[33]
g[35]=hin*S_toilet_room_wall

g[36]=hin*S_toilet_wall_E
g[37]=(wall['lambda'][1]/wall['w'][1])*S_toilet_wall_E
g[38]=g[37]
g[39]=(wall['lambda'][0]/wall['w'][0])*S_toilet_wall_E
g[40]=g[39]
g[41]=hout*S_toilet_wall_E

g[42]=hin*S_toilet_wall_N
g[43]=(wall['lambda'][1]/wall['w'][1])*S_toilet_wall_N
g[44]=g[43]
g[45]=(wall['lambda'][0]/wall['w'][0])*S_toilet_wall_N
g[46]=g[45]
g[47]=hout*S_toilet_wall_N

g[48]=GLW+Gv


G=np.diag(g)

#matrice C 


c=np.zeros(39)
c[1]=S_wall_W*wall['w'][1]*wall['rho*c'][1]
c[3]=S_wall_W*wall['w'][0]*wall['rho*c'][0]
c[6]=S_wall_S*wall['w'][1]*wall['rho*c'][1]
c[8]=S_wall_S*wall['w'][0]*wall['rho*c'][0]
c[10]=S_wall_N*wall['w'][1]*wall['rho*c'][1]
c[12]=S_wall_N*wall['w'][0]*wall['rho*c'][0]
c[15]=S_wall_E*wall['w'][1]*wall['rho*c'][1]
c[17]=S_wall_E*wall['w'][0]*wall['rho*c'][0]
c[22]=S_toilet_room_wall*wall['w'][0]*wall['rho*c'][0]
c[25]=S_toilet_room_wall*wall['w'][0]*wall['rho*c'][0]
c[28]=S_toilet_wall_E*wall['w'][0]*wall['rho*c'][0]
#c[29]=S_toilet_wall_E*wall['w'][1]*wall['rho*c'][1]
c[33]=S_toilet_wall_N*wall['w'][0]*wall['rho*c'][0]
#c[35]=S_toilet_wall_N*wall['w'][1]*wall['rho*c'][1]
c[37]=50000   
c[38]=50000   


C=np.diag(c)


#matrice A
A=np.zeros((49,39))

A[ 0, 0]=1
A[ 1, 0] = - 1
A[ 1, 1] = 1
A[ 2, 1] = - 1
A[ 2, 2] = 1
A[ 3, 2] = - 1 
A[ 3, 3] = 1
A[ 4, 3] = - 1
A[ 4, 4] = 1
A[ 5, 38] = 1


A[ 6, 5] = 1
A[ 7, 5] = - 1
A[ 7, 6] = 1
A[ 8, 6] = - 1
A[ 8, 7] = 1
A[ 9, 7] = - 1
A[ 9, 8] = 1
A[ 10, 8] = - 1
A[ 10, 9] = 1
A[ 11, 9] = - 1
A[ 11, 38] = 1


A[ 12, 10] = 1
A[ 13, 10] = - 1
A[ 13, 11] = 1 
A[ 14, 11] = - 1
A[ 14, 12] = 1
A[ 15, 12] = - 1
A[ 15, 13] = 1 
A[ 16, 13] = - 1
A[ 16, 38] = 1


A[ 17, 14] = 1
A[ 18, 14] = - 1
A[ 18, 15] = 1
A[ 19, 15] = - 1
A[ 19, 16] = 1
A[ 20, 16] = - 1
A[ 20, 17] = 1
A[ 21, 17] = - 1
A[ 21, 18] = 1
A[ 22, 18] = - 1
A[ 22, 38] = 1


A[ 23, 38] = 1


A[ 24, 38] = 1


A[ 25, 19] = 1
A[ 26, 19] = - 1 
A[ 26, ] = 1


A[ 27, 20] = 1 
A[ 28, 20] = - 1
A[ 28, 38] = 1


A[ 29, 37] = - 1
A[ 29, 38] = 1


A[ 30, 38] = 1
A[ 30, 21] = - 1
A[ 31, 21] = 1 
A[ 31, 22] = - 1
A[ 32, 22] = 1
A[ 32, 23] = - 1
A[ 33, 23] = 1
A[ 33, 24] = - 1
A[ 34, 24] = 1 
A[ 34, 25] = - 1
A[ 35, 25] = 1 
A[ 35, 37] = - 1


A[ 48, 37] = - 1 


A[ 36, 37] = 1
A[ 36, 27] = - 1
A[ 37, 27] = 1 
A[ 37, 28] = - 1
A[ 38, 28] = 1 
A[ 38, 29] = - 1
A[ 39, 29] = 1 
A[ 39, 30] = - 1
A[ 40, 30] = 1 
A[ 40, 31] = - 1
A[ 41, 31] = 1 


A[ 42, 37] = 1
A[ 42, 32] = - 1 
A[ 43, 32] = 1
A[ 43, 33] = - 1 
A[ 44, 33] = 1
A[ 44, 34] = - 1 
A[ 45, 34] = 1
A[ 45, 35] = - 1
A[ 46, 35] = 1
A[ 46, 36] = - 1
A[ 47, 36] = 1

#vecteur b

b=np.zeros((49,1))
b[0]=20
b[6]=20
b[12]=20
b[17]=20
b[23]=20
b[24]=20
b[25]=20
b[27]=20
b[41]=20
b[43]=20

#vecteur f

f=np.zeros((39,1))


#state-space model

ytc = np.linalg.inv((A.transpose()).dot(G).dot(A)).dot((A.transpose()).dot(G).dot(b)+ f)

#steady state

b=np.zeros(49)
b[[0, 5, 12, 17, 23, 24, 25, 27, 41, 47, 48]] = 10 + np.array([0, 50, 120, 170, 230, 240, 250, 270, 410, 470, 480])
f = np.zeros(39)
f[[0, 4, 5, 9, 13, 14, 18, 19, 20, 21, 27, 31, 32, 36, 37, 38]] = 100 + np.array([0, 400, 500, 900, 1300, 1400, 1800, 1900, 2000, 2100, 2700, 3100, 3200, 3600, 3700, 3800])
y = np.ones(39)


ytc = np.linalg.inv(A.T @ G @ A) @ (A.T @ G @ b + f)

u = np.hstack([b[np.nonzero(b)], f[np.nonzero(f)]])
[As, Bs, Cs, Ds] = dm4bem.tc2ss(A, G, b, C, f, y)


yss = (-Cs.dot(np.linalg.inv(As)).dot(Bs) + Ds).dot( u)

print(np.array_str(yss, precision=3, suppress_small=True))
print(np.array_str(ytc, precision=3, suppress_small=True))
print(f'Max error in steady-state between thermal circuit and state-space:\
 {max(abs(yss - ytc)):.2e}')

#dynamic simulation

b=np.zeros(49)
b[[0, 5, 12, 17, 23, 24, 25, 27, 41, 47, 48]] =1
f = np.zeros(39)
f[[0, 4, 5, 9, 13, 14, 18, 19, 20, 21, 27, 31, 32, 36, 37, 38]]=1

y=np.zeros(39)
y[[38]]=1

[As, Bs, Cs, Ds] = dm4bem.tc2ss(A, G, b, C, f, y)


#time step

dtmax = min(-2. / np.linalg.eig(As)[0])
print(f'Maximum time step: {dtmax:.2f} s')

# dt = 5
#dt = 360
dt=150
duration = 3600 * 24 * 2 
n = int(np.floor(duration / dt))
t = np.arange(0, n * dt, dt)

# Vectors of state and input (in time)
n_tC = As.shape[0]              # no of state variables (temps with capacity)
# u = [To To To Tsp Phio Phii Qaux Phia]
u = np.zeros([27, n])
u[0:10, :] = np.ones([10, n])

# u=np.ones((len(np.hstack([b[np.nonzero(b)],f[np.nonzero(f)]])),n))
# for i in range (n):
#     u[:,i]=np.hstack([b[np.nonzero(b)],f[np.nonzero(f)]])

temp_exp = np.zeros([n_tC, t.shape[0]])
temp_imp = np.zeros([n_tC, t.shape[0]])

I = np.eye(n_tC)
for k in range(n - 1):
    temp_exp[:, k + 1] = (I + dt * As) @\
        temp_exp[:, k] + dt * Bs @ u[:, k]
    temp_imp[:, k + 1] = np.linalg.inv(I - dt * As) @\
        (temp_imp[:, k] + dt * Bs @ u[:, k])

y_exp = Cs @ temp_exp + Ds @  u
y_imp = Cs @ temp_imp + Ds @  u

fig, ax = plt.subplots()
ax.plot(t / 3600, y_exp.T,'r', t / 3600, y_imp.T, 'g')
ax.set(xlabel='Time [h]',
        ylabel='$T_i$ [°C]',
        title='Step input: To = 1°C')
plt.show()


b = np.zeros(49)
b[[0, 5, 12, 17, 23, 24, 25, 27, 41, 47, 48]] = 1
f = np.zeros(39)


ytc = np.linalg.inv(A.T @ G @ A) @ (A.T @ G @ b + f)
print('Steady-state indoor temperature obtained with:')
print(f'- DAE model: {ytc[6]:.4f} °C')
print(f'- response to step input:{float(y_exp[:, -2]):.4f} °C')

#simulation with weather data

filename = 'FRA_Lyon.074810_IWEC.epw'
start_date = '2000-01-01 12:00:00'
end_date = '2000-01-30 18:00:00'

# Read weather data from Energyplus .epw file
[data, meta] = dm4bem.read_epw(filename, coerce_year=None)
weather = data[["temp_air", "dir_n_rad", "dif_h_rad"]]
del data
weather.index = weather.index.map(lambda t: t.replace(year=2000))
weather = weather[(weather.index >= start_date) & (
    weather.index < end_date)]

surface_orientation = {'slope': 90,
                       'azimuth': 0,
                       'latitude': 45}
albedo = 0.2
rad_surf1 = dm4bem.sol_rad_tilt_surf(weather, surface_orientation, albedo)
rad_surf1['Φt1'] = rad_surf1.sum(axis=1)

data = pd.concat([weather['temp_air'], rad_surf1['Φt1']], axis=1)
data = data.resample(str(dt) + 'S').interpolate(method='linear')
data = data.rename(columns={'temp_air': 'To'})

data['Ti'] = 20 * np.ones(data.shape[0])
data['Qa'] = 0 * np.ones(data.shape[0])

t = dt * np.arange(data.shape[0])

u = pd.concat([data['To'], data['To'], data['Ti'],
               α_wSW * S_c * data['Φt1'],
               τ_gSW * α_wSW * S_g * data['Φt1'],
               data['Qa'],
               α_gSW * S_g * data['Φt1']], axis=1)

temp_exp = 20 * np.ones([As.shape[0], u.shape[0]])

for k in range(u.shape[0] - 1):
    temp_exp[:, k + 1] = (I + dt * As) @ temp_exp[:, k]\
        + dt * Bs @ u.iloc[k, :]
        
y_exp = Cs @ temp_exp + Ds @ u.to_numpy().T
q_HVAC = Kp * (data['Ti'] - y_exp[0, :])

fig, axs = plt.subplots(2, 1)
# plot indoor and outdoor temperature
axs[0].plot(t / 3600, y_exp[0, :], label='$T_{indoor}$')
axs[0].plot(t / 3600, data['To'], label='$T_{outdoor}$')
axs[0].set(xlabel='Time [h]',
           ylabel='Temperatures [°C]',
           title='Simulation for weather')
axs[0].legend(loc='upper right')

# plot total solar radiation and HVAC heat flow
axs[1].plot(t / 3600,  q_HVAC, label='$q_{HVAC}$')
axs[1].plot(t / 3600, data['Φt1'], label='$Φ_{total}$')
axs[1].set(xlabel='Time [h]',
           ylabel='Heat flows [W]')
axs[1].legend(loc='upper right')

fig.tight_layout()