"""
Created on Sat Nov 26 15:07:24 2022

@author: Hiroki
"""
#%%
'''
GENERATE ANIMATION:
# Ball radius = 1.0
# Number of balls = 20
# Radius of container = 20.0
# Number of collisions = 1000
# Initial positins and velocities randomised.
# mean velocity = 0.0
# Velocity width = 2.0
'''
from Simulation import Simulation

sim = Simulation(
                  n_balls = 20, \
                  random_positions = True, \
                  random_velocities = True, \
                 )

sim.run(
         num_frames = 1000, \
         animate = True 
       )


#%%
'''
# Task 9 - Plot histograms for distances of balls from the origin and \
    distances between pairs of balls:
#===============================================================================
# Ball radius = 1.0
# Number of balls = 20
# Radius of container = 20.0
# Number of collisions = 1000
# Initial positins and velocities randomised.
# mean velocity = 0.0
# Velocity width = 2.0
# Animation turned off.
'''
from Simulation import Simulation

sim1 = Simulation(
                  n_balls = 20, \
                  random_positions = True, \
                  random_velocities = True, \
                 )
sim1.run(
         num_frames = 1000, \
         animate = False, \
         dist_fig = True, \
         ke_sys_fig = False, \
         mom_sys_fig = False, \
         temp_fig = False, \
         pressure_fig = False
        )


#%%
'''
# Tasks 10, 11 Part 1 - Plots for checking conservation laws.
#===============================================================================
# Ball radius = 1.0
# Container radius = 20.0
# number of balls = 20
# 1000 collisions. 
# Initial positins and velocities randomised.
# Velocity width = 2.0
# Animation turned off.

# Quantities conserved in the simulation are:
    - system's total kinetic energy.
    - pressure exerted by the gas particles on the container.
    - magnitude of the system's total momentum.
    - temperature of the gas.

# Total kinetic energy vs. time graph outputs a straight horizontal line.
  It suggests that this quantity is conserved in the simulation
  (Since the gas particles collide elastically).
# The pressure vs. time graph starts off with a large fluctuation that \
  quickly converges to a mean equilibrium pressure value. This quanitty is \
  also conserved but take more time than the total kinetic energy to settle to \
  an equilibrium value.
# Temperature vs. time graph outputs a straight horizontal line, which is \
  expected since it is proportional to the gas particles' total KE. This \
  quantity is therefore conserved in the simulation.
# As for the magnitude of the total momentum, the graph fluctuates slightly \
  over time but seemingly around a mean equilibrium value. However, unlike \
  with pressure, the fluctuation seems to persist for as long as the \
  simulation is run. \
# This fluctuation is most likely due to neglecting the momentum of the \
  container (which is negligibly small due to it being practically stationary \
  with a mass of 1.0e99).
# Fluctuations in the momentum vs. time graph may smoothen-out if the number \
  of balls and the number of collisions increases.
'''
from Simulation import Simulation

sim3 = Simulation(
                  n_balls = 20, \
                  random_positions = True, \
                  random_velocities = True, \
                  velocity_width = 2.0
                 )
sim3.run(
         num_frames = 1000, \
         animate = False, \
         dist_fig = False, \
         ke_sys_fig = False, \
         mom_sys_fig = False, \
         temp_fig = False, \
         pressure_fig = True
        )

'''
Outputs:
========
mean total KE =  66.23950257391547 +/- 3.040923684938172e-17
mean total momentum =  45.771551217524575 +/- 0.0016791766435533828
mean temperature =  3.3119751286957726 +/- 1.5376007716398207e-18
mean pressure =  0.08748175895979828 +/- 4.615775074098776e-06

'''

#%%
'''
# Task 11 Part 2- Plot Pressure vs. Temperature graph with a linear fit.
#===============================================================================
# Ball radius = 1.0
# Container radius = 50.0
# number of balls = 100
# 20 data points.
# 1000 collisions per data point.
# Animation turned off.
'''
import Figure as fig
radius_ball = [1.0]
radius_container = 50.0
n_balls = 100
velocity_widths = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, \
                   6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]
n_collisions = 1000

fig.pressure_vs_temperature(
                            radius_ball_list = radius_ball, \
                            radius_cont = radius_container, \
                            num = n_balls, \
                            v_width_list = velocity_widths, \
                            collisions = n_collisions
                            )
'''
Outputs:
=========
slope =  0.01282779755077852 +/- 0.00037579202742948095
y-int =  0.01462805583266944 +/- 0.018283960809822572
'''

#%%
'''
# Task 11 Part 3 - Plot Pressure and Temperature vs. Container Area. 
#===============================================================================
# Ball radius = 1.0
# Number of balls = 50
# Velocity width = 5.0
# 9 different container radii.
# 1000 collisions per data point. 
# Animation turned off.

# An inverse relationship is obtained between pressure and area of container.
# Curve fitting to determine the mean temperature.

# A rouhgly constant relationship is obtained between temperature and area of \
  # container.
# This suggests that when the initial velocity distribution of the gas is kept \
  # constant, the temperature of the gas remins rougly constant as the area of \
  # the container increases. 
# Temperature remains constant as it is an intensive quantity.
'''
import Figure as fig

radius_container = [20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]
n_balls = 50
velocity_width = 5.0
n_collisions = 1000

fig.pressure_temperature_vs_area(
                     radius_container_list = radius_container, \
                     num = n_balls, \
                     v_width = velocity_width, \
                     collisions = n_collisions, \
                     pressure_area = True, \
                     temp_area = True
                     )
'''
Output:
=======
Estimated Mean Temp =  32.45396873827356 +/- 0.7774495532499253
Estimated Pressure =  -0.00015137410768252174 +/- 26986.930208080805
temperature =  [23.249584584097313, 23.244220165654582, 22.95346769654953, \
                29.981102496205093, 19.10006581493321, 23.578025129923567, \
                23.6426271775741, 29.6603059582236, 32.33923261426634, \
                26.46392293783292, 22.43681215050276, 23.613995264609663, \
                27.007636969246555, 18.181181244597948]
pressure =  [1.261199964374068, 0.48019247629559164, 0.2553706873190667, \
             0.1971146444763413, 0.08844198367156378, 0.08104853915708397, \
             0.06007344392393825, 0.06049887663654496, 0.05443891632117944, \
             0.037125330150968264, 0.02528041236308012, 0.022320685166069047, \
             0.022993070280420667, 0.013001812425767067]
'''
'''
# Task 11 Part 4 - Plot Pressure and Temperature vs Number of balls. 
#===============================================================================
# Ball radius = 1.0
# Container radius = 50.0
# Number of balls = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
# Velocity width = 5.0
# 1000 collisions per data point. 
# Animation turned off.

# Pressure is proportional to the number of particles. 
# Temperature remains constant as it is an intensive quantity.
'''
import Figure as fig

radius_container = 50.0
n_balls = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
velocity_width = 5.0
n_collisions = 1000

fig.pressure_temperature_vs_n(
                              radius_cont = radius_container, \
                              num_list = n_balls, \
                              v_width = velocity_width, \
                              collisions = n_collisions, \
                              pressure_n = True, \
                              temp_n = True
                              )
'''
Output:
=======
Pressure vs. n:
m =  0.003318996854439371 +/- 0.0002504875187563542
c =  0.0016040231974533414 +/- 0.015542340888180767
             
Temperature vs. n:
m =  -0.13204189046124637 +/- 0.022034615106956375
c =  34.385027081366886 +/- 1.3672119128335951
'''

#%%
'''
# Task 12 - Plot Pressure vs. Temperature graph for different ball radii. 
#===============================================================================
# Container radius = 50.0
# Number of balls = 20
# 6 different ball radii
# 10 data points each
# 1000 collisions per data point. 
# Animation turned off.

# We find that the slope of the graph (pressure vs. temperature) increases as \
  the radius of the balls/particles increases.
'''
import Figure as fig

radii_ball = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
radius_container = 50.0
n_balls = 50
velocity_widths = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5]
n_collisions = 1000

fig.pressure_vs_temperature(
                            radius_ball_list = radii_ball, \
                            radius_cont = radius_container, \
                            num = n_balls, \
                            v_width_list = velocity_widths, \
                            collisions = n_collisions
                           )
'''
Outputs:
=========

ball radius = 0.5:
------------------
slope =  0.006453274497615831 +/- 6.542702317231762e-05
y-int =  0.0011064320470820196 +/- 0.001452720476666322

ball radius = 1.0:
------------------
slope =  0.006776440390582686 +/- 0.0001487014703384431
y-int =  -0.0012589103314321103 +/- 0.003128956721973173

ball radius = 1.5:
------------------
slope =  0.00724337409932803 +/- 0.0001118717185029321
y-int =  0.0006692025059140723 +/- 0.0027276215789328166

ball radius = 2.0:
------------------
slope =  0.007499139776687339 +/- 0.00018475932074092183
y-int =  0.00419688546474 +/- 0.004342297266532854


ball radius = 2.5:
------------------
slope =  0.008812735127118956 +/- 0.00017492760549630084
y-int =  -0.003117649595327962 +/- 0.003966760296085439

ball radius = 3.0:
------------------
slope =  0.010534406218689342 +/- 0.00025636663457909723
y-int =  -0.0027330793618638785 +/- 0.005472072227116614
'''
#%%
'''
Plot of the slope of Pressure vs. Temperature graph against Ball radius:
#===============================================================================
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

gradient = [0.006453274497615831, 0.006776440390582686, 0.00724337409932803, \
            0.007499139776687339, 0.008812735127118956, 0.010534406218689342]
unc_gradient = [6.542702317231762e-05, 0.0001487014703384431, \
                0.0001118717185029321, 0.00018475932074092183, \
                0.00017492760549630084, 0.00025636663457909723]
radii_ball = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]


def linear_fit(x, m, c):
  return m * x + c

params_1, params_cov_1 = curve_fit(
                                  linear_fit, radii_ball[:4], gradient[:4], \
                                  sigma = unc_gradient[:4], \
                                  absolute_sigma=True
                                  )
params_2, params_cov_2 = curve_fit(
                                  linear_fit, radii_ball[3:], gradient[3:], \
                                  sigma = unc_gradient[3:], \
                                  absolute_sigma=True
                                  )
x_1 = np.linspace(radii_ball[0], radii_ball[3], 100)
x_2 = np.linspace(radii_ball[3], radii_ball[5], 100)
plt.errorbar(
            radii_ball, gradient, yerr = unc_gradient, fmt = 'o', \
            ls='', mew=0.6, ms=0.5, capsize=3, label = 'Simulation Data'
            )
plt.plot(x_1, linear_fit(x_1, *params_1), label = 'Shallower Region')
plt.plot(x_2, linear_fit(x_2, *params_2), label = 'Steeper Region')
plt.title('Slope of PT Graph vs. Ball Radius')
plt.xlabel('Ball Radius')
plt.ylabel('Slope of PT Graph')
plt.legend()
plt.savefig('Slope of PT graph vs. ball radius.png', dpi = 1000)
plt.show()
print('shallower slope = ', params_1[0], '+/-', np.sqrt(params_cov_1[0,0]))
print('steeper slope = ', params_2[0], '+/-', np.sqrt(params_cov_2[0,0]))

'''
Output:
=======
shallower slope =  0.0007412126807328149 +/- 9.909100223525292e-05
steeper slope =  0.00297726989402224 +/- 0.0003088212158046258
'''


#%%
'''
# Task 13 - Compare histograms of particle velocities and speeds against \
  Maxwell-Boltzmann Distribution. 
#===============================================================================
# Container radius = 50.0
# Number of balls = 100
# 20 bins
# 2000 collisions. 
'''
import Figure as fig

radius_container = 50.0
n_balls = 100
velocity_width = 10.0
n_collisions = 2000

fig.speed_distribution(
                       radius_cont = radius_container, \
                       num = n_balls, \
                       v_width = velocity_width, \
                       collisions = n_collisions
                      )
'''
Outputs:
========
scale factor =  4235.1148090037 +/- 130.20567909322938
Estimated temperature =  85.33117413230555 +/- 130.20567909322938
Actual temperature =  86.79032286762434 +/- 
'''
# %%
'''
# Task 14 - Van der Waals' Law (deducing contants a and b).
#========================================================================
# Ball radius = 1.0
# Container radius = 50.0
# number of balls = 100
# 20 data points.
# 1000 collisions per data point.
# Animation turned off.
'''
import Figure as fig

radius_container = 50.0
n_balls = 100
velocity_widths = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, \
                   6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]
n_collisions = 1000

fig.van_der_waals(
                  radius_cont = radius_container, \
                  num = n_balls, \
                  v_width_list = velocity_widths, \
                  collisions = n_collisions
                  )
'''
Outputs:
========
b from van der Waals =  9.8063662950968 +/- 0.5696149395122787
'''
