"""
Created on Sun Nov 13 11:41:19 2022

@author: hk621
"""
# File for Unit Testing and Integration Testing of pices of code:

#%%
'''
# Task 2 Test 1 - Initialise Ball class and test time_to_collision, move, and \
    collide methods for single ball with container.
#===============================================================================
'''
from Ball import *

# Ball object
b1 = Ball(
          mass = 1.0, \
          radius = 1.0, \
          position = np.array([0.0, 0.0]), \
          velocity = np.array([1.0, 0.0]), \
          iscontainer = False
         )
# Container object
b2 = Ball(
          mass = 1.0e99, \
          radius = 10.0, \
          position = np.array([0.0, 0.0]), \
          velocity = np.array([0.0, 0.0]), \
          iscontainer = True
         )

# Calculate time to next collision of ball (b2) and container (b2).
dt_1 = b1.time_to_collision(b2) 
 # Move ball to the new position at time dt.
b1.move(dt_1) 
# Calculate new velocities of ball and container. Velocities of container \
    # should be unchanged (zero).
b1.collide(b2)  
# Check that collide method worked:
b1, b2 
print(b1, b2)

'''
# Output for b1:
Ball Class: mass = 1, radius = 1, position = [9. 0.], velocity = [-1.  0.], \
            is a container = False, collision = 1
# Output for b2:
Ball Class: mass = 1e+99, radius = 10, position = [0. 0.], \
            velocity = [0. 0.], is a container = True, collision = 1
# SUCCESS
'''
#%%
'''
# Task 2 Test 2 - test time_to_collision, move, and collide methodsball for \
    ball with arbitrary initial velocity and container.
#===============================================================================
'''
b3 = Ball(
          mass = 1.0, \
          radius = 1.0, \
          position = np.array([0.0, 0.0]), \
          velocity = np.array([1.0, 2.0]), \
          iscontainer = False
         )

dt_3 = b3.time_to_collision(b2) 
b3.move(dt_3) 
b3.collide(b2)
b3
print(b3)

'''
# Output for b2:
Ball Class: mass = 1.0, radius = 1.0, position = [4.02492236 8.04984472], \
            velocity = [-1. -2.], is a container = False, collision = 1
# SUCCESS
'''
#%%
'''
# Task 2 Test 3 - test time_to_collision, move, and collide methodsball for \
    ball with another ball:
#===============================================================================
'''
b4 = Ball(
          mass = 1.0, \
          radius = 1.0, \
          position = np.array([5.0, 5.0]), \
          velocity = np.array([-1.0, -1.0]), \
          iscontainer = False
         )

b5 = Ball(
          mass = 1.0, \
          radius = 1.0, \
          position = np.array([-5.0, -5.0]), \
          velocity = np.array([1.0, 1.0]), \
          iscontainer = False
         )

dt_4 = b4.time_to_collision(b5) 
b4.move(dt_4) 
b5.move(dt_4)
b4.collide(b5)
b4, b5
print(b4, b5)

'''
# Output for b4:
Ball Class: mass = 1.0, radius = 1.0, position = [0.70710678 0.70710678], \
            velocity = [1. 1.], is a container = False, collision = 1 

# Output for b5:
Ball Class: mass = 1.0, radius = 1.0, position = [-0.70710678 -0.70710678], \
            velocity = [-1. -1.], is a container = False, collision = 1
# SUCCESS
'''



#%%
'''
# Task 4 Test 1- testing the simulation class made in task 3:
#===============================================================================
'''
'''
# Test 1:
# Initialise the Simulation object with container radius == 10, ball \
  radius == 1, ball mass == 1, initial position of ball == [-5, 0], and \
  initial velocity of ball == [1, 0].
'''
from Simulation import *

sim1 = Simulation(
                  radius_container = 10.0, \
                  random_positions = False, \
                  random_velocities = False, \
                  set_position = [-5.0, 0.0], \
                  set_velocity = [1.0, 0.0], \
                 )

# Execute next_collision method:
sim1.next_collision()
sim1._balls[0]
print(sim1._balls[0])

# Run method again:
sim1.next_collision()
sim1._balls[0]
print(sim1._balls[0])

'''
# Output after 1st collision:
Ball Class: mass = 1.0, radius = 1.0, position = [9. 0.], \
            velocity = [-1.  0.], is a container = False, collision = 1
            
# Output after 2nd collision:
Ball Class: mass = 1.0, radius = 1.0, position = [-9.  0.], 
            velocity = [1. 0.], is a container = False, collision = 2

# SUCCESS
'''
#%%
'''
# Task 4 Test 2 - Initialise the Simulation object with initial velocity of \
    ball == [1.0, 2.0] (velocity having x and y components):
#===============================================================================
'''
sim2 = Simulation(
                  radius_container = 10.0, \
                  random_positions = False, \
                  random_velocities = False, \
                  set_position = [-5.0, 0.0], \
                  set_velocity = [1.0, 2.0], \
                 )

#Execute next_collision method:
sim2.next_collision()
sim2._balls[0]
print(sim2._balls[0])

#Run method again:
sim2.next_collision()
sim2._balls[0]
print(sim2._balls[0])

'''
# Output after 1st collision:
Ball Class: mass = 1.0, radius = 1.0, position = [-0.50715016  8.98569968], \
            velocity = [ 1.21869128 -1.87477774], is a container = False, \
            collision = 1

# Output after 2nd collision:
Ball Class: mass = 1.0, radius = 1.0, position = [ 8.00626111 -4.11093456], \
            velocity = [-2.23373685 -0.10207686], is a container = False, \
            collision = 2

# SUCCESS
'''

#%%
'''
# Task 5 - testing animation for a single ball inside a container:
#===============================================================================
# Radius of container = 20.0
# initial position of ball = [-5.0, 0.0]
# initial velocity of ball = [1.0, 0.0]
# SUCCESS

# Animation looks correct. Every frame, the ball's position switches between \
  the maximum and minimum x coordinates within the container (bounces from \
  side to side) as expected.
'''
from Simulation import *

sim3 = Simulation(
                  n_balls = 1, \
                  random_positions = False, \
                  random_velocities = False, \
                  set_position = np.array([-5.0, 0.0]), \
                  set_velocity = np.array([1.0, 0.0])
                 )

sim3.run(
         num_frames = 100, \
         animate = True
        )

#%%
'''
# Task 5 (Test 2) - testing animation for a single ball with arbitrary initial \
  velocity with x and y components.
#===============================================================================
# Radius of container = 20.0
# initial position of ball = [-5.0, 0.0]
# initial velocity of ball = [1.0, 2.0]
# SUCCESS

# Animation looks correct. Every frame, the ball moves to a different \
  posision along the wall of the container as expected (every frame cuts to \
  a new collision with the container).
'''
from Simulation import Simulation

sim4 = Simulation(
                  n_balls = 1, \
                  random_positions = False, \
                  random_velocities = False, \
                  set_position = np.array([-5.0, 0.0]), \
                  set_velocity = np.array([1.0, 2.0])
                 )

sim4.run(
         num_frames = 100, \
         animate = True
        )

# %%
'''
# Tasks 6, 7 & 8 Test - testing animation for multiple balls and plotting \
    kinetic energy and pressure of the gas over time:
#===============================================================================

1) Test set_random_positions() and set_random_velocities() methods:

Note: uncomment the print statement at the end of each method in the \
    Simulation.py file to print the list of positions and velocities.

# Number of balls = 20
# Radius of container = 20.0
# Velocity mean = 0.0
# Velocity width = 2.0
'''
import numpy as np
from Simulation import *

sim5 = Simulation(
                  n_balls = 20, \
                  random_positions = True, \
                  random_velocities = True, \
                  velocity_width = 2.0
                 )

'''
Outputs:

random initial positions 
=  [array([2.05459691, 9.88710293]), array([12.29135141, -13.67679863]), \
    array([5.41486051, -16.77654056]), array([-13.7295733 ,   7.64576421]), \
    array([7.16679803, 0.61974027]), array([-9.28363731, 12.21463304]), \
    array([-3.700038  ,  3.42916338]), array([-15.44339521,   1.06864524]), \
    array([-1.59559559, -9.27992363]), array([-15.85240116,   8.12107746]), \
    array([12.41380625, 11.68058996]), array([ 2.9295606 , -9.73519712]), \
    array([-4.9719224, -1.6629199]), array([-4.25476259,  7.35326151]), \
    array([12.26859618, -5.35011504]), array([-7.33781879,  7.1545335 ]), \
    array([4.79392505, 13.54496187]), array([2.29246027, -18.08626515]), \
    array([0.96665938, 18.65090841]), array([-0.84700919,  2.87528294])]

random initial velocities 
=  [array([-1.7368967 , -0.41984076]), array([4.15060283, -0.03408295]), \
    array([0.50976529, 3.08959349]), array([1.22817213, -3.56180875]), \
    array([-0.14632486,  1.38105415]), array([-0.17801719,  0.26756733]), \
    array([-0.36918568, -0.11197521]), array([3.01382798, 3.37427777]), \
    array([-0.81353526, -1.85653726]), array([0.25123701, 1.02516729]), \
    array([0.58871846, 2.49434748]), array([1.76002727, -0.99980284]), \
    array([0.23062711, -0.92645225]), array([-3.01549511,  1.4886416 ]), \
    array([-0.38034176,  0.03332307]), array([2.64702877, -0.01111107]), \
    array([3.56235275, -2.95326902]), array([-0.78173145, -0.194219  ]), \
    array([2.12673476, 1.03557115]), array([0.779944  , 0.01155013])]

# SUCCESS
'''

'''
2) Test animation for multiple balls + kinetic energy and pressure plots:

# Number of collisions = 1000
# Initial positins and velocities randomised.
# SUCCESS

# Animation looks correct. No balls overlap with one another and all balls \
  stay inside the container.
# As expected, the two figures generated below suggest that both the total \
  kinetic energy of the gas particles and the total pressure exerted on the \
  container are conserved. 
# Kinetic energy vs. time graph outputs a straight horizontal line, while the \
  pressure vs. time graph starts off with a large fluctuation that quickly \
  converges to a mean equilibrium pressure value.
# Total kinetic energy is conserved since the gas particles collide elastically. 
'''
sim5.run(
         num_frames = 1000, \
         animate = True, \
         ke_sys_fig = True, \
         pressure_fig = True, \
        )

# %%
