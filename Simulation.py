"""
Created on Fri Nov 18 09:14:54 2022

@author: Hiroki
"""

import numpy as np
import pylab as pl
import itertools as iter
import matplotlib.pyplot as plt
from Ball import Ball

# To prevent 'RuntimeWarning: invalid value encountered in subtract'
#np.seterr(all='raise')

# To prevent 'FloatingPointError: invalid value encountered in subtract'
#np.seterr(invalid='warn')

class Simulation(object):
    '''
    This class defines the simulation of balls bouncing and elastically \
    colliding inside a container as a proxy for a system of ideal gas particles.
    
    The class encompasses methods for initialising and updating the animation \
    of particle collisions with random initial positions and velocities, a \
    method for selecting and exectuing the next collision and updating the \
    simulation by minimum time until next collision, and a run method for \
    iterating the above methods for a given number of collisions/frames, \
    collecting data, and generating relevant plots.

    Attributes:
    ===========
        self._mass_ball (float): Mass of the ball.
        
        self._radius_ball (float): Radius of the ball.
        
        self._radius_container (float): Radius of the container.
        
        self._n_ball (int): Number of balls.
        
        self._random_positions (bool): True for initialising balls at random \
            positions inside the container. False for initialising a single \
            ball with specified position.
            
        self._random_velocities (bool): True for initialising balls at random \
            velocities within a specified range. False for initialising a \
            single ball with specified velocity.
        
        self._set_position (list or 1 x 2 numpy.ndarray of floats):
            Specific initial position for a single-ball system. 
            
        self._set_velocity (list or 1 x 2 numpy.ndarray of floats):
            Specific initial velocity of a single-ball system. 
        
        self._velocity_width (float): Speed corresponding to the width of the \
        normal distribution for the randomly generated velocities of balls. 
        
        self._container (class): Container object inside which balls collide.
        
        self._balls (dict of classes): Dictionary containing ball objects as \
            values and integers beginning from 0 as keys.
        
        self._dt_min (float): Minimum time to next collision. 
        
        self._mom_transferred_container (float): Magnitude of the total \
            momentum transfered to the container during the entire simulation.
        
        self._global_time (float): Duration of the simulation. 
        
        self._time_list (list of floats): List of minimum dt at each stage of \
            the simulation.
        
        self._distance_centre_list (list of floats): List containing distances \
            of balls from the centre of the container.  
        
        self._separation_relative_list (list of floats): List of distances \
            between pairs of balls.
        
        self._ke_system_list (list of floats): List of the system's total \
            kinetic energy at each stage of the simulation.
        
        self._momentum_system_list (list of floats): List of the magnitudes of \
            the total momentum of the system at each stage of the simulation.
    
        self._temperature_list (list of floats): List of the gas' temperature \
            at each stage of the simulation.
        
        self._pressure_list (list of floats): List of the pressure exerted on \
            the container by the gas particles at each stage of the simulation. 
        
        self._speed_list (list of floats): List of balls' speeds. 
    
        self._kinetic_energy_eq (float): Total kinetic energy of gas particles \
            at termodynamic equilibrium.
    
        self._momentum_system_eq (float): Magnitude of the total momentum of \
            the system at thermodynamic equilibrium.
        
        self._temperature_eq (float): Temperature of the gas at thermodynamic \
            equilibrium. 
        
        self._pressure_eq (float): Pressure exerted on the container by the \
            gas particles at thermodynamic equilibrium.
        
        
    Methods:
    ========
        
        set_random_positions(): Generates random positions for the ball \
            objects according to a uniform distribution.
        
        set_random_velocities(): Generates random velocities for the ball \
            objects according to a normal distribution.
    
        generate_patches(): Generates patches for the ball and container \
            objects to be used in the animation.
        
        next_collision(): Makes use of the time_to_collision, move, and \
            collide methods of the Ball class. Determines the shortest time \
            until the next collision (self._dt_min), moves the simulation by \
            that time, then executes collision by updating the balls' \
            velocities. It also increments self._global_time by self._dt_min. 
        
        run(num_frames = 100, animate = False, dist_fig = True, \
            ke_sys_fig = True, mom_sys_fig = True, temp_fig = True, \
            pressure_fig = True, speeds_fig = True):
            Generates the animation, executes the next_collision method, \
            collects and returns relevant data, updates the animation, and \
            plots the optional figures.
        
        distance_centre(): Returns a list of the distances of balls from \
            the centre of the container at an instance in the simulation.
        
        separation_relative(): Returns a list of the distances between all \
            the pairs of balls at an instance in the simulation.
        
        ke_system(): Returns the system's total kinetic energy at an \
            instance in the simulation.     
        
        pressure(): Returns the magnitude of the pressure exerted on the \
            container by the gas particles over the duration of the simulation \
            (self._global_time).
        
        temperature(): Returns the temperature of the gas particles at an \
            instance in the simulation.
        
        momentum_system(): Retunrs the magnitude of the system's total \
            momentum at an instance in the simulation.
        
        particle_speed(): Returns a list of the speeds of the gas particles at \
            an instance in the simulation.     
    '''
    
    def __init__(
                 self, \
                 mass_ball = 1.0, \
                 radius_ball = 1.0, \
                 radius_container = 20.0, \
                 n_balls = 1, \
                 random_positions = True, \
                 random_velocities = True, \
                 set_position = [-5.0, 0.0], \
                 set_velocity = [1.0, 0.0], \
                 velocity_mean = 0.0, \
                 velocity_width = 2.0
                 ):
        '''
        This method initialises the Simulation class. 
        
        Args:
        =====
            mass_balls (float): Mass of ball object. Default is 1.0.
            
            radius_ball (float): Radius of ball object. Default is 1.0.
            
            radius_container (float): Radius of the container. Default is 20.0.
            
            n_balls (int): Number of balls in the simulation. Default is 1.
            
            random_positions (bool): If True, the balls are initialised with \
                random positions generated from a uniform distribution centred \
                at the origin [0.0, 0.0]. If False for a single ball, the ball \
                is initialised with a specific position given as one of the \
                arguments. Default is True.
            
            random_velocities (bool): If True, the balls are initialised with \
                random velocities generated from a normal distribution \
                centred at [0.0, 0.0] and with a speed corresponding to the \
                distribution's width (1 standard deviation from the mean) \
                given as one of the arguments. If False for a single ball, the \
                ball is initialised with a specific velocity given as one of \
                the arguments. Default is True.
        
            set_position (list of 1 x 2 numpy.ndarray): If random_positions is \
                set to False for a single ball, the ball is initialised with \
                position given here. Default is [-5.0, 0.0].
            
            set_velocity (list of 1 x 2 numpy.ndarray): If random_velocity is \
                set to False for a single ball, the ball is initialised with \
                velocity given here. Default is [1.0, 0.0].
            
            self._velocity_mean (float): Mean value of the normal distribution \
                for randomly generated velocities of the balls. Default is 0.0.
            
            velocity_width (float): Speed corresponding to the width of the \
                normal distribution for the randomly generated velocities of \
                the balls. Default is 2.0.
    
        Returns:
        ========
            N/A
        '''
        self._mass_ball = mass_ball
        self._radius_ball = radius_ball
        self._radius_container = radius_container
        self._n_ball = n_balls
        self._random_positions = random_positions
        self._random_velocities = random_velocities
        self._set_position = set_position
        self._set_velocity = set_velocity
        self._velocity_mean = velocity_mean
        self._velocity_width = velocity_width
        
        # Container object:
        self._container = Ball(
                               mass = 1.0e99, \
                               radius = self._radius_container, \
                               iscontainer = True
                               )
        
        # Dictionary of ball objects:
        self._balls = dict()
        for n in range(self._n_ball):
            self._balls[n] = Ball(
                                  mass = self._mass_ball, \
                                  radius = self._radius_ball, \
                                  iscontainer = False
                                  )
        
        if self._random_positions == True:
            self.set_random_positions()
        elif self._random_positions == False:
            self._balls[0].set_pos(self._set_position)
            
        if self._random_velocities == True:
            self.set_random_velocities()
        elif self._random_velocities == False:
            self._balls[0].set_vel(self._set_velocity)
        
        self._dt_min = 0.0
        self._mom_transferred_container = 0.0
        self._global_time = 0.0 
        self._time_list = []
        self._distance_centre_list = []
        self._separation_relative_list = [] 
        self._ke_system_list = []
        self._momentum_system_list = []
        self._temperature_list = []
        self._pressure_list = []
        self._speed_list = []
        self._kinetic_energy_eq = 0.0
        self._momentum_system_eq = 0.0
        self._temperature_eq = 0.0
        self._pressure_eq = 0.0


    def __repr__(self):
        '''
        Returns the ball's mass, radius, number of balls, range of velocities, \
        and the container's radius when the Simulation object is called.
        '''
        return f'Simulation Class: ball mass = {self._mass_ball}, \
                 ball radius = {self._radius_ball}, \
                 number of balls = {self._n_ball}, \
                 range of ball velocities = +/- {self._velocity_width}, \
                 container radius = {self._container._radius}'
    
    
    def __str__(self):
        '''
        Returns the same status of the simulation as __repr__ when the ball \
        object is called within a print statement.
        '''
        return f'Simulation Class: ball mass = {self._mass_ball}, \
                 ball radius = {self._radius_ball}, \
                 number of balls = {self._n_ball}, \
                 range of ball velocities = +/- {self._velocity_width}, \
                 container radius = {self._container._radius}'
        
    
    def set_random_positions(self):
        '''
        This method is used to generate random positions for the ball objects \
        according to a uniform distribution centred at the origin of the \
        container [X = 0.0, y = 0.0].
        
        A while-loop is implemented to iterate the process until the newly \
        generated position does not overlap with any existing balls. 
        
        The positions of the balls are systematically updated to the randomly \
        generated positions by calling the set_pos(new_position) method inside 
        a for-loop.
        
        Returns:
        ========
            N/A
        
        Raises:
        =======
            Exception: If the while loop for generating random positions has \
                been iterated for an arbitrary large number of times (1e7), it \
                indicates that the container area is too small to fit all the \
                balls. An exception is raised to inform the user to either \
                increase the container radius or to decrease the number of \
                balls.
        '''
        initial_position_list = []
        for n in range(self._n_ball):
            position = np.zeros(2)
            x_y_range = self._container._radius - self._radius_ball - 0.2
            # Subtracted 0.2 to ensure that the balls are initially \
            # sufficiently spaced from the container wall (i.e. prevent balls \
            # touching wall initially). 
            
            # Counter for the number of loops.
            counter = 0 
            
            # Reference for nested while loops:
            # https://realpython.com/python-while-loop/
            while True:         
                x = np.random.uniform(-x_y_range, x_y_range)
                y = np.random.uniform(-x_y_range, x_y_range)
                
                # Repeat until the random x and y coordinates of each ball are \
                # within the container. 
                while np.sqrt(x * x + y * y) >= self._container._radius \
                                                - self._radius_ball - 0.2:
                    x = np.random.uniform(-x_y_range, x_y_range)
                    y = np.random.uniform(-x_y_range, x_y_range)

                    if np.sqrt(x * x + y * y) < self._container._radius \
                                                - self._radius_ball - 0.2:
                        break
                    
                position = np.array([x, y])
                # Indicator for overlaps.
                overlap = True
                 
                # Asess whether generated random position causes ball to \
                # overlap with existing balls. If ball overlaps with any other \
                # ball, break for-loop.
                for m in range(0, n):
                    distance_between_balls = np.sqrt(np.dot
                                            (
                                             self._balls[m].pos() \
                                             - position, self._balls[m].pos() \
                                             - position
                                             )       )
                    
                    if distance_between_balls <= 2 * self._radius_ball:
                        overlap = True
                        counter += 1
                        break

                    else:
                        overlap = False
                
                if counter > 1e7:
                    raise Exception('Container area is too small to fit all \
                                     the balls. Increase container area or \
                                     reduce the number of balls.')
                
                # If n is 0 or the ball does not overlap with any other balls, \
                # break while-loop.
                if n == 0 or overlap == False:
                    break
            
            self._balls[n].set_pos(position)
            initial_position_list.append(position)
        #print('random initial positions = ', initial_position_list)
    
    
    def set_random_velocities(self):
        '''
        This method is used to generate random velocities for the ball objects \
        according to a normal distribution with its mean located at the value \
        of self._velocity_mean and its spread (standard devaition) equal to \
        self._velocity_width.
        
        The velocity of the balls are systematically updated to the randomly \
        generated velocities by calling the set_vel(new_velocity) method inside 
        a for-loop.
        
        Returns:
        ========
            N/A        
        '''
        initial_velocity_list = []
        for n in range(self._n_ball):            
            velocity = np.random.normal(
                                        loc = self._velocity_mean, \
                                        scale = self._velocity_width, \
                                        size = 2
                                         )
            self._balls[n].set_vel(velocity)
            initial_velocity_list.append(velocity)
        #print('random initial velocities = ', initial_velocity_list)
        
        
    def generate_patches(self):
        '''
        This method initialises and draws out the patches for the ball and \
        the container objects (inside a for-loop) to be used in the animation. \
        Patches are drawn as matplotlib.pylab.Circle objects.
        
        Returns:
        ========
            N/A
        '''
        patch_ball_list = []
        for n in range(self._n_ball):
            patch_ball_list.append(self._balls[n].get_patch())
        
        self._ball_patch = patch_ball_list
        self._container_patch = self._container.get_patch()
        
        pl.figure(num = 'Animation')
        ax = pl.axes(
                     xlim = (
                             -self._container._radius, \
                             self._container._radius
                             ), \
                     ylim = (
                             -self._container._radius, \
                             self._container._radius
                             ), \
                     aspect = '1'
                     )
                
        ax.add_patch(self._container_patch)
        for n in range(0, len(self._ball_patch)):
            ax.add_patch(self._ball_patch[n])
        
    
    def update_patches(self):
        '''
        This method systematically updates the patch positions for ball \
        objects in the animation by calling the "center" attribute of the ball \
        patch and assinging to it the updated position using the pos() method \
        called on the corresponding ball object.
        
        Returns:
        ========
            N/A
        '''    
        for n in range(0, len(self._ball_patch)):
            self._ball_patch[n].center = self._balls[n].pos()
    
    
    def next_collision(self):
        '''
        This Method:
        - Calculates the time-to-collision (dt) for all possible pairs of ball \
          objects with eachother and with the container. This is done by \
          creating a dictionary containing all the ball and container objects \
          (balls_and_container), using the itertools.combinations function to \
          find all the unique pairs of objects (object_pairs), and calling the \
          time_to_collision(other) method on every pair inside a for-loop to \
          store the dt values in a list (dt_list).
        - Finds the shortest time-to-collision, i.e, the time to the next \
          collision (dt_min) by applying the sorted() function on a dictionary \
          (dt_dict) containing dt values and their corresponding object pairs \
          (pair_dt_min) to reorder them in an ascending order of dt values.
        - Moves all the balls by dt_min by calling the move(dt) method on the \
          ball objects inside a for-loop.
        - Executes collision for the pair of objects corresponding to dt_min \
          (update their velocities) by calling the collide(other) method.
        - Increments the total momentum transfered to the container \
          (self._mom_transferred_container) if the collision is between a ball \
          and the container.
        - Increments the duration of the simulation (self._global_time) by \
          dt_min.
                
        Returns:
        ========
            N/A
        '''
        balls_and_container = dict()
        for n in range(self._n_ball):
            balls_and_container[n] = list(self._balls.values())[n]
        balls_and_container[self._n_ball] = self._container
        # By inserting the container object at the end of the dictionary, we \
        # ensure that the container is always treated as "other" or object_2 \
        # when a ball collides with it (due to the pairing algorith employed \
        # in itertools.combinations below). 
        
        object_pairs = list(iter.combinations(balls_and_container.values(), 2))
        
        dt_list = []
        for i in range(len(object_pairs)):
            dt_list.append(
                           object_pairs[i][0].time_to_collision
                                                               (
                           object_pairs[i][1]
                                                                )
                           )
        
        # Create a dictionary with elements {dt : object_pair}.
        dt_dict = dict(zip(dt_list, object_pairs))
        # Reference: 
        # https://stackoverflow.com/questions/209840/how-can-i-make-a-dictionary-dict-from-separate-lists-of-keys-and-values
        
        # (dt_min, pair_dt_min) is the first tuple in the dictionary.
        dt_min, pair_dt_min = sorted(dt_dict.items(), key = lambda x: x[0])[0]
        # Source:
        # https://www.freecodecamp.org/news/sort-dictionary-by-value-in-python/#:~:text=To%20correctly%20sort%20a%20dictionary,with%20the%20item()%20method
        
        self._dt_min = dt_min
        object_1, object_2 = pair_dt_min[0], pair_dt_min[1]     
        
        for n in range(self._n_ball):
            self._balls[n].move(dt_min)
        
        # Record velocities of object pairs BEFORE collision.
        vi_object_1 = object_1.vel()
        vi_object_2 = object_2.vel()
     
        object_1.collide(object_2)
        
        # Record velocities of object pairs AFTER collision.
        vf_object_1 = object_1.vel()
        vf_object_2 = object_2.vel()
        
        dv_object_1 = vf_object_1 - vi_object_1
        dv_object_2 = vf_object_2 - vi_object_2
        
        if object_1._iscontainer == True:
            self._mom_transferred_container \
            += self._mass_ball \
            * np.sqrt(np.dot(dv_object_2, dv_object_2))
            
        if object_2._iscontainer == True:
            self._mom_transferred_container \
            += self._mass_ball \
            * np.sqrt(np.dot(dv_object_1, dv_object_1))
    
        self._global_time += self._dt_min
        

    def run(
            self, \
            num_frames = 1000, \
            animate = True, \
            dist_fig = False, \
            ke_sys_fig = False, \
            mom_sys_fig = False, \
            temp_fig = False, \
            pressure_fig = False, \
           ):
        '''
        This method:
        - Draws the animation including the ball patches and the container \
          patch by calling the generate_patch() method 
          (Optional).
        - Runs the simulation for a specified number of frames/collisions \
          (num_frames) by iterating the next_collision method inside a for-loop.
        - Collects relevant data (time, distance of balls from the origin, \
          distance between ball pairs, system's total kinetic energy, \
          magnitude of total momentum, temperature, pressure, and particle \
          speeds) for every collision and appends them onto designated lists \
          for every collision.
        - Updates the animation after each collision by calling the \
          update_patches() method. 
          (Optional)
        - Plots the following figures (optional):
            - histograms of the balls' distances from the origin and distances \
              between ball pairs.
            - graphs of the system's total kinetic energy, magnitude of total \
              momentum, temperature, and pressure against time.
        - Returns the equilibrium temperature and its standard error of the \
          mean, the pressure of the gas, as well as a list of ball speeds.
        
        Args:
        =====
            num_frames (int): Number of frames in the animation (i.e., total \
                number of collisions in the simulation). Default is 1000.
            
            animate (bool): If True, it generates an animation of the \
                simulation. Defaulst is True.
            
            dist_fig (bool): If True, it plots histograms of the balls' \
                distances from the origin and the distances between ball pairs. 
                Default is False.
                
            ke_sys_fig (bool): If True, it plots a graph of the system's \
                total kinetic energy vs. time. Default is False.
                
            mom_sys_fig (bool): If True, it plots a graph of the magnitude \
                of the system's total momentum vs. time. Default is False.
                
            temp_fig (bool): If True, it plots a graph of the temperature of \
                the gas vs. time. Default is False.
                
            pressure_fig (bool): If True, it plots a graph of the pressure \
                exerted on the container by the gas vs. time. Default is False.
        
        Returns:
        ========
            Equilibrium temperature of the gas (float).
            
            Standard error of the mean of the equilibrium temperature (float).
            
            Equilibrium pressure exerted on the container by the gas (float).
            
            List of speeds of the gas particles over the entire simulation \
                (list of floats)
        '''    
        
        if animate:
            
            self.generate_patches()
            pl.pause(0.001)
    
        for _ in range(num_frames):
            self.next_collision()
            
            # Collect data from the simulation and append to the lists \
                # defined in the init() method. 
            self._time_list.append(self._global_time)
            self._distance_centre_list.extend(self.distance_centre())
            self._separation_relative_list.extend(self.separation_relative())
            self._ke_system_list.append(self.ke_system())
            self._momentum_system_list.append(self.momentum_system())
            self._temperature_list.append(self.temperature())
            self._pressure_list.append(self.pressure())
            self._speed_list.extend(self.particle_speed())
            
            if animate:
                self.update_patches()
                pl.pause(0.001)
        
        # Calculate mean kinetic energy, mean total momentum of the balls, \
            # mean temperature, mean pressure exerted on the container, and \
            # their standard errors of the means. 
        self._kinetic_energy_eq = np.sum(self._ke_system_list) / num_frames
        SEM_ke = np.std(self._ke_system_list) / num_frames
        self._momentum_system_eq = np.sum(self._momentum_system_list) \
                                   / num_frames
        SEM_p = np.std(self._momentum_system_list) / num_frames
        self._temperature_eq = np.sum(self._temperature_list) / num_frames 
        SEM_T = np.std(self._temperature_list) / num_frames
        self._pressure_eq = np.sum(self._pressure_list) / num_frames
        SEM_P = np.std(self._pressure_list) / num_frames
        
        if animate:
            pl.show()
        
        if dist_fig == True:

            # Histogram of the particles' distances from centre.
            n, bins, patches = plt.hist(
                                    x = np.array(self._distance_centre_list), \
                                    bins = 100
                                        )
            plt.xlabel('Distance')
            plt.ylabel('Number of Counts')
            plt.title('Histogram of Distances of Particles from Centre, 20 Balls')
            plt.savefig('20 Balls, Distance from Centre Histogram.png', \
                         dpi=1000)
            plt.show()
            
            # Histogram of separation of each pair of particles. 
            n, bins, patches = plt.hist(
                                x = np.array(self._separation_relative_list), \
                                bins = 200
                                        )            
            plt.xlabel('Distance')
            plt.ylabel('Number of Counts')
            plt.title('Histogram of Separations Between Particles, 20 Balls')
            plt.savefig('20 Balls, Distance between ball pairs Histogram.png', \
                         dpi=1000)
            plt.show()
        
        if ke_sys_fig == True:
            
            # Plot for conservation of kinetic energy of the system.
            plt.plot(
                     np.array(self._time_list), \
                     np.round(np.array(self._ke_system_list), 7)
                     )
            plt.xlabel('Time')
            plt.ylabel('Kinetic Energy')
            plt.title('Kinetic Energy of the System vs. Time, 20 Balls')
            plt.savefig('20 Balls, KE vs Time.png', dpi=1000)
            plt.show()
            print('mean total KE = ', self._kinetic_energy_eq), '+/-', SEM_ke
        
        if mom_sys_fig == True:
            
            # Plot for conservation of momentum of the system.
            plt.plot(
                     np.array(self._time_list), \
                     np.array(self._momentum_system_list), \
                     label = 'Simulation Data'
                     )
            plt.axhline(
                        y = self._momentum_system_eq, \
                        color = 'orange', \
                        linestyle = '--', \
                        label = 'Equilibrium Momentum'
                        )
            plt.xlabel('Time')
            plt.ylabel('Momentum')
            plt.ylim(5, 100)
            plt.legend(loc = "upper right")
            plt.title('Magnitude of Total Momentum of Particles vs. Time, 20 \
                Balls')
            plt.savefig('20 Balls, Momentum vs Time.png', dpi=1000)
            plt.show()
            print('mean total momentum = ', self._momentum_system_eq, '+/-', \
                   SEM_p)
            
        if temp_fig == True:
            
            # Plot of the temperature of the gas against time.
            plt.plot(
                     np.array(self._time_list), \
                     np.round(np.array(self._temperature_list), 7)
                     )
            plt.xlabel('Time')
            plt.ylabel('Temperature')
            plt.title('Temperature of Gas vs. Time, 20 Balls')
            #plt.text()
            plt.savefig('20 Balls, Temperature vs Time.png', dpi=1000)
            plt.show()
            print('mean temperature = ', self._temperature_eq, '+/-', SEM_T)
        
        if pressure_fig == True:
            
            # Plot of the pressure exerted on the container by the gas \
                # particles against time.
            plt.plot(
                     np.array(self._time_list), \
                     np.array(self._pressure_list), \
                     label = 'Simulation Data'
                     )
            plt.axhline(
                        y = self._pressure_eq, \
                        color = 'orange', \
                        linestyle = '--', \
                        label = 'Equilibrium Pressure'
                        )
            plt.xlabel('Time')
            plt.ylabel('Pressure')
            plt.legend(loc = "upper right")
            plt.title('Pressure Exterted on the Container vs. Time, 20 Balls')
            #plt.text()
            plt.savefig('20 Balls, Pressure vs Time.png', dpi=1000)
            plt.show()
            print('mean pressure = ', self._pressure_eq, '+/-', SEM_P)
    
        return self._temperature_eq, SEM_T, self._pressure_eq, self._speed_list
        
    
    def distance_centre(self):
        '''        
        Returns:
        ========
            list of the distances (magnitudes of positions) of balls from the \
            centre of the container at an instance in the simulation. (list)
        '''
        distance_centre = []
        for n in range(self._n_ball):
            distance_centre.append(
                     np.sqrt(np.dot(self._balls[n].pos(), self._balls[n].pos()))
                                   )
        
        return distance_centre
        
        
    def separation_relative(self):
        '''
        Returns:
        ========
            list of the distances (magnitudes of positions) between the pairs \
            of balls at an instance in the simulation. (list)
        '''
        separation_relative = []
        ball_pairs = list(iter.combinations(self._balls.values(), 2))
        for n in range(len(ball_pairs)):
            position_relative = ball_pairs[n][0].pos() - ball_pairs[n][1].pos()
            separation_relative.append(
                        np.sqrt(np.dot(position_relative, position_relative))
                                       )
        return separation_relative
 
        
    def ke_system(self):
        '''            
        Returns:
        ========
            System's total kinetic energy. (float). The kinetic energy of the \
            container is effectively zero. 
        '''       
        
        ke_list = []
        for n in range(self._n_ball):
            ke_list.append(self._balls[n].ke())
        
        return np.sum(ke_list)
    
    
    def pressure(self):
        '''
        - Calculates the pressure exerted on the container by the gas \
        particles over the duration of the simulaiton (self._global_time).
        - Uses pressure == (momentum transferred / time) / circumference of \
        container.        
        
        Returns:
        ========
            Pressure exerted on the container by the gas particles. (float)
        
        '''            
        return (self._mom_transferred_container / self._global_time) \
                / (2 * np.pi * self._container._radius)


    def temperature(self):
        '''
        - Calculates the temperature of the gas particles at a point in time. 
        - Use mean KE of a particle == (f / 2) * k_B * T
        - Since the simulation is in 2D, f == 2. 
        - Adjust units such that Boltzmann constant == 1. 
        
        Returns:
        ========
        Temperature of the gas (float)
        '''
        return self.ke_system() / self._n_ball
    
    
    def momentum_system(self):
        '''
        Returns:
        ========
        Magnitude of the system's total momentum. (float) The momentum of the \
        container is effectively zero.
        '''
        total_momentum = []
        for n in range(self._n_ball):
            total_momentum.append(self._balls[n].magnitude_momentum())
            
        return np.sum(total_momentum)


    def particle_speed(self):
        '''
        Returns:
        ========
        List of particles' speeds at an instance in the simulation. \
        (list of floats)
        '''
        particle_speed_list = []
        for n in range(self._n_ball):
                particle_speed_list.append(
                    np.sqrt(np.dot(self._balls[n].vel(), self._balls[n].vel()))
                                           )
        return particle_speed_list
