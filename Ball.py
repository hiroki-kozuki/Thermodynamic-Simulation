"""
Created on Fri Nov 11 18:25:33 2022

@author: Hiroki
"""
import numpy as np
import pylab as pl

# To prevent 'RuntimeWarning: invalid value encountered in subtract'
#np.seterr(all='raise')
# To prevent 'FloatingPointError: invalid value encountered in subtract'
#np.seterr(invalid='warn')

class Ball(object):
    '''
    This class defines the ball objects used to represent circular gas \
    particles that will undergo elastic collisions inside a circular \
    container. The same class will also be used to define the container object.
    
    Attributes:
    ===========
        self._mass (float): Mass of the object.
        
        self._radius (float): Radius of the object.
        
        self._position (1 x 2 numpy.ndarray of floats): 
            Current position of the object.
        
        self._velocity (1 x 2 numpy.ndarray of floats): 
            Current velocity of the object.
        
        self._iscontainer (bool, True or False): Identifies whether the object \
            is a ball or a container.     
             
        self._collision_counter (int): Counts the number of times the object \
            collides.
            
        self._patch (matplotlib.pylab.Circle object): Patches used to indicate \
            the moving and colliding objects in the animation. 
    
    Methods:
    ========
        
        pos(): Returns current position of the object (ball or container).
    
        set_pos(position): Updates the positions of the object.
    
        vel(): Retunrs the current position of the object.
        
        set_vel(velocity): Updates the velocity of the object.
    
        time_to_collision(other): Calculates the time until self collides with \
            other object.
    
        move(dt): Moves the object to a new position corresponding to dt.
    
        collide(other): Updates the post-collision velocities of self and \
            other object.
        
        get_patch(): Returns the object's patch to be used in animation. 
    
        ke(): Returns the object's kinetic energy. 
        
        magnitude_momentum(): Returns the magnitude of object's momentum. 
    '''
    
    def __init__(self, \
                 mass = 1.0, \
                 radius = 1.0, \
                 position = np.array([0.0, 0.0]), \
                 velocity = np.array([0.0, 0.0]), \
                 iscontainer = False
                ):
        '''
        This method initialises the Ball class.  
        Ball and container objects are created here.
        
        Args: 
        =====
            mass (float): Mass of the object. Default is 1.0.
            
            radius (float): Radius of the object. Default is 1.0.
            
            position (list or 1 x 2 npmpy.ndarray containing floats): 
                Current position of the object. Default is [0.0, 0.0].
                
            velocity (list or 1 x 2 npmpy.ndarray containing floats): Current \
                velocity of the object. Default is [0.0, 0.0].
         
            iscontainer (bool): If true, object is a container and if False, \
                object is a ball. Default is False.       
        
        Returns:
        ========
            N/A
        
        Raises:
        =======
            Exception: If self._mass, self._radius or elements of \
                self._position and self._velocity are not stored as floats.
        '''
        self._mass = mass
        if self._mass != float(self._mass):
            raise Exception('Mass must be given as a float.')
        
        self._radius = radius
        if self._radius != float(self._radius):
            raise Exception('Radius must be given as a float.')
        
        self._position = position
        if self._position[0] != float(self._position[0]) \
        or self._position[1] != float(self._position[1]):
            raise Exception('Elements of position must be given as a float.')
        
        self._velocity = velocity
        if self._velocity[0] != float(self._velocity[0]) \
        or self._velocity[1] != float(self._velocity[1]):
            raise Exception('Elements of velocity must be given as a float.')
        
        self._iscontainer = iscontainer    
        self._collision_counter = 0
        
        if self._iscontainer == True:
            self._patch = pl.Circle(self._position, self._radius, ec = 'blue', \
                                    fill = False, ls = 'solid')
        if self._iscontainer == False:
            self._patch = pl.Circle(self._position, self._radius, \
                                    ec = 'black', fc = 'red')
    
    
    def __repr__(self):
        '''
        Returns the object's mass, radius, current position, current velocity, \
        identity (ball or container), and number of collisions when the object \
        is called.
        '''
        return f'Ball Class: mass = {self._mass}, radius = {self._radius}, \
                 position = {self._position}, velocity = {self._velocity}, \
                 is a container = {self._iscontainer}, \
                 collision = {self._collision_counter}'
    
    
    def __str__(self):
        '''
        Returns the same information as __repr__ when the object is called \
        within a print statement.
        '''
        return f'Ball Class: mass = {self._mass}, radius = {self._radius}, \
                 position = {self._position}, velocity = {self._velocity}, \
                 is a container = {self._iscontainer}, \
                 collision = {self._collision_counter}'
    
    
    def pos(self):
        '''
        This method is used to access the current position of the object. 
        
        Returns:
        ========
            The position of the object (1 x 2 numpy.ndarray).
        '''
        return np.array(self._position)
    
    
    def set_pos(self, new_position):
        '''
        This method is used to update the position of the ball object.
        
        Args:
        =====
            New_position (list or 1 x 2 numpy.ndarray of floats):
                New position of the ball. 
                
        Returns:
        ========
            N/A
                
        Raises:
        =======
            Exception: If the elements of new self._position are not stored \
                as floats.
        '''
        self._position = np.array(new_position)
        if self._position[0] != float(self._position[0]) \
        or self._position[1] != float(self._position[1]):
            raise Exception('Elements of position must be given as a float.')
    
    
    def vel(self):
        '''
        This method is used to access the current velocity of the object. 
        
        Returns:
        ========
            The velocity of the object (1 x 2 numpy.ndarray).
        '''
        return np.array(self._velocity)
    
    
    def set_vel(self, new_velocity):
        '''
        This method is used to update the velocity of the ball object.
        
        Args:
        ==========
            new_velocity (list or 1 x 2 numpy.ndarray of floats):
                New velocity of the ball.
        
        Returns:
        ========
            N/A
                
        Raises:
        =======
            Exception: If the elements of new self._velocity are not stored \
                as floats.
        '''
        self._velocity = np.array(new_velocity)
        if self._velocity[0] != float(self._velocity[0]) \
        or self._velocity[1] != float(self._velocity[1]):
            raise Exception('Elements of velocity must be given as a float.')
    
    
    def time_to_collision(self, other):
        '''
        This method is used to calculate the time it takes until the self \
        object collides with the other object given as the argument. 
        
        A quadratic equation is solved for time (dt) in terms of the relative \
        position and relaive velocity of the self and other objects. The \
        solution(s) to the quadratic equation is then classified based on the \
        sign of the determinant and smallest real positive solution is \
        selected in every case. If no real positive solution is found, an \
        arbitrarily large dt value (1.0e15) is returned. 
        
        Args:
        =====
            other (class): A second ball or container class/object.
        
        Returns:
        ========
            The time until self and other obejcts collide (float). 
        
        Classifying dt to find time-to-collision:
        ==============================
            dt is determined by solving the quadratic equation of the form \
            (r + v * dt) * (r + v * dt) == R ** 2 where r == r1 - r2 \
            (relative position), v == v1 - v2 (relative velocity) and \
            R == R1 +/- R2 (collision distance). Rearranging the quadratic \
            equation yields: dt = (-b +/- sqrt(b * b - 4 * a * c) / (2 * a). \
            dt is classified by first the sign of a, then the sign of the \
            determinant, the sign of c, then the sign of b. 
            
            - When "a" is zero, the two objects are stationary with respect to \
              each other and hence they will never collide. 
            - When the discriminant is negative, dt is complex. When the \
              discriminant is zero, dt exists if the soltion is positive. 
            - When the decriminent is positive, there exist two solutions \
              (pos_sqrt_soln and neg_sqrt_soln) and the signs of c and b are \
              used to select the smallest positive solution. 
            - Note that in a ball-container collision, dt has a valid solution \
              if c <= 0 (when ball is inside or toucing the container) and in \
              a ball-ball collision, dt has a valid solution if c > 0 (balls \
              are not overlapping).
            - In each case, dt has a valid solution when b < 0. 
        '''
        r = self._position - other._position
        v = self._velocity - other._velocity
        a = np.dot(v, v)
        b = 2 * np.dot(r, v)
                                                 
        if a == 0:
            return 1.0e15   
        
        elif a > 0:                   
            if other._iscontainer == True:                     
                R = other._radius - self._radius
                
            else:
                R = self._radius + other._radius
            
            c = np.dot(r, r) - R * R            
            discriminant = b * b - 4 * a * c   
            
            if discriminant < 0:
                    return 1.0e15
                
            elif discriminant == 0:
                if -b / (2 * a) <= 0:
                    # dt is 0 or negative.
                    return 1.0e15
                
                else:
                    return -b / (2 * a)
            
            elif discriminant > 0:  
                # When c < b * b / (4 * a).
                pos_sqrt_soln = (-b + np.sqrt(discriminant)) / (2 * a)
                neg_sqrt_soln = (-b - np.sqrt(discriminant)) / (2 * a)
                
                # Collision with the container.
                if other._iscontainer == True:
                    if c == 0: 
                        # Ball touching the container.
                        if b < 0:
                            # neg_sqrt_soln > 0 while neg_sqrt_soln == 0.
                            return pos_sqrt_soln
                              
                        elif b > 0:
                            # Ball is overlapping with the container or has \
                                # escaped.
                            return 1.0e15
                        
                        else:
                            # Ball is stuck to the container.
                            return 1.0e15
                    
                    elif c < 0:
                        # pos_sqrt_soln > 0 while neg_sqrt_soln < 0 for all b.
                        return pos_sqrt_soln 
                    
                    else:
                        # Ball is overlapping or with the container or has \
                            # escaped.
                        return 1.0e15
                    
                # Collision with another ball.
                if other._iscontainer == False:
                    if c == 0:
                        if b > 0:
                            # Balls are moving away from eachtoher.
                            return 1.0e15
                        
                        elif b < 0:
                            # Balls are overlapping.
                            return 1.0e15
                        
                        else:
                            # Balls are stuck together.
                            return 1.0e15
                    
                    if c > 0:
                        if b > 0:
                            # Both solutions are negative.
                            return 1.0e15
                        
                        elif b < 0:
                            # neg_sqrt_soln < pos_sqrt_soln while both \
                                # solutions are positive.
                            return neg_sqrt_soln 

                        else:
                            # dt is complex.
                            # Or the balls are moving in parallel / \
                                # anti-parallel direction with r not parallel \
                                # to v (not always true but a statistically \
                                # reasonable assumption).
                            return 1.0e15
                            
                    else:
                        # Balls are overlapping.
                        return 1.0e15
  
                                              
    def move(self, dt):
        '''
        This method calculates the new position of the ball based on its \
        current position, current velocity, and time dt given as the argument. \
        It then updates the position of the ball to this new position by \
        calling the set_pos(new_position) method on self. It is used to move \
        the ball to a position where it will next collide with the other object. 
        
        Args:
        =====
            dt (float): Time until the next collision.
        
        Returns:
        ========
            N/A
        '''
        new_position = self._position + self._velocity * dt
        self.set_pos(new_position)
        
        
    def collide(self, other):
        '''
        This methodn is used to calculate the new velocities of objects after \
        they collide elastically. 
        
        For the ball-container collision, the relative velocity vector v \
        (which is just the velocity of the ball - self._velocity) is resolved \
        into parallel (vi_self_para) and perpendicular (vi_self_perp) \
        components w.r.t the relative position vector r (which is also just \
        the position vector of the ball - self._position). The sign of \
        vi_self_para changes upon collision to yield vf_self_para while the \
        perpendicular component of v remains the same. 
        
        For the ball-ball collision, combining the conservation laws for \
        momentum and kinetic energy in 2D yields a comprehensive equation, \
        which expressed the new velocities of the two balls (vf_self and \
        vf_other) in terms of their masses (self._mass and other._mass), \
        initial velocities (self._velocity and other._velocity), r, and v. 
        
        In each case, the velocities of objects are updated by calling \
        set_vel(new_velocity) method on self and other. Note that the velocity \
        of the container is kept as zero throughout the simulation. 
        
        Args:
        =====
            other (class): A second ball or container class/object.
        
        Returns:
        ========
            N/A    
        '''
        r = self._position - other._position
        v = self._velocity - other._velocity 
        magnitude_r_squared = np.dot(r, r)

        # Ball - container collision:
        if other._iscontainer == True:                      
            vi_self_para = r * np.abs(np.dot(self._velocity, r)) / \
                           magnitude_r_squared
            vi_self_perp = self._velocity - vi_self_para
            
            vf_self_para = -vi_self_para 
            # Since the component of self._velocity parallel to r is radial. 
            vf_self_perp = vi_self_perp  
            self.set_vel(vf_self_para + vf_self_perp)
            
            self._collision_counter += 1
            other._collision_counter += 1
        
        # Ball - ball collision:   
        elif other._iscontainer == False:
            
            vf_self = self._velocity - (
                                        2 * other._mass / \
                                        (self._mass + other._mass)
                                        ) \
                                        * (np.dot(v, r) / np.dot(r, r)) * r
            
            vf_other = other._velocity - (
                                          2 * self._mass / \
                                          (self._mass + other._mass)
                                          ) \
                                          * (np.dot(-v, -r) / np.dot(-r, -r)) \
                                          * -r
# Reference:
# https://stackoverflow.com/questions/35211114/2d-elastic-ball-collision-physics
            
            self.set_vel(vf_self)
            other.set_vel(vf_other)

            self._collision_counter += 1
            other._collision_counter += 1
            

    def get_patch(self):
        '''        
        This method is used to access the patch for the object to be used in \
        animation. 
        
        Returns:
        ========
            Patches used to indicate the moving and colliding objects in the \
            animation (matplotlib.pylab.Circle object).
        '''
        return self._patch
    
    
    def ke(self):
        '''        
        This method is used to calculate the current kinetic energy of the ball.
        
        Returns:
        ========
            Current kinetic energy of the ball (float).
        '''
        return self._mass * np.dot(self._velocity, self._velocity) / 2
    
    
    def magnitude_momentum(self):
        '''
        This method is used to calculate the magnitude of the current momentum \
        of the ball.
        
        Returns:
        ========
            Magnitude of the current momentum of the ball (float)
        '''
        return self._mass * np.sqrt(np.dot(self._velocity, self._velocity))
    
