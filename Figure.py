"""
Created on Fri Dec 02 20:11:47 2022

@author: hk621
"""
import numpy as np
from scipy.optimize import curve_fit
from Simulation import Simulation
import matplotlib.pyplot as plt
#%%
def pressure_vs_temperature(
                            radius_ball_list, \
                            radius_cont, \
                            num, \
                            v_width_list, \
                            collisions
                            ):
    '''
    This function generates a graph of Pressure vs. Temperature for a range of \
    ball radii suplemented with a line of best fit.
    
    Args:
    =====
        radius_ball_list (list of floats): List containing a range of ball \
            radii. To keep the balls radii constant, input a list with a \
            single element equal to the desired radius (e.g., [1.0]). 
            
        radius_cont (float): Radius of the container.
        
        num (int): Number of balls.
        
        v_width_list (list of floats): List containing a range of \
            velocity width (standard deviations) for the random initial \
            velcities of the balls. These values affect the pressure and \
            temperature of the gas.
        
        collisions (int): Number of collisions / number of frames in the \
            simulation.
    
    Returns:
    =======
        
        Fit coefficients (floats): Coefficients of the best fit linear \
            functions determined from least sqaures regression.
        
        Absolute uncertainties in the fit coefficients (floats): 
            Square roots of the diagonal elements in the covariance matrix \
            corresponding to the fit coefficients.
    '''  
    for radius in radius_ball_list:
        temp_list = []
        pressure_list = []
        for width in v_width_list:
            sim = Simulation(
                             radius_ball = radius, \
                             radius_container = radius_cont, \
                             n_balls = num, \
                             velocity_width = width
                             )
            
            temp_eq, unc_temp_eq, pressure_eq, speeds = sim.run(
                                                    num_frames = collisions, \
                                                    animate = False
                                                    )
            
            temp_list.append(temp_eq)
            pressure_list.append(pressure_eq)
            
        def linear_fit(T, m, c):
            return m * T + c
            
        # Generate coefficients for linear fit.
        params, cov_params = curve_fit(linear_fit, temp_list, pressure_list)

        # Plot simulation data and fit. 
        plt.plot(np.array(temp_list), np.array(pressure_list), 'o')
        x = np.linspace(temp_list[0], temp_list[-1])
        plt.plot(
                    x, \
                    linear_fit(x, *params), \
                    label = 'ball radius = %.1f' % radius
                    )   
        print('slope = ', params[0], '+/-', np.sqrt(cov_params[0,0]))
        print('y-int = ', params[1], '+/-', np.sqrt(cov_params[1,1]))
    
    plt.xlabel('Temperature')
    plt.ylabel('Pressure')
    plt.title('Pressure vs. Temperature')
    plt.legend()
    plt.savefig('Pressure vs. Temperature at different ball radii with 50 \
        balls 50 radius 200 collisions.png', dpi=1000)
    plt.show()

    return params[0], np.sqrt(cov_params[0,0]), params[1], \
           np.sqrt(cov_params[1,1])

#%%
def pressure_temperature_vs_area(
                                radius_container_list, \
                                num, \
                                v_width, \
                                collisions, \
                                pressure_area = True, \
                                temp_area = True
                                ):
    '''
    This function generates graphs of Pressure and Temperature against the \
        Area of the container. 
    
    Args:
    =====
        radius_container_list (list of floats): List containing a range of \
            container radii. 
            
        num (int): Number of balls.
        
        v_width (float): Width (standard deviation) of the normal distributiom \
            from which random initial velocities of balls are generated.
        
        collisions (int): Number of collisions / number of frames in the \
            simulation.
        
        pressure_area (bool): If True, it generates the graph of Pressure \
            vs. Area of container. Default is True.
        
        temp_area (bool): If True, it generates the graph of Temperature \
            vs. Area of container. Default is True.
    
    Returns:
    ========
        Estimated temperature of the gas for the Pressure vs. area plot \
            determined via least squares regression. (float)
        
        Absolute uncertainties in the estimated temperature (float): 
            Square roots of the first diagonal element in the covariance matrix.
        
        Estimated average pressure of the gas for the Temperature vs. area \
            plot determined via least squares regression. (float)
        
        Absolute uncertainties in the estimated average pressure (float): 
            Square roots of the first diagonal element in the covariance matrix.
    '''
    temp_list = []
    pressure_list = []
    for radius in radius_container_list:
        sim = Simulation(
                         radius_container = radius, \
                         n_balls = num, \
                         velocity_width = v_width
                         )
        
        temp_eq, unc_temp_eq, pressure_eq, speeds = sim.run(
                                               num_frames = collisions, \
                                               animate = False
                                               )
        pressure_list.append(pressure_eq)
        temp_list.append(temp_eq)
    
    container_area = np.pi * np.square(np.array(radius_container_list))
    
    if pressure_area == True:
        N = num
        # k_B == 1.0
        def inverse_fit(V, T, c):
            return N * T / V + c
            
        # Generate coefficients for linear fit.
        params_P, cov_params_P = curve_fit(
                                    inverse_fit, \
                                    container_area, \
                                    pressure_list
                                    )

        # Plot simulation data and fit. 
        plt.plot(
                container_area, \
                np.array(pressure_list), \
                'o', \
                label = 'Simulation Data'
                )
        x = np.linspace(container_area[0], container_area[-1])
        plt.plot(
                x, \
                inverse_fit(x, *params_P), \
                label = 'Inverse Fit'
                )
        
        print('Estimated Temp = ', params_P[0], '+/-', \
               np.sqrt(cov_params_P[0,0]))
        plt.xlabel('Container Area')
        plt.ylabel('Pressure')
        plt.title('Pressure vs. Container Area')
        plt.legend()
        plt.savefig('Pressure vs. Area with 50 balls 1000 collisions 5.0 v \
            width.png', dpi=1000)
        plt.show()
        
    
    if temp_area == True:
        # k_B == 1.0
        def linear_fit(V, P, N, c):
            return P * V / N + c
            
        # Generate coefficients for linear fit.
        params_T, cov_params_T = curve_fit(
                                    linear_fit, \
                                    container_area, \
                                    temp_list
                                    )

        # Plot simulation data and fit. 
        plt.plot(
                 container_area, \
                 np.array(temp_list), \
                 'o', \
                 label = 'Simulation Data'
                 )
        x = np.linspace(container_area[0], container_area[-1])
        plt.plot(
                 x, \
                 linear_fit(x, *params_T), \
                 label = 'Linear Fit'
                 )
        
        print('Estimated Pressure = ', params_T[0], '+/-', \
               np.sqrt(cov_params_T[0,0]))
        plt.xlabel('Container Area')
        plt.ylim(0, 40)
        plt.ylabel('Temperature')
        plt.title('Temperature vs. Container Area')
        plt.legend()
        plt.savefig('Temperature vs. Area with 50 balls 1000 collisions 5.0 v \
            width.png', dpi=1000)
        plt.show()
    
    return params_P[0], np.sqrt(cov_params_P[0,0]), params_T[0], \
           np.sqrt(cov_params_T[0,0])


#%%
def pressure_temperature_vs_n(
                              radius_cont, \
                              num_list, \
                              v_width, \
                              collisions, \
                              pressure_n = True, \
                              temp_n = True
                              ):
    '''
    This function generates graphs of Pressure and Temperature against the \
        Number of balls. 
    
    Args:
    =====
        radius_cont (float): Radius of the container.
            
        num_list (list of integers): List of number of balls.
        
        v_width (float): Width (standard deviation) of the normal distributiom \
            from which random initial velocities of balls are generated.
        
        collisions (int): Number of collisions / number of frames in the \
            simulation.
        
        pressure_n (bool): If True, it generates the graph of Pressure \
            vs. Number of balls. Default is True.
        
        temp_area (bool): If True, it generates the graph of Temperature \
            vs. Number of balls. Default is True.
    
    Returns:
    ========
        Estimated temperature of the gas for the Pressure vs. Number of balls \
            plot determined via least squares regression. (float)
        
        Absolute uncertainties in the estimated temperature (float): 
            Square roots of the first diagonal element in the covariance matrix.
        
        Estimated average pressure of the gas for the Temperature vs. Number \
            of balls plot determined via least squares regression. (float)
        
        Absolute uncertainties in the estimated average pressure (float): 
            Square roots of the first diagonal element in the covariance matrix.
    '''
    temp_list = []
    pressure_list = []
    for n in num_list:
        sim = Simulation(
                         radius_container = radius_cont, \
                         n_balls = n, \
                         velocity_width = v_width
                         )
        
        temp_eq, unc_temp_eq, pressure_eq, speeds = sim.run(
                                               num_frames = collisions, \
                                               animate = False
                                               )
        pressure_list.append(pressure_eq)
        temp_list.append(temp_eq)
    
    def linear_fit(N, m, c):
            return m * N + c
    
    if pressure_n == True:
        # Generate coefficients for linear fit.
        params_P, cov_params_P = curve_fit(
                                    linear_fit, \
                                    num_list, \
                                    pressure_list
                                    )

        # Plot simulation data and fit. 
        plt.plot(
                 np.array(num_list), \
                 np.array(pressure_list), \
                 'o', \
                 label = 'Simulation Data'
                 )
        x = np.linspace(num_list[0], num_list[-1])
        plt.plot(
                x, \
                linear_fit(x, *params_P), \
                label = 'Linear Fit'
                )
        
        plt.xlabel('Number of Balls')
        plt.ylabel('Pressure')
        plt.title('Pressure vs. Number of Balls')
        plt.legend()
        plt.savefig('Pressure vs. Number of Balls with 50 cont radius 1000 \
            collisions 5.0 v width.png', dpi=1000)
        print('m = ', params_P[0], '+/-', np.sqrt(cov_params_P[0,0]))
        print('c = ', params_P[1], '+/-', np.sqrt(cov_params_P[1,1]))
        plt.show()
        
    
    if temp_n == True:            
        # Generate coefficients for linear fit.
        params_T, cov_params_T = curve_fit(
                                    linear_fit, \
                                    num_list, \
                                    temp_list
                                    )

        # Plot simulation data and fit. 
        plt.plot(
                np.array(num_list), \
                np.array(temp_list), \
                'o', \
                label = 'Simulation Data'
                )
        x = np.linspace(num_list[0], num_list[-1])
        plt.plot(
                x, \
                linear_fit(x, *params_T), \
                label = 'Linear Fit'
                )
        
        plt.xlabel('Number of Balls')
        plt.ylabel('Temperature')
        plt.title('Temperature vs. Number of Balls')
        plt.legend()
        plt.savefig('Temperature vs. Number of Balls with 50 cont rad 1000 \
            collisions 5.0 v width.png', dpi=1000)
        print('m = ', params_T[0], '+/-', np.sqrt(cov_params_T[0,0]))
        print('c = ', params_T[1], '+/-', np.sqrt(cov_params_T[1,1]))
        plt.show()
    
    return params_P[0], np.sqrt(cov_params_P[0,0]), params_T[0], \
           np.sqrt(cov_params_T[0,0])



#%%
def speed_distribution(radius_cont, num, v_width, collisions):
    '''
    This function generates a histogram of gas pasticles' speeds for a given \
    initial width for particle velocities and fits a Maxwell-Boltzmann \
    distribution fit. It outputs a scale factor and the estimated temperature \
    of the gas according to the fit.
    
    Args:
    =====
        radius_cont (float): Radius of the container. 
            
        num (int): Number of balls.
        
        v_width (float): width (standard deviation) of the normal distributiom \
            from which random initial velocities of balls are generated.
        
        collisions (int): Number of collisions / number of frames in the \
            simulation.
    
    Returns:
    ========
        
        Scale factor (float): Scale factor / fit coefficient in the \
            Maxwell-Boltzmann fit determined from least sqaures regression.
        
        Absolute uncertainty in the scale factor (float): 
            Square root of the first diagonal element in the covariance matrix \
            corresponding to the fit coefficient.
        
        Estimated temperature (float): Temperature of the gas estimated from \
            the Maxwell-Boltzmann fit.
        
        Absolute uncertainty in the estimated temperature (float): 
            Square root of the second diagonal element in the covariance \
            matrix corresponding to the fit coefficient.
        
        Actual tempretaure calculated using the simulation (float).
        
        Standard error of the mean of the actual temperature (float).
    '''
    sim = Simulation(
                     radius_container = radius_cont, \
                     n_balls = num, \
                     velocity_width = v_width
                     )
            
    temp_eq, unc_temp_eq, pressure_eq, speeds = sim.run(
                                           num_frames = collisions, \
                                           animate = False
                                           )
    # Histogram of distances from centre.
    n, bins, patches = plt.hist(x = np.array(speeds), bins = 20)
    
    # Let k_B (Boltzmann constant) be equal to 1.
    mass = sim._mass_ball
    def maxwell_boltzmann(v, A, T):
        return A * v * np.exp(-(mass * v * v / (2 * T)))
    
    v = (bins[:-1] + bins[1:]) / 2
    # Source https://stackoverflow.com/questions/72688853/get-center-of-bins-histograms-python
    
    plt.xlabel('Speed of the Gas Particles')
    plt.ylabel('Number of Counts')
    plt.title('Histogram of Gas Particle Speeds')
    
    # Generate coefficients for Maxwell-Boltzmann fit.
    params, cov_params = curve_fit(maxwell_boltzmann, v, n)
    
    x = np.linspace(bins[0], bins[-1], 1000)
    plt.plot(x, maxwell_boltzmann(x, *params), label = 'Maxwell-Boltzmann Fit')
    plt.legend()
    plt.savefig('Maxwell Boltzmann with 1 ball rad, 50 container rad, 100 \
                 balls, 10 v width.png', dpi=1000) 
    plt.show()
    print('scale factor = ', params[0], '+/-', np.sqrt(cov_params[0,0]))
    print('Estimated temperature = ', params[1], '+/-', \
            np.sqrt(cov_params[0,0]))
    print('Actual temperature = ', temp_eq, '+/-', unc_temp_eq)
    
    return params[0], np.sqrt(cov_params[0,0]), params[1], \
           np.sqrt(cov_params[0,0]), temp_eq, unc_temp_eq


#%%
def van_der_waals(radius_cont, num, v_width_list, collisions):
    '''
    This function generates a best fit trendline on the graph of Pressure vs. \
    Temperature over a range of inital velocity widths accordning to the van \
    der Waal's law. 
     
    Note:
    Since the gas particles are colliding elastically, the 'a' term in the \
    van der Waals' equation corresponding to the average attraction between \
    particles is zero.
    
    Args:
    =====
        radius_cont (float): Radius of the container. 
            
        num (int): Number of balls.
        
        v_width (float): width (standard deviation) of the normal distributiom \
            from which random initial velocities of balls are generated.
        
        collisions (int): Number of collisions / number of frames in the \
            simulation.
    
    Returns:
    ========
        
        Constant b (float): Volume excluded by a single gas particle.
        
        Absolute uncertainty in b (float):
            Square root of the second diagonal element in the covariance \
            matrix corresponding to the fit coefficient.
    '''
     
    temp_list = []
    pressure_list = []
    
    for width in v_width_list:
        sim = Simulation(
                         radius_container = radius_cont, \
                         n_balls = num, \
                         velocity_width = width
                         )
        
        temp_eq, unc_temp_eq, pressure_eq, speeds = sim.run(
                                               num_frames = collisions, \
                                               animate = False
                                               )
        temp_list.append(temp_eq)
        pressure_list.append(pressure_eq)
        
    N = num
    # Contianer area as a proxy for volume.
    V = np.pi * radius_cont * radius_cont 
    def VDW_law(T, b):
        return T * N / (V - N * b)
    
    # Generate coefficients for van der Waals.
    params, cov_params = curve_fit(VDW_law, temp_list, pressure_list)

    # Plot simulation data and fit. 
    plt.plot(
             np.array(temp_list), \
             np.array(pressure_list), \
             'o', \
             label = 'Simulation Data'
             )
    x = np.linspace(temp_list[0], temp_list[-1])
    plt.plot(
             x, \
             VDW_law(x, *params), \
             label = "van der Waals' Fit"
             )
    
    print('b from van der Waals = ', params[0], '+/-', np.sqrt(cov_params[0,0]))
    plt.xlabel('Temperature')
    plt.ylabel('Pressure')
    plt.title('Van der Waals Law')
    plt.legend()
    plt.savefig('van der Waals fit with 100 balls 50 cont radius 1000 \
        collisions per data point.png', dpi=1000)
    plt.show()
    
    #return params[0], np.sqrt(cov_params[0,0]), temp_list, pressure_list