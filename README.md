# Thermodynamic-Simulation
A 2nd year computing project for the undergraduate physics course at Imperial College London.

Investigating the Properties of Ideal Gases Using a 2D Thermodynamic Simulation:
================================================================================

Note to Users:
--------------
- The codes within these files were created using Visual Studio Code which 
  displays graphics in a separate window by default. To generate animation on 
  Spyder IDE, choose the inline or windowed graphics option by entering the 
  IPython magic command:

	%matplotlib inline
	or 
	%matplotlib auto

  in the console. 

- If the animation still does not work, try changing the graphics backend:

	- Windows 10: 
	  Tools -> Preferences -> IPython console -> graphics -> graphics backend
	- macOS: 
	  python -> Preferences -> IPython console -> graphics -> graphics backend

- Change:
	- UG PCs (Windows 10 + Spyder 3.3.1): use Tkinter.
	- Windows 10 + Spyder 4.1.4: use Automatic or Tkinter.
	- macOS + Spyder 4.1.4: use Automatic or OSX.

- To update changes, restart python kernel.


ABOUT THE PROJECT:
------------------
This repository contains Python module files which can be used to simulate ideal gas 
collisions inside a rigid container in 2D and generate an animation of collision 
events to provide visual instuition and numerical results for various 
relationships. This folder contains .py files with the following modules:


Modules:
--------
- Ball.py:
	Contains the Ball class which initialises the ball and container objects. 
	Ball class includes methods that can access and update the position an 
	velocities of the ball objects, as well as a method that calculate the 
	time until two objects collide, a method which moves objects by specified 
	time, and a method that can update the ball obejcts' velocities after 
	collision. Full list of methods can be seen in the docstrings.

- Simulation.py:
	Contains the Simulation class which initialises the simulation object by 
	initialising a container object and multiple ball objects with randomly 
	generated initial positions and velocities. Simulation class includes 
	methods that can generate patch objects to animate the collision of balls 
	inside the container, a method which finds and executes the next collision 
	in the simualtion, a method which iterates collision for a given number 
	of frames and outputs animation as well as figures. Full list of methods 
	can be seen in the docstrings.
	
- Figure.py:
	Contains functions which initialise the Simulation class and plot simulation 
	data with least squares fits to characterise the relationships of ideal 
	gases. It can be used to generate a plot of Pressure vs. Temperature for a 
	range of particle radii, a histogram of particle speeds with a 
	Maxwell-Boltzmann fit, and generate a fit based on van der Waals' equation 
	on a Pressure vs. Temperature graph. Full list of methods can be seen in 
	the docstrings.

- Investigations Script.py:
	Contains Scripts for running the animation and plotting figures for the 
	investigation section of the project.

- Tests.py:
	Contains Scrips for unit testing and integrating testing of various methods 
	and classes, making sure that they work both individually and together in a 
	sequence.


Documentation:
--------------
Run the help() command to read the docstrings for each class/method/function 
(which will provide you with more detailed breakdowns of their functionalities). 


Installation and Usage:
-----------------------
To conduct the investigation from scratch, import relevant modules (Ball.py, 
Simulation.py, and Figure.py) into a .py file, initialise the Simulation class, 
and call the relevant methods/functions.

