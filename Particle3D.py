"""
 CMod Ex3: Particle3D, a class to describe 3D particles
"""

import numpy as np

class Particle3D(object):
    """
    Class to describe 1D particles.

    Properties:
    position(numpy array) - position vector
    velocity(numpy array) - velocity vector
    mass(float) - particle mass

    Methods:
    * formatted output
    * kinetic energy
    * first-order velocity update
    * first- and second order position updates
    """

    def __init__(self, pos, vel, mass,label):
        """
        Initialise a Particle3D instance

        :param pos: position as numpy array
        :param vel: velocity as numpy array
        :param mass: mass as float
        """

        self.position = pos
        self.velocity = vel
        self.mass = mass
        self.label=label

    def __str__(self,label):
        """
        Define output format.
        For particle p=([1.0,2.0,3.0], [4.0,5.0,6.0], 7.0,barry) this will print as
        "label=barry, x = 1.0, y=2.0, z=3.0, m = 7.0"
        """

        return "x = " + str(self.position[0]) + "y = " + str(self.position[1]) + "z = " + str(self.position[2]) + " m = " + str(self.mass)

    def kinetic_energy(self):
        """
        Return kinetic energy as
        1/2*mass*vel^2
        """

        return 0.5*self.mass*np.inner(self.velocity,self.velocity)

    # Time integration methods

    def leap_velocity(self, dt, force):
        """
        First-order velocity update,
        v(t+dt) = v(t) + dt*F(t)

        :param dt: timestep as float
        :param force: force on particle as float
        """

        self.velocity = self.velocity + (dt*force)/self.mass

    def leap_pos1st(self, dt):
        """
        First-order position update,
        x(t+dt) = x(t) + dt*v(t)

        :param dt: timestep as float
        """

        self.position = self.position + dt*self.velocity

    def leap_pos2nd(self, dt, force):
        """
        Second-order position update,
        x(t+dt) = x(t) + dt*v(t) + 1/2*dt^2*F(t)

        :param dt: timestep#read in first line from input file
    constants_1=infile.readline()
    #split first line up
    tokens1=constants.split(",")

    #set up initial parameters for first particle
    D=tokens1[0]
    alpha=tokens1[2]
    r_e=tokens1[1] as float
        :param force: current force as float
        """

        self.position = self.position + dt*self.velocity + 0.5*dt**2*force/self.mass

    #Static Method that takes info from a file and creates a particle

    def file_read(file_handle):
        #read in one line from the file
        particle_info=file_handle.readline()
        #split line up
        tokens=particle_info.split(",")
        #assign properties of the particle to one of the tokens
        label=str(tokens[7])
        pos=np.array([float(tokens[0]),float(tokens[1]),float(tokens[2])])
        vel=np.array([float(tokens[3]),float(tokens[4]),float(tokens[5])])
        mass=float(tokens[6])

        return Particle3D(pos,vel,mass,label)




    #Static method to return vectorseparation of two particles

    def separation(p1,p2):
        """
        :param p1: first particle
        :param p2: second particle
        """
        return p1.position-p2.position
