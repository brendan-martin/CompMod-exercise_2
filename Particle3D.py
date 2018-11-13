"""Particle3D Class"""

Class Particle3D(object):
    """
    Class to describe 3D particles

    Properties:
    x(float)-position along the x direction
    y(float)-position along the y direction
    z_position(float)-position along the z direction

    v_x(float)-velocity along the x direction
    v_y(float)-velocity along the y direction
    v_z(float)-velocity along the z direction

    mass(float)-particle mass

    Methods:
    *formatted output
    *kinetic energy
    *first order velocity update
    *first and second order position updates
    """

    def __init__(self,x_position,y_position,z_position,x_velocity,y_velocity,z_velocity,mass):
        """
        Initialise a Particle3D instance

        :param x_position: x position as float
        :param y_position: x position as float
        :param z_position: x position as float
        :param x_velocity: x velocity as float
        :param y_velocity: y velocity as float
        :param z_velocity: z velocity as float
        :param mass: mass as float
        """

        self.x=x_position
        self.y=y_position
        self.z=z_position
        self.v_x=x_velocity
        self.v_y=y_velocity
        self.v_z=z_velocity
        self.m=mass

    def __str__(self):
        """
        Define output format
        """

        return "x="+str(self.x)+"y="+str(self.y)+"x="+str(self.y)+"v_x="+str(self.v_x)+"v_y="+str(self.v_y)+"v_z="+str(self.v_z)+"m="+str(self.m)

    def kinetic_energy(self):
        """
        Return kinetic energy as
        1/2*mass*vel^2
        """

        return 0.5*self.m*((self.v_x**2)+(self.v_y**2)+(self.v_z**2))

    #Time integration Methods

     def leap_velocity(self, dt, force_x,force_y,force_z):
        """
        First-order velocity update,
        v(t+dt) = v(t) + dt*F(t)

        :param dt: timestep as float
        :param force_x: force on particle in the x direction as float
        :param force_y: force on particle in the y direction as float
        :param force_z: force on particle in the z direction as float
        """

        self.v_x = self.v_x + dt*force_x/self.m
        self.v_y = self.v_y + dt*force_y/self.m
        self.v_z = self.v_z + dt*force_z/self.m

     def leap_pos1st(self, dt):
        """
        First-order position update,
        x(t+dt) = x(t) + dt*v(t)

        :param dt: timestep as float
        """

        self.x = self.x + dt*self.v_x
        self.y = self.y + dt*self.v_y
        self.z = self.z + dt*self.v_z

    def leap_pos2nd(self, dt, force_x,force_y,force_z):
        """
        Second-order position update,
        x(t+dt) = x(t) + dt*v(t) + 1/2*dt^2*F(t)

        :param dt: timestep as float
        :param force_x: current force in x direction as float
        :param force_y: current force in y direction as float
        :param force_z: current force in z direction as float
        """

        self.x = self.x + dt*self.v_x + 0.5*dt**2*force_x/self.m
        self.y = self.y + dt*self.v_y + 0.5*dt**2*force_y/self.m
        self.z = self.z + dt*self.v_z + 0.5*dt**2*force_z/self.m
