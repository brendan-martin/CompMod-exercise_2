"""
CompMod Exercise_2

Symplectic Euler time integration of a particle moving in a 3D Morse Potential

Produces plots of the position of the particle
and its energy, both as function of time. Also
saves both to file.

"""

import sys
import math as m
import numpy as np
import matplotlib.pyplot as pyplot
from Particle3D import Particle3D

def force_morse(p1,p2, D, alpha,r_e):
    """
    Method to return the force on a particle
    in a morse potential.
    Force is given by
    F(p1,p2) = (-2*alpha*D*(1-m.exp(-alpha*(r_12-r_e)))*(m.exp(-alpha*(r_12-r_e)))*Particle3D.separation(p1,p2))/r_12

    :param p1: Particle3D instance 1
    :param p2: Particle3D instance 2
    :param D: parameter D from morse potential
    :param alpha: parameter alpha from morse potential
    :param r_e: parameter r_e from morse potential
    :return: force acting on particle p1 due to p2 as Numpy array
    """

    #the magnitude of the particles' separation
    r_12=np.linalg.norm(Particle3D.separation(p1,p2))
    #the force on p1 due to p2
    coeff = (-2.0*float(alpha)*float(D)*(1.0-m.exp(-1.0*float(alpha)*(float(r_12)-float(r_e))))*(m.exp(-1.0*float(alpha)*(float(r_12)-float(r_e)))))*(1.0/(float(r_12)))
    force=float(coeff)*Particle3D.separation(p1,p2)
    return force

def pot_energy_morse(r_12, D, alpha,r_e):
    """
    Method to return the potential energy of particle
    in a morse potential.
    Potential is given by
    U(p1,p2) = D*(((1-m.exp(-alpha*(r_12-r_e)))**2)-1)

    :param r_12: the distance between the two particles for which we are finding the pair potential
    :param D: parameter D from morse potential
    :param alpha: parameter alpha from morse potential
    :param r_e: parameter r_e from morse potential
    :return: potential energy of particle p1 due to morse potential as float

    """

    potential=float(D)*(((float(1.0)-m.exp(-float(alpha)*(float(r_12)-float(r_e))))**float(2.0))-float(1.0))
    return potential


def symp_euler(particles,dt,D,alpha,r_e):

    #calculate new position of each particle in the list "particles"
    for i in range(len(particles)):
        particles[i].leap_pos1st(dt)

    #calculate new separations of particles (pairwise) and return as a list, called separations
    #separations=map(lambda x: x.separation(particles[i]),filter(lambda x: x != particles[i],particles))
    separations=[np.linalg.norm(i.separation(j)) for i in particles for j in particles if i!=j]

    for i in range(len(particles)):
        force=0
        #for each ith updated particle, calculate the total force on that particle by summing pairwise over all other particles
        for j in range(len(particles)):
            #don't calculate the force on the particle due to itself
            if i!=j:
                force=force+force_morse(particles[i],particles[j], D, alpha,r_e)

        particles[i].leap_velocity(dt,force)


    for i in range(len(particles)):
        force=sum(map(lambda x: force_morse(particles[i],x, D, alpha,r_e), filter(lambda x: x != particles[i], particles)))
        particles[i].leap_velocity(dt,force)

    new_pot=0

    #loop over all particles calculating distance between each pair of particles and their pairwise potential, making sure not to double count
    for i in range(len(particles)):

        for j in range(len(particles)):

            if i<j:
                r=np.linalg.norm(particles[i].separation(particles[j]))

                new_pot=new_pot+pot_energy_morse(r, D, alpha,r_e)

    #calculate the total energy after the time step
    total_energy=new_pot+sum(map(lambda x: x.kinetic_energy(), particles))

    #return the list separations and the total enery as a list
    return [separations,total_energy]







def main():

    if len(sys.argv)!=3:
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] + "<input file>" + " <output file>")
        quit()
    else:
        outfile_name = sys.argv[2]
        infile_name=sys.argv[1]

    # Open output file
    outfile = open(outfile_name, "w")

    #read in input file
    infile=open(infile_name, "r")



    #Get the specific constants for the particle from the input file
    #read in first line from input file
    constants_1=infile.readline()
    #split first line up
    tokens1=constants_1.split(",")

    #set up initial parameters for first particle
    D=tokens1[0]
    alpha=tokens1[2]
    r_e=tokens1[1]

    #read in second line from input file
    constants_2=infile.readline()
    #split second line up
    tokens2=constants_2.split(",")

    #set up initial parameters for second particle
    D_2=tokens2[0]
    alpha_2=tokens2[2]
    r_e2=tokens2[1]

    dt = 0.01
    numstep = 5000
    time = 0.0

    #set up the initial state of the particles

    p1=Particle3D.file_read(infile)
    p2=Particle3D.file_read(infile)

    #create a list of particles p1,p2
    particles=[p1,p2]

    sep_vector=p1.separation(p2)

    #calculate the magnitude of the separation between p1 and p2
    distance=np.linalg.norm(sep_vector)

    #get the total initial energy of the system
    energy=Particle3D.kinetic_energy(p1)+Particle3D.kinetic_energy(p2)+pot_energy_morse(distance, D, alpha,r_e)

    #write the initial time, initial separation and initial total energy to the output file
    outfile.write("{0:f} {1:f} {2:12.8f}\n".format(time,distance,energy))

    # Initialise data lists for plotting later
    time_list = []
    distance_list = []
    energy_list = []



    # Start the time integration loop
    for i in range(numstep):
        #use symplectic euler to get the new separations and new total energy after the time step
        new_data=symp_euler(particles,dt,D,alpha,r_e)

        #get list of new separations
        sep_list=new_data[0]

        #get first element of sep_list
        sep=float(sep_list[0])

        #get  new total energy
        new_energy=new_data[1]


        # Output particle information
        outfile.write("{0:f} {1:f} {2:12.8f}\n".format(time,sep,new_energy))

        # Append information to data lists
        time_list.append(time)
        distance_list.append(sep)
        energy_list.append(new_energy)

        # Increase time
        time = time + dt

    # Post-simulation:

    # Close output file
    outfile.close()

    # Plot particle trajectory to screen
    pyplot.title('Symplectic Euler: separation vs time')
    pyplot.xlabel('Time')
    pyplot.ylabel('Separation')
    pyplot.plot(time_list, distance_list)
    pyplot.show()

    # Plot particle energy to screen
    pyplot.title('Symplectic Euler: total energy vs time')
    pyplot.xlabel('Time')
    pyplot.ylabel('Energy')
    pyplot.plot(time_list, energy_list)
    pyplot.show()


# Execute main method:
main()
