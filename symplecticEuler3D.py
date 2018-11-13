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
    r_12=m.sqrt(np.inner(Particle3D.separation(p1,p2),Particle3D.separation(p1,p2)))
    #the force on p1 due to p2
    force = (float(-2.0)*float(alpha)*float(D)*(float(1.0)-m.exp(-float(alpha)*(float(r_12)-float(r_e))))*(m.exp(-float(alpha)*(float(r_12)-float(r_e))))*Particle3D.separation(p1,p2))/float(r_12)
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
    #create a list to hold the updated positions
    updated_position=[]

    #create a list to hold the updated velocities
    updated_velocity=[]

    #calculate new position of each particle in the list "particles"
    for i in range(len(particles)):
        new_pos=particles[i].leap_pos1st(dt)
        updated_position.append(new_pos)

    for i in range(len(particles)):
        #for each ith updated particle, calculate the total force on that particle by summing pairwise over all other particles
        for j in range(len(particles)):
            #don't calculate the force on the particle due to itself
            if i!=j:
                force=force_morse(updated_position[i],particles[j], D, alpha,r_e)
        new_vel=particles[i].leap_velocity(dt,force)
        updated_velocity.append(new_vel)

    #create a list of updated separations between particles
    separation_list=[]
    #create a list to hold the pairwise potential energies
    potential_list=[]

    #loop over all particles calculating distance between each pair of particles and their pairwise potential, making sure not to double count
    for i in range(len(particles)):
        for j in range(len(particles)):
            if i<j:
                dist=updated_position[i].append(updated_position[j])
                potential_energy= pot_energy_morse(dist, D, alpha,r_e)
                separation_list.append(dist)
                potential_list.append(potential_energy)


    return updated_position
    return updated_velocity
    return separation_list
    return potential_list






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
    D2=tokens2[0]
    alpha2=tokens2[2]
    r_e2=tokens2[1]

    dt = 0.01
    numstep = 2000
    time = 0.0

    #set up the initial state of the particles

    p1=Particle3D.file_read(infile)
    p2=Particle3D.file_read(infile)

    #create a list of particles p1,p2
    particles=[p1,p2]

    #calculate the magnitude of the separation between p1 and p2
    distance=m.sqrt(np.inner(Particle3D.separation(p1,p2),Particle3D.separation(p1,p2)))

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
        updated_info=[symp_euler(particles,dt,D,alpha,r_e)]
        updated_position=updated_info[0]
        updated_velocity=updated_info[1]
        updated_separation=updated_info[2]
        updated_potential=updated_info[3]

        #calculate new total potential energy
        new_total_pot=sum(i for i in updated_potential)

        #calculate new total kinetic energy
        new_tot_kin=sum(0.5*particles[i].mass*np.inner(updated_velocity[i],updated_velocity[i]) for i in len(particles))

        #calculate new total energy
        new_tot_energy=new_total_pot+new_tot_kin

        # Output particle information
        outfile.write("{0:f} {1:f} {2:12.8f}\n".format(time,(updated_separation[i] for i in len(particles)),new_tot_energy))

        # Append information to data lists
        time_list.append(time)
        distance_list.append(updated_separation[i] for i in len(particles))
        energy_list.append(new_tot_energy)

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
