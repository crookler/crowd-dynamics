# RUN ONLY AFTER SIMULATION GSD HAS BEEN PLACED IN OUTPUTS FOLDER
import hoomd
import gsd.hoomd
import numpy as np
import matplotlib.pyplot as plt

number_escaped = []
frames = []

# Open hoomd file where inFile is path to your gsd file
with gsd.hoomd.open(name='outputs/wipSimulation.gsd', mode='r') as f:
    start = 0 #first frame to process
    dumps = int(f.__len__()) #final frame to process
    number_escaped = np.zeros(dumps-start)
    frames = np.arange(start, dumps, 1)
                                      
    box_data = f[0].configuration.box #get box dimensions
    lx_box = box_data[0] #box length
    ly_box = box_data[1] #box width 

    #determine types (assumes types remain constant)
    types = f[0].particles.typeid #particle type
    a_particles = np.transpose(np.where(f[0].particles.typeid==0)) #calculate which particles are type 0 (corresponds to 'A' or the person particle)
    w_particles = np.transpose(np.where(f[0].particles.typeid==1)) #calculate which particles are type 1 (corresponds to walls)                         
    wall_x_coordinate = f[0].particles.position[w_particles[0][0]][0] #assuming wall verticle straight line

    #loop over frames
    for frame in range(start, dumps):
        snapshot = f[frame] #current frame data
        current_escapes = 0

        #arrays of particle data
        positions = snapshot.particles.position #current positions
        positions = positions[:,:-1] #splice off z dimension

        for particle in a_particles:
            if (positions[particle[0]][0] > wall_x_coordinate):
                current_escapes += 1
        
        number_escaped[frame] = current_escapes

plt.plot(number_escaped, frames)
plt.show()




        