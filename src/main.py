import hoomd
import gsd.hoomd
import numpy as np
import datetime, time
from align import *

#helper functions to display simulation status
class printStatus(hoomd.custom.Action):
    def act(self, timestep):
        global init_time, last_output
        try:
            sec_remaining = int((self._state._simulation.final_timestep - timestep) / self._state._simulation.tps)
        except ZeroDivisionError:
            sec_remaining = 0
        print(
            "Time", str(datetime.timedelta(seconds=int(time.time() - init_time))),
            "| Step:", timestep,
            "| TPS:", str(round(float(self._state._simulation.tps),3)),
            "| ETR:", str(datetime.timedelta(seconds=sec_remaining))
        )
        last_output = time.time()

class trigger_every_n_sec(hoomd.trigger.Trigger):
    def __init__(self):
        hoomd.trigger.Trigger.__init__(self)

    def compute(self, timestep):
        global last_output, update_every_n
        return ((time.time() - last_output) >=update_every_n)

#system variable definitions
diameter = 1 #diameter of particles
dt = 1e-6 #time step size
N = 100 #number of particles
W = 14 #number of walls

#molecular definitions (come back to this later if needed)
eps = 1 # epsilon is a constant affecting Lennard Jones function steepness/energy
kT = 1 #system thermal energy

particlePositions = [] #coordinates for each particle
doorPosition = (0, 0, 0) #coordinates for door

#set up p particles 
for i in range(10):        
  for j in range(10):
    particlePositions.append([i-12+0.5, j-5+0.5, 0])

#set of w particles
for i in range(16):
  if i < 7 or i > 8:
    particlePositions.append([0, i-15/2, 0])

start = gsd.hoomd.Frame() #initial frame
start.particles.position = list(particlePositions) 
start.particles.N = N + W
start.particles.typeid = ([0] * N + [1] * W) #particles types (related to above indices)
start.particles.types = ["A", "W"] #A is particle and W is wall and B is blackhole
start.configuration.box = [30,15,0,0,0,0]
with gsd.hoomd.open(name='outputs/wipIC.gsd', mode='w') as f:
    f.append(start)

#Setting up HOOMD simulation object
simulation = hoomd.Simulation(device=hoomd.device.auto_select(), seed=1) #create simulation object
simulation.create_state_from_gsd(filename='outputs/wipIC.gsd') # load initial state

#Custom movement
fire = hoomd.md.force.Active(filter = hoomd.filter.All())
fire.active_force['A'] = (30, 0, 0) #accerlates back to this maximum (after collision)

#Setting up particle interactions
collision = hoomd.md.pair.LJ(
    nlist=hoomd.md.nlist.Cell(buffer=0.5),
    default_r_cut=0.75 #stop applying force at 0.75
)
collision.params[('A', 'A')] = dict(epsilon=eps, sigma=diameter) #collision between particles [A,A]
collision.params[('A', 'W')] = dict(epsilon=eps, sigma=diameter)  #collision between particles and wall [A,W]
collision.params[('W', 'W')] = dict(epsilon=0, sigma=diameter) #collision between walls (set to zero by epsilon)

randomness = hoomd.md.methods.Brownian(filter=hoomd.filter.Type(['A']), kT = 1*kT)

integrator = hoomd.md.Integrator(dt = dt) #define an integrator
integrator.methods = [randomness] #add random Brownian motion to particle movements
integrator.forces = [fire,collision] #add forces to integrator
simulation.operations.integrator = integrator  #put integrator in sim object

update_every_n = 5
simulation.operations.writers.append(hoomd.write.CustomWriter(action = printStatus(),trigger = trigger_every_n_sec()))

gsd_writer = hoomd.write.GSD(trigger = hoomd.trigger.Periodic(int(1e4)), filename = "outputs/wipSimulation.gsd", mode = 'wb', filter = hoomd.filter.All(), dynamic=['property', 'momentum', 'attribute'])
simulation.operations.writers.append(gsd_writer)

alignmentUpdate = hoomd.update.CustomUpdater(action=align(numParticles=N, doorX=doorPosition[0], doorY=doorPosition[1]),trigger=hoomd.trigger.Periodic(1000)) #trigger alignment every n timesteps
simulation.operations.updaters.append(alignmentUpdate)

#Run simulation
init_time = time.time()                              
last_output = init_time                            
simulation.run(1_000_000) #simulation length
gsd_writer.flush()





