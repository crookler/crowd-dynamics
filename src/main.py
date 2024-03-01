import hoomd
import gsd.hoomd
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
P = 100 #number of particles
W = 14 #number of walls

#system environment definitions
length = 30
width = 15
main_wall_x = 0 #assuming straight vertical wall
external_walls = [hoomd.wall.Plane(origin=(-length/2,0,0), normal=(1,0,0)),
                  hoomd.wall.Plane(origin=(length/2,0,0), normal=(-1,0,0)),
                  hoomd.wall.Plane(origin=(0,-width/2,0), normal=(0,1,0)),
                  hoomd.wall.Plane(origin=(0,width/2,0), normal=(0,-1,0))]

particlePositions = [] #coordinates for each particle
doorPosition = (0, 0, 0) #coordinates for door

#molecular definitions (come back to this later if needed)
eps = 1 # epsilon is a constant affecting Lennard Jones function steepness/energy
kT = 1 #system thermal energy

#set up p particles 
for i in range(10):        
  for j in range(10):
    particlePositions.append([i-12+0.5, j-5+0.5, 0]) #place manually for now

#set of w particles
for i in range(16):
  if i < 7 or i > 8:
    particlePositions.append([main_wall_x, i-width/2, 0])

start = gsd.hoomd.Frame() #initial frame
start.particles.position = particlePositions
start.particles.N = P + W
start.particles.typeid = ([0] * P + [1] * W) #particles types (related to above indices)
start.particles.types = ["A", "W"] #A is particle and W is wall
start.configuration.box = [length,width,0,0,0,0] #length, width, height, 3D origin coordinates
start.configuration.dimensions = 2
with gsd.hoomd.open(name='outputs/wipIC.gsd', mode='w') as f:
    f.append(start)

#Setting up HOOMD simulation object
simulation = hoomd.Simulation(device=hoomd.device.auto_select(), seed=1) #create simulation object
simulation.create_state_from_gsd(filename='outputs/wipIC.gsd') # load initial state

#Setting up boundary conditions
boundaries = hoomd.md.external.wall.LJ(walls=external_walls)
boundaries.params[('A')] = dict(epsilon=eps, sigma=diameter, r_cut=0.75)
boundaries.params[('W')] = dict(epsilon=0, sigma=diameter, r_cut = 0.75)

#Setting up particle interactions
collision = hoomd.md.pair.LJ(
    nlist=hoomd.md.nlist.Cell(buffer=0.5),
    default_r_cut=0.75 #stop applying force at 0.75
)
collision.params[('A', 'A')] = dict(epsilon=eps, sigma=diameter) #collision between particles [A,A]
collision.params[('A', 'W')] = dict(epsilon=eps, sigma=diameter)  #collision between particles and wall [A,W]
collision.params[('W', 'W')] = dict(epsilon=0, sigma=diameter) #collision between walls (set to zero by epsilon)

#Custom movement
fire = hoomd.md.force.Active(filter = hoomd.filter.All())
fire.active_force['A'] = (30, 0, 0) #accerlates back to this maximum (after collision)

randomness = hoomd.md.methods.Brownian(filter=hoomd.filter.Type(['A']), kT = 1*kT)

integrator = hoomd.md.Integrator(dt = dt) #define an integrator
integrator.methods = [randomness] #add random Brownian motion to particle movements
integrator.forces = [boundaries,collision,fire] #add forces to integrator
simulation.operations.integrator = integrator  #put integrator in sim object

update_every_n = 5
simulation.operations.writers.append(hoomd.write.CustomWriter(action = printStatus(),trigger = trigger_every_n_sec()))

gsd_writer = hoomd.write.GSD(trigger = hoomd.trigger.Periodic(int(1e4)), filename = "outputs/wipSimulation.gsd", mode = 'wb', filter = hoomd.filter.All(), dynamic=['property', 'momentum', 'attribute'])
simulation.operations.writers.append(gsd_writer)

alignmentUpdate = hoomd.update.CustomUpdater(action=align(numParticles=P, doorX=doorPosition[0], doorY=doorPosition[1]),trigger=hoomd.trigger.Periodic(1000)) #trigger alignment every n timesteps
simulation.operations.updaters.append(alignmentUpdate)

#Run simulation
init_time = time.time()                              
last_output = init_time                            
simulation.run(1_000_000) #simulation length
gsd_writer.flush()





