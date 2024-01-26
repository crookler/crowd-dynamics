##### import python packages
import hoomd
import gsd.hoomd
import numpy as np
import datetime, time

##### helper functions to display simulation status
class print_sim_state(hoomd.custom.Action):
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

#System variable definitions
sigma = 1     # sigma is diameter of particles
eps = 1       # epsilon is a constant affecting Lennard Jones function steepness/energy
kT = 1        # kT is the system thermal energy
dt = 1e-6     # dt is time step size
N = 100       # N is number of particles

#Set up wall and particles
particlePositions = []     
for i in range(10):        
  for j in range(10):
    particlePositions.append([i-12+0.5, j-5+0.5, 0])

wall_particles = []
for i in range(16):
  if i < 7 or i > 8:
    wall_particles.append([0, i-15/2, 0])

s = gsd.hoomd.Frame() #initial frame
s.particles.position = list(particlePositions) + wall_particles
s.particles.N = N + len(wall_particles) 
s.particles.typeid = ([0] * N + [1] * len(wall_particles)) #particles types (related to above indices)
s.particles.types = ["A", "W"] #A is particle and W is wall
s.configuration.box = [30,15,0,0,0,0]
with gsd.hoomd.open(name='wipIC.gsd', mode='w') as f:
    f.append(s)


#Setting up HOOMD simulation object
sim = hoomd.Simulation( # define simulation object and specify device and seed
    device=hoomd.device.auto_select(),
    seed=1
)
sim.create_state_from_gsd(filename='wipIC.gsd') # load initial state

#Custom movement
fire = hoomd.md.force.Active(
    filter = hoomd.filter.Type(['A'])
)
fire.active_force['A'] = (30, 0, 0) #A interacting with itself question mark?

#Setting up particle interactions
lj = hoomd.md.pair.LJ(
    nlist=hoomd.md.nlist.Cell(buffer=0.5),
    default_r_cut=0.75
)
lj.params[('A', 'A')] = dict(epsilon=eps, sigma=sigma) #lj between particles [A,A]
lj.params[('A', 'W')] = dict(epsilon=eps, sigma=sigma)  #lj between particles and wall [A,W]
lj.params[('W', 'W')] = dict(epsilon=0, sigma=sigma) #lj between walls (set to zero by epsilon)

sim.operations.integrator = hoomd.md.Integrator(
    dt = dt,
    methods=[
        hoomd.md.methods.Brownian(
            filter = hoomd.filter.Type(['A']), #just integrate particles (walls are fixed)
            kT = 1 * kT, #system thermal energy (influences particle speed)                                                      
        )
    ],
    forces = [lj, fire]
)
update_every_n = 5
sim.operations.writers.append(
    hoomd.write.CustomWriter(
        action = print_sim_state(),
        trigger = trigger_every_n_sec()
    )
)
gsd_writer = hoomd.write.GSD(
    trigger = hoomd.trigger.Periodic(int(1e4)),
    filename = "wipSimulation.gsd",
    mode = 'wb',
    filter = hoomd.filter.All(),
    dynamic=['property', 'momentum', 'attribute']
)
sim.operations.writers.append(gsd_writer)

#Run simulation
init_time = time.time()                              
last_output = init_time                            
sim.run(500_000) #simulation length
gsd_writer.flush()




