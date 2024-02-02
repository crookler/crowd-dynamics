import hoomd
import numpy as np

class align(hoomd.custom.Action):
    def __init__(self, numParticles, doorX, doorY):
       self.numParticles = numParticles
       self.doorX = doorX
       self.doorY = doorY

    def act(self, timestep):
        snapshot = self._state.get_snapshot()

        #print(snapshot.particles.orientation[0])
        #for particle in range(self.numParticles):
        #    snapshot.particles.orientation[particle][0] = 0

        self._state.set_snapshot(snapshot)