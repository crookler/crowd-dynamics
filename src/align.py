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
        for particle in range(self.numParticles):
            coordinates = (snapshot.particles.position[particle][0], snapshot.particles.position[particle][1]) #create tuple of x and y coordinate
            
            xToDoor = self.doorX-coordinates[0]
            yToDoor = self.doorY-coordinates[1]
            newAngle = np.arctan(yToDoor/xToDoor)
            snapshot.particles.orientation[particle] = [np.cos(newAngle/2), 0, 0, np.sin(newAngle/2)]
            
        self._state.set_snapshot(snapshot)