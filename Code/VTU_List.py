import numpy

class VTU_List:
    def __init__(self, timeSteps, numberOfParticles):
        self.positions = numpy.zeros([numberOfParticles,timeSteps,3])
        self.velocities = numpy.zeros([numberOfParticles,timeSteps,3])
        self.rad = numpy.zeros([numberOfParticles])
    def UpdateValues():
        pass