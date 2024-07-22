import numpy
import math

import settings


class Standaard_Body:
    def __init__(self, mass, pos, vel, acc, force, index, radius = 1):
        self.mass = mass
        self.radius = radius

        self.pos = pos
        self.posOld = ([0,0,0])
        self.vel = vel
        self.velOld = ([0,0,0])
        self.acc = acc
        self.accOld = ([0,0,0])
        self.force = force

        self.box = []

        self.index = index
        self.FirstTimeStep = True

    def Update_acc(self):
        self.accOld = numpy.copy(self.acc)
        self.acc = self.force/self.mass

    def Leap_Frog(self):
        self.pos += self.vel * settings.DELTA_TIME + 0.5 * self.accOld * settings.DELTA_TIME ** 2
        self.vel += 0.5 * (self.accOld + self.acc) * settings.DELTA_TIME

    def LF_K_and_D(self):
        self.vel += self.accOld * 0.5 * settings.DELTA_TIME
        self.pos += self.vel * settings.DELTA_TIME
        self.vel += self.acc * 0.5 * settings.DELTA_TIME
        
    def Euler(self):
        self.pos += self.vel * settings.DELTA_TIME
        self.vel += self.acc * settings.DELTA_TIME

    def Semi_Impl_Euler(self):
        self.vel += self.acc * settings.DELTA_TIME
        self.pos += self.vel * settings.DELTA_TIME

    def Verlet(self):
        tempPos = numpy.copy(self.pos)
        self.pos = - self.oldPos + 2 * self.pos + self.acc * settings.DELTA_TIME ** 2
        self.vel = (self.pos - self.oldPos) / (2 * settings.DELTA_TIME)
        self.oldPos = tempPos

    def Velocity_Verlet(self):
        self.vel += ((self.accOld + self.acc) / 2) * settings.DELTA_TIME
        self.pos += self.vel * settings.DELTA_TIME + 0.5 * self.acc * settings.DELTA_TIME ** 2

    def Direct_Velocity_verlet(self):
        if self.FirstTimeStep == False:
            self.vel = self.pos - self.oldPos
        else:
            self.vel *= settings.DELTA_TIME
            self.FirstTimeStep = False
        self.oldPos = numpy.copy(self.pos)
        self.pos += self.vel + self.acc * settings.DELTA_TIME ** 2
            
    def Reset_force(self):
        self.force = numpy.double([0,0,0])

    def Update_Position(self,method):
        self.Update_acc()
        if method == settings.INTEGRATION_OPTIONS[0]:
            self.Leap_Frog()
        elif method == settings.INTEGRATION_OPTIONS[1]:
            self.LF_K_and_D()
        elif method == settings.INTEGRATION_OPTIONS[2]:
            self.Euler()
        elif method == settings.INTEGRATION_OPTIONS[3]:
            self.Semi_Impl_Euler()
        elif method == settings.INTEGRATION_OPTIONS[4]:
            self.Verlet()
        elif method == settings.INTEGRATION_OPTIONS[5]:
            self.Velocity_Verlet()
        elif method == settings.INTEGRATION_OPTIONS[6]:
            self.Direct_Velocity_verlet()
        else:
            exit(-1)
        self.Reset_force()

class Astroid(Standaard_Body):
    def __init__(self, mass, pos, vel, acc, force, index, radius = 1):
        super().__init__(mass, pos, vel, acc, force, index, radius = 1)
        self.color = 'brown'
        self.radius = radius

class Planet(Standaard_Body):
    def __init__(self, mass, pos, vel, acc, force, index, color, radius = 1):
        super().__init__(mass, pos, vel, acc, force, index, radius = 1)
        self.color = color
        self.radius = radius

class Sun(Standaard_Body):
    def __init__(self, mass, pos, vel, acc, force, index, radius = 1):
        super().__init__(mass, pos, vel, acc, force, index, radius = 1)
        self.color = 'gold'
        self.radius = radius