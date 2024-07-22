import math
import numpy
import random
import matplotlib.pyplot as plt

import FunctionLibrary
import Bodies
import VTU_List
import OctTree_V2

import settings


def LoadPreset(method):
    if method == settings.PRESETS[0]:
        List = ThreeBody()
    elif method == settings.PRESETS[1]:
        List = AstroidProblem()
    else:
        exit(-1)
    return List

def ThreeBody():
    List = []
    for i in range(1):
        List.append(Bodies.Planet( mass = settings.MASS_MAGNIFIER * 1, 
                                    pos = numpy.longdouble([-10,0,0]), 
                                    vel = settings.SPEED_MAGNIFIER * numpy.longdouble([0,-1,-2]), 
                                    acc = numpy.longdouble([0,0,0]), 
                                    force = numpy.double([0,0,0]),
                                    index = 1,
                                    color = 'royalblue',
                                    radius = 10))
        List.append(Bodies.Planet( mass = settings.MASS_MAGNIFIER * 2, 
                                    pos = numpy.longdouble([0,0,0]), 
                                    vel = settings.SPEED_MAGNIFIER * numpy.longdouble([0,0,2]), 
                                    acc = numpy.longdouble([0,0,0]), 
                                    force = numpy.double([0,0,0]),
                                    index = 2,
                                    color = 'springgreen',
                                    radius = 20))
        List.append(Bodies.Planet( mass = settings.MASS_MAGNIFIER * 1, 
                                    pos = numpy.longdouble([10,0,0]), 
                                    vel = settings.SPEED_MAGNIFIER * numpy.longdouble([0,1,-2]), 
                                    acc = numpy.longdouble([0,0,0]), 
                                    force = numpy.double([0,0,0]),
                                    index = 3,
                                    color = 'steelblue',
                                    radius = 10))
    return List

def AstroidProblem():
    List = []
    List.append(Bodies.Sun(     mass = settings.MASS_MAGNIFIER * 100, 
                                pos = numpy.longdouble([0,0,0]), 
                                vel = settings.SPEED_MAGNIFIER * numpy.longdouble([0,0,0]), 
                                acc = numpy.longdouble([0,0,0]), 
                                force = numpy.double([0,0,0]),
                                index = 1,
                                radius = 20))
    List.append(Bodies.Planet(  mass = settings.MASS_MAGNIFIER * 3, 
                                pos = numpy.longdouble([30,0,0]), 
                                vel = settings.SPEED_MAGNIFIER * numpy.longdouble([0,-10,0]), 
                                acc = numpy.longdouble([0,0,0]), 
                                force = numpy.double([0,0,0]),
                                index = 1,
                                color = 'steelblue',
                                radius = 10))
    List.append(Bodies.Planet(  mass = settings.MASS_MAGNIFIER * 3, 
                                pos = numpy.longdouble([-30,0,0]), 
                                vel = settings.SPEED_MAGNIFIER * numpy.longdouble([0,10,0]), 
                                acc = numpy.longdouble([0,0,0]), 
                                force = numpy.double([0,0,0]),
                                index = 1,
                                color = 'springgreen',
                                radius = 10))
    for i in range(100):
        randomRotation = random.random() * 2 * math.pi
        randomSkew = (random.random() * 0.2 - 0.05) + (math.pi/2)
        randomMag = (random.random() + 0.5) * 40
        randomVel = (random.random() + 0.5) / 1.5

        pos = PolarToCarth(randomRotation,randomSkew,randomMag)
        velDir = RotateVector(pos)
        velDir *= randomVel

        List.append(Bodies.Astroid( mass = settings.MASS_MAGNIFIER * 0.00001, 
                                    pos = numpy.copy(pos), 
                                    vel = settings.SPEED_MAGNIFIER * velDir * (5000 / randomMag) ** 0.5, 
                                    acc = numpy.longdouble([0,0,0]), 
                                    force = numpy.double([0,0,0]),
                                    index = 1))
    return List

def PolarToCarth(theta,rho,mag):
    x = mag * math.cos(theta) * math.sin(rho)
    y = mag * math.sin(theta) * math.sin(rho)
    z = mag * math.cos(rho)
    return numpy.longdouble([x,y,z])

def RotateVector(pos):
    rotVector = numpy.double([0, 0, 1])
    dir = numpy.cross(pos,rotVector)
    return dir / numpy.sqrt(numpy.sum(dir**2))
