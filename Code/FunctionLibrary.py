import math
import numpy
import random
import matplotlib.pyplot as plt

import Bodies
import VTU_List
import OctTree_V2

import settings


def ResetPlot(ax):
    ax.clear()
    ax.set_xlim((-settings.SIZE / 2, settings.SIZE / 2))
    ax.set_ylim((-settings.SIZE / 2, settings.SIZE / 2))
    ax.set_zlim((-settings.SIZE / 2, settings.SIZE / 2))

def Update(List,ax,positions,velocities,rad,timestep):
    ResetPlot(ax)
    OriginTree = OctTree_V2.Octtree(OctTree_V2.Box(pos = numpy.double([-settings.SIZE/2,-settings.SIZE/2,-settings.SIZE/2]), size = numpy.double([settings.SIZE,settings.SIZE,settings.SIZE])))
    for point in List:
        ax.plot(*point.pos,marker="o",markersize=point.radius,color = point.color)
        OriginTree.insertPoint(point)
    for point in List:
        OriginTree.CalculateBodyForces4(point)
    for point in List:
        point.Update_Position(settings.INTEGRATION_CHOICE)
    if settings.WRITE_TO_VTU:
        i = 0
        for point in List:
            positions[i,timestep] = numpy.copy(point.pos)
            velocities[i,timestep] = numpy.copy(point.vel)
            rad[i] = numpy.copy(point.radius)
            i += 1
    if settings.DRAW_OCTTREE_BOXES:
        OriginTree.DrawBoxes(ax)
    if settings.DRAW_BOX_MASS:
        OriginTree.DrawBoxMass(ax)
    plt.pause(settings.TIME_BETWEEN_FRAMES)
    del OriginTree
    return