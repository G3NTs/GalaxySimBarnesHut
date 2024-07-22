import numpy
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import random

from functools import partial
from IPython.display import HTML

from Bodies import body_v4

def Start(List):
    for i in range(10):
        List.append(body_v4(    mass = 10, 
                                radius = 0.1, 
                                pos = numpy.double([random.random()-0.5,random.random()-0.5,random.random()-0.5]), 
                                vel = numpy.double([0,0,0]), 
                                acc = numpy.double([0,0,0]), 
                                force = numpy.double([0,0,0]), 
                                index = i))
    return List

