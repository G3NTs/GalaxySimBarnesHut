import matplotlib.pyplot as plt
import numpy

import SimulationPresets
import FunctionLibrary
import my_lib

import settings

#=================INITIALIZE=LIST=================#

fig = plt.figure(figsize=(10, 10))
ax = plt.subplot(1,1,1, projection="3d")

fig.tight_layout()
ax.view_init(45, -45)

Astroid_List = SimulationPresets.LoadPreset(settings.PRESET_CHOICE)

if settings.WRITE_TO_VTU:
    positions = numpy.empty([len(Astroid_List),settings.TIME_STEPS],dtype=object)
    velocities = numpy.empty([len(Astroid_List),settings.TIME_STEPS],dtype=object)
    rad = numpy.empty([len(Astroid_List)],dtype=float)
else:
    positions = []
    velocities = []
    rad = []

for i in range(settings.TIME_STEPS):
    FunctionLibrary.Update(Astroid_List,ax,positions,velocities,rad,i)



plt.pause(10)


my_lib.save_to_vtu(positions, velocities, rad)
