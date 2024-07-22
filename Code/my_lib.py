# My library of python functions


def main():
    print("Hello World!")

if __name__ == "__main__":
    main()

	
print("Started reading the script")

def fact(N):
    ret = 1
    for i in range (N):
        ret = ret*(i+1)
    return ret


def sin(x):
    _num_of_terms = 4
    ret = 0
    for i in range(_num_of_terms):
        ret = ret + (-1)**i * (x**(2*i+1)/fact(2*i+1))
    return ret

	
from numba import jit

@jit
def jfact(N):
    ret = 1
    for i in range (N):
        ret = ret*(i+1)
    return ret

@jit
def jsin(x):
    _num_of_terms = 4
    ret = 0
    for i in range(_num_of_terms):
        ret = ret + (-1)**i * (x**(2*i+1)/jfact(2*i+1))
    return ret

    
def save_to_vtu(positions, velocities, rad):
    import numpy as np
    import os
    path = "./vtu/"
    if ( not os.path.exists( path )):
        os.mkdir( path )   
    Num_part = np.shape(positions)[0]
    Num_timesteps = np.shape(positions)[1]
    

    for timestep in range(Num_timesteps):
        pos_string = ""
        vel_string = ""
        rad_string = ""
        for num in range(Num_part):
            point_pos = positions[num,timestep]
            point_vel = velocities[num,timestep]
            pos_string = pos_string + "%.4f " % point_pos[0]
            pos_string = pos_string + "%.4f " % point_pos[1]
            pos_string = pos_string + "%.4f " % point_pos[2]
            vel_string = vel_string + "%.4f " % point_vel[0]
            vel_string = vel_string + "%.4f " % point_vel[1]
            vel_string = vel_string + "%.4f " % point_vel[2]
            rad_string = rad_string + "%.4f " % rad[num]
            
            
        pos_string = pos_string + "\n"
        vel_string = vel_string + "\n"
        rad_string = rad_string + "\n"
        
        textfile = open("./vtu/timestep_"+str(timestep)+".vtu", "w") 
        textfile. write('<?xml version="1.0"?>\n')
        textfile. write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        textfile. write(' <UnstructuredGrid>\n')
        textfile. write('  <Piece NumberOfPoints="'+str(Num_part)+'" NumberOfCells="0">\n')
        textfile. write('   <Cells>\n')
        textfile. write('    <DataArray type="Int32" name="connectivity" format="ascii">\n')
        textfile. write('       0\n')
        textfile. write('    </DataArray>\n')
        textfile. write('    <DataArray type="Float32" name="offset" format="ascii">\n')
        textfile. write('       0\n')
        textfile. write('    </DataArray>\n')
        textfile. write('    <DataArray type="UInt8" name="types" format="ascii">\n')
        textfile. write('       1\n')
        textfile. write('    </DataArray>\n')
        textfile. write('   </Cells>\n')
        textfile. write('   <Points>\n')
        textfile. write('<DataArray type="Float32" NumberOfComponents="3" format="ascii">\n')
        textfile. write(pos_string)
        textfile. write('</DataArray>\n')
        textfile. write('    </Points>\n')
        textfile. write('    <PointData>\n')
        textfile. write('<DataArray type="Float32" Name="Position" NumberOfComponents="3" format="ascii">\n')
        textfile. write(pos_string)
        textfile. write('</DataArray>\n')
        textfile. write('<DataArray type="Float32" Name="Velocity" NumberOfComponents="3" format="ascii">\n')
        textfile. write(vel_string)
        textfile. write('</DataArray>\n')
        textfile. write('<DataArray type="Float32" Name="Scale" NumberOfComponents="1" format="ascii">\n')
        textfile. write(rad_string)
        textfile. write('</DataArray>\n')
        textfile. write('</PointData>\n')
        textfile. write('<CellData/>\n')
        textfile. write('</Piece>\n')
        textfile. write('</UnstructuredGrid>\n')
        textfile. write('</VTKFile>\n')
        textfile. close()
    
    return 

print("Ended reading the script")