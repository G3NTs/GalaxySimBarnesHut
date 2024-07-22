import numpy
import matplotlib.pyplot as plt
import math

if 'GRAV_CONST' and 'DELTA_TIME' and 'SIZE' in globals():
    pass
else:
    GRAV_CONST = 6.674e-11
    DELTA_TIME = 0.1
    SIZE = 10


# Box Class, which is used as a boundary to check if a given point falls within the area.

class Box:
    def __init__(self,pos = numpy.double([0,0,0]), size = numpy.double([0,0,0])):
        self.pos = pos  
        self.size = size
        self.centre = pos + size/2

    # checks if a point is within the box
    def contains(self, point):
        if (all(point.pos >= self.pos) and all(point.pos < (self.pos + self.size))):
            return True
        return False
    
    # checks if an other range(box) intersects with the box.
    def Intersects(self,checkRange):
        if(all(checkRange.pos < (self.pos + self.size)) and all((checkRange.pos + checkRange.size) >= self.pos )):
            return True
        return False

class Octtree:
    def __init__(self,box,max_points=1,depth=0):
        self.box = box
        self.max_points = max_points
        self.points = []
        self.depth = depth
        self.divided = False
        self.childMass = 0
        self.massPos = box.centre
        self.subTree = []
        self.force = numpy.double([0,0,0])
        self.acc = numpy.double([0,0,0])
        self.vel = numpy.double([0,0,0])
        self.deltaPos = 0

    def insertPoint(self, point):
        if (not self.box.contains(point)):
            #print(point.pos," does NOT fall within box: ",self.box.pos,",",self.box.size)
            return False
        #print(point.pos," does fall within box: ",self.box.pos,",",self.box.size)
        
        # edit by adding child mass and velocity to the box parameters
        self.vel = (point.vel*point.mass + self.massPos*self.childMass)/(self.childMass+point.mass)
        self.massPos = ((point.pos * point.mass) + (self.childMass * self.massPos))/(self.childMass+point.mass)
        self.childMass += point.mass

        if (len(self.points) < self.max_points and self.divided == False): # check if room in box
            self.points.append(point)
            #print("assigned point")
            return True # end loop
        elif (self.divided == False): # if no room, subdivide
            self.Subdivide()
            #print("subdividing")
            for oldPoint in self.points: # then remove old points and reassign them to the new subdivided boxes
                if (self.subTree[0].insertPoint(oldPoint) or
                    self.subTree[1].insertPoint(oldPoint) or 
                    self.subTree[2].insertPoint(oldPoint) or 
                    self.subTree[3].insertPoint(oldPoint) or 
                    self.subTree[4].insertPoint(oldPoint) or 
                    self.subTree[5].insertPoint(oldPoint) or 
                    self.subTree[6].insertPoint(oldPoint) or 
                    self.subTree[7].insertPoint(oldPoint)):
                    pass
                else:
                    print("ERROR: old point wasn't reassigned!")
        return (self.subTree[0].insertPoint(point) or
                self.subTree[1].insertPoint(point) or 
                self.subTree[2].insertPoint(point) or 
                self.subTree[3].insertPoint(point) or 
                self.subTree[4].insertPoint(point) or 
                self.subTree[5].insertPoint(point) or 
                self.subTree[6].insertPoint(point) or 
                self.subTree[7].insertPoint(point)) # then add the newest point
    
    def Subdivide(self):
        #======Index=representation=====#
        #           3                   #
        #   1       |                   #
        #   |       |   7      z        #
        #   |   5   2   |      |        #
        #   0   |       |      O --- y  #
        #       |       6       \       #
        #       4                x      #
        #===============================#

        i = 0
        for x in range(2):
            for y in range(2):
                for z in range(2):
                    self.subTree.append(Octtree(Box( numpy.double([ self.box.pos[0] + x*self.box.size[0]/2,
                                                                    self.box.pos[1] + y*self.box.size[1]/2,
                                                                    self.box.pos[2] + z*self.box.size[2]/2]),
                                                    numpy.double([  self.box.size[0]/2,
                                                                    self.box.size[1]/2,
                                                                    self.box.size[2]/2])),
                                                    self.max_points,
                                                    self.depth + 1))
                    #print(self.box.pos[0] + x*self.box.size[0]/2)
                    #print(self.box.pos[1] + y*self.box.size[1]/2)
                    #print(self.box.pos[2] + z*self.box.size[2]/2)
                    #print("=====")
                    #print(self.subTree[i].box.pos)
                    #print("==========================")
                    i+=1
        self.divided = True
    
    def Find(self,checkRange,found_points = []):
        if (self.box.Intersects(checkRange) == False):
            return False
        if (self.divided == True):
            self.subTree[0].Find(checkRange,found_points)
            self.subTree[1].Find(checkRange,found_points)
            self.subTree[2].Find(checkRange,found_points)
            self.subTree[3].Find(checkRange,found_points)
            self.subTree[4].Find(checkRange,found_points)
            self.subTree[5].Find(checkRange,found_points)
            self.subTree[6].Find(checkRange,found_points)
            self.subTree[7].Find(checkRange,found_points)
        else:
            for point in self.points:
                if (checkRange.contains(point)):
                    found_points.append(point)
        return found_points
    
    def CalculateBodyForces(self):
        if self.subTree:
            for index, tree1 in enumerate(self.subTree):
                tree1.CalculateBodyForces()
                for tree2 in self.subTree[index+1:]:
                    #print("forces are calculated for tree point")
                    distance = tree1.massPos - tree2.massPos
                    magnitude = numpy.sqrt(distance.dot(distance))
                    bodyForces = (GRAV_CONST*tree1.childMass*tree2.childMass*distance)/(numpy.power(magnitude,3))
                    tree1.force += bodyForces
                    tree2.force -= bodyForces
        elif self.points:
            #print("point exists")
            #print(self.points)
            for index, point1 in enumerate(self.points):
                for point2 in self.points[index+1:]:
                    #print("forces are calculated for point")
                    distance = point1.pos - point2.pos
                    magnitude = numpy.sqrt(distance.dot(distance))
                    bodyForces = (GRAV_CONST*point1.mass*point2.mass*distance)/(numpy.power(magnitude,3))
                    point1.force -= bodyForces
                    point2.force += bodyForces
        #print("point doesn't exist")
        return
    

    def CalculateBodyForces2(self):
        if self.subTree:
            for index1, tree1 in enumerate(self.subTree):
                for tree2 in self.subTree[index1+1:]:
                    distance = tree2.massPos - tree1.massPos
                    magnitude = numpy.sqrt(distance.dot(distance))
                    bodyForces = (GRAV_CONST*tree1.childMass*tree2.childMass*distance)/(numpy.power(magnitude,3))
                    tree1.force += bodyForces
                    tree2.force -= bodyForces
            for tree in self.subTree:
                tree.force += self.force
                tree.CalculateBodyForces2()
        elif self.points:
            for index2, point1 in enumerate(self.points):
                for point2 in self.points[index2+1:]:
                    distance = point2.pos - point1.pos
                    magnitude = numpy.sqrt(distance.dot(distance))
                    bodyForces = (GRAV_CONST*point1.mass*point2.mass*distance)/(numpy.power(magnitude,3))
                    point1.force += bodyForces
                    point2.force -= bodyForces
            for point in self.points:
                point.force += self.force
        return

    def CalculateBodyForces3(self):
        if not self.subTree:
            for index1, point1 in enumerate(self.points):
                point1.force += self.force
                for point2 in self.points[index1+1:]:
                    distance = point2.pos - point1.pos
                    magnitude = numpy.sqrt(distance.dot(distance))
                    bodyForces = (GRAV_CONST*point1.mass*point2.mass*distance)/(numpy.power(magnitude,3))
                    point1.force += bodyForces
                    point2.force -= bodyForces
        else:
            for index2, tree1 in enumerate(self.subTree):
                if tree1.points:
                    tree1.force += self.force
                    print("why is this running?")
                    for tree2 in self.subTree[index2+1:]:
                        if tree2.points:
                            distance = tree2.massPos - tree1.massPos
                            magnitude = numpy.sqrt(distance.dot(distance))
                            bodyForces = (GRAV_CONST*tree1.childMass*tree2.childMass*distance)/(numpy.power(magnitude,3))
                            tree1.force += bodyForces
                            tree2.force -= bodyForces
                    tree1.CalculateBodyForces3()
    
    def UpdatePositions(self):
        for point in self.points:
            point.acc = point.force/point.mass
            point.vel += point.acc * DELTA_TIME
            point.pos += point.vel * DELTA_TIME
        return
    
    def DrawBoxes(self,ax,fig):
        if self.divided == True:
            self.subTree[0].DrawBoxes(ax,fig)
            self.subTree[1].DrawBoxes(ax,fig)
            self.subTree[2].DrawBoxes(ax,fig)
            self.subTree[3].DrawBoxes(ax,fig)
            self.subTree[4].DrawBoxes(ax,fig)
            self.subTree[5].DrawBoxes(ax,fig)
            self.subTree[6].DrawBoxes(ax,fig)
            self.subTree[7].DrawBoxes(ax,fig)

        x1 = self.box.pos
        x2 = self.box.pos + [self.box.size[0],0,0]
        x3 = self.box.pos + [self.box.size[0],self.box.size[1],0]
        x4 = self.box.pos + [0,self.box.size[1],0]

        x5 = self.box.pos + [0,0,self.box.size[2]]
        x6 = self.box.pos + [self.box.size[0],0,self.box.size[2]]
        x7 = self.box.pos + [self.box.size[0],self.box.size[1],self.box.size[2]]
        x8 = self.box.pos + [0,self.box.size[1],self.box.size[2]]

        ax.plot((x1[0],x2[0]),(x1[1],x2[1]),(x1[2],x2[2]), marker = 'o', color = 'r')
        ax.plot((x2[0],x3[0]),(x2[1],x3[1]),(x2[2],x3[2]), marker = 'o', color = 'r')
        ax.plot((x3[0],x4[0]),(x3[1],x4[1]),(x3[2],x4[2]), marker = 'o', color = 'r')
        ax.plot((x4[0],x1[0]),(x4[1],x1[1]),(x4[2],x1[2]), marker = 'o', color = 'r')

        ax.plot((x1[0],x5[0]),(x1[1],x5[1]),(x1[2],x5[2]), marker = 'o', color = 'r')
        ax.plot((x2[0],x6[0]),(x2[1],x6[1]),(x2[2],x6[2]), marker = 'o', color = 'r')
        ax.plot((x3[0],x7[0]),(x3[1],x7[1]),(x3[2],x7[2]), marker = 'o', color = 'r')
        ax.plot((x4[0],x8[0]),(x4[1],x8[1]),(x4[2],x8[2]), marker = 'o', color = 'r')

        ax.plot((x6[0],x5[0]),(x6[1],x5[1]),(x6[2],x5[2]), marker = 'o', color = 'r')
        ax.plot((x7[0],x6[0]),(x7[1],x6[1]),(x7[2],x6[2]), marker = 'o', color = 'r')
        ax.plot((x8[0],x7[0]),(x8[1],x7[1]),(x8[2],x7[2]), marker = 'o', color = 'r')
        ax.plot((x5[0],x8[0]),(x5[1],x8[1]),(x5[2],x8[2]), marker = 'o', color = 'r')