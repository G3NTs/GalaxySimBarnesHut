import numpy
import matplotlib.pyplot as plt
import math

import settings

#==========================================================================#

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
    
#==========================================================================#

class Octtree:
    def __init__(self,box,max_points=settings.CAPACITY,depth=0):
        self.box        = box
        self.depth      = depth
        self.divided    = False
        self.max_points = max_points

        self.subTree    = []

        self.points     = []
        self.childMass  = 0
        self.massPos    = box.centre

        self.force      = numpy.double([0,0,0])
        self.acc        = numpy.double([0,0,0])
        self.vel        = numpy.double([0,0,0])

    def insertPoint(self, point):
        if (not self.box.contains(point)):
            return False
        
        self.massPos = ((point.pos * point.mass) + (self.childMass * self.massPos))/(self.childMass+point.mass)
        self.childMass += point.mass

        if (len(self.points) < self.max_points and self.divided == False) or self.depth == settings.MAX_DEPTH: # check if room in box
            self.points.append(point)
            point.box = self.box
            return True # end loop
        elif (self.divided == False): # if no room, subdivide
            self.Subdivide()
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
    
    def GetChildren(self):
        if self.divided == True:
            return [self.subTree[0], self.subTree[1], self.subTree[2], self.subTree[3], self.subTree[4], self.subTree[5], self.subTree[6], self.subTree[7]]
        else:
            return self.points
    
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

    def CalculateBodyForces4(self,point):
        distance = self.massPos - point.pos
        magnitude = numpy.sqrt(distance.dot(distance))
        if (((self.box.size[0] + self.box.size[1] + self.box.size[2])/3)/magnitude) >= settings.THRESHOLD:
            for child in self.subTree:
                child.CalculateBodyForces4(point)
            if not self.subTree:
                for treePoint in self.points:
                    if point != treePoint:
                        distance = treePoint.pos - point.pos
                        magnitude = numpy.sqrt(distance.dot(distance))
                        magnitude = max(min(magnitude,50),5)
                        bodyForces = (settings.GRAV_CONST*point.mass*treePoint.mass*distance)/(numpy.power(magnitude,3))
                        point.force += bodyForces
        else:
            bodyForces = (settings.GRAV_CONST*point.mass*self.childMass*distance)/(numpy.power(magnitude,3))
            point.force += bodyForces

    def UpdateOctTree(self):
        for point in self.points:
            if not point.box.contains(point):
                del self
                return
    
    def DrawBoxes(self,ax):
        if self.divided == True:
            self.subTree[0].DrawBoxes(ax)
            self.subTree[1].DrawBoxes(ax)
            self.subTree[2].DrawBoxes(ax)
            self.subTree[3].DrawBoxes(ax)
            self.subTree[4].DrawBoxes(ax)
            self.subTree[5].DrawBoxes(ax)
            self.subTree[6].DrawBoxes(ax)
            self.subTree[7].DrawBoxes(ax)

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
    
    def DrawBoxMass(self,ax):
        if self.divided == True:
            self.subTree[0].DrawBoxMass(ax)
            self.subTree[1].DrawBoxMass(ax)
            self.subTree[2].DrawBoxMass(ax)
            self.subTree[3].DrawBoxMass(ax)
            self.subTree[4].DrawBoxMass(ax)
            self.subTree[5].DrawBoxMass(ax)
            self.subTree[6].DrawBoxMass(ax)
            self.subTree[7].DrawBoxMass(ax)
        ax.plot(*self.massPos, marker = 'o', color = 'b')