#!/usr/bin/env python 

# Some suitable functions and data structures for drawing a map and particles

import time
import random
import math
from version3 import *

# Functions to generate some dummy particles data:
def calcX():
    return random.gauss(80,3) + 70*(math.sin(t)); # in cm

def calcY():
    return random.gauss(70,3) + 60*(math.sin(2*t)); # in cm

def calcW():
    return random.random();

def calcTheta():
    return random.randint(0,360);

# A Canvas class for drawing a map and particles:
#     - it takes care of a proper scaling and coordinate transformation between
#      the map frame of reference (in cm) and the display (in pixels)
class Canvas:
    def __init__(self,map_size=210):
        self.map_size    = map_size;    # in cm;
        self.canvas_size = 768;         # in pixels;
        self.margin      = 0.05*map_size;
        self.scale       = self.canvas_size/(map_size+2*self.margin);
        
    # line is a tuple that contains the x,y for the start point and end point
    def drawLine(self,line):
        x1 = self.__screenX(line[0]);
        y1 = self.__screenY(line[1]);
        x2 = self.__screenX(line[2]);
        y2 = self.__screenY(line[3]);
        print "drawLine:" + str((x1,y1,x2,y2))

    # data is the cloud of particles
    def drawParticles(self,data):
        display = [(self.__screenX(d[0]),self.__screenY(d[1])) + d[2:] for d in data];
        print "drawParticles:" + str(display);
    
    # given the x,y display it on the screen
    def __screenX(self,x):
        return (x + self.margin)*self.scale

    def __screenY(self,y):
        return (self.map_size + self.margin - y)*self.scale

# A Map class containing walls
class Map:
    def __init__(self):
        self.walls = [];

    def add_wall(self,wall):
        self.walls.append(wall);

    def clear(self):
        self.walls = [];

    def draw(self):
        for wall in self.walls:
            canvas.drawLine(wall);

# Simple Particles set
class Particles:
    def __init__(self):
        self.n = 10;     # how many particles we use
        self.data = [];

    def update(self):
        self.data = [(calcX(), calcY(), calcTheta(), calcW()) for i in range(self.n)];
    
    def draw(self):
        canvas.drawParticles(self.data);

canvas = Canvas();    # global canvas we are going to draw on

mymap = Map();
# Definitions of walls
# a: O to A
# b: A to B
# c: C to D
# d: D to E
# e: E to F
# f: F to G
# g: G to H
# h: H to O
mymap.add_wall((0,0,0,168));        # a
mymap.add_wall((0,168,84,168));     # b
mymap.add_wall((84,126,84,210));    # c
mymap.add_wall((84,210,168,210));   # d
mymap.add_wall((168,210,168,84));   # e
mymap.add_wall((168,84,210,84));    # f
mymap.add_wall((210,84,210,0));     # g
mymap.add_wall((210,0,0,0));        # h
mymap.draw();

particles = Particles();

'''
t = 0;
while True:
    particles.update();
    particles.draw();
    t += 0.05;
    time.sleep(0.05);
'''
#webSimulation()
numOfWall = 8
def calculate_likelihood(x,y,theta,z):
    m = 9999 # predicted distance to the wall
    for wall_index in range(numOfWall):
        Ax = mymap.walls[wall_index][0]
        Ay = mymap.walls[wall_index][1]
        Bx = mymap.walls[wall_index][2]
        By = mymap.walls[wall_index][3]
        
        dist = ((By-Ay)*(Ax-x)-(Bx-Ax)*(Ay-y)) / ((By-Ay)*math.cos(theta)-(Bx-Ax)*math.sin(theta)+0.00000000001)
        
        if dist >= 0:
            x_hypo = x + dist * math.cos(theta)
            y_hypo = y + dist * math.sin(theta)
            if (x_hypo-Ax)*(x_hypo-Bx) + (y_hypo-Ay)*(y_hypo-By) <= 0 and dist<m:
                m = dist
        sigma_s = 1
    return math.exp(-(z-m)**2/(2*sigma_s**2)) #un normalised
        
    
print calculate_likelihood(10,10,0.1,200)

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    