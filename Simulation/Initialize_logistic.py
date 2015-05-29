from dolfin import *
import numpy, ast
import time
import os

start = time.clock()

# Set parameters
memTickness = 0.009 #membrane thickness in um


# Create mesh and finite element
directory = "Mesh"
mesh = Mesh(os.path.join(directory, "mesh.xml"))


# Create function spaces
W1 = FunctionSpace(mesh, "CG", 1)
V = MixedFunctionSpace([W1,W1])
W0 = FunctionSpace(mesh, "DG", 0)
V0 = MixedFunctionSpace([W0,W0])



# Define subdomains
subdomains = MeshFunction("size_t", mesh, os.path.join(directory, "subdomains.xml"))

# Read Chan-Vese polygons from file
init_innerMem = numpy.loadtxt(os.path.join(directory, "init_membrane_in.txt"), delimiter = ',')
innerMem = numpy.loadtxt(os.path.join(directory, "membrane_in.txt"), delimiter = ',')

init_outerMem = numpy.loadtxt(os.path.join(directory, "init_membrane_out.txt"), delimiter = ',')
outerMem = numpy.loadtxt(os.path.join(directory, "membrane_out.txt"), delimiter = ',')

Mem_center = numpy.loadtxt(os.path.join(directory, "nucleus.txt"), delimiter = ',')


# Determine if a point is inside a given polygon or not
# Algorithm name: "Ray Casting Method".
def point_in_poly(x,y,polyx,polyy):
    n = len(polyx)
    inside = False
    x = int(round(x*10**4,0))
    y = int(round(y*10**4,0))
    p1x = int(round(polyx[0]*10**4,0))
    p1y = int(round(polyy[0]*10**4,0))
    for i in range(1,n+1):
        p2x = int(round(polyx[i % n]*10**4,0))
        p2y = int(round(polyy[i % n]*10**4,0))
        if y > min(p1y,p2y) and y <= max(p1y,p2y) and x <= max(p1x,p2x):
            if p1y != p2y:
                xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
            if p1x == p2x or x <= xints:
                inside = not inside
        p1x,p1y = p2x,p2y
    return inside


# Constants for reaction parameter "1/eps * 1/(1+exp(+-k*|x-xm|))"
eps = 0.0001
k = 10000


# Claas for creation of km
class ReactionParameterKm(Expression):
    def eval_cell(self, value, x, ufc_cell):
        if subdomains[ufc_cell.index] == 0 or subdomains[ufc_cell.index] == 1:
            dist = []
            # euclidean distances from the other points
            dist = ((Mem_center-numpy.array([x[0],x[1]]))**2).sum(axis=1)
            # find the nearest neighbours indexed in a sorted list
            idx = numpy.argsort(dist)
            inReactangle = False
            i = 1
            while inReactangle == False:
                Pidx = idx[0]
                Qidx = idx[i]
                xs = [init_innerMem[Pidx][0],init_outerMem[Pidx][0],init_outerMem[Qidx][0],init_innerMem[Qidx][0]]
                ys = [init_innerMem[Pidx][1],init_outerMem[Pidx][1],init_outerMem[Qidx][1],init_innerMem[Qidx][1]]
                inReactangle = point_in_poly(x[0],x[1],xs,ys)
                i = i+1
            # Vector direction between the two points
            alpha = Mem_center[Pidx][0] - Mem_center[Qidx][0] + DOLFIN_EPS
            beta = Mem_center[Pidx][1] - Mem_center[Qidx][1] + DOLFIN_EPS
            # Distance from point to membrane center
            a= 1/alpha
            b= - 1/beta
            c = -Mem_center[Pidx][0]/alpha + Mem_center[Pidx][1]/beta
            dist = abs(a*x[0] + b*x[1] + c)/sqrt(a**2 + b**2)
            # Set reaction parameter to be a function of the distance from membrane egde
            if subdomains[ufc_cell.index] == 0:
                value[0] = 1/eps * 1/(1+exp(k*(-dist))) #Outer-membrane
            else:
                value[0] = 1/eps * 1/(1+exp(k*dist))    #Inner-membrane
        elif subdomains[ufc_cell.index] == 2 or subdomains[ufc_cell.index] == 3:
            value[0] = 1/eps
        else:
            value[0] = 0.0

#Class for creation of kp
class ReactionParameterKp(Expression):
    def eval_cell(self, value, x, ufc_cell):
        if subdomains[ufc_cell.index] == 0 or subdomains[ufc_cell.index] == 1:
            dist = []
            # euclidean distances from the other points
            dist = ((outerMem-numpy.array([x[0],x[1]]))**2).sum(axis=1)
            # find the nearest neighbours indexed in a sorted list
            idx = numpy.argsort(dist)
            inReactangle = False
            i = 1
            while inReactangle == False:
                Pidx = idx[0]
                Qidx = idx[i]
                xs = [init_outerMem[Pidx][0],init_innerMem[Pidx][0],init_innerMem[Qidx][0],init_outerMem[Qidx][0]]
                ys = [init_outerMem[Pidx][1],init_innerMem[Pidx][1],init_innerMem[Qidx][1],init_outerMem[Qidx][1]]
                inReactangle = point_in_poly(x[0],x[1],xs,ys)
                i = i+1
            # Vector direction between the two points
            alpha = Mem_center[Pidx][0] - Mem_center[Qidx][0]+DOLFIN_EPS
            beta = Mem_center[Pidx][1] - Mem_center[Qidx][1]+DOLFIN_EPS
            # Distance from point to membrane egde
            a= 1/alpha
            b= - 1/beta
            c = -Mem_center[Pidx][0]/alpha + Mem_center[Pidx][1]/beta
            dist = abs(a*x[0] + b*x[1] + c)/sqrt(a**2 + b**2)
            # Set reaction parameter to be a function of the distance from membrane egde
            if subdomains[ufc_cell.index] == 0:
                value[0] = 1/eps * 1/(1+exp(-k*(-dist))) #Outer-membrane
            else:
                    value[0] = 1/eps * 1/(1+exp(-k*dist))    #Inner-membrane
        elif subdomains[ufc_cell.index] == 2 or subdomains[ufc_cell.index] == 3:
            value[0] = 0.0
        else:
            value[0] = 1/eps


# Set km
km = Function(W1)
temp = ReactionParameterKm()
km.interpolate(temp)
#plot(km)
#interactive()
File("initialization_logistic/km.xml") << km
#File("initialization_logistic/km.pvd") << km   # For visualization in Paraview
print "km done"

# Set kp
kp = Function(W1)
temp1 = ReactionParameterKp()
kp.interpolate(temp1)
#plot(km)
#interactive()
File("initialization_logistic/kp.xml") << kp
#File("initialization_logistic/kp.pvd") << kp   # For visualization in Paraview
print "kp done"



# Constants for initial conditions "1/eps * 1/(1+exp(+-k*|x-xm|))"
umax = 53.9   #Scale for cytoplasm
vmax = 119.4  #Scale for nucleus
k = 10000

# Claas for creation of initial conditions
class InitialConditions(Expression):
    def eval_cell(self, value, x, ufc_cell):
        #u0
        if subdomains[ufc_cell.index] == 0 or subdomains[ufc_cell.index] == 1:
            dist = []
            # euclidean distances from the other points
            dist = ((Mem_center-numpy.array([x[0],x[1]]))**2).sum(axis=1)
            # find the nearest neighbours indexed in a sorted list
            idx = numpy.argsort(dist)
            inReactangle = False
            i = 1
            while inReactangle == False:
                Pidx = idx[0]
                Qidx = idx[i]
                xs = [init_innerMem[Pidx][0],init_outerMem[Pidx][0],init_outerMem[Qidx][0],init_innerMem[Qidx][0]]
                ys = [init_innerMem[Pidx][1],init_outerMem[Pidx][1],init_outerMem[Qidx][1],init_innerMem[Qidx][1]]
                inReactangle = point_in_poly(x[0],x[1],xs,ys)
                i = i+1
            # Vector direction between the two points
            alpha = Mem_center[Pidx][0] - Mem_center[Qidx][0] + DOLFIN_EPS
            beta = Mem_center[Pidx][1] - Mem_center[Qidx][1] + DOLFIN_EPS
            # Distance from point to membrane center
            a= 1/alpha
            b= - 1/beta
            c = -Mem_center[Pidx][0]/alpha + Mem_center[Pidx][1]/beta
            dist = abs(a*x[0] + b*x[1] + c)/sqrt(a**2 + b**2)
            # Set reaction parameter to be a function of the distance from membrane egde
            if subdomains[ufc_cell.index] == 0:
                value[0] = umax * 1/(1+exp(k*(-dist))) #Outer-membrane
            else:
                value[0] = umax * 1/(1+exp(k*dist))    #Inner-membrane
        elif subdomains[ufc_cell.index] == 2 or subdomains[ufc_cell.index] == 3:
            value[0] = umax
        else:
            value[0] = 0.0
        #v0
        if subdomains[ufc_cell.index] == 0 or subdomains[ufc_cell.index] == 1:
            dist = []
            # euclidean distances from the other points
            dist = ((outerMem-numpy.array([x[0],x[1]]))**2).sum(axis=1)
            # find the nearest neighbours indexed in a sorted list
            idx = numpy.argsort(dist)
            inReactangle = False
            i = 1
            while inReactangle == False:
                Pidx = idx[0]
                Qidx = idx[i]
                xs = [init_outerMem[Pidx][0],init_innerMem[Pidx][0],init_innerMem[Qidx][0],init_outerMem[Qidx][0]]
                ys = [init_outerMem[Pidx][1],init_innerMem[Pidx][1],init_innerMem[Qidx][1],init_outerMem[Qidx][1]]
                inReactangle = point_in_poly(x[0],x[1],xs,ys)
                i = i+1
            # Vector direction between the two points
            alpha = Mem_center[Pidx][0] - Mem_center[Qidx][0]+DOLFIN_EPS
            beta = Mem_center[Pidx][1] - Mem_center[Qidx][1]+DOLFIN_EPS
            # Distance from point to membrane egde
            a= 1/alpha
            b= - 1/beta
            c = -Mem_center[Pidx][0]/alpha + Mem_center[Pidx][1]/beta
            dist = abs(a*x[0] + b*x[1] + c)/sqrt(a**2 + b**2)
            # Set reaction parameter to be a function of the distance from membrane egde
            if subdomains[ufc_cell.index] == 0:
                value[1] = vmax * 1/(1+exp(-k*(-dist))) #Outer-membrane
            else:
                value[1] = vmax * 1/(1+exp(-k*dist))    #Inner-membrane
        elif subdomains[ufc_cell.index] == 2 or subdomains[ufc_cell.index] == 3:
            value[1] = 0.0
        else:
            value[1] = vmax
    def value_shape(self):
        return (2,)



# Create intial conditions and interpolate
c0 = Function(V)
c_init = InitialConditions()
c0.interpolate(c_init)


#(u0, v0) = as_vector((c0[0], c0[1]))
#plot(c0)
#time.sleep(2)
#plot(v0)
#interactive()

File("initialization_logistic/c0.xml") << c0


print "Initialization done in %g seconds" %(time.clock() - start)

