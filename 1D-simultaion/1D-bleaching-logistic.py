from dolfin import *
import numpy, ast
import time
import os


# set parameters
t = 0.0;            #Start time
dt = 0.1;           #Step size
T = 5;             #Total simulation time
memTickness = 0.009 #membrane thickness in um

# define piecewise const diffusion coefficients
alpha_C = Constant(25)       # Diffusioncoefficient for u in cytoplasma
beta_N = Constant(25)       # Diffusioncoefficient v in nucleus
alpha_M = Constant(0.0005)          # Diffusioncoefficient for u in membrane
beta_M = Constant(0.00025)           # Diffusioncoefficient for v in membrane
alpha_N = Constant(1.0) #Constant(0.001)       # Small diff. coef.
beta_C = Constant(1.0) #Constant(0.001)        # Small diff. coef.



# define mesh
mesh = IntervalMesh(1000, 0.0 , 0.029)


# Create function spaces
W1 = FunctionSpace(mesh, "CG", 1)
V = MixedFunctionSpace([W1,W1])
W0 = FunctionSpace(mesh, "DG", 0)
V0 = MixedFunctionSpace([W0,W0])



#
# set symmetric initial data
#

# left part
class Left(SubDomain):
    def inside(self, x, on_boundary):
        return ( (x[0] < 0.01) )
left = Left()

# right part
class Right(SubDomain):
    def inside(self, x, on_boundary):
        return ( (x[0] > 0.019) )
right = Right()

# Initialize mesh function for interior domains
subdomains = CellFunction("size_t", mesh)
subdomains.set_all(0)
right.mark(subdomains, 2)
left.mark(subdomains, 1)






# Define new measures associated with the interior domains
dxx = Measure("dx")[subdomains]      
#plot(subdomains)
#interactive()

#Coefficients for heaviside function
#works
#k=8000
#epsilon=0.00000000001
k=10000
epsilon=0.0001


# Class for creation of kp
class ReactionParameterKp(Expression):
    def eval_cell(self, value, x, ufc_cell):
        value[0] = 1/epsilon*(1/(1+exp(-k*(x[0]-0.0145))))


#Class for creation of km
class ReactionParameterKm(Expression):
    def eval_cell(self, value, x, ufc_cell):
        value[0] = 1/epsilon*(1/(1+exp(k*(x[0]-0.0145))))


kp = Function(W1)
temp1 = ReactionParameterKp()
kp.interpolate(temp1)
#plot(kp)
#time.sleep(3)
File("initialization/kp.xml") << kp
File("initialization/kp.pvd") << kp

km = Function(W1)
temp = ReactionParameterKm()
km.interpolate(temp)
plot(km)
#time.sleep(3)
File("initialization/km.xml") << km
File("initialization/km.pvd") << km






#
# Create intial conditions and functions
#

## Class representing the intial conditions
class InitialConditions(Expression):
    def eval_cell(self, value, x, ufc_cell):
        # u
        if subdomains[ufc_cell.index] == 0:
            value[0] = 0.5
        elif subdomains[ufc_cell.index] == 1:
            value[0] = 1.0
        else:
            value[0] = 0.0
        # v
        if subdomains[ufc_cell.index] == 1:
            value[1] = 0.0
        elif subdomains[ufc_cell.index] == 2:
            value[1] = 1.0
        else:
            value[1] = 0.5
    def value_shape(self):
        return (2,)


(u, v) = TrialFunctions(V)
(phi, psi) = TestFunctions(V)

# Create intial conditions and interpolate
c0 = Function(V)
c_init = InitialConditions()
c0.interpolate(c_init)

(u0, v0) = as_vector((c0[0], c0[1]))

c1 = Function(V)
#(u1, v1) = as_vector((c1[0], c1[1]))

#plot(u0)
#time.sleep(2)
#plot(v0)
#interactive()


# Boundary
#cell_boundary = AutoSubDomain(lambda x, on_boundary: on_boundary and not near(x[0], 0.029))
#
#bc = DirichletBC(V.sub(0), Constant(0.0), cell_boundary)



Fu = (1/dt)*(u - u0)*phi *dx \
    + alpha_C*inner(grad(u), grad(phi)) *dxx(1) \
    + alpha_M*inner(grad(u), grad(phi)) *dxx(0) \
    + alpha_N*inner(grad(u), grad(phi)) *dxx(2) \
    + inner(kp*u - km*v, phi) *dxx(0)


Fv = (1/dt)*(v - v0)*psi *dx \
    + beta_N*inner(grad(v), grad(psi)) *dxx(2) \
    + beta_M*inner(grad(v), grad(psi)) *dxx(0) \
    + beta_C*inner(grad(v), grad(psi)) *dxx(1) \
    + inner(km*v - kp*u, psi) *dxx(0)


F = Fu + Fv

a = lhs(F); L = rhs(F)



# preassembly
A = assemble(a)
b = None



file1 = File("solutions/u.pvd")
file2 = File("solutions/v.pvd")
file3 = File("solutions/c.pvd")
t += dt

while t < T:
    # Print initial data to file
    if t == dt:
        (u0, v0) = c0.split()
        file1 << u0
        file2 << v0
        
    
    b = assemble(L, tensor=b)
#    bc.apply(A, b)
    solve(A, c1.vector(), b)

    c0.assign(c1)
    (u1, v1) = c1.split()
    
    # Continuous solution
    c = TrialFunction(W1)
    rho = TestFunction(W1)
    L1 = u1*rho*dx + v1*rho*dx
    a1 = c*rho*dx
    c = Function(W1)
    solve(a1 == L1, c)

#    # Discontinuous solution
#    c = TrialFunction(W0)
#    v = TestFunction(W0)
#    L1 = c1m*v*dx + c1p*v*dx
#    a1 = c*v*dx
#    c = Function(W0)
#    solve(a1 == L1, c)

#    plot(u1 , title = "u"); time.sleep(1)
#    plot(v1, title = "v"); time.sleep(1)
#    plot(c, title = "c"); time.sleep(1)
    print " time = %g, |c_trial| = %g" %(t, c.vector().norm("l2"))

    # Move to next time step
    
    t += dt


#    plot(c1m)
    file1 << u1
    file2 << v1
    file3 << c

#plot(u1)
#interactive()


