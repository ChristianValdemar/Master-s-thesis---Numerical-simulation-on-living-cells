#
# Logistic reaction functions and logistic initial data
#
from dolfin import *
import numpy, ast, math, csv
import time
import os

start = time.clock()
print "Reading mesh, kp, km, c0 and subdomains from files"

# set parameters
t = 0.0;            #Start time
dt = 0.1;          #Step size
T = 80             #Total simulation time
t_frame = 1.6       #time between each frame
t_bleach = 0.4      #Bleaching time
t_to_bleach = 0.1    #t_to_bleach*t_frame = first bleaching time

# define piecewise const diffusion coefficients
alpha_C = Constant(25)          # Diffusioncoefficient for u in cytoplasma
beta_N = Constant(25)           # Diffusioncoefficient for v in nucleus
alpha_M = Constant(0.0005)          # Diffusioncoefficient for u in membrane
beta_M = Constant(0.00025)           # Diffusioncoefficient for v in membrane
alpha_N = Constant(0.001)        # Small diff. coef.
beta_C = Constant(0.001)         # Small diff. coef.

# Define reaction constants
q = 2 #Constant depending on laser intensitet
k = 240.0*q/(1+q);  #Bleaching effect constant
Kp = 1    #Constant for scaling the kp reaction parameter
Km = 1   #Constant for scaling the km reaction parameter
nucleus = []
cyto = []


# create mesh and finite element
directory = "Mesh"
mesh = Mesh(os.path.join(directory, "mesh.xml"))


# Create function spaces
W1 = FunctionSpace(mesh, "CG", 1)
V = MixedFunctionSpace([W1,W1])
W0 = FunctionSpace(mesh, "DG", 0)
V0 = MixedFunctionSpace([W0,W0])



# Define subdomains
subdomains = MeshFunction("size_t", mesh, os.path.join(directory, "subdomains.xml"))



# Define new measures associated with the interior domains
dxx = Measure("dx")[subdomains]
# dx(0) is outer membrane, dx(1) is inner membrane, dx(2) is bleaching, dx(3) is cytoplasma dx(4) is nucleus
#plot(subdomains)
#interactive()


# Read coefficients for the reaction part
kp = Function(W1,"initialization_logistic/kp.xml")
km = Function(W1,"initialization_logistic/km.xml")

# Print kp and km to files, for visualization in Paraview
#filek = File(os.path.join(directory, "reactionconst.pvd"))
#filek << kp
#filek << km

# Show kp and km in FEniCS
#plot(kp)
#time.sleep(2)
#plot(km)
#interactive()





# Create intial conditions and functions
(u, v) = TrialFunctions(V)
(phi, psi) = TestFunctions(V)

c0 = Function(V,"initialization_logistic/c0.xml")
(u0, v0) = as_vector((c0[0], c0[1]))


c1 = Function(V)
(u1, v1) = as_vector((c1[0], c1[1]))

# Plot u0 and v0 in FEniCS
#plot(u0)
#time.sleep(2)
#plot(v0)
#interactive()



print "Mesh and function imported. In time: %g seconds" %(time.clock() - start)


# # dx(0) is outer membrane, dx(1) is inner membrane, dx(2) is bleaching, dx(3) is cytoplasma dx(4) is nucleus       


Fv = (1/dt)*(v - v0)*psi *dx \
    + beta_N*inner(grad(v), grad(psi)) *dxx(4) \
    + beta_C*inner(grad(v), grad(psi)) *(dxx(2)+dxx(3)) \
    + beta_M*inner(grad(v), grad(psi)) *(dxx(0)+dxx(1)) \
    + inner(Km*km*v - Kp*kp*u, psi) *(dxx(0)+dxx(1))

# Fu_b with bleaching
Fu_b = (1/dt)*(u - u0)*phi *dx \
    + alpha_C*inner(grad(u), grad(phi)) *(dxx(2)+dxx(3)) \
    + alpha_M*inner(grad(u), grad(phi)) *(dxx(0)+dxx(1)) \
    + alpha_N*inner(grad(u), grad(phi)) *dxx(4) \
    + inner(Kp*kp*u - Km*km*v, phi) *(dxx(0)+dxx(1)) \
    + k*u*phi *dxx(2)

# Fp without bleaching
Fu = (1/dt)*(u - u0)*phi *dx \
    + alpha_C*inner(grad(u), grad(phi)) *(dxx(2)+dxx(3)) \
    + alpha_M*inner(grad(u), grad(phi)) *(dxx(0)+dxx(1)) \
    + alpha_N*inner(grad(u), grad(phi)) *dxx(4) \
    + inner(Kp*kp*u - Km*km*v, phi) *(dxx(0)+dxx(1)) 

F = Fu + Fv                     #System without bleaching
F_b = Fu_b + Fv                 #System with bleaching

a = lhs(F); L = rhs(F)          #System without bleaching
a_b = lhs(F_b); L_b = rhs(F_b)  #System with bleaching



# preassembly
A = assemble(a)
A_b = assemble(a_b)
b = None



file1 = File(os.path.join(directory, "solutions_logistic/u.pvd"))
file2 = File(os.path.join(directory, "solutions_logistic/v.pvd"))
file3 = File(os.path.join(directory, "solutions_logistic/c.pvd"))
t += dt

while t < T:
    # Print initial data to file
    if t == dt:
        (u0, v0) = c0.split()
        file1 << u0
        file2 << v0

    # Bleach only if we have been simulating for more than t_to_bleach*t_frame and only if it is bleaching time.
    if t%t_frame <= t_bleach and t >= t_to_bleach*t_frame:
        b = assemble(L_b, tensor=b)
        solve(A_b, c1.vector(), b)
    else:
        b = assemble(L, tensor=b)
        solve(A, c1.vector(), b)

    c0.assign(c1)
    (u1, v1) = c1.split()
    
    # Continuous solution c=u+v
    c = TrialFunction(W1)
    rho = TestFunction(W1)
    L1 = v1*rho*(dxx(0)+dxx(1)+dxx(4)) + u1*rho*(dxx(0)+dxx(1)+dxx(2)+dxx(3))
    a1 = c*rho*dx
    c = Function(W1)
    solve(a1 == L1, c)


    print " time = %g, |c_trial| = %g" %(t, c.vector().norm("l2"))
    
    

#    plot(c1m)
#   Write solutions to files, which can be displayed in Paraview
    file1 << u1
    file2 << v1
    file3 << c

#   Find average values
    if "{0:.1f}".format(t%t_frame) == "{0:.1f}".format((t_frame+t_bleach)/2):
        
        g = TrialFunctions(V)
        (gu, gv) = as_vector((g[0], g[1]))
        Fg = gu*phi*dx - u1*phi*(dxx(2)+dxx(3)) + gv*psi*dx - v1*psi*dxx(4)
        ag = lhs(Fg); Lg = rhs(Fg)
        g = Function(V)
        solve(ag == Lg, g)
        (gu, gv) = g.split()
        
        g_idx = TrialFunctions(V)
        (gu_idx, gv_idx) = as_vector((g_idx[0], g_idx[1]))
        Fg = gu_idx*phi*dx - Constant(1)*phi*(dxx(2)+dxx(3)) + gv_idx*psi*dx - Constant(1)*psi*dxx(4)
        ag = lhs(Fg); Lg = rhs(Fg)
        g_idx = Function(V)
        solve(ag == Lg, g_idx)
        (gu_idx, gv_idx) = g_idx.split()
        
        gu = project(gu,W1)
        gv = project(gv,W1)
        gu_idx = project(gu_idx,W1)
        gv_idx = project(gv_idx,W1)
        
        a_u = gu.vector().sum()/gu_idx.vector().sum()
        a_v = gv.vector().sum()/gv_idx.vector().sum()
        cyto.append(a_u)
        nucleus.append(a_v)
        print a_u
        print a_v

    # Move to next time step
    t += dt


# Write average data to file
# open a file for writing.
csv_out = open('simulation-results-between-bleach.csv', 'wb')
mywriter = csv.writer(csv_out)

# writerow - one row of data at a time.
for row in zip(cyto, nucleus):
    mywriter.writerow(row)

csv_out.close()

