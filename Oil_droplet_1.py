#---------------------------------------------------------------------- 
# Import Libraries
#----------------------------------------------------------------------
from fenics import *
from mshr   import *
from ufl import *
from dolfin import *
from fenics import *
import random as random
import numpy as np
import time
import csv
#import matplotlib.pyplot as py
from mshr import *

# Set compiler optimizations
parameters["form_compiler"]["cpp_optimize"] = True
parameters["allow_extrapolation"]           = True


#----------------------------------------------------------------------
# Procedures and general Fenics options
#----------------------------------------------------------------------

# Create intial conditions and interpolate phase-field
class SphericalInitialConditions(UserExpression):
	def __init__(self, xxc, yyc, r, we, **kwargs):
		random.seed(5 + MPI.rank(MPI.comm_world))
		self.r   = r
		self.we = we
		self.xc  = xxc
		self.yc  = yyc
		super().__init__(**kwargs)
	def eval(self, values, x):
		if sqrt((x[0]-self.xc)*(x[0]-self.xc)+(x[1]-self.yc)*(x[1]-self.yc)) < self.r :
			values[0] = 0.
			values[1] = self.we
			values[2] = 0.
			values[3] = 0.
			values[4] = 0.
		else:
			values[0] = 0.
			values[1] = 0.
			values[2] = 0.
			values[3] = 0.
			values[4] = 0.
	def value_shape(self):
		return(5,)

class LenghtInitialConditions(UserExpression):
	def __init__(self, xxc, yyc, r, we, **kwargs):
		random.seed(5 + MPI.rank(MPI.comm_world))
		self.r   = r
		self.we = we
		self.xc  = xxc
		self.yc  = yyc
		super().__init__(**kwargs)
	def eval(self, values, x):
		if (sqrt((x[0]-self.xc)*(x[0]-self.xc)+(x[1]-self.yc)*(x[1]-self.yc)) < self.r) and (x[1] < 0.88e-3) :
			values[0] = 0.
			values[1] = self.we
			values[2] = 0.
			values[3] = 0.
			values[4] = 0.
		else:
			values[0] = 0.
			values[1] = 0.
			values[2] = 0.
			values[3] = 0.
			values[4] = 0.
	def value_shape(self):
		return(5,)

def Max_height(ch_, mesh):
    chsub1 = ch_.split()[1]
    F = Function(FunctionSpace(mesh, 'CG', 1))
    F.interpolate(Expression('c>0.5?sqrt(x[0]*x[0]):0', degree = 1, c = chsub1))
    ymax = np.max(F.vector().get_local())
    return ymax

#----------------------------------------------------------------------
# Experimental data (20091030)
#----------------------------------------------------------------------

# Pipette radius
R_PIPE = 0.88e-3

# Initial radius of the droplet
RIC = 4.5e-3

# Imposed pressure on the right
p_right    = Constant(2.)

#----------------------------------------------------------------------
# General setting
#----------------------------------------------------------------------

# Mesh_type (1 = half mesh; 2 = full mesh)
mesh_type = 1

# Mesh_ref (1 coarse, 2 intermediate, 3 fine)
mesh_ref = 3

# Case (1, 2, 3)
case  = 3

#----------------------------------------------------------------------
# Model parameters
#----------------------------------------------------------------------

# Cells species and interface parameters
cequ     = 1.
sigma_cm = 0.02  

# Epsilon
epsilon  = R_PIPE/30.

# Diff_c to be identified
#Mc = 2.e-8    
#Mc = 1.e-8
#Mc = 8.e-9
Mc = 2.e-9

#----------------------------------------------------------------------

# if case == 1:
	# eta_c  = 1.

# if case == 2:
	# eta_c  = 3.

# if case == 3:
	# eta_c  = 5. 

eta_c  = 1.
    
if case == 1:
	Mc = 2.e-9

if case == 2:
	Mc = 3.e-9

if case == 3:
	Mc = 4.e-9 


eta_m  = 0.001    

#----------------------------------------------------------------------


#----------------------------------------------------------------------
# Other geometrical parameters
#----------------------------------------------------------------------
T_PIPE = 0.30e-3

L1 = 90.e-3
#L1 = 45.e-3
H1 = 3. * R_PIPE

L2 = L1 - (4.* H1)
#L2 = 45.e-3 - (4.* 3. * R_PIPE)

H2 = T_PIPE  

L3 = L2 - (0.5*H1)
H3 = H1 - R_PIPE - T_PIPE
RC = 0.5*T_PIPE 

FACT = 1.10
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# Relevant undimensional numbers
#----------------------------------------------------------------------
# The caracteristic lenght is the typical radius of the tube
l0 = R_PIPE
# The caracteristic velocity is a typical mean advancement 
# velocity in the experiment 100 microns / 200 sec.
v0 = 1.e-8
# Density cell
rhoc = 1000.
# Reynolds number
Re = (rhoc * v0 * l0)/eta_c
# Capillary number
Ca = (v0 * eta_c)/sigma_cm
# Peclet number
Pe = (v0*l0*l0)/(Mc*sigma_cm)
# Viscosity ratio
Vr = eta_m / eta_c
# Delta
De = l0/epsilon
# Theta
Th = 180

# Printing of the adimensional numbers
print('Re =', Re)
print('Ca =', Ca)
print('Pe =', Pe)
print('Vr =', Vr)
print('De =', De)
print('Th =', Th)
#----------------------------------------------------------------------



#----------------------------------------------------------------------
# Construction of the FE mesh
#----------------------------------------------------------------------
R1 = Rectangle(Point(0., 0.), Point(L1, H1))
R2 = Rectangle(Point(2*H1, R_PIPE), Point((2*H1 + L2), (R_PIPE+H2)))
R3 = Rectangle(Point((2.25*H1), (R_PIPE+T_PIPE)), Point((2.25*H1 + L3), (R_PIPE+T_PIPE+H3)))
R4 = Rectangle(Point(0., 0.), Point(L1, R_PIPE)) 
R5 = Rectangle(Point(0., (R_PIPE + T_PIPE)), Point((2.1*H1), (R_PIPE + 1.1*T_PIPE))) 
R6 = Rectangle(Point((L1 - 2.1*H1), (R_PIPE + T_PIPE)), Point(L1, (R_PIPE + 1.1*T_PIPE)))

C1 = Circle(Point(2*H1,(R_PIPE + 0.5*T_PIPE)), (FACT*RC), 16)
C2 = Circle(Point((2*H1 + L2),(R_PIPE + 0.5*T_PIPE)), (FACT*RC), 16)
domain1 = R1 - R2 - R3 - C1 - C2 - R4 - R5 - R6;
domain2 = domain1 + R4 + R5 + R6;

if mesh_type == 1:
    R7 = Rectangle(Point(0.5*L1, 0.), Point(L1, H1)) 
    domain2 = domain1 + R4 + R5 + R6 - R7

    
if mesh_ref == 1:
	mesh   = generate_mesh(domain2, 160)

if mesh_ref == 2:
	mesh   = generate_mesh(domain2, 200)

if mesh_ref == 3:
	mesh   = generate_mesh(domain2, 240)


markers = MeshFunction("bool", mesh, 2)
markers.set_all(False)
for c in cells(mesh):
    # Mark cells with facet midpoints near y == 1.0
    for f in facets(c):
        if ((f.midpoint()[0] < 0.017) and (f.midpoint()[0] > 0.0015)):
            markers[c] = True

mesh = refine(mesh, markers, redistribute=False)


# Mesh writting in file
file_mesh = File('mesh.pvd')
file_mesh << mesh
#----------------------------------------------------------------------





#----------------------------------------------------------------------
# Time discretization of the problem
#----------------------------------------------------------------------

TFIN0   = 60. 
TPHASE1 = 30.

dt0 = 1. 
num_steps = TFIN0 / dt0

TFIN1 = TFIN0 + TPHASE1
dt1   = 0.1
num_steps = (TFIN1 - TFIN0) / dt1
#----------------------------------------------------------------------



#----------------------------------------------------------------------
# Definition of functional space, trial & test functions
#----------------------------------------------------------------------
#___________________________________________________________________________
# Define function space
CP1  = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
CV2  = VectorElement("Lagrange", mesh.ufl_cell(), 2)
CHS  = FunctionSpace(mesh, MixedElement(CP1, CP1, CP1, CV2))
TENSOR = FunctionSpace(mesh, MixedElement(CP1, CP1, CP1, CP1, CP1, CP1, CP1, CP1, CP1))

#___________________________________________________________________________
# Define function & parameters
dsol              = TrialFunction(CHS)
v_m, v_c, q, w    = TestFunctions(CHS)
# Solution vector
sol                  = Function(CHS)     # Current solution
sol0                 = Function(CHS)     # Previous solution
m,  c,  p,  u      = split(sol)
m0, c0, p0, u0     = split(sol0)

# Definition of the normal vector
n      = FacetNormal(mesh)
x = SpatialCoordinate(mesh)
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# Definition of relevant bounds
# This time with boundary markers to extract subdomains
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# Definition of relevant bounds
# This time with boundary markers to extract subdomains
#----------------------------------------------------------------------
boundary_markers = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
boundary_markers.set_all(0)
tol  = 1.0e-5
tolc = 1.0e-5

#----------------------------------------------------------------------
# Bounds needed for the fluid flow problem
#----------------------------------------------------------------------
# Boundary_1(Left bound)
class Boundary_1(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[0], 0, tol)

bound_1 = Boundary_1()
bound_1.mark(boundary_markers, 1)

# Boundary_2U(upper bound big channel)
class Boundary_2U(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[1], H1 , tol)

bound_2U = Boundary_2U()
bound_2U.mark(boundary_markers, 2)


# Boundary_3U(Left and right mid vertical bound )
class Boundary_3U(SubDomain):
	def inside(self, x, on_boundary):
		return (near(x[0], ((L1 - L3)/2), tol) and (x[1] >= (R_PIPE + T_PIPE ))) \
                or (near(x[0], (L1 - L3)/2 + L3, tol) and (x[1] >= (R_PIPE + T_PIPE )))

bound_3U = Boundary_3U()
bound_3U.mark(boundary_markers, 3)


# Boundary_4U(Left and right upper bound small channel)
class Boundary_4U(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[1], (R_PIPE + T_PIPE ) , tol) and (x[0] >= ((0.5*(L1 - L2))-T_PIPE)) and (((0.5*(L1 - L2) + L2)+T_PIPE) >= x[0])

bound_4U = Boundary_4U()
bound_4U.mark(boundary_markers, 4)


# Boundary_5U(Left and right upper circular bound)

xc1 = 2*H1
yc1 = R_PIPE + 0.5*T_PIPE
xc2 = 2*H1 + L2
yc2 = R_PIPE + 0.5*T_PIPE

class Boundary_5U(SubDomain):
	def inside(self, x, on_boundary):
		return (abs(sqrt((x[1]-yc1)*(x[1]-yc1)+(x[0]-xc1)*(x[0]-xc1))-(FACT*RC)) < tolc and on_boundary) \
                or (abs(sqrt((x[1]-yc2)*(x[1]-yc2)+(x[0]-xc2)*(x[0]-xc2))-(FACT*RC)) < tolc and on_boundary)

bound_5U = Boundary_5U()
bound_5U.mark(boundary_markers, 5)



# Boundary_6U(upper bound small channel)
class Boundary_6U(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[1], R_PIPE, tol) and (x[0] > ((0.5*(L1 - L2))-T_PIPE)) and (x[0] < ((0.5*(L1 - L2) + L2)+T_PIPE))

bound_6U = Boundary_6U()
bound_6U.mark(boundary_markers, 6)


# Boundary_7(Right bound)

x_maxx = L1

if mesh_type == 1:
    x_maxx = 0.5*L1 

class Boundary_7(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[0], x_maxx, tol)

bound_7 = Boundary_7()
bound_7.mark(boundary_markers, 7)

# Boundary_8 (symmetry plan)
class Boundary_8(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[1], 0., tol)

bound_8 = Boundary_8()
bound_8.mark(boundary_markers, 8)

# Boundaries writting in file
file_boundary = File('boundary.pvd')
file_boundary << boundary_markers

# Redefiniction dx and ds
dx = Measure("dx", domain = mesh)
ds = Measure("ds", domain = mesh, subdomain_data = boundary_markers)
#----------------------------------------------------------------------




#----------------------------------------------------------------------
# Definition of dirichlet boundary conditions
#----------------------------------------------------------------------
v_null = Constant((0.0, 0.0))
zero   = Constant(0)

# No-slip boundary condition for velocity
bcNS2    = DirichletBC(CHS.sub(3), v_null, boundary_markers, 2)
bcNS3    = DirichletBC(CHS.sub(3), v_null, boundary_markers, 3)
bcNS4    = DirichletBC(CHS.sub(3), v_null, boundary_markers, 4)
bcNS5    = DirichletBC(CHS.sub(3), v_null, boundary_markers, 5)
bcNS6    = DirichletBC(CHS.sub(3), v_null, boundary_markers, 6)

# outflow boundary condition for velocity
bcNS7y = DirichletBC(CHS.sub(3).sub(1), zero,   boundary_markers, 7)

# inflow vy = 0
bcNS1y = DirichletBC(CHS.sub(3).sub(1), zero, boundary_markers, 1)

# Symmetry plan vy = 0
bcNS8y = DirichletBC(CHS.sub(3).sub(1), zero, boundary_markers, 8)

# Collect boundary conditions
bc_tot = [bcNS2, bcNS3, bcNS4, bcNS5, bcNS6, bcNS7y, bcNS1y, bcNS8y]

# No-wet condition
bc_ch4  = DirichletBC(CHS.sub(1), zero, boundary_markers, 4)
bc_ch5  = DirichletBC(CHS.sub(1), zero, boundary_markers, 5)
bc_ch6  = DirichletBC(CHS.sub(1), zero, boundary_markers, 6)

bc_tot = [bcNS2, bcNS3, bcNS4, bcNS5, bcNS6, bcNS7y, bcNS1y, bcNS8y, bc_ch4, bc_ch5, bc_ch6]
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# Definition of Initial condition
#----------------------------------------------------------------------
XC_0 = 10.e-3
YC_0 = 0.

u_init = LenghtInitialConditions(XC_0, YC_0, RIC, cequ, degree = 1)
sol.interpolate(u_init)
sol0.interpolate(u_init)
#----------------------------------------------------------------------


# Define variational problem
import dolfin as df
I = df.Identity(3)

beta = 1./4.
dfdc = 2*beta*c*(cequ-c)**2 - 2*(cequ-c)*beta*(c**2)

#----------------------------------------------------------------------
# Weak form
#----------------------------------------------------------------------

# Defivition of divrgence and gradient in cylindrycal coordinates
r = abs(x[1])

def grad_cyl(u):
	return as_tensor([[u[0].dx(0), u[0].dx(1), 0.], [u[1].dx(0), u[1].dx(1), 0.], [0., 0., u[1]/r]])

def div_cyl(u):
	return u[1]/r + u[0].dx(0) + u[1].dx(1)

def strain_rate(u): 
	return 0.5*(grad_cyl(u) + grad_cyl(u).T)

def sigma(u,p,eta_i): 
	return -1.*p*I + 2*eta_i*strain_rate(u)

T     = TensorFunctionSpace(mesh, 'CG', 1, (3,3))
stress_out = Function(TENSOR)
strain_out = Function(TENSOR)

# Wilson-theta
theta   = 0.5

# Interface regularization parameter
alfa     = (6*sqrt(2))/(cequ**3.)

# Theta-consistent quantities (to ajust no consistently implemented)
m_mid  = (1.0-theta)*m0 + theta*m
u_mid  = (1.0-theta)*u0 + theta*u
c_mid  = (1.0-theta)*c0 + theta*c
eta    = eta_m + c_mid*(eta_c - eta_m)

# Equation (28)
L0 = c*v_m*r*dx - c0*v_m*r*dx + dt0*Mc*dot(grad(m_mid), grad(v_m))*r*dx + dt0*dot(u_mid, grad(c_mid))*v_m*r*dx

# Equation (29)
L1 = (1.0/(epsilon*epsilon))*m*(epsilon/(alfa*sigma_cm))*v_c*r*dx - (1.0/(epsilon*epsilon))*dfdc*v_c*r*dx - dot(grad(c), grad(v_c))*r*dx 

# NS equations Divergence form
L2 =  1000.*dot(u, w)*(1./dt0)*r*dx - 1000.*dot(u0, w)*(1./dt0)*r*dx + eta*inner(grad_cyl(u), grad_cyl(w))*r*dx +  eta*inner(grad_cyl(u).T, grad_cyl(w))*r*dx -  p*div_cyl(w)*r*dx - m_mid*dot(grad(c_mid), w)*r*dx 

L3 = q*div_cyl(u)*r*dx 

# Assembling of the system of eqs (4) + (5) + (6) in weak form
L   = L0 + L1 + L2 + L3
#----------------------------------------------------------------------


#----------------------------------------------------------------------
#----------------------------------------------------------------------
# Solution of the nonlinear problem (PHASE 0 - for initial solution)
#----------------------------------------------------------------------
#----------------------------------------------------------------------

# Compute directional derivative about u in the direction of du (Jacobian)
a = derivative(L, sol, dsol)

# Create nonlinear problem and Newton solver
problem = NonlinearVariationalProblem(L, sol, bcs=bc_tot, J=a)
solver = NonlinearVariationalSolver(problem)
solver.parameters["newton_solver"]["linear_solver"] = "lu"
solver.parameters["newton_solver"]["convergence_criterion"] = "incremental"
solver.parameters["newton_solver"]["relative_tolerance"] = 1e-6
solver.parameters["newton_solver"]['maximum_iterations'] = 50

# Output file
tps = np.array([]) ; radi = np.array([])

# Saving results
if case  == 1:
    file0 = File("pipette_case1/potential.pvd",    "compressed")
    file1 = File("pipette_case1/cells.pvd",        "compressed")
    file2 = File("pipette_case1/pressure.pvd",     "compressed")
    file3 = File("pipette_case1/velocity.pvd",     "compressed")
    file4 = File("pipette_case1/strain_rate.pvd",  "compressed")
    file5 = File("pipette_case1/stress.pvd",       "compressed")

if case  == 2:
    file0 = File("pipette_case2/potential.pvd",    "compressed")
    file1 = File("pipette_case2/cells.pvd",        "compressed")
    file2 = File("pipette_case2/pressure.pvd",     "compressed")
    file3 = File("pipette_case2/velocity.pvd",     "compressed")
    file4 = File("pipette_case2/strain_rate.pvd",  "compressed")
    file5 = File("pipette_case2/stress.pvd",       "compressed")

if case  == 3:
    file0 = File("pipette_case3/potential.pvd",    "compressed")
    file1 = File("pipette_case3/cells.pvd",        "compressed")
    file2 = File("pipette_case3/pressure.pvd",     "compressed")
    file3 = File("pipette_case3/velocity.pvd",     "compressed")
    file4 = File("pipette_case3/strain_rate.pvd",  "compressed")
    file5 = File("pipette_case3/stress.pvd",       "compressed")

# Step in time
t = 0.
i = 1

# Calculation of strain rate and stress tensors
piii = sol.split()[2] 
uiii = sol.split()[3]
ciii = sol.split()[1]
eta_i = eta_m + ciii*(eta_c - eta_m)

dddi = strain_rate(uiii)
sssi = sigma(uiii,piii,eta_i)
strain_ratei = project(dddi,T)
stressi      = project(sssi,T)
stress_out.assign(stressi) 
strain_out.assign(strain_ratei) 

file0 << (sol.split()[0], t)
file1 << (sol.split()[1], t)
file2 << (sol.split()[2], t)
file3 << (sol.split()[3], t)
file4 << (strain_out,     t)
file5 << (stress_out,     t)

while (t < (TFIN0  - 0.0001)):

	t += dt0
	print('t =', t, 'dt = ', dt0, 'iteration numero', i)

	sol0.vector()[:] = sol.vector()
	solver.solve()

		# File outputs each 3 time steps (if i%3 == 0:)
	if i%1 == 0:
		piii = sol.split()[2] 
		uiii = sol.split()[3]
		ciii = sol.split()[1]
		eta_i = eta_m + ciii*(eta_c - eta_m)

		dddi = strain_rate(uiii)
		sssi = sigma(uiii,piii,eta_i)
		strain_ratei = project(dddi,T)
		stressi      = project(sssi,T)
		stress_out.assign(stressi) 
		strain_out.assign(strain_ratei) 

		file0 << (sol.split()[0], t)
		file1 << (sol.split()[1], t)
		file2 << (sol.split()[2], t)
		file3 << (sol.split()[3], t)
		file4 << (strain_out,     t)
		file5 << (stress_out,     t)
	
	i += 1;


#----------------------------------------------------------------------
#----------------------------------------------------------------------
# Solution of the nonlinear problem (PHASE 1 - aspiration)
#----------------------------------------------------------------------
#----------------------------------------------------------------------
m0, c0, p0, u0  = sol0.split()

# Equation (28)
L0 = c*v_m*r*dx - c0*v_m*r*dx + dt1*Mc*dot(grad(m_mid), grad(v_m))*r*dx + dt1*dot(u_mid, grad(c_mid))*v_m*r*dx

# Update of the boundary condition
L2 =  1000.*dot(u, w)*(1./dt1)*r*dx - 1000.*dot(u0, w)*(1./dt1)*r*dx + eta*inner(grad_cyl(u), grad_cyl(w))*r*dx +  eta*inner(grad_cyl(u).T, grad_cyl(w))*r*dx -  p*div_cyl(w)*r*dx - m_mid*dot(grad(c_mid), w)*r*dx + p_right*dot(n, w)*r*ds(7)

L   = L0 + L1 + L2 + L3

# Compute directional derivative about u in the direction of du (Jacobian)
a = derivative(L, sol, dsol)

# Create nonlinear problem and Newton solver
problem = NonlinearVariationalProblem(L, sol, bcs=bc_tot, J=a)
solver = NonlinearVariationalSolver(problem)
solver.parameters["newton_solver"]["linear_solver"] = "lu"
solver.parameters["newton_solver"]["convergence_criterion"] = "incremental"
solver.parameters["newton_solver"]["relative_tolerance"] = 1e-6
solver.parameters["newton_solver"]['maximum_iterations'] = 50


# Calculation of strain rate and stress tensors
piii = sol.split()[2] 
uiii = sol.split()[3]
ciii = sol.split()[1]
eta_i = eta_m + ciii*(eta_c - eta_m)

dddi = strain_rate(uiii)
sssi = sigma(uiii,piii,eta_i)
strain_ratei = project(dddi,T)
stressi      = project(sssi,T)
stress_out.assign(stressi) 
strain_out.assign(strain_ratei) 

# Step in time
i = 1
file0 << (sol.split()[0], t)
file1 << (sol.split()[1], t)
file2 << (sol.split()[2], t)
file3 << (sol.split()[3], t)
file4 << (strain_out,     t)
file5 << (stress_out,     t)


while (t < TFIN1):
	t += dt1
	print('t =', t, 'dt = ', dt1, 'iteration numero', i)
	
	sol0.vector()[:] = sol.vector()
	solver.solve()

		# File outputs each 3 time steps (if i%3 == 0:)
	if i%10 == 0:
		piii = sol.split()[2] 
		uiii = sol.split()[3]
		ciii = sol.split()[1]
		eta_i = eta_m + ciii*(eta_c - eta_m)

		dddi = strain_rate(uiii)
		sssi = sigma(uiii,piii,eta_i)
		strain_ratei = project(dddi,T)
		stressi      = project(sssi,T)
		stress_out.assign(stressi) 
		strain_out.assign(strain_ratei) 

		file0 << (sol.split()[0], t)
		file1 << (sol.split()[1], t)
		file2 << (sol.split()[2], t)
		file3 << (sol.split()[3], t)
		file4 << (strain_out,     t)
		file5 << (stress_out,     t)

	# At each time step, store time, aggregate radius and pressure
	tps = np.append(tps, t); radi = np.append(radi, (Max_height(sol, mesh)*1.E3))
	if case  == 1: np.savetxt("pipette_case1/sauvegarde.csv", radi, delimiter=",")
	if case  == 2: np.savetxt("pipette_case2/sauvegarde.csv", radi, delimiter=",")
	if case  == 3: np.savetxt("pipette_case3/sauvegarde.csv", radi, delimiter=",")
	
	i += 1;


if case  == 1:
    with open('pipette_case1/aspiration_phase1.csv', 'w') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        spamwriter.writerow(['tps', 'height'])
        for i in range (np.size(tps)):
            spamwriter.writerow([tps[i], radi[i]])

if case  == 2:
    with open('pipette_case2/aspiration_phase1.csv', 'w') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        spamwriter.writerow(['tps', 'height'])
        for i in range (np.size(tps)):
            spamwriter.writerow([tps[i], radi[i]])

if case  == 3:
    with open('pipette_case3/aspiration_phase1.csv', 'w') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        spamwriter.writerow(['tps', 'height'])
        for i in range (np.size(tps)):
            spamwriter.writerow([tps[i], radi[i]])
