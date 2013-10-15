#! /usr/bin/python2

dim = 2

#tin = 0.
#tout = 0.25
#dt = 2e-3

lx = 3.
ly = 0.2
lz = 0.2
nx = 15
ny = 1
nz = 1

E = 1e8
ni = 0.3

from dolfin import *
 
E = 1e8
ni = 0.3
mu = Constant( E / (2*(1+ni)) )
lmbda = Constant( E*ni / ((1+ni)*(1-2*ni)) )

def structure():
  
  #parameters["form_compiler"]["optimize"]     = True
  #parameters["form_compiler"]["cpp_optimize"] = True

  #n = 4
  #mesh = UnitSquareMesh(n,n,"crossed")
  #mesh = UnitCubeMesh(4,4,4)
  #mesh = Mesh("beam2d.xml")
  if(dim == 2):
    mesh = RectangleMesh(0., 0., lx, ly, nx, ny)
  elif(dim == 3):
    mesh = BoxMesh(0., 0., 0., lx, ly, lz, nx, ny, nz)

  V = VectorFunctionSpace( mesh, 'CG', 1);

  if dim == 2:
    zero = Constant(( 0., 0.))
  elif dim == 3:
    zero = Constant(( 0., 0., 0.))

  def Bleft(x, on_boundary):
    return on_boundary and x[0] == 0.

  bc = DirichletBC(V, zero, Bleft)

  d = TrialFunction(V)
  b = TestFunction(V)
  
  I = Identity(V.cell().d)

  if dim == 2:
    f = Constant(( 0., 1.))
  elif dim == 3:
    f = Constant(( 0., 0., 1.))

  def eps(v):
    return sym(grad(v))

  def sigma(v):
    return 2.0*mu*eps(v) + lmbda*tr(eps(v))*I

  #dt = 200.

  a = inner(sigma(d),eps(b))*dx
  A = assemble(a)

  datafile = File('structure.pvd')

  L = dot(f,b)*dx
  rhs = assemble(L)
  bc.apply(A,rhs)
  disp = Function(V)
  solve( A, disp.vector(), rhs)
  datafile << disp

if __name__ == '__main__':
  structure()
