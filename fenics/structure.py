#! /usr/bin/python2

dim = 3

tin = 0.
tout = 0.25
dt = 2e-3

lx = 3.
ly = 0.2
lz = 0.2
nx = 15
ny = 1
nz = 1

E = 1e8
ni = 0.3

from dolfin import *

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


  V = VectorFunctionSpace( mesh, 'CG', 2);
  Q = FunctionSpace( mesh, 'CG', 1);
  W = MixedFunctionSpace([V, V, Q])

  if dim == 2:
    zero = Constant(( 0., 0.))
  elif dim == 3:
    zero = Constant(( 0., 0., 0.))

  def Bleft(x, on_boundary):
    return on_boundary and abs(x[0] - 0.) < DOLFIN_EPS

  bc_d = DirichletBC(W.sub(0), zero, Bleft)
  bc_u = DirichletBC(W.sub(1), zero, Bleft)

  bcs = [bc_d, bc_u]

  ( d, u, p ) = TrialFunctions(W)
  ( b, v, q ) = TestFunctions(W)
  d_old = Function(V)
  u_old = Function(V)
  p_old = Function(Q)

  d_old = interpolate(zero,V)
  u_old = interpolate(zero,V)

  I = Identity(V.cell().d)
  F = I + grad(d_old)
  C = F.T*F

  Ic = tr(C)
  J = det(F)

  if dim == 2:
    f = Constant(( 0., 1.))
  elif dim == 3:
    f = Constant(( 0., 0., 1.))

  mu = Constant( E / 2*(1+ni) )
  lmbda = Constant( E*ni / (1+ni)*(1-2*ni) )

  def eps(v):
    return sym(grad(v))

  #def sigma(v):
  #  return 2.0*mu*eps(v) + lmbda*tr(eps(v))*I

  #  a = dot(d,b)*dx \
  #    - dt*dot(u,b)*dx \
  #    + dot(u,v)*dx \
  #    + dt*inner(sigma(d),eps(v))*dx \
  #    + dt*p*q*dx
  a = dot(d,b)*dx \
    - dt*dot(u,b)*dx \
    + dot(u,v)*dx \
    + dt*2.*mu*inner(eps(d),eps(v))*dx \
    + dt*p*div(v)*dx \
    + dt*q*div(d)*dx \
    - dt*(1./lmbda)*p*q*dx

  A = assemble(a)
  for bc in bcs:
    bc.apply(A)
  
  #prec = dot(d,b)*dx + dot(u,v)*dx + p*q*dx
  #P = assemble(prec)
  #for bc in bcs:
  #  bc.apply(P)
  #solver = KrylovSolver("gmres", "ilu")
  #solver.set_operator(A)
  #solver.parameters['monitor_convergence']=True

  dfile = File('d.pvd')
  ufile = File('u.pvd')

  sol = Function(W)
  t = tin
  dfile << (d_old,t)
  ufile << (u_old,t)
  while t + dt < tout+DOLFIN_EPS:
    t += dt
    print "time = %f" % t

    L = dot(d_old,b)*dx \
      + dot(u_old+dt*f,v)*dx

    rhs = assemble(L)
    for bc in bcs:
      bc.apply(rhs)

    #solver.solve(sol.vector(), rhs)
    solve( A, sol.vector(), rhs)

    d_old, u_old, p_old = sol.split()
    dfile << (d_old,t)
    ufile << (u_old,t)

if __name__ == '__main__':
  structure()
