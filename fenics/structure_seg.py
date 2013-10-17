#! /usr/bin/python2

output_dir = 'output/seg'

#mesh_file = 'hex2.xml'
mesh_file = None
dim = 2

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

import os
from dolfin import *

def structure():
  
  #parameters["form_compiler"]["optimize"]     = True
  #parameters["form_compiler"]["cpp_optimize"] = True

  #n = 4
  #mesh = UnitSquareMesh(n,n,"crossed")
  #mesh = UnitCubeMesh(4,4,4)
  #mesh = Mesh("beam2d.xml")
  if mesh_file is not None:
   mesh = Mesh(mesh_file)
  if dim == 2:
    mesh = RectangleMesh(0., 0., lx, ly, nx, ny)
  elif dim == 3:
    mesh = BoxMesh(0., 0., 0., lx, ly, lz, nx, ny, nz)

  V = VectorFunctionSpace( mesh, 'CG', 2);
  Q = FunctionSpace( mesh, 'CG', 1);
  W = MixedFunctionSpace([V, Q])

  if(dim == 2):
    zero = Constant(( 0., 0.))
  elif(dim == 3):
    zero = Constant(( 0., 0., 0.))

  def Bleft(x, on_boundary):
    return on_boundary and abs(x[0] - 0.) < DOLFIN_EPS

  bc_d = DirichletBC(V, zero, Bleft)
  bc_u = DirichletBC(W.sub(0), zero, Bleft)

  d = TrialFunction(V)
  ( u, p ) = TrialFunctions(W)
  b = TestFunction(V)
  ( v, q ) = TestFunctions(W)
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

  if (dim == 2):
    f = Constant(( 0., 1.))
  elif (dim == 3):
    f = Constant(( 0., 0., 1.))

  mu = Constant( E / 2*(1+ni) )
  lmbda = Constant( E*ni / (1+ni)*(1-2*ni) )
  ilambda = 1. / lmbda

  def eps(v):
    return sym(grad(v))

  def sigma(v):
    return 2.0*mu*eps(v) + lmbda*tr(eps(v))*I

  a_d = dot(d,b)*dx
  #a_u = dot(u,v)*dx \
      #+ dt*inner(sigma(dt*u),eps(v))*dx \
      #+ dt*p*q*dx
  a_u = dot(u,v)*dx \
      + dt*2.*mu*inner(eps(dt*u),eps(v))*dx \
      + dt*p*div(v)*dx \
      + dt*q*div(u)*dx \
      - ilambda*p*q*dx

  A_d = assemble(a_d)
  A_u = assemble(a_u)

  bc_d.apply(A_d)
  bc_u.apply(A_u)
  
  #prec = dot(d,b)*dx + dot(u,v)*dx + p*q*dx
  #P = assemble(prec)
  #for bc in bcs:
  #  bc.apply(P)
  #solver_d = KrylovSolver("minres", "mg")
  #solver_d.set_operator(A_d)
  #solver_d.parameters['monitor_convergence']=True

  if not os.path.exists(output_dir):
        os.makedirs(output_dir)
  file_d = File(output_dir + '/d.pvd')
  file_u = File(output_dir + '/u.pvd')
  file_p = File(output_dir + '/p.pvd')

  sol = Function(W)
  t = tin
  file_d << (d_old,t)
  file_u << (u_old,t)
  file_p << (p_old,t)
  while t + dt < tout+DOLFIN_EPS:
    t += dt
    print "time = %f" % t

    #L_u = dot(u_old+dt*f,v)*dx \
        #- dt*inner(sigma(d_old),eps(v))*dx
    L_u = dot(u_old+dt*f,v)*dx \
        - dt*2.*mu*inner(eps(d_old),eps(v))*dx \
        - div(d_old)*q*dx

    rhs_u = assemble(L_u)
    bc_u.apply(rhs_u)

    #solver.solve(sol.vector(), rhs)
    solve( A_u, sol.vector(), rhs_u)
    u_old, p_old = sol.split()
    
    L_d = dot(d_old+dt*u_old,b)*dx #+ dt*dot(u_old,b)*dx
    rhs_d = assemble(L_d)
    bc_d.apply(rhs_d)
    solve( A_d, d_old.vector(), rhs_d)

    file_d << (d_old,t)
    file_u << (u_old,t)
    file_p << (p_old,t)

if __name__ == '__main__':
  structure()

