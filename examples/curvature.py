"""
Example written by Niklaus Johner (niklaus.johner@a3.epfl.ch)

This is an example of curvature calculation done on a cylinder capped with half spheres.
It should be run from within the "example" directory, otherwise the path
to the python modules will have to be set.
"""
sys.path.append("../")
import surface_alg
import math

#First we build an Entity with atoms on a cylinder with spherical caps
eh_cyl=mol.CreateEntity()
edi=eh_cyl.EditXCS()
c=edi.InsertChain('A')
r=edi.AppendResidue(c,'b')
for i in range(30):
  for j in range(20):
    #edi.InsertAtom(r,'CA',geom.Vec3(i,4.*math.sin(j*math.pi/10.)+0.05*random.random(),4.*math.cos(j*math.pi/10.)+0.05*random.random()),'C')
    edi.InsertAtom(r,'CA',geom.Vec3(i,4.*math.sin(j*math.pi/10.),4.*math.cos(j*math.pi/10.)),'C')
i=29
for j in range(20):
  for k in range(1,5):
    edi.InsertAtom(r,'CA',geom.Vec3(i+4*math.sin(k*math.pi/10.),math.cos(k*math.pi/10.)*4.*math.sin(j*math.pi/10.),math.cos(k*math.pi/10.)*4.*math.cos(j*math.pi/10.)),'C')
k=5
edi.InsertAtom(r,'CA',geom.Vec3(i+4*math.sin(k*math.pi/10.),math.cos(k*math.pi/10.)*4.*math.sin(j*math.pi/10.),math.cos(k*math.pi/10.)*4.*math.cos(j*math.pi/10.)),'C')

#Now we display the capped cylinder in the scene.
go_cyl=gfx.Entity('cylinder',eh_cyl)
scene.Add(go_cyl)

#We set the normals as property to the atoms of the cylinder
for a in eh_cyl.atoms:
  n=geom.Normalize(geom.Vec3(0.0,a.pos.y,a.pos.z))
  a.SetFloatProp('nx',n.x)
  a.SetFloatProp('ny',n.y)
  a.SetFloatProp('nz',n.z)
for a in eh_cyl.Select('x>29').atoms:
  n=geom.Normalize(a.pos-geom.Vec3(29,0,0))
  a.SetFloatProp('nx',n.x)
  a.SetFloatProp('ny',n.y)
  a.SetFloatProp('nz',n.z)

#We display the normals in green in the scene
go_n=gfx.PrimList('cyl_normals')
for a in eh_cyl.atoms:
  n=geom.Vec3(a.GetFloatProp('nx'),a.GetFloatProp('ny'),a.GetFloatProp('nz'))
  go_n.AddLine(a.pos,a.pos+n)
scene.Add(go_n)
go_n.SetColor(gfx.GREEN)

#We calculate the curvature
surface_alg.CalculateCurvature(eh_cyl)

#We display the principal directions in red and blue
go_e1=gfx.PrimList('e1_cyl')
go_e2=gfx.PrimList('e2_cyl')
for a in eh_cyl.atoms:
  c=a.pos
  e1=geom.Vec3(a.GetFloatProp('e1x'),a.GetFloatProp('e1y'),a.GetFloatProp('e1z'))
  e2=geom.Vec3(a.GetFloatProp('e2x'),a.GetFloatProp('e2y'),a.GetFloatProp('e2z'))
  k1=a.GetFloatProp('k1')
  k2=a.GetFloatProp('k2')
  #if k1<0.1:
  go_e1.AddLine(c,c+2*k1*e1,gfx.RED)
  #if k2<0.1:
  go_e2.AddLine(c,c+2*k2*e2,gfx.BLUE)
scene.Add(go_e1)
scene.Add(go_e2)

#Now we print some average values of curvature
m1=sum([a.GetFloatProp('k1') for a in eh_cyl.Select('x<20').atoms])/float(eh_cyl.Select('x<20').GetAtomCount())
s1=sum([a.GetFloatProp('k1')**2.0-m1**2. for a in eh_cyl.Select('x<20').atoms])/float(eh_cyl.Select('x<20').GetAtomCount())**0.5
m2=sum([a.GetFloatProp('k2') for a in eh_cyl.Select('x<20').atoms])/float(eh_cyl.Select('x<20').GetAtomCount())
s2=sum([a.GetFloatProp('k2')**2.0-m2**2. for a in eh_cyl.Select('x<20').atoms])/float(eh_cyl.Select('x<20').GetAtomCount())**0.5
print "Cylinder and spherical cap have radius 5."
print 'The mean k1 on cylinder is {0} +/- {1}'.format(round(m1,2),round(s1,6))
print 'The mean k2 on cylinder is {0} +/- {1}'.format(round(m2,2),round(s2,6))
m1=sum([a.GetFloatProp('k1') for a in eh_cyl.Select('x>29').atoms])/float(eh_cyl.Select('x>29').GetAtomCount())
s1=sum([a.GetFloatProp('k1')**2.0-m1**2. for a in eh_cyl.Select('x>29').atoms])/float(eh_cyl.Select('x>29').GetAtomCount())**0.5
m2=sum([a.GetFloatProp('k2') for a in eh_cyl.Select('x>29').atoms])/float(eh_cyl.Select('x>29').GetAtomCount())
s2=sum([a.GetFloatProp('k2')**2.0-m2**2. for a in eh_cyl.Select('x>29').atoms])/float(eh_cyl.Select('x>29').GetAtomCount())**0.5
print 'The mean k1 on the sphere is {0} +/- {1}'.format(round(m1,2),round(s1,6))
print 'The mean k2 on the sphere is {0} +/- {1}'.format(round(m2,2),round(s2,6))




