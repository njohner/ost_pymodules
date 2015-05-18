"""
Example written by Niklaus Johner (niklaus.johner@a3.epfl.ch)

This is an example of curvature calculation done on a cylinder capped with half spheres.
It should be run from within the directory.
"""

#Test case: cylinder
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
go_cyl=gfx.Entity('cylinder',eh_cyl)
scene.Add(go_cyl)
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
go_n=gfx.PrimList('cyl_normals')
for a in eh_cyl.atoms:
  n=geom.Vec3(a.GetFloatProp('nx'),a.GetFloatProp('ny'),a.GetFloatProp('nz'))
  go_n.AddLine(a.pos,a.pos+n)
scene.Add(go_n)
CalculateCurvature(eh_cyl)
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

testdir='/Work/cubic_phase/test_curvature'
scene.Export(os.path.join(testdir,'cylinder_principal_curvature.png'))
m1=sum([a.GetFloatProp('k1') for a in eh_cyl.Select('x<20').atoms])/float(eh_cyl.Select('x<20').GetAtomCount())
s1=sum([a.GetFloatProp('k1')**2.0-m1**2. for a in eh_cyl.Select('x<20').atoms])/float(eh_cyl.Select('x<20').GetAtomCount())**0.5
m2=sum([a.GetFloatProp('k2') for a in eh_cyl.Select('x<20').atoms])/float(eh_cyl.Select('x<20').GetAtomCount())
s2=sum([a.GetFloatProp('k2')**2.0-m2**2. for a in eh_cyl.Select('x<20').atoms])/float(eh_cyl.Select('x<20').GetAtomCount())**0.5
print 'mean k1 on cylinder',m1,s1
print 'mean k2 on cylinder',m2,s2
m1=sum([a.GetFloatProp('k1') for a in eh_cyl.Select('x>29').atoms])/float(eh_cyl.Select('x>29').GetAtomCount())
s1=sum([a.GetFloatProp('k1')**2.0-m1**2. for a in eh_cyl.Select('x>29').atoms])/float(eh_cyl.Select('x>29').GetAtomCount())**0.5
m2=sum([a.GetFloatProp('k2') for a in eh_cyl.Select('x>29').atoms])/float(eh_cyl.Select('x>29').GetAtomCount())
s2=sum([a.GetFloatProp('k2')**2.0-m2**2. for a in eh_cyl.Select('x>29').atoms])/float(eh_cyl.Select('x>29').GetAtomCount())**0.5
print 'mean k1 on sphere',m1,s1
print 'mean k2 on sphere',m2,s2



