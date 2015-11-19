"""
.. codeauthor:: Niklaus Johner <niklaus.johner@a3.epfl.ch>

This module is used to mainly to calculate curvatures of a point set surface
Such surfaces can be obtained from densities using functions form the density_alg module
"""
from ost import *
import math
import numpy as npy
import scipy
from scipy import optimize
import time
import sys
sys.path.append('/Work/python_modules/')
import pbc_utilities

__all__=('CalculateSurface','CalculateSurface2','CalculateNormals','CalculateNormalsFast2',\
        'CalculateNormalsFast','OrientNormalsAlongDensity','OrientNormalsFast','OrientNormals'\
        ,'CleanMaxDensitySurface','CleanMaxDensitySurface2','CalculateCurvature')

def CalculateSurface(eh,sampling):
  """
  This function estimates the surface of a point set surface
  by simply summing up the number of points times the infinitesimal surface
  given by sampling**2.0
  """
  ds=sampling*sampling
  return ds*eh.GetAtomCount()

def CalculateSurface2(eh,within_size=5,PBC=False,cell_center=None,cell_size=None,float_prop_key=None):
  """
  This function estimates the surface of a point set surface
  by simply summing up the number of points times the infinitesimal surface
  given by sampling**2.0
  """ 
  if PBC:
    if not cell_center:
      print 'Entity bounds center used as cell_center'
      cell_center=eh.bounds.center
    if not cell_size:
      print 'Entity bounds size used as cell_size'
      cell_size=eh.bounds.size
  ds=math.pi*within_size*within_size
  s=0.
  for a in eh.atoms:
    if PBC:within=pbc_utilities.FindWithinWithPBC(eh,a.pos,within_size,cell_center,cell_size)
    else:within=mol.CreateViewFromAtoms([a2 for a2 in eh.FindWithin(a.pos,within_size)])
    s+=ds/float(within.GetAtomCount())
    if float_prop_key:a.SetFloatProp(float_prop_key,ds/float(within.GetAtomCount()))
  return s


def CalculateNormals(eh,within_size=5,PBC=False,cell_center=False,cell_size=False):
  """
  This function assigns normals to each atom in eh
  No specific orientation of the normals will come out
  """
  if PBC:
    if not cell_center:
      print 'Entity bounds center used as cell_center'
      cell_center=eh.bounds.center
    if not cell_size:
      print 'Entity bounds size used as cell_size'
      cell_size=eh.bounds.size
  count=0
  count_tot=0
  n_tot=eh.GetAtomCount()
  t1=time.time()
  for a in eh.atoms:
    count+=1
    count_tot+=1
    if PBC:within=pbc_utilities.FindWithinWithPBC(eh,a.pos,within_size,cell_center,cell_size)
    else:within=mol.CreateViewFromAtoms([a2 for a2 in eh.FindWithin(a.pos,within_size)]).Select('')
    vl=geom.Vec3List()
    for a2 in within.atoms:
      vl.append(a2.pos)
    if PBC:vl=geom.WrapVec3List(vl,a.pos,cell_size)
    n=vl.principal_axes.GetRow(0)
    a.SetVec3Prop('n',n)
    if count==5000:
      count=0
      print count_tot,'normals out of',n_tot,'in',time.time()-t1,'seconds'
  return

def CalculateNormalsFast2(eh,within_size=5,PBC=False,cell_center=False,cell_size=False):
  """
  This function assigns normals to each atom in eh
  No specific orientation of the normals will come out
  """
  if PBC:
    if not cell_center:
      print 'Entity bounds center used as cell_center'
      cell_center=eh.bounds.center
    if not cell_size:
      print 'Entity bounds size used as cell_size'
      cell_size=eh.bounds.size
  count=0
  count_tot=0
  n_tot=eh.GetAtomCount()
  t1=time.time()
  for a in eh.atoms:
    a.SetBoolProp('done',False)
  print 'set done to False',time.time()-t1
  for ref_a in eh.atoms:
    if ref_a.GetBoolProp('done'):continue
    if PBC:within_zone=pbc_utilities.FindWithinWithPBC(eh,ref_a.pos,2.*within_size,cell_center,cell_size)
    else:within_zone=mol.CreateViewFromAtoms([a2 for a2 in eh.FindWithin(ref_a.pos,2.*within_size)]).Select('')
    if PBC:within2=pbc_utilities.FindWithinWithPBC(within_zone,ref_a.pos,within_size,cell_center,cell_size)
    else:within2=mol.CreateViewFromAtoms([a2 for a2 in within_zone.FindWithin(ref_a.pos,within_size)]).Select('')    
    for a in within2.atoms:
      if a.GetBoolProp('done'):continue
      count+=1
      count_tot+=1
      if PBC:within=pbc_utilities.FindWithinWithPBC(within_zone,a.pos,within_size,cell_center,cell_size)
      else:within=mol.CreateViewFromAtoms([a2 for a2 in within_zone.FindWithin(a.pos,within_size)]).Select('')
      vl=geom.Vec3List()
      for a2 in within.atoms:
        vl.append(a2.pos)
      if PBC:vl=geom.WrapVec3List(vl,a.pos,cell_size)
      n=vl.principal_axes.GetRow(0)
      a.SetVec3Prop('n',n)
      a.SetBoolProp('done',True)
      if count==5000:
        count=0
        print count_tot,'normals out of',n_tot,'in',time.time()-t1,'seconds'
  return

def CalculateNormalsFast(eh,within_size=15,within_size2=7,PBC=False,cell_center=False,cell_size=False):
  """
  This function assigns normals to each atom in eh
  No specific orientation of the normals will come out
  """
  if PBC:
    if not cell_center:
      print 'Entity bounds center used as cell_center'
      cell_center=eh.bounds.center
    if not cell_size:
      print 'Entity bounds size used as cell_size'
      cell_size=eh.bounds.size
  count=0
  count_tot=0
  n_tot=eh.GetAtomCount()
  t1=time.time()
  for a in eh.atoms:
    a.SetBoolProp('done',False)
  print 'set done to False',time.time()-t1
  
  for a in eh.atoms:
    if a.GetBoolProp('done'):continue
    if PBC:within=pbc_utilities.FindWithinWithPBC(eh,a.pos,within_size,cell_center,cell_size)
    else:within=mol.CreateViewFromAtoms([a2 for a2 in eh.FindWithin(a.pos,within_size)]).Select('')
    vl=geom.Vec3List()
    for a2 in within.atoms:
      vl.append(a2.pos)
    if PBC:vl=geom.WrapVec3List(vl,a.pos,cell_size)
    n=vl.principal_axes.GetRow(0)
    if PBC:within=pbc_utilities.FindWithinWithPBC(eh,a.pos,within_size2,cell_center,cell_size)
    else:within=mol.CreateViewFromAtoms([a2 for a2 in eh.FindWithin(a.pos,within_size2)]).Select('')
    a.SetVec3Prop('n',n)
    a.SetBoolProp('done',True)
    for a2 in within.atoms:
      if a2.GetBoolProp('done'):continue
      count+=1
      count_tot+=1
      a2.SetVec3Prop('n',n)
      a2.SetBoolProp('done',True)
    if count>=5000:
      count=0
      print count_tot,'normals out of',n_tot,'in',time.time()-t1,'seconds'
  return

def OrientNormalsAlongDensity(eh,density):
  for a in eh.atoms:
    n=a.GetVec3Prop('n')
    d1=density.GetReal(img.Point(density.CoordToIndex(a.pos)))
    d2=density.GetReal(img.Point(density.CoordToIndex(a.pos+4*n)))
    if d2>d1:
      a.SetVec3Prop('n',-n)
  return

def OrientNormalsFast(eh,within_size=15,PBC=False,cell_center=None,cell_size=None):
  if PBC:
    if not cell_center:
      print 'Entity bounds center used as cell_center'
      cell_center=eh.bounds.center
    if not cell_size:
      print 'Entity bounds size used as cell_size'
      cell_size=eh.bounds.size
  t1=time.time()
  n_atoms=eh.GetAtomCount()
  print 'create ref view'
  for a in eh.atoms:
    a.SetIntProp('done',0)
  print 'prop set',time.time()-t1
  ref_view=eh.CreateEmptyView()
  grid_size=within_size/(2.0)
  for i in range(-int(cell_size[0]/(2*grid_size)),int(cell_size[0]/(2*grid_size))+1):
    xmin=cell_center[0]+grid_size*(i-0.5)
    xmax=cell_center[0]+grid_size*(i+0.5)
    v1=eh.Select('x>'+str(xmin)+' and x<='+str(xmax))
    for j in range(-int(cell_size[1]/(2*grid_size)),int(cell_size[1]/(2*grid_size))+1):
      xmin=cell_center[1]+grid_size*(j-0.5)
      xmax=cell_center[1]+grid_size*(j+0.5)
      v2=v1.Select('y>'+str(xmin)+' and y<='+str(xmax))
      for k in range(-int(cell_size[2]/(2*grid_size)),int(cell_size[2]/(2*grid_size))+1):
        xmin=cell_center[2]+grid_size*(k-0.5)
        xmax=cell_center[2]+grid_size*(k+0.5)
        within=v2.Select('z>'+str(xmin)+' and z<='+str(xmax))
        if len(within.atoms)>0:
          ref_view.AddAtom(within.atoms[0],mol.ViewAddFlag.CHECK_DUPLICATES)
          within.atoms[0].SetIntProp('done',1)
  print grid_size,'number of atoms in ref view',ref_view.GetAtomCount()
  print 'start orientation correction for ref view',time.time()-t1 
  starting_point=ref_view.bounds.min+cell_size/4.
  within=pbc_utilities.FindWithinWithPBC(ref_view,starting_point,within_size,cell_center,cell_size)
  if not within.IsValid():within=ref_view.CreateEmptyView()
  if within.GetAtomCount()==0:
    print 'start from random point'
    starting_point=ref_view.atoms[0].pos
  OrientNormals(ref_view,starting_point,within_size,PBC,cell_center,cell_size)
  total_atoms=ref_view.GetAtomCount()
  print 'start orientation correction for the other atoms',time.time()-t1 
  for ref_a in ref_view.atoms:
    ref_n=ref_a.GetVec3Prop('n')
    if PBC:within=pbc_utilities.FindWithinWithPBC(eh.Select('gadone!=1'),ref_a.pos,within_size,cell_center,cell_size)
    else:within=mol.CreateViewFromAtoms([a for a in eh.Select('gadone!=1').FindWithin(ref_a.pos,within_size)]).Select('')
    if not within.IsValid():within=eh.CreateEmptyView()
    for a in within.atoms:
      n=a.GetVec3Prop('n')
      if geom.Dot(ref_n,n)<0.0:
        a.SetVec3Prop('n',-n)
      a.SetIntProp('done',1)
      total_atoms+=1
      if total_atoms%5000==0:print total_atoms,'out of',n_atoms,'in',time.time()-t1,'seconds'
  print 'total number of atoms treated',total_atoms,'out of',n_atoms,'in',time.time()-t1,'seconds'
  return

  
def OrientNormals(eh,starting_point,within_size=10,PBC=False,cell_center=None,cell_size=None):
  if PBC:
    if not cell_center:
      print 'Entity bounds center used as cell_center'
      cell_center=eh.bounds.center
    if not cell_size:
      print 'Entity bounds size used as cell_size'
      cell_size=eh.bounds.size
  t1=time.time()
  count=0
  count2=0
  count3=0
  n_atoms=eh.GetAtomCount()
  total_atoms=0
  print 'start orientation correction'  
  if PBC:within=pbc_utilities.FindWithinWithPBC(eh,starting_point,within_size,cell_center,cell_size)
  else:within=mol.CreateViewFromAtoms([a for a in eh.FindWithin(starting_point,within_size)])
  if not within.IsValid():within=eh.CreateEmptyView()
  if within.GetAtomCount()==0:
    print 'no atoms close to the starting point'
    return
  total_atoms+=within.GetAtomCount()
  ref_a=within.atoms[0]
  ref_n=ref_a.GetVec3Prop('n')
  for a in within.atoms:
    count+=1
    if count%5000==0:print count,'out of',n_atoms
    n=a.GetVec3Prop('n')
    if geom.Dot(ref_n,n)<0.0:
      a.SetVec3Prop('n',-n)
  l=max(eh.bounds.size.data)
  ref_within=within.Copy()
  i=1
  flag=False
  while flag==False:
    i+=1
    radius=within_size*i
    ref_radius=within_size*(i-1)
    if PBC:
      within=pbc_utilities.FindWithinWithPBC(eh,starting_point,radius,cell_center,cell_size)
      ref_within=pbc_utilities.FindWithinWithPBC(eh,starting_point,ref_radius,cell_center,cell_size)
    else:
      within=mol.CreateViewFromAtoms([a for a in eh.FindWithin(starting_point,radius)]).Select('')
      ref_within=mol.CreateViewFromAtoms([a for a in eh.FindWithin(starting_point,ref_radius)]).Select('')
    if within.GetAtomCount()==eh.GetAtomCount():flag=True
    complement=mol.Difference(within,ref_within).Select('')
    total_atoms+=complement.GetAtomCount()
    print within.GetAtomCount(),complement.GetAtomCount(),total_atoms,'out of',n_atoms
    for a in complement.atoms:
      count+=1
      if count%5000==0:print count,'out of',n_atoms
      if PBC:within2=pbc_utilities.FindWithinWithPBC(ref_within,a.pos,within_size,cell_center,cell_size)
      else:within2=mol.CreateViewFromAtoms([a2 for a2 in ref_within.FindWithin(a.pos,within_size)])
      if not within2.IsValid():within2=ref_within.CreateEmptyView()
      if within2.GetAtomCount()==0:
        count2+=1
        if PBC:within2=pbc_utilities.FindWithinWithPBC(ref_within,a.pos,radius,cell_center,cell_size)
        else:within2=mol.CreateViewFromAtoms([a2 for a2 in ref_within.FindWithin(a.pos,radius)]).Select('')
        if not within2.IsValid():within2=ref_within.CreateEmptyView()
      try:a2=within2.atoms[0]
      except:
        count3+=1
        continue
      n=a.GetVec3Prop('n')
      ref_n=a2.GetVec3Prop('n')
      if geom.Dot(ref_n,n)<0.0:
        a.SetVec3Prop('n',-n)
  print count2,'times no atoms in within out of',n_atoms,'and ',count3,'atoms were not oriented'
  print 'total number of atoms treated',total_atoms
  return

def CleanMaxDensitySurface(eh,den_map,n_steps=2,step_size=1,PBC=False,cell_center=None,cell_size=None):
  """
  We keep only points that are at a maximum of the density in the direction of the normal
  """
  #The step size has to be at least of the same size as the spatial sampling of the density map
  if max(den_map.spatial_sampling.data)>step_size:
    step_size=max(den_map.spatial_sampling.data)
    print 'step_size is reset to the max(den_map.spatial_sampling)=',step_size
  if PBC:
    if not cell_center:
      print 'Entity bounds center used as cell_center'
      cell_center=eh.bounds.center
    if not cell_size:
      print 'Entity bounds size used as cell_size'
      cell_size=eh.bounds.size
  n_steps=int(n_steps)
  step_size=max(int(step_size),1)
  edi=eh.handle.EditXCS(mol.BUFFERED_EDIT)
  counter=0
  for j,a in enumerate(eh.atoms):
    if j%10000==0:print j,'atoms done and deleted',counter
    n=a.GetVec3Prop('n')
    d1=den_map.GetReal(img.Point(den_map.CoordToIndex(a.pos)))
    for i in range(step_size,n_steps*step_size+1):
      if PBC:
        x1=geom.WrapVec3(a.pos+i*n,cell_center,cell_size)
        x2=geom.WrapVec3(a.pos-i*n,cell_center,cell_size)
      else:
        x1=a.pos+i*n
        x2=a.pos-i*n      
      d2=den_map.GetReal(img.Point(den_map.CoordToIndex(x1)))
      d3=den_map.GetReal(img.Point(den_map.CoordToIndex(x2)))
      if d1<d2:
        counter+=1
        edi.DeleteAtom(a.handle)
        break
      if d1<d3:
        counter+=1
        edi.DeleteAtom(a.handle)
        break
  print 'deleted',counter,'atoms'
  return

def CleanMaxDensitySurface2(eh,den_map,clean_length=10,step_size=1,r=0.7,PBC=False,cell_center=None,cell_size=None):
  """
  We keep only points that have a larger density than the other points along the normal
  """
  t1=time.time()
  n_steps=int(clean_length/step_size)+1
  if PBC:
    if not cell_center:
      print 'Entity bounds center used as cell_center'
      cell_center=eh.bounds.center
    if not cell_size:
      print 'Entity bounds size used as cell_size'
      cell_size=eh.bounds.size
  edi=eh.handle.EditXCS(mol.BUFFERED_EDIT)
  counter=0
  for a in eh.atoms:
    a.SetBoolProp('done',False)
    a.SetFloatProp('delete',0)
  print 'set done to False',time.time()-t1
  for ref_a in eh.atoms:
    if ref_a.GetBoolProp('done') or ref_a.GetFloatProp('delete'):continue
    if PBC:within_zone=pbc_utilities.FindWithinWithPBC(eh,ref_a.pos,2.*(clean_length+r),cell_center,cell_size)
    else:within_zone=mol.CreateViewFromAtoms([a2 for a2 in eh.FindWithin(ref_a.pos,2.*(clean_length+r))]).Select('')
    if PBC:within2=pbc_utilities.FindWithinWithPBC(within_zone,ref_a.pos,(clean_length+r),cell_center,cell_size)
    else:within2=mol.CreateViewFromAtoms([a2 for a2 in within_zone.FindWithin(ref_a.pos,(clean_length+r))]).Select('')    
    #print ref_a, ref_a.pos,within_zone.GetAtomCount()
    for a in within2.atoms:
      if a.GetBoolProp('done') or a.GetFloatProp('delete'):continue
      #print a, within_zone.GetAtomCount()
      n=a.GetVec3Prop('n')
      d1=den_map.GetReal(img.Point(den_map.CoordToIndex(a.pos)))
      p=a.pos
      for i in range(1,n_steps):
        if PBC:
          x1=geom.WrapVec3(p+i*step_size*n,cell_center,cell_size)
          x2=geom.WrapVec3(p-i*step_size*n,cell_center,cell_size)
        else:
          x1=p+i*step_size*n
          x2=p-i*step_size*n  
        if PBC:w2=pbc_utilities.FindWithinWithPBC(within_zone,x1,r,cell_center,cell_size).Select('gadelete=0')
        else:w2=mol.CreateViewFromAtoms([a2 for a2 in within_zone.FindWithin(x1,r)]).Select('gadelete=0')    
        #return (i,within_zone,x1,r,cell_center,cell_size)
        for a2 in w2.atoms:
          d2=den_map.GetReal(img.Point(den_map.CoordToIndex(a2.pos)))
          if d2>d1:a.SetFloatProp('delete',1)
          elif d2<d1:a2.SetFloatProp('delete',1)
        if PBC:w2=pbc_utilities.FindWithinWithPBC(within_zone,x2,r,cell_center,cell_size).Select('gadelete=0')
        else:w2=mol.CreateViewFromAtoms([a2 for a2 in within_zone.FindWithin(x2,r)]).Select('gadelete=0')    
        for a2 in w2.atoms:
          d2=den_map.GetReal(img.Point(den_map.CoordToIndex(a2.pos)))
          if d2>d1:a.SetFloatProp('delete',1)
          elif d2<d1:a2.SetFloatProp('delete',1)
    for a in within_zone.Select('gadelete=1').atoms:
      edi.DeleteAtom(a.handle)
      counter+=1
    edi.ForceUpdate()
    print 'deleted',counter,'atoms'
  return



def CalculateCurvature(eh,within_size=5,normal_corr=False,PBC=False,cell_center=None,cell_size=None):
  """
  This function calculates the principal curvatures k1 and k2 and their
  direction (principal directions e1 and e2), as well as the gaussian
  and mean curvature for each atom in the entity and assigns them as FloatProp.
  Each atom should have a normal vector assigned to it beforehand, 
  as done by the CalculateNormals function. 
  """
  if PBC:
    if not cell_center:
      print 'Entity bounds center used as cell_center'
      cell_center=eh.bounds.center
    if not cell_size:
      print 'Entity bounds size used as cell_size'
      cell_size=eh.bounds.size
  t1=time.time()
  count=0
  count_tot=0
  n_tot=eh.GetAtomCount()
  for a in eh.atoms:
    try:
      ai=a.GetIndex()
      if PBC:within=pbc_utilities.FindWithinWithPBC(eh,a.pos,within_size,cell_center,cell_size)
      else:within=mol.CreateViewFromAtoms([a2 for a2 in eh.FindWithin(a.pos,within_size)]).Select('')
      count+=1
      count_tot+=1
      if count==250:
        count=0
        print 'curvature for',count_tot,'out of',n_tot,'in',time.time()-t1,'sec. Number of atoms in vicinity used for calculation',within.GetAtomCount()
      p=a.pos
      n=a.GetVec3Prop('n')
      #We define the local coordinate system
      try:
        psi=math.acos(n.z)
        if n.x!=0.0:phi=math.atan(n.y/n.x)
        elif n.y>0.0:phi=math.pi/2.
        elif n.y<=0.0:phi=-math.pi/2.
        x=geom.Vec3(-math.sin(phi),math.cos(phi),0)
        y=geom.Cross(n,x)
        #y=geom.Vec3(math.cos(psi)*math.cos(phi),math.cos(psi)*math.sin(phi),-math.sin(psi))
      except:
        print 'could not determine basis vec',ai,n.x,n.y,n.z
        v=geom.Vec3(1.,0,0)
        x=geom.Cross(n,v)
        y=geom.Cross(n,x)
      #Now we define the posistions and normals in the local coordinate system
      pl=geom.Vec3List()
      nl=geom.Vec3List()
      for a2 in within.atoms:
        if a2.handle==a.handle:continue
        if PBC:pi=geom.WrapVec3(a2.pos,p,cell_size)-p
        else:pi=a2.pos-p
        ni=a2.GetVec3Prop('n')
        pl.append(geom.Vec3(geom.Dot(pi,x),geom.Dot(pi,y),geom.Dot(pi,n)))
        niz=geom.Dot(ni,n)
        if normal_corr:
          if niz>=0.0:nl.append(geom.Vec3(geom.Dot(ni,x),geom.Dot(ni,y),niz))
          else:nl.append(geom.Vec3(-geom.Dot(ni,x),-geom.Dot(ni,y),-niz))
        else:nl.append(geom.Vec3(geom.Dot(ni,x),geom.Dot(ni,y),niz))
      #We calculate the normal curvature for each point
      kl=FloatList()
      for pi,ni in zip(pl,nl):
        kl.append(_CalculateNormalCurvature(pi,ni))
      #We calculate the angle to x at each point
      theta_l=FloatList()    
      for pi in pl:
        ti=geom.SignedAngle(geom.Vec2(1.0,0.0),geom.Vec2(pi.x,pi.y))
        if npy.isnan(ti):ti=0.0
        theta_l.append(ti)
      #Initial guess for k1,k2 and theta
      im=int(npy.argmax(kl))
      k1=kl[im]
      k2=min(kl)
      e1=pl[im]
      theta=geom.SignedAngle(geom.Vec2(1.0,0.0),geom.Vec2(e1.x,e1.y))
      #Now we optimize k1,k2 and theta
      xmin=optimize.fmin_powell(_CurvatureScore,[theta,k1,k2],args=(kl,theta_l),maxiter=20,full_output=True,maxfun=2000,disp=False)#,callback=_PrintX)
      [theta,k1,k2]=xmin[0]
      if npy.isnan(xmin[1]):
        [theta,k1,k2]=[npy.nan,npy.nan,npy.nan]
        print 'atom',ai,'did not converge','x',x,'y',y,'n',n,'k1',kl[im],'k2',min(kl),'theta',geom.SignedAngle(geom.Vec2(x.x,x.y),geom.Vec2(e1.x,e1.y))
      e1=geom.AxisRotation(n,theta)*x
      e2=geom.AxisRotation(n,theta+math.pi/2.)*x
    except:
      e1=e2=geom.Vec3(npy.nan,npy.nan,npy.nan)
      k1=k2=npy.nan    
    #We assign them to the atom
    a.SetFloatProp('k1',k1)
    a.SetFloatProp('k2',k2)
    a.SetFloatProp('e1x',e1.x)
    a.SetFloatProp('e1y',e1.y)
    a.SetFloatProp('e1z',e1.z)
    a.SetFloatProp('e2x',e2.x)
    a.SetFloatProp('e2y',e2.y)
    a.SetFloatProp('e2z',e2.z)
    a.SetFloatProp('Gauss',k1*k2)
    a.SetFloatProp('Mean',0.5*(k1+k2))
  print 'total time for',eh.GetAtomCount(),'atoms:',time.time()-t1
  return



def _PrintX(x):
  print 'current solution:',x

def _CalculateNormalCurvature(pi,ni):
  r=math.sqrt(pi.x*pi.x+pi.y*pi.y)
  nxy=(pi.x*ni.x+pi.y*ni.y)/r
  kni=-nxy/(math.sqrt(nxy*nxy+ni.z*ni.z)*r)
  if abs(kni)>10:print kni,pi,ni
  return kni
    
def _CurvatureScore(x,kl,theta_l):
  [theta,k1,k2]=x
  s=0.0
  for ki,ti in zip(kl,theta_l):
    s+=(k1*(math.cos(ti-theta)**2.0)+k2*(math.sin(ti-theta)**2.0)-ki)**2.0
  return s



#Test case: cylinder
"""
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
  a.SetVec3Prop('n',n)
for a in eh_cyl.Select('x>29').atoms:
  n=geom.Normalize(a.pos-geom.Vec3(29,0,0))
  a.SetVec3Prop('n',n)
go_n=gfx.PrimList('cyl_normals')
for a in eh_cyl.atoms:
  n=a.GetVec3Prop('n')
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
"""



