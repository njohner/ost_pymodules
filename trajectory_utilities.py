"""
Module written by Niklaus Johner (niklaus.johner@a3.epfl.ch) 01.2013
This module contains functions for filtering and building trajectories
"""
from ost import *
import time
import numpy as npy
import os,math

__all__=('CalculateInterfaceFromTraj','CalculateGeneralizedCorrelations','CalculateMutualInformation',\
        'CalculateCovariance','AverageMembraneThickness','WrapFloatList','RepresentPrincipalComponentOnStruccture',\
        'ProjectOnPrincipalComponentsAtomWise','ProjectOnPrincipalComponent','ReconstructTrajFromPrincipalComponents',\
        'CalculatePrincipalComponents','GetCellVectorsList','GetCellAnglesList','GetCellSizeList',\
        'CreatePerResidueCMTrajectory','ExtendTrajectoryToNeighboringUnitCells','TranslateFrames',\
        'WrapTrajectoryInPeriodicCell')

def WrapTrajectoryInPeriodicCell(t,centers,cell_sizes=None,cell_angles=None,group_res=False,follow_bonds=False):
  n=t.GetFrameCount()
  eh=t.GetEntity()
  if cell_sizes==None:
    cell_sizes=geom.Vec3List()
    for i in range(t.GetFrameCount()):
      f=t.GetFrame(i)
      if f.GetCellSize()==geom.Vec3(0,0,0):print 'cell size is 0 at step', i
      cell_sizes.append(f.GetCellSize())
  if cell_angles==None:
    cell_angles=geom.Vec3List()
    for i in range(t.GetFrameCount()):
      f=t.GetFrame(i)
      if f.GetCellAngles()==geom.Vec3(0,0,0):print 'cell angles are 0 at step', i, 'will be considered orthogonal'
      cell_angles.append(f.GetCellAngles())
  for i in range(t.GetFrameCount()):
    t.CopyFrame(i)
    mol.alg.WrapEntityInPeriodicCell(eh,centers[i],cell_sizes[i],cell_angles[i],group_res,follow_bonds)
    t.Capture(i)
  return

def TranslateFrames(t,trans_list):
  """
  This function translates each frame in the trajectory t by the corresponding Vec3
  in the trans_list. The transformation is done in place, so that the provided trajectory
  gets transformed.
  """
  eh=t.GetEntity()
  for i in range(t.GetFrameCount()):
    T=geom.Transform()
    eh.SetTransform(T)
    t.CopyFrame(i)
    T.SetTrans(trans_list[i])
    eh.SetTransform(T)
    t.Capture(i)
  return

def ExtendTrajectoryToNeighboringUnitCells(t,vecs_to_neighbor_ucells_list,cell_size_mult_factor=1):
  import entity_alg
  if not len(vecs_to_neighbor_ucells_list)==t.GetFrameCount():
    print 'you should provide a list of vectors for each frame'
  nreplicas=len(vecs_to_neighbor_ucells_list[0])
  if not all([len(vl)==nreplicas for vl in vecs_to_neighbor_ucells_list]):
    print 'Needs the same number of vectors for each frame'
  eh=t.GetEntity()
  extended_eh=entity_alg.ExtendEntityToNeighboringUnitCells(eh,vecs_to_neighbor_ucells_list[0])
  extended_t=mol.CreateCoordGroup(extended_eh.atoms)
  for i in range(t.GetFrameCount()):
    T=geom.Transform()
    eh.SetTransform(T)
    t.CopyFrame(i)
    vl=eh.GetAtomPosList()
    for v in vecs_to_neighbor_ucells_list[i]:
      T.SetTrans(v)
      eh.SetTransform(T)
      vl.extend(eh.GetAtomPosList())
    extended_t.AddFrame(vl,cell_size_mult_factor*t.GetFrame(i).GetCellSize(),t.GetFrame(i).GetCellAngles())
  return extended_t

def CreatePerResidueCMTrajectory(eh,t,cm_sele_list=[],atom_name_list=[],superposition_view=None):
  """
  This function will generate a trajectory containing for each residue the center of masse
  for each selection in the selection list. The atoms corresponding to these COMs will be named
  according to the atom_name_list.
  If provided the trajectory will first be aligned on the superposition_view (has to be from the same entity).
  As default it gives one CM for the backbone heavy atoms and one for the sidechain heavy atoms,
  except for GLY and PRO for which only the backbone is used.
  """
  time1=time.time()
  if not cm_sele_list:cm_sele_list=['aname=N,CA,C,O','aname!=N,CA,C,O and ele!=H and rname!=GLY,PRO']
  if not atom_name_list:atom_name_list=['BB','SC']
  if len(atom_name_list)!=len(cm_sele_list):atom_name_list=['A'+str(i) for i in range(len(cm_sele_list))]
  if superposition_view:
    t=mol.alg.SuperposeFrames(t,superposition_view,superposition_view)
    eh=t.GetEntity()
  print 'trajectory superposed in',time.time()-time1,'seconds'
  time2=time.time()
  eh2=mol.CreateEntity()
  edi=eh2.EditXCS()
  vl_list=[]
  for chain in eh.chains:
    c=edi.InsertChain(chain.name)
    for res in chain.residues:
      r=edi.AppendResidue(c,res.name,res.number.num)
      for aname,sele in zip(atom_name_list,cm_sele_list):
        v=res.Select(sele)
        if v.GetAtomCount()==0:continue
        edi.InsertAtom(r,aname,geom.Vec3())
        vl_list.append(mol.alg.AnalyzeCenterOfMassPos(t,res.Select(sele)))
  print 'cm positions determined in',time.time()-time2
  time2=time.time()
  t2=mol.CreateCoordGroup(eh2.atoms)
  for i in range(t.GetFrameCount()):
    v=geom.Vec3List()
    for vl in vl_list:
      v.append(vl[i])
    t2.AddFrame(v)
  t2.CopyFrame(0)
  print 'finished trajectory generation in',time.time()-time2,'. Total time',time.time()-time1
  return t2

def GetCellSizeList(t):
  cell_sizes=geom.Vec3List()
  for i in range(t.GetFrameCount()):
    f=t.GetFrame(i)
    if f.GetCellSize()==geom.Vec3(0,0,0):print 'cell size is 0 at step', i
    cell_sizes.append(f.GetCellSize())
  return cell_sizes
    
def GetCellAnglesList(t):
  cell_angles=geom.Vec3List()
  for i in range(t.GetFrameCount()):
    f=t.GetFrame(i)
    if f.GetCellAngles()==geom.Vec3(0,0,0):print 'cell angles are 0,0,0 at step', i
    cell_angles.append(f.GetCellAngles())
  return cell_angles
 
def GetCellVectorsList(t):
  cell_vectors=[t.GetFrame(i).GetCellVectors() for i in range(t.GetFrameCount())]
  return cell_vectors

def CalculatePrincipalComponents(t,pc_sele,superposition_sele=None,first=0,last=-1,stride=1):
  """
  This function calculates the principal components for the positions
  of the N atoms in ev, i.e the 3N eigenvalues and eigenvectors.
  Specifically it performs an svd decomposition A=U*S*V
  Return:
    -The unitary matrix U containing the eigenvectors in its columns
    -The singular values S, so that S*S are the eigenvalues
    -The unitary matrix V
  """
  if not _import_numpy:return False
  if last==-1:last=t.GetFrameCount()
  nframes=(last-first)/stride
  if superposition_sele!=None:t2=mol.alg.SuperposeFrames(t,t.GetEntity().Select(superposition_sele))
  else:t2=t
  ev=t2.GetEntity().Select(pc_sele)
  natoms=ev.GetAtomCount()
  atom_pos_list=npy.zeros([3*natoms,nframes])
  mean_atom_pos=npy.zeros(3*natoms)
  for i,a in enumerate(ev.atoms):
    vl=mol.alg.AnalyzeAtomPos(t2,a.handle,stride)[first:last]
    for j in range(3):
      vj=[v[j] for v in vl]
      vjm=npy.mean(vj)
      mean_atom_pos[3*i+j]=vjm
      atom_pos_list[3*i+j]=vj-vjm
  (U,S,VT)=npy.linalg.svd(atom_pos_list,full_matrices=False)
  return (U,S,VT.T,mean_atom_pos,atom_pos_list)

def ReconstructTrajFromPrincipalComponents(ev,U,S,V,mean_atom_pos,pc_indices=[0,1]):
  natoms=ev.GetAtomCount()
  nframes=V.shape[0]
  atom_pos_list=npy.dot(U[:,pc_indices],npy.dot(npy.diag(S)[npy.ix_(pc_indices,pc_indices)],V[:,pc_indices].T))
  atom_pos_list=npy.array([atom_pos_list[i]+mean_atom_pos[i] for i in range(3*natoms)])
  eh=mol.CreateEntityFromView(ev,1)
  t=mol.CreateCoordGroup(eh.atoms)
  for j in range(nframes):t.AddFrame(geom.Vec3List([geom.Vec3(atom_pos_list[3*i,j],atom_pos_list[3*i+1,j],atom_pos_list[3*i+2,j]) for i in range(natoms)]))
  return t

def ProjectOnPrincipalComponent(atom_pos_list,U,pc_index=0,first=0,last=-1):
  if last==-1:last=atom_pos_list.shape[1]
  return npy.dot(U[:,pc_index],atom_pos_list[:,first:last])

def ProjectOnPrincipalComponentsAtomWise(t,ev,first=0,last=-1,stride=1):
  """
  Calculates the principal components for each atom and projects the trajectory
  on them.
  It returns:
    - The list of principal axes in which the ith element is a 3x3 matrix
    with in its columns the principal axes for the ith atom in ev.
    - The positions projected on the principal axes
  """
  if last==-1:last=t.GetFrameCount()
  pvl=[]
  pcl=geom.Vec3List()
  for a in ev.atoms:
    vl=mol.alg.AnalyzeAtomPos(t,a.handle,stride)[first:last]
    pc=vl.principal_axes
    c=vl.center
    vl=geom.Vec3List([el-c for el in vl])
    pc1=pc.GetRow(0)
    pc2=pc.GetRow(1)
    pc3=pc.GetRow(2)
    pvl.append(geom.Vec3List([geom.Vec3(geom.Dot(pc1,v),geom.Dot(pc1,v),geom.Dot(pc1,v)) for v in vl]))
    pcl.append(pc)
  return (pcl,pvl)

def RepresentPrincipalComponentOnStruccture(ev,U,pc_index=0,go_name='pc0',color=gfx.RED):
  go=gfx.PrimList(go_name)
  for i in range(ev.GetAtomCount()):
    d=geom.Vec3(U[3*i,pc_index],U[3*i+1,pc_index],U[3*i+2,pc_index])
    p=ev.atoms[i].pos
    go.AddLine(p,p+50*d,color)
  return go
  
def WrapFloatList(fl,center,period):
  for i in range(len(fl)):
    fl[i]=fl[i]-period*round((fl[i]-center)/period)
  return fl
  
def AverageMembraneThickness(t,upper,lower,step_size=2):
  if not _import_numpy:return False
  import entity_alg
  step_size=float(step_size)
  nf=t.GetFrameCount()
  nu=upper.GetAtomCount()
  nl=lower.GetAtomCount()
  xu=npy.zeros([nf*nu,3])
  xl=npy.zeros([nf*nl,3])
  for i in range(nf):
    t.CopyFrame(i)
    for j,a in enumerate(upper.atoms):
      p=a.pos
      xu[i*nu+j]=(p[0],p[1],p[2])
  for i in range(nf):
    t.CopyFrame(i)
    for j,a in enumerate(lower.atoms):
      p=a.pos
      xl[i*nl+j]=(p[0],p[1],p[2])
  xmin=int(min(xu[0:].min(),xl[0:].min()))
  xmax=int(max(xu[0:].max(),xl[0:].max()))+1
  ymin=int(min(xu[1:].min(),xl[1:].min()))
  ymax=int(max(xu[1:].max(),xl[1:].max()))+1
  n1=npy.ceil((xmax-xmin)/step_size)
  n2=npy.ceil((ymax-ymin)/step_size)
  zu=npy.zeros([n1,n2])
  cu=npy.zeros([n1,n2])
  zl=npy.zeros([n1,n2])
  cl=npy.zeros([n1,n2])
  for p in xu:
    zu[(p[0]-xmin)/step_size,(p[1]-ymin)/step_size]+=p[2]
    cu[(p[0]-xmin)/step_size,(p[1]-ymin)/step_size]+=1
  for p in xl:
    zl[(p[0]-xmin)/step_size,(p[1]-ymin)/step_size]+=p[2]
    cl[(p[0]-xmin)/step_size,(p[1]-ymin)/step_size]+=1
  zu=zu/cu
  zl=zl/cl
  thickness=zu-zl
  #Now we create Vec3Lists with the positions of the beads for the lower and upper leaflets
  pu=geom.Vec3List()
  pl=geom.Vec3List()
  for i in range(zu.shape[0]):
    for j in range(zu.shape[1]):
      pu.append(geom.Vec3(xmin+(i+0.5)*step_size,ymin+(j+0.5)*step_size,zu[i,j]))
      pl.append(geom.Vec3(xmin+(i+0.5)*step_size,ymin+(j+0.5)*step_size,zl[i,j]))
  eh_upper=entity_alg.CreateEntityFromVec3List(pu,chain_name='MEMU')
  eh_lower=entity_alg.CreateEntityFromVec3List(pl,chain_name='MEML')
  n1=zu.shape[0]
  for i in range(zu.shape[0]):
    for j in range(zu.shape[1]):
      eh_upper.atoms[n1*i+j].SetFloatProp('thickness',thickness[i,j])
      eh_upper.atoms[n1*i+j].SetFloatProp('count',cu[i,j])
      eh_lower.atoms[n1*i+j].SetFloatProp('thickness',thickness[i,j])
      eh_lower.atoms[n1*i+j].SetFloatProp('count',cl[i,j])
  #we clean the entities:
  edi=eh_upper.EditXCS(mol.BUFFERED_EDIT)
  v=mol.Difference(eh_upper.Select(''),eh_upper.Select('(z>=0 or z<0)'))
  for a in v.atoms:
    edi.DeleteAtom(a.handle)
  edi.ForceUpdate()
  edi=eh_lower.EditXCS(mol.BUFFERED_EDIT)
  v=mol.Difference(eh_lower.Select(''),eh_lower.Select('(z>=0 or z<0)'))
  for a in v.atoms:
    edi.DeleteAtom(a.handle)
  edi.ForceUpdate()
  return (eh_lower,eh_upper)

def _import_numpy():
  try:
    import numpy as npy
  except: 
    raise ImportError("This function uses numpy, which could not be imported")
    return False
  return True

def CalculateCovariance(t,view,superposition_view=None):
  """
  Calculates the covariance matrix C(3Nx3N) for the N atoms in the view.
  If superposition_view is defined, the trajectory is first aligned.
  """
  if not _import_numpy:return False
  if superposition_view:t2=mol.alg.SuperposeFrames(t,superposition_view)
  else:t2=t
  natoms=view.GetAtomCount()
  pos_list=[mol.alg.AnalyzeAtomPos(t2,a.handle) for a in view.atoms]
  pos_list=npy.array([npy.array([[p[0] for p in pl],[p[1] for p in pl],[p[2] for p in pl]]) for pl in pos_list])
  pos_list=pos_list.reshape(pos_list.shape[0]*pos_list.shape[1],pos_list.shape[2])
  pos_list=npy.array([(pl-npy.average(pl))/npy.std(pl) for pl in pos_list])
  C=npy.cov(pos_list)
  return C

def CalculateMutualInformation(t,view,superposition_view=None):
  """
  Calculates the mutual information between the N atoms in the view.
  If superposition_view is defined, the trajectory is first aligned.
  """
  if not _import_numpy:return False
  natoms=view.GetAtomCount()
  d=3.
  C=CalculateCovariance(t,view,superposition_view)
  D=npy.zeros([natoms,natoms])
  for i in range(natoms):
    M=C[3*i:3*i+3,3*i:3*i+3]
    D[i,i]=npy.linalg.det(M)
    for j in range(i+1,natoms):
      M=npy.vstack((npy.hstack((C[3*i:3*i+3,3*i:3*i+3],C[3*i:3*i+3,3*j:3*j+3])),npy.hstack((C[3*j:3*j+3,3*i:3*i+3],C[3*j:3*j+3,3*j:3*j+3]))))
      D[i,j]=npy.linalg.det(M)
      D[j,i]=D[i,j]
  I=npy.zeros([natoms,natoms])
  for i in range(natoms):
    I[i,i]=npy.inf
    for j in range(i+1,natoms):
      I[i,j]=0.5*(math.log(D[i,i])+math.log(D[j,j])-math.log(D[i,j]))
      I[j,i]=I[i,j]
  return I

def CalculateGeneralizedCorrelations(t,view, superposition_view):
  """
  Calculates the generalized correlation coefficient between the N atoms in the view,
  as defined in: "Generalized Correlation for biomolecular dynamics", O. F. Lange and H. Grubmuller,
  Protein: Structure, Function and Bioinformatics 62:1053-1061 (2006)
  If superposition_view is defined, the trajectory is first aligned.
  """
  if not _import_numpy:return False
  natoms=view.GetAtomCount()
  d=3.
  I=CalculateMutualInformation(t,view)
  R=npy.zeros([natoms,natoms])
  for i in range(natoms):
    R[i,i]=1.0
    for j in range(i+1,natoms):
      R[i,j]=(1-math.exp(-2*I[i,j]/d))**0.5
      R[j,i]=R[i,j]
  return R

  
def _SetPBC(sele_view,PBC,cell_center,cell_size):
  if PBC:
    if not cell_center:
      print 'EntityView bounds center used as cell_center'
      cell_center=sele_view.bounds.center
    if not cell_size:
      print 'EntityView bounds size used as cell_size'
      cell_size=sele_view.bounds.size
  return (PBC,cell_center,cell_size)

def CalculateInterfaceFromTraj(t,waters,lipids,PBC=False,cell_center=None,cell_size=None,low_pass_filter_length=10,\
                      density_cutoff=None,outdir='',outname='',stride=1):
  import density_alg
  (PBC,cell_center,cell_size)=_SetPBC(waters,PBC,cell_center,cell_size)
  if PBC:
    water_den=density_alg.CreateDensityMapFromTrajectoryWithPBC(t,waters,cell_center=cell_center,cell_size=cell_size,stride=stride)
    lipid_den=density_alg.CreateDensityMapFromTrajectoryWithPBC(t,lipids,cell_center=cell_center,cell_size=cell_size,stride=stride)
  else:
    water_den=density_alg.CreateDensityMapFromTrajectory(t,waters,stride=stride)
    lipid_den=density_alg.CreateDensityMapFromTrajectory(t,lipids,stride=stride)
  water_filtered=water_den.Apply(img.alg.LowPassFilter(low_pass_filter_length))
  lipid_filtered=lipid_den.Apply(img.alg.LowPassFilter(low_pass_filter_length))
  if outdir:
    io.SaveMap(water_den,os.path.join(outdir,outname+'water_density.mrc'))
    io.SaveMap(water_filtered,os.path.join(outdir,outname+'water_density_filtered.mrc'))
    io.SaveMap(lipid_den,os.path.join(outdir,outname+'lipid_density.mrc'))
    io.SaveMap(lipid_filtered,os.path.join(outdir,outname+'lipid_density_filtered.mrc'))
  if not density_cutoff:
    dml=0;dmw=0
    for el in lipid_filtered:
      if lipid_filtered.GetReal(el)>dml:dml=lipid_filtered.GetReal(el)
    for el in water_filtered:
      if water_filtered.GetReal(el)>dmw:dmw=water_filtered.GetReal(el)
    density_cutoff=(dml+dmw)/20.
  boundary=density_alg.GetBoundaryBetweenDensities(water_filtered,lipid_filtered,density_cutoff,PBC,cell_center,cell_size)
  return (water_filtered,lipid_filtered,boundary)

  


