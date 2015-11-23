"""
.. codeauthor: Niklaus Johner <niklaus.johner@a3.epfl.ch>

This module contains functions for filtering and building trajectories
"""
from ost import *
import time
import numpy as npy
import os,math

__all__=('CalculateInterfaceFromTraj','CalculateGeneralizedCorrelations','CalculateMutualInformation',\
        'CalculateCovariance','AverageMembraneThickness','SmoothAverageMembrane','WrapFloatList',\
        'GetCellVectorsList','GetCellAnglesList','GetCellSizeList',\
        'CreatePerResidueCMTrajectory','ExtendTrajectoryToNeighboringUnitCells','ExtendTrajectory','TranslateFrames',\
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

def ExtendTrajectoryToNeighboringUnitCells(t,extension_directions,cell_size_mult_factors=(1,1,1)):
  """
  This function is used to extend a trajectory to its neighboring unit cells.
  Specifically, it will copy and translate the simulation box for each frame *i* in each of the directions
  given by **extension_directions**. For example if **extension_directions=[[1,0,0],[0,1,0],[1,1,0]]**
  it will extend the simulation in the x, y and x+y directions.
  It will also reset the cell size for each frame by multiplying it by **cell_size_mult_factors**.
  Specifically, the unit cell is specified by the length of the 3 unit cell vectors and 3 angles.
  So for each frame *i*, the size of each of the unit cell vectors will be multiplied by the corresponding
  component of **cell_size_mult_factors[i]**.

  :param t: the trajectory
  :param extension_directions: The list of lists used to define the directions in which the 
    unit cell should be copied.
  :param cell_size_mult_factors:

  :type t: :class:`~ost.mol.CoordGroupHandle`
  :type extension_directions: :class:`list` (:class:`list`)
  :type cell_size_mult_factors: :class:`~ost.geom.Vec3`
  """
  vecs_to_neighbor_ucells_list=[]
  for i in range(t.GetFrameCount()):
    ucell_vec=t.GetFrame(i).GetCellVectors()
    ext_vecs=[geom.Vec3(el[0]*ucell_vec[0]+el[1]*ucell_vec[1]+el[2]*ucell_vec[2]) for el in extension_directions]
    vecs_to_neighbor_ucells_list.append(geom.Vec3List(ext_vecs))
  return ExtendTrajectory(t,vecs_to_neighbor_ucells_list,cell_size_mult_factors)

def ExtendTrajectory(t,vecs_to_neighbor_ucells_list,cell_size_mult_factors=(1,1,1)):
  """
  This function is used to extend a trajectory to its neighboring unit cells.
  Specifically, it will copy and translate the simulation box for each frame *i* by each
  vector in **vecs_to_neighbor_ucells_list[i]**. 
  It will also reset the cell size for each frame by multiplying it by **cell_size_mult_factors**.
  Specifically, the unit cell is specified by the length of the 3 unit cell vectors and 3 angles.
  So for each frame *i*, the size of each of the unit cell vectors will be multiplied by the corresponding
  component of **cell_size_mult_factors[i]**.

  :param t: the trajectory
  :param vecs_to_neighbor_ucells_list: The list of :class:`Vec3List` used to translate the 
    unit cell. This list should contain one :class:`Vec3List` for each frame in the trajectory.
  :param cell_size_mult_factors:

  :type t: :class:`~ost.mol.CoordGroupHandle`
  :type vecs_to_neighbor_ucells_list: :class:`list` (:class:`~ost.geom.Vec3List`)
  :type cell_size_mult_factors: :class:`~ost.geom.Vec3`
  """
  import pbc_utilities
  if not len(vecs_to_neighbor_ucells_list)==t.GetFrameCount():
    print 'you should provide a list of vectors for each frame'
  nreplicas=len(vecs_to_neighbor_ucells_list[0])
  if not all([len(vl)==nreplicas for vl in vecs_to_neighbor_ucells_list]):
    print 'Needs the same number of vectors for each frame'
  cell_size_mult_factors=geom.Vec3(*cell_size_mult_factors)
  eh=t.GetEntity()
  extended_eh=pbc_utilities.ExtendEntityToNeighboringUnitCells(eh,vecs_to_neighbor_ucells_list[0])
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
    cell_size=t.GetFrame(i).GetCellSize()
    new_cell_size=geom.CompMultiply(cell_size_mult_factors,cell_size)
    extended_t.AddFrame(vl,new_cell_size,t.GetFrame(i).GetCellAngles())
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

  
def WrapFloatList(fl,center,period):
  for i in range(len(fl)):
    fl[i]=fl[i]-period*round((fl[i]-center)/period)
  return fl
  
def AverageMembraneThickness(t,upper,lower,step_size=2,z_midplane=None):
  """
  This function creates a discrete mesh in the xy plane and calculates the average
  z position of the atoms in the upper and lower EntityViews over the trajectory t.
  Typically upper and lower contain the phosphate atoms of the lipids in the upper and
  lower leaflets of the membrane.
  It returns two Entities, representing these average positions in the upper and lower leaflets,
  as well as two matricies containing the same information.
  Each atom of the output entities is assigned a FloatProp 'thickness', 
  representing the local thickness of the membrane, another FloatProp ''
  'monolayer_thickness' representing the local thickness of the monolayer,
  an FloatProp 'count' which is the number times that particular bin was sampled and
  a FloatProp 'count_thickness' which is the min between the counts for this
  bin and the corresponding bin on the other leaflet.
  Relevant atoms for which the bin was sampled at least i times 
  can therefore be selected with eh_upper.Select('gacount>=i')
  whereas relevant atoms for which the thickness was sampled at least i times 
  can therefore be selected with eh_upper.Select('gacount_thickness>=i')
  Input:   -t          : CoordGroup, the trajectory
           -upper      : EntityView, the upper membrane leaflet
           -lower      : EntityView, the lower membrane leaflet 
           -step_size  : size of the bins in xy
           -z_midplane : average height of the midplane
           
  Output:  -eh_upper   : EntityHandle, average membrane upper leaflet
           -eh_lower   : EntityHandle, average membrane lower leaflet
           -zu         : numpy matrix, binned average z positions upper leaflet
           -zl         : numpy matrix, binned average z positions lower leaflet
  """
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
  zu=npy.where(npy.isnan(zu),0.0,zu)
  zl=npy.where(npy.isnan(zl),0.0,zl)
  thickness=zu-zl
  if not z_midplane:z_midplane=0.5*(npy.average(zu,weights=cu)+npy.average(zl,weights=cl))
  thickness_upper=zu-z_midplane
  thickness_lower=z_midplane-zl
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
      eh_upper.atoms[n1*i+j].SetFloatProp('count_thickness',min(cu[i,j],cl[i,j]))
      eh_lower.atoms[n1*i+j].SetFloatProp('count_thickness',min(cu[i,j],cl[i,j]))
      eh_upper.atoms[n1*i+j].SetFloatProp('monolayer_thickness',thickness_upper[i,j])
      eh_lower.atoms[n1*i+j].SetFloatProp('monolayer_thickness',thickness_lower[i,j])
  return (eh_upper,eh_lower,zu,zl)

def _gaussxy(p,refp,sigma):
  return npy.exp(geom.Length(geom.Vec2(p.x,p.y)-geom.Vec2(refp.x,refp.y))**2.0/(2.0*sigma))

def SmoothAverageMembrane(eh_upper,eh_lower,zu,zl,smooth_sigma=2.0,smooth_sampling_cutoff=5,z_midplane=None):
  """
  This function smoothes a membane surface obtained from the
  AverageMembraneThickness function.
  Input:   -eh_upper   : EntityHandle, average membrane upper leaflet
           -eh_lower   : EntityHandle, average membrane lower leaflet
           -zu         : numpy matrix, binned average z positions upper leaflet
           -zl         : numpy matrix, binned average z positions lower leaflet
           -smooth_sigma: standard deviation of the gaussian used for smoothing
           -smooth_sampling_cutoff: atoms with FloatProp 'count' smaller than this
                                    do not participate in smoothing of their neighboring atoms
           -z_midplane : average height of the midplane
                      
  Output:  -eh_upper_smooth : EntityHandle, smoothed average membrane upper leaflet
           -eh_lower_smooth : EntityHandle, smoothed average membrane lower leaflet
  """
  if z_midplane==None:
    z1=npy.average(npy.array([a.pos.z for a in eh_upper.atoms]),weights=npy.array([a.GetFloatProp('count') for a in eh_upper.atoms]))
    z2=npy.average(npy.array([a.pos.z for a in eh_lower.atoms]),weights=npy.array([a.GetFloatProp('count') for a in eh_lower.atoms]))
    z_midplane=0.5*(z1+z2)
  smooth_cutoff=2.0*smooth_sigma
  eh_upper_smooth=eh_upper.Copy()
  eh_lower_smooth=eh_lower.Copy()
  edi_upper=eh_upper_smooth.EditXCS(mol.BUFFERED_EDIT)
  edi_lower=eh_lower_smooth.EditXCS(mol.BUFFERED_EDIT)
  n1=zu.shape[0]
  vu=eh_upper.Select('gacount>=1')
  vl=eh_lower.Select('gacount>=1')
  vua=mol.AtomHandleList([a.handle for a in vu.atoms])
  vla=mol.AtomHandleList([a.handle for a in vl.atoms])
  vus=eh_upper.Select('gacount>={0}'.format(smooth_sampling_cutoff))
  vls=eh_lower.Select('gacount>={0}'.format(smooth_sampling_cutoff))
  for i in range(zu.shape[0]):
    for j in range(zu.shape[1]):
      au=eh_upper.atoms[n1*i+j]
      al=eh_lower.atoms[n1*i+j]
      if au in vua:
        within=mol.AtomHandleList([a.handle for a in vus.FindWithin(au.pos,smooth_cutoff)])
        if not au.handle in within:within.append(au)
        pu=npy.average([a.pos.z for a in within],weights=[_gaussxy(a.pos,au.pos,smooth_sigma) for a in within])
      else:pu=npy.nan
      if al in vla:
        within=mol.AtomHandleList([a.handle for a in vls.FindWithin(al.pos,smooth_cutoff)])
        if not al.handle in within:within.append(al)
        pl=npy.average([a.pos.z for a in within],weights=[_gaussxy(a.pos,al.pos,smooth_sigma) for a in within])
      else:pl=npy.nan
      if not (pu==npy.nan or pl==npy.nan):thickness=pu-pl
      else:thickness=npy.nan
      au2=eh_upper_smooth.atoms[n1*i+j]
      al2=eh_lower_smooth.atoms[n1*i+j]
      edi_upper.SetAtomPos(au2,geom.Vec3(au2.pos.x,au2.pos.y,pu))
      edi_lower.SetAtomPos(al2,geom.Vec3(al2.pos.x,al2.pos.y,pl))
      au2.SetFloatProp('thickness',thickness)
      al2.SetFloatProp('thickness',thickness)
      au2.SetFloatProp('monolayer_thickness',au2.pos.z-z_midplane)
      al2.SetFloatProp('monolayer_thickness',z_midplane-al2.pos.z)
  return (eh_upper_smooth,eh_lower_smooth)

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
    print "Generating water density"
    water_den=density_alg.CreateDensityMapFromTrajectoryWithPBC(t,waters,cell_center=cell_center,cell_size=cell_size,stride=stride)
    print "Generating lipid density"
    lipid_den=density_alg.CreateDensityMapFromTrajectoryWithPBC(t,lipids,cell_center=cell_center,cell_size=cell_size,stride=stride)
  else:
    print "Generating water density"
    water_den=density_alg.CreateDensityMapFromTrajectory(t,waters,stride=stride)
    print "Generating lipid density"
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
  print "Calculating the interface between the two densities"
  boundary=density_alg.GetBoundaryBetweenDensities(water_filtered,lipid_filtered,density_cutoff,PBC,cell_center,cell_size)
  return (water_filtered,lipid_filtered,boundary)

  


