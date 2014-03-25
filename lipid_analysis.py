"""
Module written by Niklaus Johner (niklaus.johner@a3.epfl.ch) 01.2013
This module contains functions to determine lipid tilt and splay angles and
Calculate the elastic properties of the membrane from there
"""
try:
  from ost import *
  import time
  import numpy as npy
  import os,math
  import entity_alg,trajectory_utilities,surface_alg,file_utilities
except:
  print 'could not import at least one of the modules nedded: ost, time, numpy, os, math, entity_alg,trajectory_utilities,surface_alg,file_utilities'


__all__=('AnalyzeSplayAngle','AssignNormalsFromSurfaceToResidues','CalculateTilts',\
        'GetBoundaryBetweenViews','AssignNormalsToLipids','AnalyzeLipidTilts',\
        'AnalyzeLipidSplay','AnalyzeLipidTiltAndSplayGeneral','AnalyzeLipidTiltAndSplayBilayer',\
        'AnalyzeLipidTiltAndSplayFixedNormal','WriteTiltDict','WriteSplayDict','OrientNormalsAlongVector',\
        'FitSplayOrTiltDistribution','AnalyzeAreaPerLipid','FitCompressionDistribution',\
        'AnalyzeLipidCompression','AnalyzeLipidCompression2','WriteCompressionDict')

def AnalyzeSplayAngle(t,v11,v12,v21,v22,v1d,v2d,n1,n2,an1,an2,distance_cutoff,tilt_angle_cutoff,dist_prot1=None,dist_prot2=None):
  h1=mol.alg.AnalyzeCenterOfMassPos(t,v11)
  h2=mol.alg.AnalyzeCenterOfMassPos(t,v21)
  t1=mol.alg.AnalyzeCenterOfMassPos(t,v12)
  t2=mol.alg.AnalyzeCenterOfMassPos(t,v22)
  v1=h1-t1
  v2=h2-t2
  cmd=npy.array(mol.alg.AnalyzeDistanceBetwCenterOfMass(t,v1d,v2d))
  if dist_prot1 and dist_prot2:
    prot_dist=[min(dp1,dp2) for dp1,dp2 in zip(dist_prot1,dist_prot2)]
    return [[di,geom.Angle(v1i,geom.AxisRotation(geom.Cross(n2i,n1i),geom.Angle(n1i,n2i))*v2i),an1i,an2i,pd] 
            for v1i,v2i,di,an1i,an2i,n1i,n2i,pd in zip(v1,v2,cmd,an1,an2,n1,n2,prot_dist) if (di<distance_cutoff 
            and (an1i<tilt_angle_cutoff or an2i<tilt_angle_cutoff))]
  else:return [[di,geom.Angle(v1i,geom.AxisRotation(geom.Cross(n2i,n1i),geom.Angle(n1i,n2i))*v2i),an1i,an2i] 
              for v1i,v2i,di,an1i,an2i,n1i,n2i in zip(v1,v2,cmd,an1,an2,n1,n2) if (di<distance_cutoff 
              and (an1i<tilt_angle_cutoff or an2i<tilt_angle_cutoff))]
  
  
CalculateInterfaceFromTraj=trajectory_utilities.CalculateInterfaceFromTraj

def AssignNormalsFromSurfaceToResidues(t,sele,surface,within_size=10):
  surface=surface.Select('')#Make sure we don't have to make that selection each time
  ln=[geom.Vec3List() for r in sele.residues]
  res_view_list=[r.Select('') for r in sele.residues]
  for i in range(t.GetFrameCount()):
    t.CopyFrame(i)
    for j,r in enumerate(res_view_list):
      within=surface.FindWithin(r.center_of_atoms,within_size)
      if not len(within)==0:a=entity_alg.FindClosestAtom(r,mol.CreateViewFromAtoms(within))
      else:a=entity_alg.FindClosestAtom(r,surface)
      ln[j].append(geom.Vec3(a.GetFloatProp('nx'),a.GetFloatProp('ny'),a.GetFloatProp('nz')))
  return ln

def CalculateTilts(t,lipids,normals,head_sele,tail_sele,prot_cm=None,bool_prop=''):
  tilts=[]
  prot_dist=[]
  for i,r in enumerate(lipids.residues):
    if bool_prop and not r.GetBoolProp(bool_prop):continue
    p1=mol.alg.AnalyzeCenterOfMassPos(t,r.Select(head_sele))
    p2=mol.alg.AnalyzeCenterOfMassPos(t,r.Select(tail_sele))
    n=normals[i]
    tilts.append([geom.Angle(ni,p1i-p2i) for ni,p1i,p2i in zip(n,p1,p2)])
    if prot_cm:
      prot_dist.append([min(geom.Length(geom.Vec2(cm-el1)),geom.Length(geom.Vec2(cm-el2))) for cm,el1,el2 in zip(prot_cm,p1,p2)])
  return (tilts,prot_dist)

def GetBoundaryBetweenViews(t,lipid_names,water_name='TIP3',outdir='',PBC=False,cell_center=None,cell_size=None,
                            den_cutoff=None,density_stride=1,within_size_normals=5.0,filename_basis=''):
  import time
  t0=time.time()
  eh=t.GetEntity()
  t2=t.Filter(eh.Select(''),stride=density_stride)
  #first we build the interafce between water and membrane
  lipid_sele='rname='+','.join(lipid_names)
  water_sele='rname='+water_name
  lipids=eh.Select(lipid_sele)
  waters=eh.Select(water_sele)
  print 'Calculating densities on',t2.GetFrameCount(),'frames'
  (water_filtered,lipid_filtered,boundary_filtered)=CalculateInterfaceFromTraj(t2,waters,lipids,PBC,cell_center,cell_size,10,den_cutoff,outdir,filename_basis)
  b_eh=entity_alg.CreateEntityFromVec3List(boundary_filtered)
  print 'calculating the normals',time.time()-t0
  surface_alg.CalculateNormals(b_eh,within_size_normals,PBC,cell_center,cell_size)
  surface_alg.OrientNormalsAlongDensity(b_eh,lipid_filtered)
  if outdir:
    io.SavePDB(b_eh,os.path.join(outdir,filename_basis+'boundary.pdb'))
    entity_alg.WriteFloatPropList(b_eh,['nx','ny','nz'],os.path.join(outdir,filename_basis+'boundary_normals.txt'),index=True)
  return (water_filtered,lipid_filtered,b_eh)
      
def AssignNormalsToLipids(t,eh,b_eh,lipid_names,head_group_dict):
  import time
  t0=time.time()
  #Now we assign a normal for each lipid in each frame
  print 'assigning normals to lipids',lipid_names
  lipid_normal_dict={}
  for lipid_name in lipid_names:
    print 'rname='+lipid_name+' and '+head_group_dict[lipid_name]
    t0=time.time()
    lipids=eh.Select('rname='+lipid_name+' and '+head_group_dict[lipid_name])
    print lipids.GetResidueCount()
    lipid_normal_dict[lipid_name]=AssignNormalsFromSurfaceToResidues(t,lipids,b_eh)
    print 'done for',lipid_name,'in',time.time()-t0,'seconds'
  return lipid_normal_dict

def AnalyzeLipidTilts(t,eh,lipid_names,lipid_normal_dict,head_group_dict,tail_dict,prot_cm=None,bool_prop=''):
  lipid_tilt_dict={}
  for ln in lipid_names:
    lipids=eh.Select('rname='+ln)
    lipid_tilt_dict[ln]=CalculateTilts(t,lipids,lipid_normal_dict[ln],head_group_dict[ln],tail_dict[ln],prot_cm,bool_prop)
  return lipid_tilt_dict
  
def AnalyzeLipidSplay(t,eh,lipid_names,head_group_dict,tail_dict,lipid_normal_dict,lipid_tilt_dict,distance_sele_dict={},distance_cutoff=10,angle_cutoff=0.175,bool_prop=''):
  import time
  t0=time.time()
  splay_dict={}
  if not distance_sele_dict:
    for el in lipid_names:distance_sele_dict[el]=''
  #first we calculate the splays for pairs of identical lipids
  print bool_prop,lipid_names
  for lipid_name in lipid_names:
    lipids=eh.Select('rname='+lipid_name)
    head_group_sele=head_group_dict[lipid_name]
    tail_sele=tail_dict[lipid_name]
    distance_sele=distance_sele_dict[lipid_name]
    al=[]
    t1l=lipid_tilt_dict[lipid_name][0]
    pd1=lipid_tilt_dict[lipid_name][1]
    n1l=lipid_normal_dict[lipid_name]
    head_group_view_list=[r.Select(head_group_sele) for r in lipids.residues]
    tail_view_list=[r.Select(tail_sele) for r in lipids.residues]
    distance_view_list=[r.Select(distance_sele) for r in lipids.residues]
    for i,r1 in enumerate(lipids.residues):
      v11=head_group_view_list[i]
      v12=tail_view_list[i]
      v13=distance_view_list[i]
      for j,r2 in enumerate(lipids.residues):
        if j>=i:continue
        if bool_prop and not (r1.GetBoolProp(bool_prop) or r2.GetBoolProp(bool_prop)):continue
        v21=head_group_view_list[j]
        v22=tail_view_list[j]
        v23=distance_view_list[j]
        if pd1:al.extend(AnalyzeSplayAngle(t,v11,v12,v21,v22,v13,v23,n1l[i],n1l[j],t1l[i],t1l[j],distance_cutoff,angle_cutoff,pd1[i],pd1[j]))
        else:al.extend(AnalyzeSplayAngle(t,v11,v12,v21,v22,v13,v23,n1l[i],n1l[j],t1l[i],t1l[j],distance_cutoff,angle_cutoff))
    splay_dict[lipid_name+'-'+lipid_name]=al
    print 'done for',lipid_name,lipid_name,time.time()-t0,len(al)
  
  #Now the splay for pairs of different lipids
  for l1,lipid_name1 in enumerate(lipid_names):
    lipids1=eh.Select('rname='+lipid_name1)
    head_group_sele1=head_group_dict[lipid_name1]
    tail_sele1=tail_dict[lipid_name1]
    distance_sele1=distance_sele_dict[lipid_name1]
    t1l=lipid_tilt_dict[lipid_name1][0]
    pd1=lipid_tilt_dict[lipid_name1][1]
    n1l=lipid_normal_dict[lipid_name1]
    head_group_view_list1=[r.Select(head_group_sele1) for r in lipids1.residues]
    tail_view_list1=[r.Select(tail_sele1) for r in lipids1.residues]
    distance_view_list1=[r.Select(distance_sele1) for r in lipids1.residues]
    for l2,lipid_name2 in enumerate(lipid_names):
      if l1>=l2:continue
      lipids2=eh.Select('rname='+lipid_name2)
      head_group_sele2=head_group_dict[lipid_name2]
      tail_sele2=tail_dict[lipid_name2]
      distance_sele2=distance_sele_dict[lipid_name2]    
      t2l=lipid_tilt_dict[lipid_name2][0]
      pd2=lipid_tilt_dict[lipid_name2][1]
      n2l=lipid_normal_dict[lipid_name2]
      head_group_view_list2=[r.Select(head_group_sele2) for r in lipids2.residues]
      tail_view_list2=[r.Select(tail_sele2) for r in lipids2.residues]
      distance_view_list2=[r.Select(distance_sele2) for r in lipids2.residues]
      al=[]
      for i,r1 in enumerate(lipids1.residues):
        if bool_prop and not r1.GetBoolProp(bool_prop):continue
        v11=head_group_view_list1[i]
        v12=tail_view_list1[i]
        v13=distance_view_list1[i]
        #v11=r1.handle.Select(head_group_sele1)
        #v12=r1.handle.Select(tail_sele1)
        #v13=r1.handle.Select(distance_sele1)
        for j,r2 in enumerate(lipids2.residues):
          #if j>=i:continue
          if bool_prop and not (r1.GetBoolProp(bool_prop) or r2.GetBoolProp(bool_prop)):continue
          v21=head_group_view_list2[j]
          v22=tail_view_list2[j]
          v23=distance_view_list2[j]
          #v21=r2.handle.Select(head_group_sele2)
          #v22=r2.handle.Select(tail_sele2)
          #v23=r2.handle.Select(distance_sele2)
          if pd1 and pd2:al.extend(AnalyzeSplayAngle(t,v11,v12,v21,v22,v13,v23,n1l[i],n2l[j],t1l[i],t2l[j],distance_cutoff,angle_cutoff,pd1[i],pd2[j]))
          else:al.extend(AnalyzeSplayAngle(t,v11,v12,v21,v22,v13,v23,n1l[i],n2l[j],t1l[i],t2l[j],distance_cutoff,angle_cutoff))
      splay_dict[lipid_name1+'-'+lipid_name2]=al
      print 'done for',lipid_name1,lipid_name2,time.time()-t0
  print 'done in',time.time()-t0
  return splay_dict

  
def AnalyzeLipidTiltAndSplayGeneral(t,lipid_names,head_group_dict,tail_dict,distance_cutoff=10.0,angle_cutoff=0.175,within_size_normals=5.0
                      ,distance_sele_dict={},water_name='TIP3',outdir='',PBC=False,cell_center=None,cell_size=None,den_cutoff=None
                      ,prot_sele=None,density_stride=10,tilt_bool_prop='',splay_bool_prop='',filename_basis='',sele_dict={}):
  import time
  t0=time.time()
  eh=t.GetEntity()
  t.CopyFrame(0)
  #first we build the interafce between water and membrane
  #and assign a normal for each lipid in each frame
  print 'Generating the water and lipid densities and boundary surface'
  (water_filtered,lipid_filtered,b_eh)=GetBoundaryBetweenViews(t,lipid_names,water_name,outdir,PBC,cell_center,cell_size,
                                                              den_cutoff,density_stride,within_size_normals,filename_basis)
  #Protein center of mass
  if prot_sele:prot_cm=mol.alg.AnalyzeCenterOfMassPos(t,eh.Select(prot_sele))
  else:prot_cm=None
  #From here we separate into the sub categories defined in sele_dict, which could be upper and lower leaflets
  s=['(rname='+ln+' and '+head_group_dict[ln]+')' for ln in lipid_names]
  lipid_sele_dict={}
  if sele_dict:
    for sele_name in sele_dict:
      lipid_sele_dict[sele_name]=' or '.join(s)+' and '+sele_dict[sele_name]
  else: 
    sele_dict['all']=''
    lipid_sele_dict['all']=' or '.join(s)
  print sele_dict,lipid_sele_dict
  lipid_normal_dict={}
  lipid_tilt_dict={}
  if tilt_bool_prop:lipid_tilt_dict_sele={}
  splay_dict={}
  for sele_name in sele_dict:
    lipid_sele=lipid_sele_dict[sele_name]
    sele=sele_dict[sele_name]
    print 'Assigning normals for',sele_name,sele,lipid_sele
    lipid_normal_dict[sele_name]=AssignNormalsToLipids(t,eh.Select(lipid_sele,mol.MATCH_RESIDUES),b_eh.Select(sele),lipid_names,head_group_dict)
    print 'Done in',time.time()-t0,'seconds'
    t0=time.time()
    print 'calculating tilts for',sele
    lipid_tilt_dict[sele_name]=AnalyzeLipidTilts(t,eh.Select(lipid_sele,mol.MATCH_RESIDUES),lipid_names,lipid_normal_dict[sele_name],
                                                  head_group_dict,tail_dict,prot_cm)
    if tilt_bool_prop:lipid_tilt_dict_sele[sele_name]=AnalyzeLipidTilts(t,eh.Select(lipid_sele,mol.MATCH_RESIDUES),lipid_names,lipid_normal_dict[sele_name],
                                                  head_group_dict,tail_dict,prot_cm,tilt_bool_prop)
    print 'Done in',time.time()-t0,'seconds'
    t0=time.time()
    print 'calculating splay for',sele
    splay_dict[sele_name]=AnalyzeLipidSplay(t,eh.Select(lipid_sele,mol.MATCH_RESIDUES),lipid_names,head_group_dict,tail_dict,
                                            lipid_normal_dict[sele_name],lipid_tilt_dict[sele_name],distance_sele_dict,distance_cutoff,angle_cutoff,splay_bool_prop)
    print 'Done in',time.time()-t0,'seconds'
    if outdir:
      try:WriteSplayDict(splay_dict[sele_name],outdir,filename_basis+sele_name+'_')
      except:print 'could not write dict'
      if tilt_bool_prop:WriteTiltDict(lipid_tilt_dict_sele[sele_name],outdir,filename_basis+sele_name+'_')
      else:WriteTiltDict(lipid_tilt_dict[sele_name],outdir,filename_basis+sele_name+'_')
  if tilt_bool_prop: return (lipid_tilt_dict_sele,lipid_normal_dict,splay_dict,b_eh)
  else: return (lipid_tilt_dict,lipid_normal_dict,splay_dict,b_eh)



def AnalyzeLipidTiltAndSplayBilayer(t,lipid_names,head_group_dict,tail_dict,distance_cutoff=10.0,angle_cutoff=0.175,within_size_normals=5.0,distance_sele_dict={},outdir=''
                      ,PBC=False,cell_center=None,cell_size=None,prot_sele=None,tilt_bool_prop='',splay_bool_prop='',filename_basis='',sele_dict={'upper':'z>0','lower':'z<0'}):
  t0=time.time()
  eh=t.GetEntity()
  t.CopyFrame(0)
  #Protein center of mass
  if prot_sele:prot_cm=mol.alg.AnalyzeCenterOfMassPos(t,eh.Select(prot_sele))
  else:prot_cm=None
  #From here on we do upper and lower leaflets separately
  s=['(rname='+ln+' and '+head_group_dict[ln]+')' for ln in lipid_names]
  lipid_sele_dict={}
  if sele_dict:
    for sele_name in sele_dict:
      lipid_sele_dict[sele_name]=' or '.join(s)+' and '+sele_dict[sele_name]
  else: 
    sele_dict['all']=''
    lipid_sele_dict['all']=' or '.join(s)
  #We determine the membrane interface from a 2D binning of the position of headgroups
  print 'Building average membrane surface'
  (eh_l,eh_u)=trajectory_utilities.AverageMembraneThickness(t,eh.Select(lipid_sele_dict['upper']),eh.Select(lipid_sele_dict['lower']),2.)
  surface_alg.CalculateNormals(eh_u,within_size_normals,PBC,cell_center,cell_size)
  surface_alg.CalculateNormals(eh_l,within_size_normals,PBC,cell_center,cell_size)
  OrientNormalsAlongVector(eh_u,geom.Vec3(0,0,1))
  OrientNormalsAlongVector(eh_l,geom.Vec3(0,0,-1))
  v=eh_l.CreateFullView()
  v.AddChain(eh_u.chains[0],mol.ViewAddFlag.INCLUDE_ALL)
  b_eh=mol.CreateEntityFromView(v,1)
  if outdir:io.SavePDB(b_eh,os.path.join(outdir,filename_basis+'boundary.pdb'),dialect='CHARMM')
  #normals
  #surface_alg.CalculateNormals(b_eh,10,PBC,cell_center,cell_size)
  print 'done with boundary'
  #sele_dict={'upper':'z>0','lower':'z<0'}
  lipid_normal_dict={}
  lipid_tilt_dict={}
  splay_dict={}
  for sele_name in sele_dict:
    lipid_sele=lipid_sele_dict[sele_name]
    sele=sele_dict[sele_name]
    print 'Assigning normals for',sele
    lipid_normal_dict[sele_name]=AssignNormalsToLipids(t,eh.Select(lipid_sele,mol.MATCH_RESIDUES),b_eh.Select(sele),lipid_names,head_group_dict)
    print 'Done in',time.time()-t0,'seconds'
    t0=time.time()
    print 'calculating tilts for',sele
    lipid_tilt_dict[sele_name]=AnalyzeLipidTilts(t,eh.Select(lipid_sele,mol.MATCH_RESIDUES),lipid_names,lipid_normal_dict[sele_name],
                                                  head_group_dict,tail_dict,prot_cm,tilt_bool_prop)
    print 'Done in',time.time()-t0,'seconds'
    t0=time.time()
    print 'calculating splay for',sele
    splay_dict[sele_name]=AnalyzeLipidSplay(t,eh.Select(lipid_sele,mol.MATCH_RESIDUES),lipid_names,head_group_dict,tail_dict,lipid_normal_dict[sele_name],
                                            lipid_tilt_dict[sele_name],distance_sele_dict,distance_cutoff,splay_bool_prop)
    print 'Done in',time.time()-t0,'seconds'
    if outdir:
      WriteTiltDict(lipid_tilt_dict[sele_name],outdir,filename_basis+sele_name+'_')
      WriteSplayDict(splay_dict[sele_name],outdir,filename_basis+sele_name+'_')
  return (lipid_tilt_dict,lipid_normal_dict,splay_dict,b_eh)

def AnalyzeLipidTiltAndSplayFixedNormal(t,lipid_names,head_group_dict,tail_dict,sele_dict,sele_normal_dict,distance_cutoff=10.0,angle_cutoff=0.175
                      ,distance_sele_dict={},outdir='',PBC=False,cell_center=None,cell_size=None,prot_sele=None,tilt_bool_prop='',splay_bool_prop='',filename_basis=''):
  import time
  t0=time.time()
  eh=t.GetEntity()
  t.CopyFrame(0)
  nframes=t.GetFrameCount()
  #Protein center of mass
  if prot_sele:prot_cm=mol.alg.AnalyzeCenterOfMassPos(t,eh.Select(prot_sele))
  else:prot_cm=None
  #From here we separate into the sub categories defined in sele_dict, which could be upper and lower leaflets
  s=['(rname='+ln+' and '+head_group_dict[ln]+')' for ln in lipid_names]
  lipid_sele_dict={}
  for sele_name in sele_dict:
    lipid_sele_dict[sele_name]=' or '.join(s)+' and '+sele_dict[sele_name]
  print sele_dict,lipid_sele_dict,sele_normal_dict
  lipid_normal_dict={}
  lipid_tilt_dict={}
  if tilt_bool_prop:lipid_tilt_dict_sele={}
  splay_dict={}
  for sele_name in sele_dict:
    lipid_sele=lipid_sele_dict[sele_name]
    sele=sele_dict[sele_name]
    n=sele_normal_dict[sele_name]
    lipid_normal_dict[sele_name]={}
    print 'Assigning normal',lipid_normal_dict[sele_name],'for',sele_name,sele,lipid_sele
    for lipid_name in lipid_names:
      print 'rname='+lipid_name+' and '+head_group_dict[lipid_name]
      lipids=eh.Select('rname='+lipid_name+' and '+head_group_dict[lipid_name])
      lipid_normal_dict[sele_name][lipid_name]=[geom.Vec3List([n for i in range(nframes)]) for r in lipids.residues]
    t0=time.time()
    print 'calculating tilts for',sele
    lipid_tilt_dict[sele_name]=AnalyzeLipidTilts(t,eh.Select(lipid_sele,mol.MATCH_RESIDUES),lipid_names,lipid_normal_dict[sele_name],
                                                  head_group_dict,tail_dict,prot_cm)
    if tilt_bool_prop:lipid_tilt_dict_sele[sele_name]=AnalyzeLipidTilts(t,eh.Select(lipid_sele,mol.MATCH_RESIDUES),lipid_names,lipid_normal_dict[sele_name],
                                                  head_group_dict,tail_dict,prot_cm,tilt_bool_prop)
    print 'Done in',time.time()-t0,'seconds'
    t0=time.time()
    print 'calculating splay for',sele
    splay_dict[sele_name]=AnalyzeLipidSplay(t,eh.Select(lipid_sele,mol.MATCH_RESIDUES),lipid_names,head_group_dict,tail_dict,
                                            lipid_normal_dict[sele_name],lipid_tilt_dict[sele_name],distance_sele_dict,distance_cutoff,angle_cutoff,splay_bool_prop)
    print 'Done in',time.time()-t0,'seconds'
    if outdir:
      try:WriteSplayDict(splay_dict[sele_name],outdir,filename_basis+sele_name+'_')
      except:print 'could not write dict'
      if tilt_bool_prop:WriteTiltDict(lipid_tilt_dict_sele[sele_name],outdir,filename_basis+sele_name+'_')
      else:WriteTiltDict(lipid_tilt_dict[sele_name],outdir,filename_basis+sele_name+'_')
  if tilt_bool_prop: return (lipid_tilt_dict_sele,lipid_normal_dict,splay_dict,b_eh)
  else: return (lipid_tilt_dict,lipid_normal_dict,splay_dict)


def WriteTiltDict(lipid_tilt_dict,outdir,filename_basis=''):
  for ln in lipid_tilt_dict:
    tl=[180.*el/math.pi  for lt in lipid_tilt_dict[ln][0] for el in lt]
    dl=[el  for lt in lipid_tilt_dict[ln][1] for el in lt]
    if len(dl)!=0:
      if len(dl)!=len(tl):print 'not the same number of tilts and distances',len(tl),len(dl)
      ll=[tl,dl]
      titles=['tilt','dist from prot']
    else:
      ll=[tl]
      titles=['tilt']
    file_utilities.WriteListOfListsInColumns(titles,ll,os.path.join(outdir,filename_basis+'tilt_'+ln+'.txt'))

def WriteSplayDict(splay_dict,outdir,filename_basis):
  if len(splay_dict.keys())==0:return
  for ln in splay_dict:
    if len(splay_dict[ln][0])==4:prot_dist_flag=False
    else:prot_dist_flag=True
    if prot_dist_flag:
      sl=[[el[0],180.*el[1]/math.pi,180.*el[2]/math.pi,180.*el[3]/math.pi,el[4]] for el in splay_dict[ln]]
      file_utilities.WriteListOfListsInLines(['dist','splay','tilt1','tilt2','dist2prot'],sl,os.path.join(outdir,filename_basis+'splay_'+ln+'.txt'))
    else:
      sl=[[el[0],180.*el[1]/math.pi,180.*el[2]/math.pi,180.*el[3]/math.pi] for el in splay_dict[ln]]
      file_utilities.WriteListOfListsInLines(['dist','splay','tilt1','tilt2'],sl,os.path.join(outdir,filename_basis+'splay_'+ln+'.txt'))
  return

def OrientNormalsAlongVector(eh,v=geom.Vec3(0,0,1)):
  for a in eh.atoms:
    n=geom.Vec3(a.GetFloatProp('nx'),a.GetFloatProp('ny'),a.GetFloatProp('nz'))
    if geom.Dot(v,n)<0.0:
      a.SetFloatProp('nx',-n.x)
      a.SetFloatProp('ny',-n.y)
      a.SetFloatProp('nz',-n.z)
  return

def gauss(x, *p):
  A, mu, sigma = p
  return A*npy.exp(-(x-mu)**2/(2.*sigma**2))

def FindIndexOfClosestValue(l,v):
  return min(enumerate(l), key=lambda x: abs(x[1]-v))[0]
  
def FitSplayOrTiltDistribution(angles_list,bins=90,outdir='',filename_basis='',gfx_plot=False,xlabel=True,ylabel=True,title=True):
  from scipy.optimize import curve_fit
  import matplotlib as mpl
  if not gfx_plot:mpl.use('Agg')
  import matplotlib.pyplot as plt
  pa=npy.histogram(angles_list,bins,range=[0,90],density=True)
  bincenters=0.5*(pa[1][1:]+pa[1][:-1])
  #We fit a gaussian to the distributions to determine the fitting range for f(alpha)
  p0 = [1, 20.,20]
  (a1,mu1,sigma1), var_matrix = curve_fit(gauss, bincenters, pa[0], p0=p0)
  range_list=[[FindIndexOfClosestValue(bincenters,max(0,mu1-i*sigma1)),FindIndexOfClosestValue(bincenters,mu1+i*sigma1)] for i in [1,0.75,0.5,1.5,2]]
  #We calculate f(alpha)
  sin_bincenters=npy.sin(bincenters/180.*math.pi)
  fa=-npy.log((180./math.pi)*pa[0]/sin_bincenters)
  x=bincenters*bincenters
  res_list=[npy.polyfit(x[r[0]:r[1]],fa[r[0]:r[1]],1) for r in range_list]
  K_list=[r[0]*2*(180.*180.)/(math.pi*math.pi) for r in res_list]
  DeltaK=npy.std(K_list)
  K=K_list[0]
  #Now we make some plots and write out a file with the distribution
  if outdir:  
    file_utilities.WriteListOfListsInColumns(['bin','dist','fa'],[bincenters,pa[0],fa],os.path.join(outdir,filename_basis+'_distribution.txt'),separator=' ')
    #We plot the distribution and fitting range
    r1=range_list[0]
    plt.figure()
    plt.plot(bincenters,pa[0])
    plt.vlines([bincenters[r1[0]],bincenters[r1[1]]],0,a1,linestyle='--',color='b')
    if xlabel:plt.xlabel('Angle')
    if ylabel:plt.ylabel('Probability distribution')
    if title:plt.title(filename_basis+' Angle distributions')
    plt.show();plt.savefig(os.path.join(outdir,filename_basis+'_distribution.png'))
    #We plot fa
    plt.figure()
    plt.plot(bincenters,fa,'o',c='r')
    plt.plot(bincenters,[npy.polyval(res_list[0],xi) for xi in x],'--',c='r',label='K=%.1f $\pm$ %.1f' %(K,DeltaK))
    if xlabel:plt.xlabel('Angle')
    if ylabel:plt.ylabel(r'$-\log\left[\frac{P(\alpha)}{\sin(\alpha)}\right]$')
    if title:plt.title(filename_basis+' Distribution of angles')
    plt.legend(loc='best')
    plt.show()
    plt.savefig(os.path.join(outdir,filename_basis+'_fa_fit.png'))
  return (K,DeltaK,K_list)

def AnalyzeAreaPerLipid(t,lipids):
  n=lipids.GetResidueCount()
  Al=FloatList()
  for i in range(t.GetFrameCount()):
    f=t.GetFrame(i)
    cell_size=f.GetCellSize()
    Al.append(cell_size[0]*cell_size[1])
  return 2.*npy.average(Al)/float(n)

def FitCompressionDistribution(compression_list,d0,area_per_lipid,bins=100,outdir='',filename_basis='',gfx_plot=False):
  from scipy.optimize import curve_fit
  import matplotlib as mpl
  if not gfx_plot:mpl.use('Agg')
  import matplotlib.pyplot as plt
  dl2=[(el)*(el)/(d0*d0) for el in compression_list]
  print max(dl2)
  pa=npy.histogram(dl2,bins,range=[0,1],density=True)
  bincenters=0.5*(pa[1][1:]+pa[1][:-1])
  plt.figure()
  plt.plot(bincenters,pa[0])
  plt.show();plt.savefig(os.path.join(outdir,filename_basis+'_distribution.png'))
  #We fit a gaussian to the distributions to determine the fitting range for f(alpha)
  p0 = [10, 0.0,0.1]
  #(a1,mu1,sigma1), var_matrix = curve_fit(gauss, bincenters, pa[0], p0=p0)
  #range_list=[[FindIndexOfClosestValue(bincenters,max(0,mu1-i*sigma1)),FindIndexOfClosestValue(bincenters,mu1+i*sigma1)] for i in [1,0.75,0.5,1.5,2]]
  range_list=[[3,20]]
  a1=1
  print range_list
  #We calculate f(alpha)
  bincenters=0.5*(pa[1][1:]+pa[1][:-1])
  fa=-npy.log(pa[0])
  x=bincenters*bincenters
  res_list=[npy.polyfit(x[r[0]:r[1]],fa[r[0]:r[1]],1) for r in range_list]
  K_list=[r[0]/area_per_lipid for r in res_list]
  DeltaK=npy.std(K_list)
  K=K_list[0]
  #Now we make some plots and write out a file with the distribution
  if outdir:  
    file_utilities.WriteListOfListsInColumns(['bin','dist','fa'],[bincenters,pa[0],fa],os.path.join(outdir,filename_basis+'_distribution.txt'),separator=' ')
    #We plot the distribution and fitting range
    r1=range_list[0]
    plt.figure()
    plt.plot(bincenters,pa[0])
    plt.vlines([bincenters[r1[0]],bincenters[r1[1]]],0,a1,linestyle='--',color='b')
    plt.xlabel('Srtain');plt.ylabel('Probability density')
    plt.title(filename_basis+' Compression distributions')
    plt.show();plt.savefig(os.path.join(outdir,filename_basis+'_distribution.png'))
    #We plot fa
    plt.figure()
    plt.plot(bincenters,fa,'o',c='r')
    plt.plot(bincenters,[npy.polyval(res_list[0],xi) for xi in x],'--',c='r',label='K=%.1f $\pm$ %.1f' %(K,DeltaK))
    plt.xlabel('Angle');plt.ylabel(r'$-\log\left[P\left(\frac{\Delta L^2}{L}\right)\right]$'),plt.title(filename_basis+' Distribution of compressions')
    plt.legend(loc='best')
    plt.show()
    plt.savefig(os.path.join(outdir,filename_basis+'_fa_fit.png'))
  return (K,DeltaK,K_list)


def AnalyzeLipidCompression(t,upper,lower,lipid_names,head_group_dict,d0_dict={},prot_cm=None,outdir='',filename_basis=''):
  compression_dict={}
  for l1,lipid_name1 in enumerate(lipid_names):
    lipid_headgroups1=upper.Select('rname='+lipid_name1+' and '+head_group_dict[lipid_name1])
    head_group_view_list1=[r.Select('') for r in lipid_headgroups1.residues]
    print lipid_name1,len(head_group_view_list1)
    for l2,lipid_name2 in enumerate(lipid_names):
      lipid_headgroups2=lower.Select('rname='+lipid_name2+' and '+head_group_dict[lipid_name2])
      head_group_view_list2=[r.Select('') for r in lipid_headgroups2.residues]
      print lipid_name2,len(head_group_view_list2)
      if lipid_name1+'-'+lipid_name2 in d0_dict:d0=d0_dict[lipid_name1+'-'+lipid_name2]
      else:
        cm1=mol.alg.AnalyzeCenterOfMassPos(t,lipid_headgroups1)
        cm2=mol.alg.AnalyzeCenterOfMassPos(t,lipid_headgroups2)
        d0=npy.average([el[2] for el in cm1])-npy.average([el[2] for el in cm2])
        d0_dict[lipid_name1+'-'+lipid_name2]=d0
        print lipid_name1,lipid_name2,d0
      dl=[]
      for ru in head_group_view_list1:
        cm_au=mol.alg.AnalyzeCenterOfMassPos(t,ru)
        if prot_cm:
          dp1=[vp-vl for vp,vl in zip(prot_cm,cm_au)]
          dp1=[geom.Length(geom.Vec2(v[0],v[1])) for v in dp1]
        for rl in head_group_view_list2:
          cm_al=mol.alg.AnalyzeCenterOfMassPos(t,rl)
          if prot_cm:
            dp2=[vp-vl for vp,vl in zip(prot_cm,cm_al)]
            dp2=[geom.Length(geom.Vec2(v[0],v[1])) for v in dp2]
            dp=[min(dp1i,dp2i) for dp1i,dp2i in zip(dp1,dp2)]
          vl=[cm_aui-cm_ali for cm_ali,cm_aui in zip(cm_al,cm_au)]
          if prot_cm:dl.extend([[(vi[2]-d0),dpi] for vi,dpi in zip(vl,dp) if geom.Length(geom.Vec2(vi[0],vi[1]))<10.0])
          else:dl.extend([[(vi[2]-d0)] for vi in vl if geom.Length(geom.Vec2(vi[0],vi[1]))<10.0])
      if lipid_name1+'-'+lipid_name2 in compression_dict:
        print 'weird'
      else:compression_dict[lipid_name1+'-'+lipid_name2]=dl
  if outdir:WriteCompressionDict(compression_dict,outdir,filename_basis)
  if outdir:file_utilities.WriteListOfListsInColumns(['pairs','d0'],[d0_dict.keys(),d0_dict.values()],os.path.join(outdir,filename_basis+'d0.txt'))
  return (compression_dict,d0_dict)


def AnalyzeLipidCompression2(t,lipids,lipid_names,head_group_dict,tail_dict):
  lipid_length_dict={}
  for ln in lipid_names:lipid_length_dict[ln]=FloatList()
  for r in lipids.residues:
    ln=r.name
    cm1=mol.alg.AnalyzeCenterOfMassPos(t,r.Select(head_group_dict[ln]))
    cm2=mol.alg.AnalyzeCenterOfMassPos(t,r.Select(tail_dict[ln]))
    #l=[abs((v1-v2)[2]) for v1,v2 in zip(cm1,cm2)]
    l=[geom.Length(v1-v2) for v1,v2 in zip(cm1,cm2)]
    #lipid_length_dict[ln].extend(mol.alg.AnalyzeDistanceBetwCenterOfMass(t,r.Select(head_group_dict[ln]),r.Select(tail_dict[ln])))
    lipid_length_dict[ln].extend(FloatList(l))
  return lipid_length_dict

def WriteCompressionDict(compression_dict,outdir,filename_basis):
  if len(compression_dict.keys())==0:return
  for ln in compression_dict:
    if len(compression_dict[ln][0])==1:titles=['compression']
    else:titles=['compression','prot_dist']
    file_utilities.WriteListOfListsInLines(titles,compression_dict[ln],os.path.join(outdir,filename_basis+'compression_'+ln+'.txt'))
  return
