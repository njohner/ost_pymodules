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
  from scipy.optimize import curve_fit
  import matplotlib as mpl
  try:gui.dng
  except:mpl.use('Agg')
  import matplotlib.pyplot as plt
  import entity_alg,trajectory_utilities,surface_alg,file_utilities
except:
  print 'could not import at least one of the modules nedded: ost, time, numpy, os, math, entity_alg,trajectory_utilities,surface_alg,file_utilities'


__all__=('CalculateSplayAngle','AssignNormalsFromSurfaceToResidues','CalculateTilts',\
        'GetBoundaryBetweenViews','AssignNormalsToLipids','AnalyzeLipidTilts',\
        'AnalyzeLipidSplay','AnalyzeLipidTiltAndSplayGeneral',\
        'WriteTiltDict','WriteSplayDict','OrientNormalsAlongVector',\
        'FitTiltDistribution','FitSplayDistribution','AnalyzeAreaPerLipid')


def CalculateSplayAngle(v11,v12,v21,v22,v1p,v2p,n1,n2,distance_cutoff):
  if geom.Dot(n1,n2)<0:return None
  x=v2p.GetCenterOfMass()-v1p.GetCenterOfMass()
  d=geom.Length(x)
  x=x-geom.Dot(x,n1)*n1
  #d=geom.Length(x)
  if d>distance_cutoff:return None
  x=geom.Normalize(x)
  v1=geom.Normalize(v11.GetCenterOfMass()-v12.GetCenterOfMass())
  v2=geom.Normalize(v21.GetCenterOfMass()-v22.GetCenterOfMass())
  #x=geom.Normalize(v2p.GetCenterOfMass()-v1p.GetCenterOfMass())
  #d=geom.Length(v2p.GetCenterOfMass()-v1p.GetCenterOfMass())
  return [(geom.Dot(v2,x)-geom.Dot(v1,x)-geom.Dot(n2,x)+geom.Dot(n1,x))/d,d]
  

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
  return (npy.array(tilts),npy.array(prot_dist))

def GetBoundaryBetweenViews(t,lipid_names,water_name='TIP3',outdir='',PBC=False,cell_center=None,cell_size=None,
                            den_cutoff=None,density_stride=1,within_size_normals=5.0,filename_basis=''):
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
  

def AnalyzeLipidSplay(t,eh,lipid_names,head_group_dict,tail_dict,lipid_normal_dict,lipid_tilt_dict,distance_sele_dict,distance_cutoff=10,bool_prop='',patch=False):
  nframes=t.GetFrameCount()
  t0=time.time()
  splay_dict={}
  mean_dist_dict={}
  area_dict={}
  n_dict={}
  
  lipid_sele_string='rname='+','.join([rname for rname in distance_sele_dict.keys()])
  distance_sele_string=' or '.join(['(rname={0} and {1})'.format(rname,distance_sele_dict[rname]) for rname in distance_sele_dict.keys()])
  print 'selecting all the lipids using',lipid_sele_string
  print 'selection for distance calculation',distance_sele_string
  lipids=eh.Select(lipid_sele_string)
  distance_sele=lipids.Select(distance_sele_string)
  print bool_prop,lipid_names
  head_group_view_list_dict={}
  tail_view_list_dict={}
  distance_view_list_dict={}
  for ln in lipid_names:
    v=lipids.Select('rname='+ln)
    head_group_view_list_dict[ln]=[r.Select(head_group_dict[ln]) for r in v.residues]
    tail_view_list_dict[ln]=[r.Select(tail_dict[ln]) for r in v.residues]
    distance_view_list_dict[ln]=[r.Select(distance_sele_dict[ln]) for r in v.residues]
  splay_dict={}#We prepare the dictionary that will contain the splays for each type of lipid pair
  for ln1 in lipid_names:
    for ln2 in lipid_names:splay_dict[ln1+'-'+ln2]=[]
  for f in range(t.GetFrameCount()):
    if f%50==0:print "frame {0} in {1} seconds".format(f,time.time()-t0)
    t.CopyFrame(f)
    for r1 in lipids.residues:
      if not r1.GetBoolProp(bool_prop):continue
      ln1=r1.GetName()
      i=r1.GetIntProp('index')
      within=distance_sele.Select("{0}<>[cname={1} and rnum={2}]".format(distance_cutoff,r1.chain.name,r1.number.num))
      within=within.Select("not (cname={0} and rnum={1})".format(r1.chain.name,r1.number.num))
      #if within.GetResidueCount()<2:continue
      v11=head_group_view_list_dict[ln1][i]
      v12=tail_view_list_dict[ln1][i]
      v13=distance_view_list_dict[ln1][i]
      n1=lipid_normal_dict[ln1][i][f]
      a1=lipid_tilt_dict[ln1][0][i][f]
      for r2 in within.residues:
        ln2=r2.name
        j=r2.GetIntProp('index')
        v21=head_group_view_list_dict[ln2][j]
        v22=tail_view_list_dict[ln2][j]
        v23=distance_view_list_dict[ln2][j]
        n2=lipid_normal_dict[ln2][j][f]
        a2=lipid_tilt_dict[ln2][0][j][f]
        if a1>0.8 or a2>0.8:continue
        s=CalculateSplayAngle(v11,v12,v21,v22,v13,v23,n1,n2,distance_cutoff)
        if s:splay_dict[ln1+'-'+ln2].append(s)
  #We finalize the dictionary by merging equivalent entries
  nsplays=0
  splay_dict_def={}
  for ln1 in lipid_names:
    for ln2 in lipid_names:
      key=ln1+'-'+ln2
      key2=ln2+'-'+ln1
      nsplays+=len(splay_dict[key])
      if key2 in splay_dict_def.keys():
        splay_dict_def[key2].extend(splay_dict[key])
      else:splay_dict_def[key]=splay_dict[key]
  print 'done in {0} seconds. Computed {1} splays'.format(time.time()-t0,nsplays)
  return splay_dict_def

def AnalyzeLipidTiltAndSplayGeneral(t,lipid_names,head_group_dict,tail_dict,distance_cutoff=10.0,within_size_normals=5.0
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
  for ln in lipid_names:
    v=eh.Select('rname={0}'.format(ln))
    for i,r in enumerate(v.residues):r.SetIntProp('index',i)
  lipid_normal_dict={}
  lipid_tilt_dict={}
  if tilt_bool_prop:lipid_tilt_dict_sele={}
  splay_dict={}
  for sele_name in sele_dict:
    lipid_sele=lipid_sele_dict[sele_name]
    sele=sele_dict[sele_name]
    print 'Assigning normals for',sele_name,sele,lipid_sele
    lipid_normal_dict[sele_name]=AssignNormalsToLipids(t,eh.Select(lipid_sele,mol.MATCH_RESIDUES),b_eh.Select(sele),lipid_names,head_group_dict)
    print 'lipid normal dict',lipid_normal_dict[sele_name].keys()
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
                                            lipid_normal_dict[sele_name],lipid_tilt_dict[sele_name],distance_sele_dict,distance_cutoff,splay_bool_prop)
    print 'Done in',time.time()-t0,'seconds'
    if outdir:
      try:WriteSplayDict(splay_dict[sele_name],outdir,filename_basis+sele_name+'_')
      except:print 'could not write dict'
      if tilt_bool_prop:WriteTiltDict(lipid_tilt_dict_sele[sele_name],outdir,filename_basis+sele_name+'_')
      else:WriteTiltDict(lipid_tilt_dict[sele_name],outdir,filename_basis+sele_name+'_')
  if tilt_bool_prop: return (lipid_tilt_dict_sele,lipid_normal_dict,splay_dict,b_eh)
  else: return (lipid_tilt_dict,lipid_normal_dict,splay_dict,b_eh)

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
    if len(splay_dict[ln])==0:continue
    if len(splay_dict[ln][0])==3:prot_dist_flag=True
    else:prot_dist_flag=False
    if prot_dist_flag:
      sl=[[el[0],el[1],el[2]] for el in splay_dict[ln]]
      file_utilities.WriteListOfListsInLines(['splay','dist','ProtDist'],sl,os.path.join(outdir,filename_basis+'splay_'+ln+'.txt'))
    else:
      sl=[[el[0],el[1]] for el in splay_dict[ln]]
      file_utilities.WriteListOfListsInLines(['splay','dist'],sl,os.path.join(outdir,filename_basis+'splay_'+ln+'.txt'))
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

def parabole(x,*p):
  a,b,x0=p
  return a+b*(x-x0)**2.0
  
def centered_parabole(x,*p):
  a,b=p
  return a+b*x**2.0

def FindIndexOfClosestValue(l,v):
  return min(enumerate(l), key=lambda x: abs(x[1]-v))[0]


def _FitGaussian(bincenters,pa):
  mu0=npy.sum(bincenters*pa)/npy.sum(pa)
  A0=max(pa)
  sigma0=npy.sqrt(npy.sum(((bincenters-mu0)**2.0)*pa)/npy.sum(pa))
  (A,mu,sigma),v=curve_fit(gauss,bincenters,pa,[A0,mu0,sigma0])
  print A0,mu0,sigma0,A,mu,abs(sigma)
  return A,mu,abs(sigma)

def _PlotGaussian(bincenters,pa,A,mu,sigma,outfile,title='',xlabel=''):
  plt.figure()
  plt.plot(bincenters,pa,'o')
  plt.plot(bincenters,gauss(bincenters,A,mu,sigma),'-',color='g',label='$\mu={0}\ ;\ \sigma={1}$'.format(round(mu,int(1/mu)+1),round(sigma,int(1/sigma)+1)))
  plt.vlines(mu,0,plt.ylim()[1],linestyle='--',color='r')
  plt.vlines([mu+sigma,mu-sigma],0,plt.ylim()[1],linestyle='--',color='c')
  plt.title(title)
  plt.xlabel(xlabel)
  plt.ylabel('Probability')
  plt.show()
  plt.savefig(outfile)
  
def _FitParabole(bincenters,fa,range,centered=False):
  first=FindIndexOfClosestValue(bincenters,range[0])
  last=FindIndexOfClosestValue(bincenters,range[1])
  mask=fa!=npy.inf
  a=min(fa)
  x0=bincenters[npy.argmin(fa)]
  xm=bincenters[mask][npy.argmax(fa[mask])]
  fm=max(fa[mask])
  b=(fm-a)/(xm-x0)**2.0
  if centered:r,v=curve_fit(centered_parabole,bincenters[first:last],fa[first:last],[a,b])
  else:r,v=curve_fit(parabole,bincenters[first:last],fa[first:last],[a,b,x0])
  return r

def _PlotParabola(bincenters,fa,a,b,x0,fitting_range,outfile,title='',xlabel='',ylabel='',fit_label=''):
  plt.figure()
  plt.plot(bincenters,fa,'o')
  ymax=plt.ylim()[1]
  ymin=plt.ylim()[0]
  plt.plot(bincenters,[parabole(xi,a,b,x0) for xi in bincenters],'--',color='r',label=fit_label)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.title(title)
  if fit_label:plt.legend(loc='best')
  plt.vlines(fitting_range,ymin,ymax,linestyle='--',color='c')
  plt.ylim(ymin,ymax)
  plt.show()
  plt.savefig(outfile)
  return

def FitSplayDistribution(splay_list,lipid_area,factor=2,nbins=100,xrange=None,outdir='',filename_basis='',title_complement=''):  
  if not xrange:
    w=npy.std(splay_list)
    m=npy.average(splay_list)
    xrange=[m-3*w,m+3*w]
  pa=npy.histogram(splay_list,nbins,range=xrange,density=True)
  bincenters=0.5*(pa[1][1:]+pa[1][:-1])
  fa=-npy.log(pa[0])
  A,mu,sigma=_FitGaussian(bincenters,pa[0])
  #ranges=[(mu-i*sigma,mu+i*sigma) for i in [1,0.75,0.5,1.5,2]]
  ranges=[(mu-i*sigma,mu+i*sigma) for i in [1,1.25,1.5,1.75,2]]
  res_list=[]
  for range in ranges:
    try:r=_FitParabole(bincenters,fa,range)
    except:r=[0,0]
    res_list.append(r)
  K_list=[factor*r[1]/lipid_area for r in res_list]
  #K_list=[2.*r[1] for r in res_list]
  DeltaK=npy.std([el for el in K_list if not el==0.0])
  K=K_list[0]
  if outdir:
    file_utilities.WriteListOfListsInColumns(['bin','dist','fa'],[bincenters,pa[0],fa],os.path.join(outdir,'_'.join([filename_basis,'splay','distribution'])+'.txt'),separator=' ')
    title='Splay distribution {0}'.format(title_complement)
    outfile=os.path.join(outdir,'_'.join([filename_basis,'splay','distribution'])+'.png')
    _PlotGaussian(bincenters,pa[0],A,mu,sigma,outfile,title=title,xlabel='Splay')
    outfile=os.path.join(outdir,'_'.join([filename_basis,'splay','fit'])+'.png')
    a,b,x0=res_list[0]
    _PlotParabola(bincenters,fa,a,b,x0,ranges[0],outfile,title,'Splay',r'$-\log\left[P(\alpha)\right]$','$K={0}\pm {1}$'.format(round(K,int(-math.log10(K))+1),round(DeltaK,int(-math.log10(DeltaK))+1)))
  return K_list,DeltaK
  
def FitTiltDistribution(tilt_list,nbins=90,range=None,outdir='',filename_basis='',title_complement='',degrees=False):
  if degrees:ac=math.pi/180.
  else:ac=1
  if None:range=[0,math.pi/(2.*ac)]
  pa=npy.histogram(tilt_list,nbins,range=range,density=True)
  bincenters=0.5*(pa[1][1:]+pa[1][:-1])
  pa=pa[0]
  pa2=pa/npy.sin(bincenters*ac)
  fa=-npy.log(pa2)
  A,mu,sigma=_FitGaussian(bincenters,pa)
  #ranges=[(max(mu-i*sigma,0),mu+i*sigma) for i in [1,0.75,0.5,1.5,2]]
  ranges=[(max(mu-i*sigma,0),mu+i*sigma) for i in [1,1.25,1.5,1.75,2.0]]
  res_list=[_FitParabole(bincenters,fa,range,centered=True) for range in ranges]
  K_list=[2.*r[1]/(ac*ac) for r in res_list]
  DeltaK=npy.std(K_list)
  K=K_list[0]
  if outdir:
    file_utilities.WriteListOfListsInColumns(['bin','dist','fa'],[bincenters,pa,fa],os.path.join(outdir,'_'.join([filename_basis,'tilt','distribution'])+'.txt'),separator=' ')
    title='Tilt distribution {0}'.format(title_complement)
    outfile=os.path.join(outdir,'_'.join([filename_basis,'tilt','distribution'])+'.png')
    _PlotGaussian(bincenters,pa,A,mu,sigma,outfile,title=title,xlabel='tilt')
    outfile=os.path.join(outdir,'_'.join([filename_basis,'tilt','fit'])+'.png')
    a,b=res_list[0]
    r=[el for el in ranges[0]]
    _PlotParabola(bincenters,fa,a,b,0.0,r,outfile,title,'tilt',r'$-\log\left[\frac{P(\alpha)}{\sin(\alpha)}\right]$','$\chi={0}\pm {1}$'.format(round(K,int(-math.log10(K))+1),round(DeltaK,int(-math.log10(DeltaK))+1)))
  return K_list,DeltaK


def AnalyzeAreaPerLipid(t,lipids):
  n=lipids.GetResidueCount()
  Al=FloatList()
  for i in range(t.GetFrameCount()):
    f=t.GetFrame(i)
    cell_size=f.GetCellSize()
    Al.append(cell_size[0]*cell_size[1])
  return 2.*npy.average(Al)/float(n)

