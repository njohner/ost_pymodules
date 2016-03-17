#------------------------------------------------------------------------------------------------
#This file is part of the ost_pymodules project (https://github.com/njohner/ost_pymodules).
#
#Copyright 2015 Niklaus Johner
#
#ost_pymodules is free software: you can redistribute it and/or modify
#it under the terms of the GNU Lesser General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#ost_pymodules is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Lesser General Public License for more details.
#
#You should have received a copy of the GNU Lesser General Public License
#along with ost_pymodules.  If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------------------------
"""
.. codeauthor:: Niklaus Johner <niklaus.johner@a3.epfl.ch>

This module uses matplotlib and provides some high level plotting functions
"""
import matplotlib.pyplot as _plt
import matplotlib.gridspec as _gridspec
import ost as _ost

__all__=('PlotList','PlotFloatProp','PlotDistanceBetwAtoms','PlotHistDistanceBetwAtoms','PlotMinDistance')

smooth=_ost.mol.alg.trajectory_analysis.smooth

def PlotList(vec_list,x_list=[],labels=[],xlabel='',ylabel='',filename=None,ncols=1,hlines=[],ylim=(),title='',smooth_level=0):
  _plt.figure()
  linestyles=['-' for i in range(5)]
  linestyles.extend(['--' for i in range(5,len(vec_list))])
  if not labels:labels=['' for el in vec_list]
  if not x_list:x_list=[range(len(v)) for v in vec_list]
  if labels:
    for x,v,l,ls in zip(x_list,vec_list,labels,linestyles):
      _plt.plot(x,smooth(v,smooth_level),label=l,linestyle=ls)
  else:
    for x,v,ls in zip(x_list,vec_list,linestyles):
      _plt.plot(x,smooth(v,smooth_level),linestyle=ls)
  _plt.xlabel(xlabel)
  _plt.ylabel(ylabel)
  if labels:_plt.legend(loc='best',ncol=ncols)
  if hlines:
    for el in hlines:_plt.hlines(el,x_list[0][0],x_list[0][-1],linestyles='dashed')
  if ylim:_plt.ylim(ylim[0],ylim[1])
  if title:_plt.title(title)
  _plt.show()
  if filename:_plt.savefig(filename)
  return _plt

def PlotFloatProp(eh,prop_name,xlabel='',ylabel='',title='',outname='',yrange=(0,0),sm=0):
  SecStrucLists={'H':[],'E':[],'C':[]}
  r_prev=eh.residues[0]
  ss_prev=r_prev.GetSecStructure()
  if ss_prev.IsHelical():ss_prev='H'
  elif ss_prev.IsExtended():ss_prev='E'
  else:ss_prev='C'
  SecStrucLists[ss_prev].append([r_prev.number.num])
  for r in eh.residues[1:]:
    ss=r.GetSecStructure()
    if ss.IsHelical():ss='H'
    elif ss.IsExtended():ss='E'
    else:ss='C'
    if ss==ss_prev:
      r_prev=r
      continue
    SecStrucLists[ss_prev][-1].append((r_prev.number.num+r.number.num)/2.)
    SecStrucLists[ss].append([(r_prev.number.num+r.number.num)/2.])
    ss_prev=ss
    r_prev=r
  SecStrucLists[ss_prev][-1].append(r_prev.number.num)
  print SecStrucLists
  #We prepare the x and y vector for plotting
  x=[]
  y=[]
  for r in eh.residues:
    x.append(r.number.num)
    y.append(r.GetFloatProp(prop_name))
  #Now we plot
  f=_plt.figure()
  gs = _gridspec.GridSpec(2, 1,height_ratios=[3,1])
  s1=_plt.subplot(gs[0])
  s2=_plt.subplot(gs[1],sharex=s1)
  s1.plot(x,smooth(y,sm),'-o')
  for b in SecStrucLists['E']:
    a=_plt.Arrow(b[0],0,b[1]-b[0],0)
    s2.add_patch(a)
  for h in SecStrucLists['H']:
    r=_plt.Rectangle((h[0],-0.2),h[1]-h[0],0.4)
    s2.add_patch(r)
  for c in SecStrucLists['C']:
    l=_plt.Rectangle((c[0],-0.01),c[1]-c[0],0.02)
    s2.add_patch(l)
  s1.set_title(title)
  s2.set_ylim(-0.5,0.5)
  s1.set_xlim(x[0],x[-1])
  s2.set_axis_off()
  _plt.subplots_adjust(hspace=0.0)
  s2.text(0.5, 0.,xlabel,horizontalalignment='center',verticalalignment='center',transform = s2.transAxes)
  s1.set_ylabel(ylabel)
  if not yrange==(0,0):s1.set_ylim(yrange[0],yrange[1])
  _plt.show()
  if not outname=='':_plt.savefig(outname)
  return f


def PlotDistanceBetwAtoms(t_list,eh_list,cname1,rnum1,aname1,cname2,rnum2,aname2,x=[],labels=[],smooth_level=0,title='',filename=None):
  d=[]
  if len(labels)==0:labels=[str(i) for i in range(len(t_list))]
  if len(x)==0:x=[range(t.GetFrameCount()) for t in t_list]
  for eh,t in zip(eh_list,t_list):
    d.append(_ost.mol.alg.AnalyzeDistanceBetwAtoms(t,eh.FindAtom(cname1,rnum1,aname1),eh.FindAtom(cname2,rnum2,aname2)))
  _plt.figure()
  for xi,di,name in zip(x,d,labels):
    _plt.plot(xi,smooth(di,smooth_level),label=name)
  _plt.xlabel('Time [ns]')
  _plt.ylabel('Distance [$\AA$]')
  _plt.legend(loc='best')
  _plt.title(title)
  _plt.show()
  if filename:_plt.savefig(filename)
  return _plt,d

def PlotHistDistanceBetwAtoms(t_list,eh_list,cname1,rnum1,aname1,cname2,rnum2,aname2,x=[],labels=[],smooth_level1=0,smooth_level2=0,title='',xlabel='',ylabel=''):
  d=[]
  if len(labels)==0:labels=[str(i) for i in range(len(t_list))]
  if len(x)==0:x=[range(t.GetFrameCount()) for t in t_list]
  for eh,t in zip(eh_list,t_list):
    d.append(_ost.mol.alg.AnalyzeDistanceBetwAtoms(t,eh.FindAtom(cname1,rnum1,aname1),eh.FindAtom(cname2,rnum2,aname2)))
  h=_plt.hist(d,normed=True,bins=40)
  _plt.close()
  x1=[]
  for i in range(len(h[1])-1):
    x1.append((h[1][i+1]+h[1][i])/2.)
  f=_plt.figure()
  gs = _gridspec.GridSpec(1, 2,width_ratios=[3,1])
  s1=_plt.subplot(gs[0])
  s2=_plt.subplot(gs[1],sharey=s1)
  for i in range(len(t_list)):
    s1.plot(x[i],smooth(d[i],smooth_level1),label=labels[i][:3])
  s1.set_xlabel(xlabel)
  s1.set_ylabel(ylabel)
  for i in range(len(t_list)):
    el=h[0][i]
    s2.plot(smooth(el,smooth_level2),x1,'-',label=labels[i])
  for label in s2.get_yticklabels():
      label.set_visible(False)
  for label in s2.get_xticklabels()[::2]:
      label.set_visible(False)
  _plt.subplots_adjust(wspace=0.05)
  s2.set_xlabel('Probability')
  t=f.text(0.5, 0.94, title, ha='center', va='center', rotation='horizontal')
  t.set_fontsize(14)
  s1.legend(loc='best')
  _plt.show()
  return _plt

def PlotMinDistance(t_list,eh_list,sele1,sele2,x=[],labels=[],smooth_level=0,title='',filename=None):
  d=[]
  if len(labels)==0:labels=[str(i) for i in range(len(t_list))]
  if len(x)==0:x=[range(t.GetFrameCount()) for t in t_list]
  for eh,t in zip(eh_list,t_list):
    d.append(_ost.mol.alg.AnalyzeMinDistance(t,eh.Select(sele1),eh.Select(sele2)))
  _plt.figure()
  for xi,di,name in zip(x,d,labels):
    _plt.plot(xi,smooth(di,smooth_level),label=name)
  _plt.xlabel('Time [ns]')
  _plt.ylabel('Distance [$\AA$]')
  _plt.legend(loc='best')
  _plt.title(title)
  _plt.show()
  if filename:_plt.savefig(filename)
  return _plt,d
  
