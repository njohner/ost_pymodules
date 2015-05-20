"""
.. codeauthor:: Niklaus Johner <niklaus.johner@a3.epfl.ch>

This module is a collection of functions used in conjunction with namd's colvar module
"""
from ost import *

__all__=('ReadNamdSMDFile','IntegrateFloatList','ReadColvarFile')

def ReadNamdSMDFile(filename):
  """
  This function opens and reads a namd ouput file
  It looks for the SMD ouput and returns it in a dictionary
  with keys ['timestep','cm_pos','force']
  """
  file=open(filename,'r')
  file.seek(0)
  l=file.next()
  s=l.split()
  fields=['timestep','cm_pos','force']
  ts=FloatList()
  cm_pos=geom.Vec3List()
  force=geom.Vec3List()
  for l in file:
    if not l.startswith('SMD '):continue
    s=l.split()
    ts.append(float(s[1]))
    cm_pos.append(geom.Vec3(float(s[2]),float(s[3]),float(s[4])))
    force.append(geom.Vec3(float(s[5]),float(s[6]),float(s[7])))
  out_vecs=[ts,cm_pos,force]
  out_dict={}
  for f,v in zip(fields,out_vecs):
    out_dict[f]=v
  return out_dict


def IntegrateFloatList(x_list,y_list,lim=-1):
  """
  Function to integrate a FloatList y over x
  from x_list[0] to x_list[lim]
  """
  if not len(x_list)==len(y_list):
    print 'x and y don\'t have the same length'
    return
  if lim==-1:lim=len(x_list)
  sum=0
  for i in range(lim-1):
    y=0.5*(y_list[i+1]+y_list[i])
    #y=y_list[i]
    dx=x_list[i+1]-x_list[i]
    sum+=dx*y
  return sum

def ReadColvarFile(filename,skip_first_n=1):
  """
  This function opens and reads a namd colvar trajectory file
  It searches for the fields in the first line of the file and then generates
  a dictionary with these fields as keys
  """
  file=open(filename,'r')
  file.seek(0)
  l=file.next()
  s=l.split()
  fields=[]
  out_vecs=[]
  #We skip the first line as it is for t=0, not in the trajectory
  for i in range(skip_first_n):
    file.next()
  for el in s[1:]:
    fields.append(el)
    out_vecs.append(FloatList())
  print 'fields in the colvar file:',fields
  for l in file:
    if l.startswith('#'):continue
    s=l.split()
    v_list=[]
    for el in s:
      try:
        v_list.append(float(el))
      except:
        print 'error with line',l
    for i,el in enumerate(v_list):  
      out_vecs[i].append(el)
  out_dict={}
  for f,v in zip(fields,out_vecs):
    out_dict[f]=v
  return out_dict
