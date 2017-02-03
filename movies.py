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

This module contains functions to render movies from OpenStructure using ffmpeg
"""
from ost import *
import time
import os

__all__=('InterpolateTrajectory','MakeMovieAndCleanImageFiles','RenderTrajSimple','Fade',\
         'Fades','ClipObject','AdvanceTraj','RenderFrames','Zoom','ZoomAndClip','ClipScene','Rotate')


def InterpolateTrajectory(t,nsteps):
  nframes=t.GetFrameCount()
  cell_sizes=[]
  for i in range(nframes):
    f=t.GetFrame(i)
    cell_sizes.append(f.GetCellSize())
  t=t.Filter(t.GetEntity().Select(''))
  eh=t.GetEntity()
  center=eh.GetCenterOfAtoms()
  step=1./nsteps
  for i in range(1,nframes):
    cell_size=cell_sizes[i-1]
    vl1=t.GetFramePositions(i-1)
    vl2=t.GetFramePositions(i)
    vl2=geom.Vec3List([geom.WrapVec3(el2,el1,cell_size) for el1,el2 in zip(vl1,vl2)])
    for j in range(nsteps):
      vl=vl1+j*step*(vl2-vl1)
      t.AddFrame(geom.WrapVec3List(vl,center,cell_size))
  t.AddFrame(vl2)
  return t.Filter(eh.Select(''),first=nframes)

def MakeMovieAndCleanImageFiles(outdir,outname,nd='4',image_rate=24,quality=30,clean=False,first=0,last=0):
  """
  Uses ffmpeg to encode a movie from png images
  -pix_fmt yuv420p is needed to replay the movie in QuickTime
  """
  filename=os.path.join(outdir,outname+'%0'+nd+'d.png')
  os.system('/opt/local/bin/ffmpeg '+'-f image2 -r '+str(image_rate)+' -i '+filename+' -vcodec libx264 -pix_fmt yuv420p -crf '+str(quality)+' '+os.path.join(outdir,outname+'.avi'))
  if clean:
    for i in range(first,last):
      os.system('rm '+os.path.join(outdir,filename % i))
  return

def RenderTrajSimple(t,go,outdir,outname,width=500,height=500,first=0,last=-1,stride=1,first_index=0,nd='4',ffmpeg_bin="/opt/local/bin/ffmpeg",r_in=24):
  import time
  scene=gfx.Scene()
  filename=os.path.join(outdir,outname+'%0'+nd+'d.png')
  if last==-1:last=t.GetFrameCount()
  j=first_index
  for i in range(first,last,stride):
    scene.StartOffscreenMode(width,height)
    t.CopyFrame(i)
    go.UpdatePositions()
    #go.UpdateView()
    scene.RequestRedraw()
    scene.Export(filename % j,width,height,False)
    scene.StopOffscreenMode()
    j+=1
  os.system(ffmpeg_bin+' -f image2 -r '+str(r_in)+' -i '+filename+' -r 24 -vb 20M '+os.path.join(outdir,outname+'.mpg'))
  j=first_index
  for i in range(first,last,stride):
    os.system('rm '+os.path.join(filename % j))
    j+=1

def Fade(go,nsteps,initial_opacity,final_opacity,outdir,outname,width=500,height=500,first_index=0,nd='4'):
  scene=gfx.Scene()
  step=(final_opacity-initial_opacity)/float(nsteps-1)
  filename=outname+'%0'+nd+'d.png'
  j=first_index
  scene.StartOffscreenMode(width,height)
  for i in range(nsteps):
    go.SetOpacity(initial_opacity+i*step)
    j=first_index+i
    scene.Export(os.path.join(outdir,filename % j),width,height)
  scene.StopOffscreenMode()
  return

def Fades(go_list,nsteps,initial_opacity_list,final_opacity_list,outdir,outname,width=500,height=500,first_index=0,nd='4'):
  scene=gfx.Scene()
  step_list=[(fop-iop)/float(nsteps-1) for fop,iop in zip(final_opacity_list,initial_opacity_list)]
  filename=outname+'%0'+nd+'d.png'
  j=first_index
  for i in range(nsteps):
    scene.StartOffscreenMode(width,height)
    for go,initial_opacity,step in zip(go_list,initial_opacity_list,step_list):
      go.SetOpacity(initial_opacity+i*step)
    j=first_index+i
    scene.Export(os.path.join(outdir,filename % j),width,height)
    scene.StopOffscreenMode()
  return  

def ClipObject(go,nsteps,initial_offset,final_offset,outdir,outname,width=500,height=500,first_index=0,nd='4'):
  scene=gfx.Scene()
  step=(final_offset-initial_offset)/float(nsteps-1)
  filename=outname+'%0'+nd+'d.png'
  j=first_index
  scene.StartOffscreenMode(width,height)
  for i in range(nsteps):
    go.clip_offset=initial_offset+i*step
    j=first_index+i
    scene.Export(os.path.join(outdir,filename % j),width,height)
  scene.StopOffscreenMode()
  return


def ChangeDensity(go_list,nsteps,initial_density_list,final_density_list,outdir,outname,width=500,height=500,first_index=0,nd='4'):
  scene=gfx.Scene()
  step_list=[(fd-id)/float(nsteps-1) for fd,id in zip(final_density_list,initial_density_list)]
  filename=outname+'%0'+nd+'d.png'
  j=first_index
  scene.StartOffscreenMode(width,height)
  for i in range(nsteps):
    for go,initial_density,final_density,step in zip(go_list,initial_density_list,final_density_list,step_list):
      go.SetLevel(initial_density+i*step)
    j=first_index+i
    scene.Export(os.path.join(outdir,filename % j),width,height)
  scene.StopOffscreenMode()
  return  

  
def AdvanceTraj(t,go,outdir,outname,width=500,height=500,first=0,last=-1,stride=1,first_index=0,nd='4'):
  scene=gfx.Scene()
  filename=os.path.join(outdir,outname+'%0'+nd+'d.png')
  if last==-1:last=t.GetFrameCount()
  j=first_index
  for i in range(first,last,stride):
    scene.StartOffscreenMode(width,height)
    t.CopyFrame(i)
    go.UpdatePositions()
    go.UpdateView()
    scene.Export(filename % j,width,height)
    scene.StopOffscreenMode()  
    j+=1
  return

def RenderFrames(nframes,outdir,outname,width,height,first_index,nd='4'):
  scene=gfx.Scene()
  filename=os.path.join(outdir,outname+'%0'+nd+'d.png')
  j=first_index
  for i in range(nframes):
    scene.StartOffscreenMode(width,height)
    scene.Export(filename % j,width,height)
    j+=1
    scene.StopOffscreenMode() 
  return

def Zoom(nsteps,outdir,outname,dz,width=500,height=500,first_index=0,nd='4'):
  import math
  step_size=(dz)/float(nsteps)
  print step_size
  scene=gfx.Scene()
  filename=os.path.join(outdir,outname+'%0'+nd+'d.png')
  scene.StartOffscreenMode(width,height)
  j=first_index
  for i in range(nsteps):
    camera=scene.transform
    camera.ApplyZAxisTranslation(step_size)
    scene.transform=camera
    scene.Export(filename % j, width, height)
    j+=1
  scene.StopOffscreenMode() 
  return

def ZoomAndClip(nsteps,outdir,outname,dz,initial_offset,final_offset,width=500,height=500,first_index=0,nd='4'):
  import math
  step_size=(dz)/float(nsteps)
  step=(final_offset-initial_offset)/float(nsteps-1)
  scene=gfx.Scene()
  filename=os.path.join(outdir,outname+'%0'+nd+'d.png')
  scene.StartOffscreenMode(width,height)
  for i in range(nsteps):
    j=first_index+i
    camera=scene.transform
    camera.ApplyZAxisTranslation(step_size)
    scene.transform=camera
    scene.near=initial_offset+i*step
    scene.Export(filename % j, width, height)
  scene.StopOffscreenMode() 
  return

def ClipScene(nsteps,initial_near,final_near,outdir,outname,width=500,height=500,first_index=0,nd='4',change_far=False,near_to_far=20):
  scene=gfx.Scene()
  step=(final_near-initial_near)/float(nsteps-1)
  filename=outname+'%0'+nd+'d.png'
  scene.StartOffscreenMode(width,height)
  for i in range(nsteps):
    j=first_index+i
    scene.near=initial_near+i*step
    if change_far:scene.far=scene.near+near_to_far
    scene.Export(os.path.join(outdir,filename % j),width,height)
  scene.StopOffscreenMode()
  return

  
def Rotate(nsteps,outdir,outname,width=500,height=500,initial_angle=0,final_angle=6.28,axis=geom.Vec3(0, 1, 0),first_index=0,nd='4'):
  import math
  step_size=(final_angle-initial_angle)/float(nsteps)
  print step_size
  scene=gfx.Scene()
  filename=os.path.join(outdir,outname+'%0'+nd+'d.png')
  scene.StartOffscreenMode(width,height)
  j=first_index
  for i in range(nsteps):
    #initial_angle+=step_size
    #if initial_angle>math.pi*2:
    #  initial_angle-=math.pi*2
    camera=scene.transform
    #rot=geom.AxisRotation(axis, initial_angle)
    #camera.SetRot(rot)
    camera.ApplyAxisRotation(step_size,axis)
    scene.transform=camera
    scene.Export(filename % j, width, height)
    j+=1
  scene.StopOffscreenMode() 
  return
  