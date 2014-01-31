#-------------------------------------------------------------------------------
# This file is part of the OpenStructure project <www.openstructure.org>
#
# Copyright (C) 2008-2011 by the OpenStructure authors
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3.0 of the License, or (at your option)
# any later version.
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#-------------------------------------------------------------------------------
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from ost import *

class TrajWidget(QWidget):
  def __init__(self, trajlist=None, golist=None, parent=None,ref_index=0):
    QWidget.__init__(self, parent, Qt.Tool)
    self.setFocusPolicy(Qt.ClickFocus)
    self.ref_index_=ref_index
    self.trajlist_=trajlist
    self.golist_=golist
    self.index_dict={}
    self.fix_dict={}
    self.modifiers=None
    for i,go in enumerate(self.golist_):
      self.index_dict[go.name]=i
    vb=QVBoxLayout()
    hb=QHBoxLayout()
    hb1=QHBoxLayout()
    hb2=QHBoxLayout()
    hb3=QHBoxLayout()
    self.callback=None   
    self._slider=QSlider(self)
    self._slider.setOrientation(Qt.Horizontal)
    self._speed_slider=QSlider(self)
    self._speed_slider.setOrientation(Qt.Horizontal)
    self._speedLabel=QLabel(self)
    self._speedLabel.setText('Speed:')
    self._speedLabel.setAlignment(Qt.AlignLeft)
    self._play=QToolButton(self)
    self._repeat=QCheckBox(self)
    self._frame=QLabel(self)
    self._frameNo=QLabel(self)
    self._frameEnd=QLabel(self)
    self._repeat.setText('Repeat')
    self._slider.setTracking(True)
    self._play.setText('Play')
    self._play.setCheckable(True)
    self._frame.setText('Frame: ')
    self._frameNo.setNum(0)
    self._frameNo.setAlignment(Qt.AlignRight)
    self._frameEnd.setAlignment(Qt.AlignLeft)
    self._timeLabel=QLabel(self)
    self._timeNo=QLabel(self)
    self._timeEnd=QLabel(self)
    self._timeUnit=QLabel(self)
    self._timeLabel.setText('Time: ')
    self._timeNo.setNum(0)
    self._timeNo.setAlignment(Qt.AlignRight)
    self._timeEnd.setAlignment(Qt.AlignLeft)
    self._timeUnit.setAlignment(Qt.AlignLeft)
    self._speed_slider.setTracking(True)
    self._speed_slider_min=-50
    self._speed_slider_max=-10
    self._speed_slider.setRange(self._speed_slider_min,self._speed_slider_max)
    self._right_arrow=QPushButton(">")
    self._right_arrow2=QPushButton(">>")
    self._right_end=QPushButton(">|")
    self._left_arrow=QPushButton("<")
    self._left_arrow2=QPushButton("<<")
    self._left_end=QPushButton("|<")
    hb.addWidget(self._play)
    hb.addWidget(self._repeat)
    hb1.addWidget(self._left_end)
    hb1.addWidget(self._left_arrow2)
    hb1.addWidget(self._left_arrow)
    hb1.addWidget(self._right_arrow)
    hb1.addWidget(self._right_arrow2)
    hb1.addWidget(self._right_end)
    hb2.addWidget(self._frame)
    hb2.addWidget(self._frameNo)
    hb2.addWidget(self._frameEnd)
    hb2.addWidget(self._timeLabel)
    hb2.addWidget(self._timeNo)
    hb2.addWidget(self._timeEnd)
    hb2.addWidget(self._timeUnit)
    hb3.addWidget(self._speedLabel)
    hb3.addWidget(self._speed_slider)
    self.setLayout(vb)
    vb.addLayout(hb)
    vb.addLayout(hb1)
    vb.addWidget(self._slider)
    vb.addLayout(hb2)
    vb.addLayout(hb3)
    self._SetSpeedSliderPos(self._speed_slider_min+0.5*(self._speed_slider_max-self._speed_slider_min))
    QObject.connect(self._play, SIGNAL('toggled(bool)'), 
                    self._TogglePlay)
    QObject.connect(self._slider, SIGNAL('valueChanged(int)'), 
                    self._SliderValueChanged)
    QObject.connect(self._speed_slider,SIGNAL('valueChanged(int)'),self._SpeedSliderValChanged)
    QObject.connect(self._right_end,SIGNAL('clicked()'),self._RightEndClicked)
    QObject.connect(self._right_arrow2,SIGNAL('clicked()'),self._Right2Clicked)
    QObject.connect(self._right_arrow,SIGNAL('clicked()'),self._RightClicked)
    QObject.connect(self._left_arrow2,SIGNAL('clicked()'),self._Left2Clicked)
    QObject.connect(self._left_arrow,SIGNAL('clicked()'),self._LeftClicked)
    QObject.connect(self._left_end,SIGNAL('clicked()'),self._LeftEndClicked)
    self._slider.setMinimum(0)
    self.frame_number_=self.traj_.GetFrameCount()-1
    self.timestep_=self.traj_.GetDelta()
    self.SetTimeUnit('ns')
    self.SetReferenceIndex(0)
    
  def p2u(self,u):
    if u=="s":
      return 1.0e-12
    elif u=="ms":
      return 1.0e-9
    elif u=="us":
      return 1.0e-6
    elif u=="ns":
      return 1.0e-3
    elif u=="ps":
      return 1.0
    elif u=="fs":
      return 1.0e3
    raise RuntimeError("expected one of s,ms,us,ns,ps or fs for unit")

  def SetTimeUnit(self,u):
    self._time_prefactor=self.p2u(u)
    self._timeEnd.setText('/ '+ str((self.frame_number_-1)*self.timestep_*self._time_prefactor))
    self._timeNo.setNum(self.current_frame*self.timestep_*self._time_prefactor)
    self._timeUnit.setText('['+u+']')
  
  def _SetSpeedSliderPos(self,pos):
    self._speed_slider.setSliderPosition(pos)
    self._SpeedSliderValChanged(pos)
    
  def _SpeedSliderValChanged(self,speed_pos):
    self.time=math.exp(-0.15*speed_pos)
    if self._play.isChecked():
      self._TogglePlay(False)
      self._TogglePlay(True)
  
  def _SetTime(self,t):
    self.time=t
    self._speed_slider.setSliderPosition(-1./0.15*math.log(t))
    if self._play.isChecked():
      self._TogglePlay(False)
      self._TogglePlay(True)
  
  
  def _SliderValueChanged(self, pos):
    self.current_frame=pos
    for traj,go in zip(self.trajlist_,self.golist_):
      if go.name in self.fix_dict:continue
      traj.CopyFrame(self.current_frame)
      go.UpdatePositions()

  def _GetCurrentFrame(self):
    return self._slider.sliderPosition()
    
  def _SetCurrentFrame(self, pos):
    if self._slider.maximum()<pos:
      if self._repeat.isChecked():
        pos=0
      else:
        pos=self._slider.maximum()
    self._slider.setSliderPosition(pos)
    self._frameNo.setNum(pos)
    self._timeNo.setNum(pos*self.timestep_*self._time_prefactor)
    
  current_frame=property(_GetCurrentFrame, _SetCurrentFrame)
  
  def _GetReferenceTraj(self):
    return self.trajlist_[self.ref_index_]
  traj_=property(_GetReferenceTraj)
  
  def timerEvent(self, event):
    #if self.callback:
    #  self.callback(self.golist_)
    self.current_frame+=1
    for go in self.golist_:
      if go.name in self.fix_dict:continue
      go.BlurSnapshot()
      go.UpdatePositions()
  
  def _TogglePlay(self, playing):
    if playing:
      self.timer_id_=self.startTimer(self.time)
    else:
      self.killTimer(self.timer_id_)
  
  def _LeftClicked(self):
    if self.current_frame>0:
      self.current_frame-=1
  
  def _RightClicked(self):
    if self.current_frame<self.frame_number_-1:
      self.current_frame+=1
  
  def _Left2Clicked(self):
    if self.current_frame>=10:
      self.current_frame-=10
  
  def _Right2Clicked(self):
    if self.current_frame<self.frame_number_-10:
      self.current_frame+=10
  
  def _LeftEndClicked(self):
    if self.current_frame>0:
      self.current_frame=0
  
  def _RightEndClicked(self):
    if self.current_frame<self.frame_number_-1:
      self.current_frame=self.frame_number_-1
        
  def _SetBlur(self, blur):
    for go in self.golist_:
      go.SetBlur(blur)

  def _GetBlur(self):
    return self.gfx_entity.GetBlur()

  blur=property(_GetBlur, _SetBlur)

  def SetReferenceIndex(self,ref_index):
    if type(ref_index)==type(''):
      self.ref_index_=self.index_dict[ref_index]
    else:
      self.ref_index_=ref_index
    self.frame_number_=self.traj_.GetFrameCount()
    self.timestep_=self.traj_.GetDelta()
    self._slider.setMaximum(self.frame_number_-1)
    self._frameEnd.setText('/ '+ str(self.frame_number_-1))
    self._timeEnd.setText('/ '+ str((self.frame_number_-1)*self.timestep_*self._time_prefactor))
  
  def keyPressEvent(self, event):
    key=event.key()
    if event.modifiers()==Qt.ControlModifier:
      self.modifiers=event.modifiers()
    if self.modifiers==Qt.ControlModifier:
      if key==Qt.Key_Left:
        self._Left2Clicked()
      elif key==Qt.Key_Right:
        self._Right2Clicked()
    else:
      if key==Qt.Key_Left:
        self._LeftClicked()
      elif key==Qt.Key_Right:
        self._RightClicked()
      elif key==Qt.Key_Space:
        self._play.setChecked(not self._play.isChecked())
      else:
        QWidget.keyPressEvent(self, event)
  
  def keyReleaseEvent(self, event):
    if event.key()==Qt.Key_Control:
      self.modifiers=None
    else:
      QWidget.keyPressEvent(self, event)
  
  def FixGfxEntity(self,index,frame_number=None):
    if not frame_number:frame_number=self.current_frame
    if type(index)==type(''):
      self.fix_dict[index]=frame_number
      i=self.index_dict[index]
    else:
      i=index
      self.fix_dict[self.golist_[i].name]=frame_number
    self.trajlist_[i].CopyFrame(frame_number)
    self.golist_[i].UpdatePositions()
  
  def ReleaseGfxEntity(self,index):
    if type(index)==type(''):self.fix_dict.pop(index)
    else:self.fix_dict.pop(self.golist_[i].name)
  
  def AddTrajectory(self,traj,go):
    self.trajlist_.append(traj)
    self.golist_.append(go)
    self.index_dict[go.name]=len(self.golist_)-1
  
  def RemoveTrajectory(self,index):
    if type(index)==type(''):
      index=self.index_dict[index]
    self.trajlist_.pop(index)
    self.golist_.pop(index)
  
  def SetSpeed(self,val):
  #Value should be between 0 and 1
    if not (val<=1. and val >=0.):
      print 'Speed should be set between 0 and 1'
      return
    else:
      val=self._speed_slider_min-val*(self._speed_slider_min-self._speed_slider_max)
      self._SetSpeedSliderPos(val)
      print val