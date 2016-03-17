#------------------------------------------------------------------------------------------------
#This file is part of the ost_pymodules project (https://github.com/njohner/ost_pymodules).
#
#Copyright 2015 Niklaus Johner, Marco Biasini
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

Module written by Marco Biasini.

This module contains functions to startup a scene from python.
"""

from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4.QtOpenGL import *
try:
  from ost import gfx
except ImportError:
  print 'can not find openstructure module. make sure it is in your python path'

import weakref

class _GLWin(gfx.GLWinBase):
  """
  We need this class to get proper redraw events when the scene changes.
  """
  def __init__(self, threed_widget):
    gfx.GLWinBase.__init__(self)
    self._threed_widget = weakref.proxy(threed_widget)
    gfx.Scene().Register(self)

  def DoRefresh(self):
    self._threed_widget.update()
  def Height(self): return self.threed_widget.height()
  def Width(self): return self.threed_widget.width()

class ThreeDWidget(QGLWidget):
  def __init__(self, format, parent=None):
    QGLWidget.__init__(self, format, parent=None)
    self._last_pos = QPoint()
    self._gl_win = _GLWin(self)
    self._key_handler = {}
  def initializeGL(self):
    gfx.Scene().InitGL(True)
    gfx.Scene().fogno=10

  def paintGL(self):
    gfx.Scene().SetShadingMode('hemilight')
    gfx.Scene().RenderGL()

  def resizeGL(self, w, h):
    gfx.Scene().Resize(w, h)

  def HandleKey(self, key, handler):
    self._key_handler[key] = handler

  def keyPressEvent(self, event):
    handler = self._key_handler.get(str(event.text()), None)
    if handler:
      handler()

  def mouseDoubleClickEvent(self,event):
    self._last_pos = QPoint(event.x(), event.y())
    gfx_obj, atom=gfx.PickAtom(gfx.Scene(), event.x(),
                               self.height()-event.y())
    if gfx_obj:
      gfx.Scene().SetCenter(atom.pos)

  def mousePressEvent(self, event):
    self._last_pos = QPoint(event.x(), event.y())
    gfx_obj, atom=gfx.PickAtom(gfx.Scene(), event.x(),
                               self.height()-event.y())
    if gfx_obj:
      print 'picked: ', gfx_obj.name, atom

  def mouseMoveEvent(self, event):
    """
    Handles translation and rotation of the scene...
    """
    delta = QPoint(event.x(), event.y()) - self._last_pos
    if event.buttons() & Qt.LeftButton:
      self._RotateDelta(delta)

    if event.buttons() & Qt.RightButton:
      self._TranslateDelta(delta)

    self._last_pos = QPoint(event.x(), event.y())
  def wheelEvent(self, event):
    """
    Scrolling zooming is...
    """
    gfx.Scene().Apply(gfx.InputEvent(gfx.INPUT_DEVICE_MOUSE,
                                     gfx.INPUT_COMMAND_TRANSZ,
                                     -0.1*event.delta()), False)
    self.update()

  def _TranslateDelta(self, delta):
    """
    Responds to a change in mouse position of delta.
    """
    if delta.y()!=0:
      gfx.Scene().Apply(gfx.InputEvent(gfx.INPUT_DEVICE_MOUSE,
                                       gfx.INPUT_COMMAND_TRANSY,
                                       -delta.y()*0.5), False)
    if delta.x()!=0:
      gfx.Scene().Apply(gfx.InputEvent(gfx.INPUT_DEVICE_MOUSE,
                                       gfx.INPUT_COMMAND_TRANSX,
                                       delta.x()*0.5), False)
    self.update()
  
  def _RotateDelta(self, delta):
    """
    Responds to a change in mouse position of delta.
    """
    if delta.y()!=0:
      gfx.Scene().Apply(gfx.InputEvent(gfx.INPUT_DEVICE_MOUSE,
                                       gfx.INPUT_COMMAND_ROTX,
                                       delta.y()*0.5), False)
    if delta.x()!=0:
      gfx.Scene().Apply(gfx.InputEvent(gfx.INPUT_DEVICE_MOUSE,
                                       gfx.INPUT_COMMAND_ROTY,
                                       delta.x()*0.5), False)
    self.update()

def App(title='OST Viewer', width=800, height=600):
  """
  """
  import sys
  app=QApplication(sys.argv)
  fmt=QGLFormat()
  fmt.setAlpha(True)
  fmt.setSampleBuffers(True)
  wnd=ThreeDWidget(fmt)
  wnd.resize(width, height)
  wnd.setWindowTitle(title)
  wnd.show()
  app.threed=wnd
  return app
if __name__ == '__main__':
  from ost import io
  # example use
  app  = App()
  structure = io.LoadPDB('1crn.pdb')
  obj = gfx.Entity('structure', gfx.HSC, structure)
  obj.SetColor(gfx.YELLOW, 'rtype=H,E')
  obj.SetDetailColor(gfx.RED, 'rtype=H,E')
  gfx.Scene().Add(obj)
  app.exec_()

