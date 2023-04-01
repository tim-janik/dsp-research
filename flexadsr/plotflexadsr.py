#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# This Source Code Form is licensed MPL-2.0: http://mozilla.org/MPL/2.0

# coding=utf-8

import sys
import subprocess
import math

from PyQt5 import QtGui, QtCore, QtWidgets

def clamp (value, minv, maxv):
  if value < minv:
    return minv
  if value > maxv:
    return maxv
  return value

def osc_cxx (params):
  process = subprocess.Popen(['testadsr', 'plot',
                              "%.5f" % params["a"].value,
                              "%.5f" % params["as"].value,
                              "%.5f" % params["d"].value,
                              "%.5f" % params["ds"].value,
                              "%.5f" % params["s"].value,
                              "%.5f" % params["r"].value,
                              "%.5f" % params["rs"].value,
                              "%.5f" % params["len"].value,
                              params["SHAPE"].value
                              ], stdout=subprocess.PIPE)
  out, err = process.communicate()
  plot_data = []
  for o in out.splitlines():
    floats = [float(x) for x in o.split()]
    f = floats[1]
    plot_data += [ f ]
  return plot_data

def line (qp, x1, y1, x2, y2):
  qp.drawLine (int (x1), int (y1), int (x2), int (y2))

class BlepWidget(QtWidgets.QWidget):
  def __init__(self):
    super(BlepWidget, self).__init__()
    self.setMinimumSize (100, 100)

  def update_params (self, params):
    try:
      self.plot_data = osc_cxx (params)
      self.repaint()
    except KeyError:
      # hacky way of avoiding to run testblep if not all params are known
      pass

  def paintEvent(self, event):
    qp = QtGui.QPainter()
    qp.begin(self)
    qp.fillRect (event.rect(), QtGui.QColor (255, 255, 255));
    xgreycolor = QtGui.QColor (240, 240, 240)
    xscale = self.width() / 48000 / 2
    xcenter = 0 #-self.width() / 5 * 3.5
    yscale = -self.height() / 1.2
    ycenter = self.height() * 0.9
    qp.setPen (QtGui.QColor (190, 190, 190))
    line (qp, 0, ycenter + 1 * yscale, self.width(), ycenter + 1 * yscale)
    line (qp, 0, ycenter + -1 * yscale, self.width(), ycenter + -1 * yscale)
    for x in range (5):
      xx = x * 48000 * xscale / 2
      line (qp, xx, 0.0, xx, self.height())
    qp.setPen (QtCore.Qt.black)
    for i in range (len (self.plot_data) - 1):
      line (qp, i * xscale + xcenter, self.plot_data[i] * yscale + ycenter, (i + 1) * xscale + xcenter, self.plot_data[i+1] * yscale + ycenter)

    qp.setPen (QtGui.QColor(168, 34, 3))
    line (qp, 0, ycenter, self.width(), ycenter)

    qp.end()

class Param:
  pass

class PlotWindow (QtWidgets.QMainWindow):
  def __init__ (self):
    QtWidgets.QMainWindow.__init__(self)

    self.init_ui()

  def init_ui (self):
    central_widget = QtWidgets.QWidget (self)
    self.grid_layout = QtWidgets.QGridLayout()
    self.blep_widget = BlepWidget ()
    self.grid_layout.addWidget (self.blep_widget, 1, 0, 3, 1)
    central_widget.setLayout (self.grid_layout)
    self.setCentralWidget (central_widget)

    self.params = dict()
    self.add_slider (1, "a", 50, 0, 100)
    self.add_slider (2, "as", 0, -100, 100)
    self.add_slider (3, "d", 50, 0, 100)
    self.add_slider (4, "ds", 0, -100, 100)
    self.add_slider (5, "s", 30, 0, 100)
    self.add_slider (6, "r", 50, 0, 100)
    self.add_slider (7, "rs", 0, -100, 100)
    self.add_slider (8, "len", 100, 0, 200)
    self.add_slider (9, "SHAPE", 0, 0, 2);

  def add_slider (self, row, label, vdef, vmin, vmax):
    slider = QtWidgets.QSlider (QtCore.Qt.Vertical)
    slider.setRange (vmin * 10, vmax * 10)
    slider.setValue (vdef * 10)
    slider.setEnabled (True)
    slider.valueChanged.connect (lambda value: self.on_param_changed (label, value))
    self.params[label] = Param()
    self.params[label].value_label = QtWidgets.QLabel()
    self.on_param_changed (label, vdef * 10)
    self.grid_layout.addWidget (slider, 3, row, 1, 1, QtCore.Qt.AlignHCenter)

    self.grid_layout.addWidget (QtWidgets.QLabel (label), 1, row, 1, 1, QtCore.Qt.AlignHCenter)
    self.grid_layout.addWidget (self.params[label].value_label, 2, row, 1, 1, QtCore.Qt.AlignHCenter)
    self.grid_layout.setColumnStretch (0, 1)

  def on_param_changed (self, label, value):
    self.params[label].value = value / 1000
    if (label == "SHAPE"):
      if (value < 6):
        S = "LIN"
      elif (value > 14):
        S = "FLEX"
      else:
        S = "EXP"
      self.params[label].value_label.setText ("%s" % S)
      self.params[label].value = S
    else:
      self.params[label].value_label.setText ("%s" % self.params[label].value)
    self.blep_widget.update_params (self.params)

  def on_source_clicked (self):
    self.blep_widget.update_params (self.params)

def main():
  app = QtWidgets.QApplication (sys.argv)
  plot = PlotWindow()
  plot.show()
  sys.exit (app.exec_())

if __name__ == "__main__":
  main()
