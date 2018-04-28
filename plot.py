# -*- coding: utf-8 -*-

import pandas as pd

from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg

app = QtGui.QApplication([])
win = pg.GraphicsWindow(title="Ramsey finger")
win.resize(1000, 600)
win.setWindowTitle('Plotting')
# Enable antialiasing for prettier plots
pg.setConfigOptions(antialias=True)
p6 = win.addPlot(title="Updating plot")
curve = p6.plot(pen='y')


def update():
    df = pd.read_csv('RamseyFinger.csv',index_col ='index')
    y= df['freq'].values - df['freq'].values[0]
    x = df.index
    curve.setData(y=y, x=x)


timer = QtCore.QTimer()
timer.timeout.connect(update)

timer.start(2000)

win.show()
app.exec_()