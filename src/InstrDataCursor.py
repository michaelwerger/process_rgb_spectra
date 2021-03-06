import matplotlib.pyplot as plt
import matplotlib.widgets as widgets
import numpy as np

class InstrDataCursor(object):
    '''
    Allows to click on a displayed plot to determine the instrument curve
    stop by closing the window
    it returns two arrays: positons and corresponding intensity values
    '''
    text_template = 'x: %0.2f\ny: %0.2f'
    x, y = 0.0, 0.0
    xoffset, yoffset = -20, 20
    text_template = 'x: %0.2f\ny: %0.2f'

    

    def __init__(self, ax):
        self.ax = ax
        self.positions = []
        self.intensities = []
        self.annotation = ax.annotate(self.text_template, 
                xy=(self.x, self.y), xytext=(self.xoffset, self.yoffset), 
                textcoords='offset points', ha='right', va='bottom',
                bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0')
                )
        self.annotation.set_visible(False)

    def __call__(self, event):
        self.event = event
        # xdata, ydata = event.artist.get_data()
        # self.x, self.y = xdata[event.ind], ydata[event.ind]
        self.x, self.y = event.mouseevent.xdata, event.mouseevent.ydata
        if self.x is not None:
            self.annotation.xy = self.x, self.y
            self.positions.append(float(self.x))
            self.intensities.append(float(self.y))
            self.annotation.set_text(self.text_template % (self.x, self.y))
            self.annotation.set_visible(True)
            event.canvas.draw()

    def getPositions(self):
        return self.positions

    def getIntensities(self):
        return self.intensities
