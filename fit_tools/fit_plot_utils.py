#!/usr/bin/env python

"""fit_plot_util.py: module is dedicated to plot fit level data to RTI and FoV."""

__author__ = "Chakraborty, S."
__copyright__ = "Copyright 2021, SuperDARN@VT"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter, num2date
from matplotlib import patches
import matplotlib.patches as mpatches

import fit_utils as utils

class RadarImages(object):
    
    def __init__(self):
        return
    
    def _add_colorbar(self, bounds, colormap, label=""):
        """
        Add a colorbar to the right of an axis.
        """
        pos = self.ax.get_position()
        cpos = [pos.x1 + 0.025, pos.y0 + 0.0125,
                0.015, pos.height * 0.9]                # this list defines (left, bottom, width, height
        cax = self.fig.add_axes(cpos)
        norm = mpl.colors.BoundaryNorm(bounds, colormap.N)
        cb2 = mpl.colorbar.ColorbarBase(cax, cmap=colormap,
                                        norm=norm,
                                        ticks=bounds,
                                        spacing="uniform",
                                        orientation="vertical")
        cb2.set_label(label)
        return
    
    def _add_key_specific_colorbar(self, colormap, label="gs", idh=1):
        """ Add a colorbar to the right of an axis. """
        pos = self.ax.get_position()
        cpos = [pos.x1 + 0.035, pos.y0 + ((idh-1)*.1+0.05)*pos.height,
                0.01, 0.02]            # this list defines (left, bottom, width, height)
        cax = self.fig.add_axes(cpos)
        cb2 = mpl.colorbar.ColorbarBase(cax, cmap=colormap,
                spacing="uniform",
                orientation="vertical")
        cb2.ax.tick_params(size=0)
        cb2.ax.text(0.5,-.3,label,fontdict={"size":6}, ha="center", va="center")
        # Remove the outer bounds in tick labels
        cb2.ax.set_yticklabels([])
        return
    
    def save(self, filepath):
        self.fig.savefig(filepath, bbox_inches="tight")
        return
    
    def close(self):
        self.fig.clf()
        plt.close()
        return

class RangeTimeIntervalPlot(RadarImages):
    """
    Create plots for velocity, width, power, elevation angle, etc.
    """
    
    def __init__(self, nrang, unique_times, fig_title="", num_subplots=3):
        self.nrang = nrang
        self.unique_gates = np.linspace(1, nrang, nrang)
        self.unique_times = unique_times
        self.num_subplots = num_subplots
        self._num_subplots_created = 0
        self.fig = plt.figure(figsize=(8, 3*num_subplots), dpi=100) # Size for website
        plt.suptitle(fig_title, x=0.075, y=0.99, ha="left", fontweight="bold", fontsize=15)
        mpl.rcParams.update({"font.size": 10})
        return
    
    def addPlot(self, df, beam, title="", p_max=100, p_min=-100, p_step=25, p_ub=9999, p_lb=-9999,
                xlabel="Time [UT]", zparam="v", label="Velocity [m/s]", cmap=mpl.pyplot.get_cmap("Spectral"),
                add_colorbar=True, colorbar_label="Velocity [m/s]", **kwargs):
        self._add_axis()
        df = df[df.bmnum==beam]
        X, Y, Z = utils.get_gridded_parameters(df, xparam="time", yparam="slist", zparam=zparam)
        bounds = list(range(p_min, p_max + 1, p_step))
        bounds.insert(0, p_lb)
        bounds.append(p_ub)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        # Configure axes
        self.ax.xaxis.set_major_formatter(DateFormatter("%H:%M"))
        hours = mdates.HourLocator(byhour=range(0, 24, 4))
        self.ax.xaxis.set_major_locator(hours)
        self.ax.set_xlabel(xlabel, fontdict={"size":12, "fontweight": "bold"})
        self.ax.set_xlim([self.unique_times[0], self.unique_times[-1]])
        self.ax.set_ylim([0, self.nrang])
        self.ax.set_ylabel("Range gate", fontdict={"size":12, "fontweight": "bold"})
        self.ax.pcolormesh(X, Y, Z.T, lw=0.01, edgecolors="None", cmap=cmap, norm=norm, vmax=p_max, vmin=p_min)
        if add_colorbar: self._add_colorbar(bounds, cmap, label=label)
        self.ax.set_title(title, loc="left", fontdict={"fontweight": "bold"})
        return
    
    def _add_axis(self):
        self._num_subplots_created += 1
        self.ax = self.fig.add_subplot(self.num_subplots, 1, self._num_subplots_created)
        return self.ax
    
    
        
        
class FoVPlot(RadarImages):
    """
    Create FoV plots for velocity, width, power, elevation angle, etc.
    """
    
    def __init__(self):
        return