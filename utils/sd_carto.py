#!/usr/bin/env python

"""sd_carto.py: utility module for Costom Carto py geoaxes to plot data on aacgmv2 coordinates."""

__author__ = "Kunduri, B.; Chakraborty, S."
__copyright__ = "Copyright 2020, SuperDARN@VT"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Kunduri, B.; Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import matplotlib
matplotlib.use("Agg")
import numpy as np

import cartopy
from cartopy.mpl.geoaxes import GeoAxes
import aacgmv2
import numpy
from shapely.geometry import  MultiLineString, mapping, LineString, Polygon
from descartes.patch import PolygonPatch
from matplotlib.projections import register_projection
import copy
import datetime as dt

import rad_fov

class SDCarto(GeoAxes):
    name = "sdcarto"

    def __init__(self, *args, **kwargs):
        self.supported_coords = [ "geo", "aacgmv2", "aacgmv2_mlt" ]
        if "coords" in kwargs:
            self.coords = kwargs.pop("coords")
            if self.coords not in self.supported_coords:
                err_str = "Coordinates not supported, choose from : "
                for _n,_sc in enumerate(self.supported_coords):
                    if _n + 1 != len(self.supported_coords): err_str += _sc + ", "
                    else: err_str += _sc
                raise TypeError(err_str)
        else: self.coords = "geo"
        if "map_projection" in kwargs: self.map_projection = kwargs.pop("map_projection")
        else: self.map_projection = cartopy.crs.PlateCarree()
        if "rad" in kwargs: self.rad = kwargs.pop("rad")
        else: self.rad = "bks"
        if "plot_date" in kwargs: self.plot_date = kwargs.pop("plot_date")
        else:
            if self.coords == "aacgmv2" or self.coords == "aacgmv2_mlt":
                raise TypeError("Need to provide a date using 'plot_date' keyword for aacgmv2 plotting")
        super().__init__(map_projection=self.map_projection,*args, **kwargs)
        return
    
    def overaly_coast_lakes(self, resolution="50m", color="black", **kwargs):
        """  Overlay AACGM coastlines and lakes  """
        kwargs["edgecolor"] = color
        kwargs["facecolor"] = "none"
        # overaly coastlines
        feature = cartopy.feature.NaturalEarthFeature("physical", "coastline",
                resolution, **kwargs)
        self.add_feature( cartopy.feature.COASTLINE, **kwargs )
        self.add_feature( cartopy.feature.LAKES, **kwargs )
        return

    def coastlines(self,resolution="50m", color="black", **kwargs):
        # details!
        kwargs["edgecolor"] = color
        kwargs["facecolor"] = "none"
        feature = cartopy.feature.NaturalEarthFeature("physical", "coastline",
                resolution, **kwargs)
        return self.add_feature(feature, **kwargs)
    
    def add_feature(self, feature, **kwargs):
        if "edgecolor" not in kwargs: kwargs["edgecolor"] = "black"
        if self.coords == "geo":
            super().add_feature(feature, **kwargs)
        else:
            aacgm_geom = self.get_aacgm_geom(feature)
            aacgm_feature = cartopy.feature.ShapelyFeature(aacgm_geom, cartopy.crs.Geodetic(), **kwargs)
            kwargs["facecolor"] = "none"
            super().add_feature(aacgm_feature, **kwargs)
        return

    def grid_on(self, tx=cartopy.crs.PlateCarree(), draw_labels=[True, True, True, True], 
                linewidth=0.5, color="gray", alpha=0.5, linestyle="--"):
        """ Adding grids to map """
        import matplotlib.ticker as mticker
        from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
        gl = self.gridlines(crs=tx, draw_labels=True,
                linewidth=linewidth, color=color, alpha=alpha, linestyle=linestyle)
        gl.top_labels = draw_labels[0]
        gl.bottom_labels = draw_labels[2]
        gl.right_labels = draw_labels[1]
        gl.left_labels = draw_labels[3]
        gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 15))
        gl.ylocator = mticker.FixedLocator(np.arange(-90, 90, 15))
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        return
    
    def add_dn_terminator(self, **kwargs):
        """ Adding day night terminator """
        from cartopy.feature.nightshade import Nightshade
        if self.plot_date:
            ns_feature = Nightshade(self.plot_date, alpha=0.2)
            super().add_feature(feature, **kwargs)
        return
    
    def get_aacgm_geom(self, feature, out_height=300. ):
        new_i = []
        # cartopy.feature.COASTLINE
        for _n,i in enumerate(feature.geometries()):
            aa = mapping(i)
            mag_list = []
            geo_coords = aa["coordinates"]
            for _ni, _list in enumerate(geo_coords):
                mlon_check_jump_list = []
                split_mag_list = None
                if len(_list) == 1:
                    _loop_list = _list[0]
                else:
                    _loop_list = _list
                for _ngc, _gc in enumerate(_loop_list):
                    _mc = aacgmv2.get_aacgm_coord(_gc[1], _gc[0], out_height, self.plot_date)
                    if numpy.isnan(_mc[0]):
                        continue 
                    mlon_check_jump_list.append( _mc[1] )
                    if self.coords == "aacgmv2":
                        mag_list.append( (_mc[1], _mc[0]) )
                    else:
                        if _mc[2]*15. > 180.:
                            mag_list.append( (_mc[2]*15.-360., _mc[0]) )
                        else:
                            mag_list.append( (_mc[2]*15., _mc[0]) )
                # check for unwanted jumps
                mlon_check_jump_list = numpy.array( mlon_check_jump_list )

                jump_arr = numpy.diff( mlon_check_jump_list )
                bad_inds = numpy.where( numpy.abs(jump_arr) > 10.)[0]
                # delete the range of bad values
                # This is further complicated because
                # in some locations mlon jumps from -177 to +178
                # and this causes jumps in the maps! To deal with 
                # this we'll split arrays of such jumps 
                # (these jumps typically have just one bad ind )
                # and make them into two seperate entities (LineStrings)
                # so that shapely will treat them as two seperate boundaries!
                if len(bad_inds) > 0:
                    if len(bad_inds) > 1:
                        mag_list = [i for j, i in enumerate(mag_list) if j-1 not in numpy.arange(bad_inds[0], bad_inds[1])]
                    else:
                        split_mag_list = mag_list[bad_inds[0]+1:]
                        mag_list = mag_list[:bad_inds[0]+1]
                mag_coords = tuple(mag_list)
                if len(mag_list) > 1:
                    new_i.append( mag_coords )
                if split_mag_list is not None:
        #             print(split_mag_list)
                    if len(split_mag_list) > 1:
                        new_i.append( tuple(split_mag_list) )

        aacgm_coast = MultiLineString( new_i )
        return aacgm_coast
    
    def mark_latitudes(self, lat_arr, lon_location=45, **kwargs):
        """
        mark the latitudes
        Write down the latitudes on the map for labeling!
        we are using this because cartopy doesn't have a 
        label by default for non-rectangular projections!
        """
        if isinstance(lat_arr, list):
            lat_arr = numpy.array(lat_arr)
        else:
            if not isinstance(lat_arr, numpy.ndarray):
                raise TypeError('lat_arr must either be a list or numpy array')
        # make an array of lon_location
        lon_location_arr = numpy.full( lat_arr.shape, lon_location )
        proj_xyz = self.projection.transform_points(\
                            cartopy.crs.Geodetic(),\
                            lon_location_arr,\
                            lat_arr
                            )
        # plot the lats now!
        out_extent_lats = False
        for _np,_pro in enumerate(proj_xyz[..., :2].tolist()):
            # check if lats are out of extent! if so ignore them
            lat_lim = self.get_extent(crs=cartopy.crs.Geodetic())[2::]
            if (lat_arr[_np] >= min(lat_lim)) and (lat_arr[_np] <= max(lat_lim)):
                self.text( _pro[0], _pro[1], str(lat_arr[_np]), **kwargs)
            else:
                out_extent_lats = True
        if out_extent_lats:
            print( "some lats were out of extent ignored them" )

    def mark_longitudes(self, lon_arr=numpy.arange(-180,180,60), **kwargs):
        """
        mark the longitudes
        Write down the longitudes on the map for labeling!
        we are using this because cartopy doesn't have a 
        label by default for non-rectangular projections!
        This is also trickier compared to latitudes!
        """
        if isinstance(lon_arr, list):
            lon_arr = numpy.array(lon_arr)
        else:
            if not isinstance(lon_arr, numpy.ndarray):
                raise TypeError('lat_arr must either be a list or numpy array')
        # get the boundaries
        [x1, y1], [x2, y2] = self.viewLim.get_points()
        bound_lim_arr = []
        right_bound = LineString(([-x1, y1], [x2, y2]))
        top_bound = LineString(([x1, -y1], [x2, y2]))
        bottom_bound = LineString(([x1, y1], [x2, -y2]))
        left_bound = LineString(([x1, y1], [-x2, y2]))
        plot_outline = MultiLineString( [\
                                        right_bound,\
                                        top_bound,\
                                        bottom_bound,\
                                        left_bound\
                                        ] )
        # get the plot extent, we'll get an intersection
        # to locate the ticks!
        plot_extent = self.get_extent(cartopy.crs.Geodetic())
        line_constructor = lambda t, n, b: numpy.vstack(\
                        (numpy.zeros(n) + t, numpy.linspace(b[2], b[3], n))\
                        ).T
        for t in lon_arr:
            xy = line_constructor(t, 30, plot_extent)
            # print(xy)
            proj_xyz = self.projection.transform_points(\
                            cartopy.crs.PlateCarree(), xy[:, 0], xy[:, 1]\
                            )
            xyt = proj_xyz[..., :2]
            ls = LineString(xyt.tolist())
            locs = plot_outline.intersection(ls)
            if not locs:
                continue
            # we need to get the alignment right
            # so get the boundary closest to the label
            # and plot it!
            closest_bound =min( [\
                            right_bound.distance(locs),\
                            top_bound.distance(locs),\
                            bottom_bound.distance(locs),\
                            left_bound.distance(locs)\
                            ] )
            if closest_bound == right_bound.distance(locs):
                ha = 'left'
                va = 'top'
            elif closest_bound == top_bound.distance(locs):
                ha = 'left'
                va = 'bottom'
            elif closest_bound == bottom_bound.distance(locs):
                ha = 'left'
                va = 'top'
            else:
                ha = 'right'
                va = 'top'
            if self.coords == "aacgmv2_mlt":
                marker_text = str(int(t/15.))
            else:
                marker_text = str(t)
            self.text( locs.bounds[0],locs.bounds[1], marker_text, ha=ha, va=va)

register_projection(SDCarto)


class RangeTimeIntervalPlot(object):
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
        plt.suptitle(fig_title, x=0.075, y=0.99, ha='left', fontweight='bold', fontsize=15)
        mpl.rcParams.update({'font.size': 10})
        return
    
    def addPlot(self, data_dict, beam, title, vel_max=200, vel_step=25, xlabel='Time UT'):
        # add new axis
        self.vel_ax = self._add_axis()
        # set up variables for plotter
        time = np.hstack(data_dict['time'])
        gate = np.hstack(data_dict['gate'])
        allbeam = np.hstack(data_dict['beam'])
        flags = np.hstack(data_dict['vel'])
        bounds = list(range(-vel_max, vel_max+1, vel_step))
        cmap = plt.cm.jet
        mask = allbeam == beam
        self._create_colormesh(self.vel_ax, time, gate, flags, mask, bounds, cmap, xlabel)
        self._tight_layout()    # need to do this before adding the colorbar, because it depends on the axis position
        self._add_colorbar(self.fig, self.vel_ax, bounds, cmap, label='Velocity [m/s]')
        self.vel_ax.set_title(title, loc='left', fontdict={'fontweight': 'bold'})
        return
    
    def _add_axis(self):
        self._num_subplots_created += 1
        ax = self.fig.add_subplot(self.num_subplots, 1, self._num_subplots_created)
        return ax
    
    def _add_colorbar(self, fig, ax, bounds, colormap, label=''):
        """
        Add a colorbar to the right of an axis.
        """
        import matplotlib as mpl
        pos = ax.get_position()
        cpos = [pos.x1 + 0.025, pos.y0 + 0.0125,
                0.015, pos.height * 0.9]                # this list defines (left, bottom, width, height
        cax = fig.add_axes(cpos)
        norm = mpl.colors.BoundaryNorm(bounds, colormap.N)
        cb2 = mpl.colorbar.ColorbarBase(cax, cmap=colormap,
                                        norm=norm,
                                        ticks=bounds,
                                        spacing='uniform',
                                        orientation='vertical')
        cb2.set_label(label)