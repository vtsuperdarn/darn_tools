#!/usr/bin/env python

"""fit_utils.py: utility module fitacf<v> level data."""

__author__ = "Chakraborty, S."
__copyright__ = "Copyright 2020, SuperDARN@VT"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import pandas as pd
import time
from netCDF4 import Dataset, date2num, num2date
import numpy as np

import pydarnio
import pydarn
import cartopy

def beams_to_pd(beams, s_params=["bmnum", "noise.sky", "tfreq", "scan", "nrang", "time"],
            v_params=["v", "w_l", "gflg", "p_l", "slist"]):
    """
        Convert beams to pandas
        Parameters:
        -----------
        beams: <list> List of beams
        s_params: Optional[<list>] Scaler parameters
        v_params: Optional[<list>] Vector parameters
    """
    _o = dict(zip(s_params+v_params, ([] for _ in s_params+v_params)))
    for b in beams:
        l = len(getattr(b, "slist"))
        for p in v_params:
            _o[p].extend(getattr(b, p))
        for p in s_params:
            _o[p].extend([getattr(b, p)]*l)
    L = len(_o["slist"])
    for p in s_params+v_params:
        if len(_o[p]) < L:
            l = len(_o[p])
            _o[p].extend([np.nan]*(L-l))
    return pd.DataFrame.from_records(_o)

def scans_to_pd(scans, s_params=["bmnum", "noise.sky", "tfreq", "scan", "nrang", "time"],
            v_params=["v", "w_l", "gflg", "p_l", "slist"]):
    """
        Convert scans to pandas
        Parameters:
        -----------
        scans: <list> List of scans
        s_params: Optional[<list>] Scaler parameters
        v_params: Optional[<list>] Vector parameters
    """
    _o = dict(zip(s_params+v_params, ([] for _ in s_params+v_params)))
    for scan in scans:
        for b in scan.beams:
            l = len(getattr(b, "slist"))
            for p in v_params:
                _o[p].extend(getattr(b, p))
            for p in s_params:
                _o[p].extend([getattr(b, p)]*l)
    L = len(_o["slist"])
    for p in s_params+v_params:
        if len(_o[p]) < L:
            l = len(_o[p])
            _o[p].extend([np.nan]*(L-l))
    return pd.DataFrame.from_records(_o)

def save_to_csv(fname, scans=None, beams=None, s_params=["bmnum", "noise.sky", "tfreq", "scan", "nrang", "time"],
            v_params=["v", "w_l", "gflg", "p_l", "slist"]):
    """
        Save beams or scans to .csv format files
        Parameters:
        -----------
        fname: <str> File name
        scans: Optional[<list>] List of scans
        beams: Optional[<list>] List of beams
        s_params: Optional[<list>] Scaler parameters
        v_params: Optional[<list>] Vector parameters
    """
    if scans is not None:
        df = scans_to_pd(scans, s_params, v_params)
        df.to_csv(fname, index=False, header=True)
    if beams is not None:
        df = beams_to_pd(beams, s_params, v_params)
        df.to_csv(fname, index=False, header=True)
    return

def save_to_netcdf(fname, scans, th=np.nan):
    """
        Save beams or scans to .nc format files
        Parameters:
        -----------
        fname: <str> File name
        scans: <list> List of scans
    """
    s_params, type_params, sdesc_params = ["bmnum", "noise.sky", "tfreq", "scan", "nrang", "intt.sc", "intt.us", "mppul"],\
                ["i1","f4","f4","i1","f4","f4","f4","i1"],\
                ["Beam numbers", "Sky Noise", "Frequency", "Scan Flag", "Max. Range Gate", "Integration sec", "Integration u.sec",
                        "Number of pulses"]
    v_params, vdesc_params = ["v", "w_l", "gflg", "p_l", "slist"],\
                ["LoS Velocity (+ve: towards the radar, -ve: away from radar)", "LoS Width",
                        "GS Flag", "LoS Power", "Gate"]
    if scans is not None:
        _u = {key: [] for key in v_params + s_params}
        slist = []
        t = []
        blen, glen = 0, 110
        for fscan in scans:
            blen += len(fscan.beams)
            for b in fscan.beams:
                t.append(b.time)
                l = len(getattr(b, "slist"))
                slist.append(getattr(b, "slist"))
                for p in v_params:
                    _u[p].extend(getattr(b, p))
                for p in s_params:
                    _u[p].append(getattr(b, p))
        
        rootgrp = Dataset(fname, "w", format="NETCDF4")
        rootgrp.description = """
                                 Fitacf++ : Boxcar filtered data.
                                 Filter parameter: weight matrix - default; threshold - {th}
                              """.format(th=th)
        rootgrp.history = "Created " + time.ctime(time.time())
        rootgrp.source = "SuperDARN - SD data processing: Median filter"
        rootgrp.createDimension("nbeam", blen)
        rootgrp.createDimension("ngate", glen)
        beam = rootgrp.createVariable("nbeam","i1",("nbeam",))
        gate = rootgrp.createVariable("ngate","i1",("ngate",))
        beam[:], gate[:] = range(blen), range(glen)
        times = rootgrp.createVariable("time", "f8", ("nbeam",))
        times.units = "hours since 1970-01-01 00:00:00.0"
        times.calendar = "julian"
        times[:] = date2num(t,units=times.units,calendar=times.calendar)
        
        for _i, k in enumerate(s_params):
            tmp = rootgrp.createVariable(k, type_params[_i],("nbeam",))
            tmp.description = sdesc_params[_i]
            tmp[:] = np.array(_u[k])
        
        for _i, k in enumerate(v_params):
            tmp = rootgrp.createVariable(k, "f4", ("nbeam", "ngate"))
            tmp.description = vdesc_params[_i]
            _m = np.empty((blen,glen))
            _m[:], x = np.nan, _u[k]
            for _j in range(blen):
                _m[_j,slist[_j]] = np.array(x[_j])
            tmp[:] = np.ma.masked_invalid(_m)
        rootgrp.close()
    return

def overlay_radar(rads, ax, tx=cartopy.crs.PlateCarree(), zorder=2, markerColor="darkblue", markerSize=15, fontSize=10, font_color="k",
            xOffset=None, yOffset=-0.1, annotate=True):
    """ 
        Adding the radar location
        Parameters:
        -----------
        rads: <list> Radar cod list
    """
    for rad in rads:
        hdw = pydarn.read_hdw_file(rad)
        print(" Radar:", self.hdw.geographic.lon, self.hdw.geographic.lat)
        self.scatter([self.hdw.geographic.lon], [self.hdw.geographic.lat], s=markerSize, marker="D",
                     color=markerColor, zorder=zorder, transform=tx)
        nearby_rad = [["adw", "kod", "cve", "fhe", "wal", "gbr", "pyk", "aze", "sys"],
                      ["ade", "ksr", "cvw", "fhw", "bks", "sch", "sto", "azw", "sye"]]
        if annotate:
            if self.rad in nearby_rad[0]: xOff, ha = 0.1 if not xOffset else xOffset, 0
            elif self.rad in nearby_rad[1]: xOff, ha = -0.1 if not xOffset else xOffset, 1
            else: xOff, ha = 0.0, 0.5
            lon, lat = self.hdw.geographic.lon, self.hdw.geographic.lat
            x, y = self.projection.transform_point(lon, lat, src_crs=tx)
    return

def overlay_fov(rads, ax, tx=cartopy.crs.PlateCarree(), maxGate=75, rangeLimits=None, beamLimits=None,
            model="IS", fov_dir="front", fovColor=None, fovAlpha=0.2,
            fovObj=None, zorder=2, lineColor="k", lineWidth=1, ls="-"):
    """ Overlay radar FoV """
    lcolor = lineColor
    from numpy import transpose, ones, concatenate, vstack, shape
    for rad in rads:
        hdw = pydarn.read_hdw_file(rad)
        sgate = 0
        egate = hdw.gates if not maxGate else maxGate
        ebeam = hdw.beams
        if beamLimits is not None: sbeam, ebeam = beamLimits[0], beamLimits[1]
        else: sbeam = 0
        self.rad_fov = rad_fov.CalcFov(hdw=hdw, ngates=egate)
        xyz = self.projection.transform_points(tx, self.rad_fov.lonFull, self.rad_fov.latFull)
        x, y = xyz[:, :, 0], xyz[:, :, 1]
        contour_x = concatenate((x[sbeam, sgate:egate], x[sbeam:ebeam, egate],
                                 x[ebeam, egate:sgate:-1],
                                 x[ebeam:sbeam:-1, sgate]))
        contour_y = concatenate((y[sbeam, sgate:egate], y[sbeam:ebeam, egate],
                                 y[ebeam, egate:sgate:-1],
                                 y[ebeam:sbeam:-1, sgate]))
        self.plot(contour_x, contour_y, color=lcolor, zorder=zorder, linewidth=lineWidth, ls=ls)
        if fovColor:
            contour = transpose(vstack((contour_x, contour_y)))
            polygon = Polygon(contour)
            patch = PolygonPatch(polygon, facecolor=fovColor, edgecolor=fovColor, alpha=fovAlpha, zorder=zorder)
            self.add_patch(patch)
    return