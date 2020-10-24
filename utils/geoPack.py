# -*- coding: utf-8 -*-
# Copyright (C) 2012  VT SuperDARN Lab
# Full license can be found in LICENSE.txt
"""Earth coordinate conversion routines
Functions
---------
geodToGeoc  : converts from geodetic to geocentric (and vice-versa)
geodToGeocAzEl : converts azimuth and elevation from geodetic to geocentric
                 (and vice-versa)
gspToGcar : converts global spherical coordinates to global cartesian
            coordinates (and vice-versa)
gcarToLcar : converts from global cartesian coordinates to local cartesian
             coordinates (and vice-versa)
lspToLcar : converts from local spherical coordinates to local cartesian
            coordinates (and vice-versa)
calcDistPnt : calculates the coordines|distance,elevation,azimuth of a point
              given a point of origin and distance, elevation, azimuth|distant
              point coordinates
greatCircleMove : Calculates the coordinates of an end point along a great
                  circle path given the original coordinates, distance, azimuth,
                  and altitude.
greatCircleAzm : Calculates the azimuth from the coordinates of a start point
                 to and end point along a great circle path.
greatCircleDist : Calculates the distance in radians along a great circle path
                  between two points.
References
----------
Based on J.M. Ruohoniemi's geopack
Based on R.J. Barnes radar.pro
Updates based on G. Chishams cnvtcoord_vhm.c

Copied from DaViTPy
"""
import logging
import numpy

def geodToGeoc(lat, lon, inverse=False):
    """Converts position from geodetic to geocentric or vice-versa.
    Based on the IAU 1964 oblate spheroid model of the Earth.
    Parameters
    ----------
    lat : float
        latitude [degree]
    lon : float
        longitude [degree]
    inverse : Optional[bool]
        inverse conversion (geocentric to geodetic).  Default is false.
    Returns
    -------
    lat_out : float
        latitude [degree] (geocentric/detic if inverse=False/True)
    lon_out : float
        longitude [degree] (geocentric/detic if inverse=False/True)
    rade : float
        Earth radius [km] (geocentric/detic if inverse=False/True)
    """
    a = 6378.16
    f = 1.0 / 298.25
    b = a * (1.0 - f)
    e2 = (a**2 / b**2) - 1.0
    
    if not inverse:
        # geodetic into geocentric
        lat_out = numpy.degrees(numpy.arctan(b**2 / a**2 * numpy.tan(numpy.radians(lat))))
        lon_out = lon
    else:
        # geocentric into geodetic
        lat_out = numpy.degrees(numpy.arctan(a**2 / b**2 * numpy.tan(numpy.radians(lat))))
        lon_out = lon
        
    rade = a / numpy.sqrt( 1. + e2 * numpy.sin(numpy.radians(lat_out))**2)
        
    return lat_out, lon_out, rade

def geodToGeocAzEl(lat, lon, az, el, inverse=False):
    """Converts pointing azimuth and elevation measured with respect to the
    local horizon to azimuth and elevation with respect to the horizon defined
    by the plane perpendicular to the Earth-centered radial vector drawn
    through a user defined point.
    Parameters
    ----------
    lat : float
        latitude [degree]
    lon : float
        longitude [degree]
    az : float
        azimuth [degree, N]
    el : float
        elevation [degree]
    inverse : Optional[bool]
        inverse conversion
    Returns
    -------
    lat : float
        latitude [degree]
    lon : float
        longitude [degree]
    Re : float
        Earth radius [km]
    az : float
        azimuth [degree, N]
    el : float
        elevation [degree]
    """
    taz = numpy.radians(az)
    tel = numpy.radians(el)
    
    # In this transformation x is east, y is north and z is up
    if not inverse:
        # Calculate deviation from vertical (in radians)
        (geocLat, geocLon, Re) = geodToGeoc(lat, lon)
        devH = numpy.radians(lat - geocLat)
        # Calculate cartesian coordinated in local system
        kxGD = numpy.cos(tel) * numpy.sin(taz)
        kyGD = numpy.cos(tel) * numpy.cos(taz)
        kzGD = numpy.sin(tel)
        # Now rotate system about the x axis to align local vertical vector
        # with Earth radial vector
        kxGC = kxGD
        kyGC = kyGD * numpy.cos(devH) + kzGD * numpy.sin(devH)
        kzGC = -kyGD * numpy.sin(devH) + kzGD * numpy.cos(devH)
        # Finally calculate the new azimuth and elevation in the geocentric
        # frame
        azOut = numpy.degrees(numpy.arctan2(kxGC, kyGC))
        elOut = numpy.degrees(numpy.arctan(kzGC / numpy.sqrt(kxGC**2 + kyGC**2)))
        latOut = geocLat
        lonOut = geocLon
    else:
        # Calculate deviation from vertical (in radians)
        (geodLat, geodLon, Re) = geodToGeoc(lat, lon, inverse=True)
        devH = numpy.radians(geodLat - lat)
        # Calculate cartesian coordinated in geocentric system
        kxGC = numpy.cos(tel) * numpy.sin(taz)
        kyGC = numpy.cos(tel) * numpy.cos(taz)
        kzGC = numpy.sin(tel)
        # Now rotate system about the x axis to align local vertical vector
        # with Earth radial vector
        kxGD = kxGC
        kyGD = kyGC * numpy.cos(-devH) + kzGC * numpy.sin(-devH)
        kzGD = -kyGC * numpy.sin(-devH) + kzGC * numpy.cos(-devH)
        # Finally calculate the new azimuth and elevation in the geocentric
        # frame
        azOut = numpy.degrees(numpy.arctan2(kxGD, kyGD))
        elOut = numpy.degrees(numpy.arctan(kzGD / numpy.sqrt(kxGD**2 + kyGD**2)))
        latOut = geodLat
        lonOut = geodLon
    
    return latOut, lonOut, Re, azOut, elOut


def gspToGcar(xin, yin, zin, inverse=False):
    """Converts a position from global spherical (geocentric) to global
    cartesian (and vice-versa).
    Parameters
    ----------
    xin : float
        latitude [degree] or global cartesian X [km]
    yin : float
        longitude [degree] or global cartesian Y [km]
    zin : float
        distance from center of the Earth [km] or global cartesian Z [km]
    inverse : Optional[bool]
        inverse conversion
    Returns
    -------
    xout : float
        global cartesian X [km] (inverse=False) or latitude [degree]
    yout : float
        global cartesian Y [km] (inverse=False) or longitude [degree]
    zout : float
        global cartesian Z [km] (inverse=False) or distance from the center of
        the Earth [km]
    Notes
    -------
    The global cartesian coordinate system is defined as:
        - origin: center of the Earth
        - x-axis in the equatorial plane and through the prime meridian.
        - z-axis in the direction of the rotational axis and through the North
          pole
    The meaning of the input (x,y,z) depends on the direction of the conversion 
    (to global cartesian or to global spherical).
    """

    if not inverse:
        # Global spherical to global cartesian
        xout = zin * numpy.cos(numpy.radians(xin)) * numpy.cos(numpy.radians(yin))
        yout = zin * numpy.cos(numpy.radians(xin)) * numpy.sin(numpy.radians(yin))
        zout = zin * numpy.sin(numpy.radians(xin))
    else:
        # Calculate latitude (xout), longitude (yout) and distance from center
        # of the Earth (zout)
        zout = numpy.sqrt(xin**2 + yin**2 + zin**2)
        xout = numpy.degrees(numpy.arcsin(zin / zout))
        yout = numpy.degrees(numpy.arctan2(yin, xin))
        
    return xout, yout, zout


def gcarToLcar(X, Y, Z, lat, lon, rho , inverse=False):
    """Converts a position from global cartesian to local cartesian
    (or vice-versa).
    Parameters
    ----------
    X : float
        global cartesian X [km] or local cartesian X [km]
    Y : flaot
        global cartesian Y [km] or local cartesian Y [km]
    Z : float
        global cartesian Z [km] or local cartesian Z [km]
    lat : float
        geocentric latitude [degree] of local cartesian system origin
    lon : float
        geocentric longitude [degree] of local cartesian system origin
    rho : float
        distance from center of the Earth [km] of local cartesian system origin
    inverse : Optional[bool]
        inverse conversion
    Returns
    -------
    X : float
        local cartesian X [km] or global cartesian X [km]
    Y : float
        local cartesian Y [km] or global cartesian Y [km]
    Z : float
        local cartesian Z [km] or global cartesian Z [km]
    Notes
    -------
    The global cartesian coordinate system is defined as:
        - origin: center of the Earth
        - Z axis in the direction of the rotational axis and through the North
          pole
        - X axis in the equatorial plane and through the prime meridian.
    The local cartesian coordinate system is defined as:
        - origin: local position
        - X: East
        - Y: North
        - Z: up
    The meaning of the input (X,Y,Z) depends on the direction of the conversion 
    (to global cartesian or to global spherical).
    """
    # First get global cartesian coordinates of local origin
    (goX, goY, goZ) = gspToGcar(lat, lon, rho)
    
    if not inverse:
        # Translate global position to local origin
        tx = X - goX
        ty = Y - goY
        tz = Z - goZ
        # Then, rotate about global-Z to get local-X pointing eastward
        rot = -numpy.radians(lon + 90.0)
        sx = tx * numpy.cos(rot) - ty * numpy.sin(rot)
        sy = tx * numpy.sin(rot) + ty * numpy.cos(rot)
        sz = tz
        # Finally, rotate about X axis to align Z with upward direction
        rot = -numpy.radians(90.0 - lat)
        xOut = sx
        yOut = sy * numpy.cos(rot) - sz * numpy.sin(rot)
        zOut = sy * numpy.sin(rot) + sz * numpy.cos(rot)
    else:
        # First rotate about X axis to align Z with Earth rotational axis
        # direction
        rot = numpy.radians(90.0 - lat)
        sx = X
        sy = Y * numpy.cos(rot) - Z * numpy.sin(rot)
        sz = Y * numpy.sin(rot) + Z * numpy.cos(rot)
        # Rotate about global-Z to get global-X pointing to the prime meridian
        rot = numpy.radians(lon + 90.)
        xOut = sx * numpy.cos(rot) - sy * numpy.sin(rot)
        yOut = sx * numpy.sin(rot) + sy * numpy.cos(rot)
        zOut = sz
        # Finally, translate local position to global origin
        xOut = xOut + goX
        yOut = yOut + goY
        zOut = zOut + goZ
    
    return xOut, yOut, zOut


def lspToLcar(X, Y, Z, inverse=False):
    """Convert a position from local spherical to local cartesian, or vice-versa
    Parameters
    ----------
    X : float
        azimuth [degree, N] or local cartesian X [km]
    Y : float
        elevation [degree] or local cartesian Y [km]
    Z : float
        distance origin [km] or local cartesian Z [km]
    inverse : Optional[bool]
        inverse conversion
    Returns
    -------
    X : float
        local cartesian X [km] or azimuth [degree, N]
    Y : float
        local cartesian Y [km] or elevation [degree]
    Z : float
        local cartesian Z [km] or distance from origin [km]
    Notes
    ------
    The local spherical coordinate system is defined as:
        - origin: local position
        - azimuth (with respect to North)
        - Elevation (with respect to horizon)
        - Altitude
    The local cartesian coordinate system is defined as:
        - origin: local position
        - X: East
        - Y: North
        - Z: up
    The meaning of the input (X,Y,Z) depends on the direction of the conversion 
    (to global cartesian or to global spherical).
    """
    if not inverse:
        # local spherical into local cartesian
        r = Z
        el = Y
        az = X
        xOut = r * numpy.cos(numpy.radians(el)) * numpy.sin(numpy.radians(az))
        yOut = r * numpy.cos(numpy.radians(el)) * numpy.cos(numpy.radians(az))
        zOut = r * numpy.sin(numpy.radians(el))
    else:
        # local cartesian into local spherical
        r = numpy.sqrt(X**2 + Y**2 + Z**2)
        el = numpy.degrees(numpy.arcsin(Z / r))
        az = numpy.degrees(numpy.arctan2(X, Y))
        xOut = az
        yOut = el
        zOut = r
    
    return xOut, yOut, zOut


# *************************************************************
def calcDistPnt(origLat, origLon, origAlt, dist=None, el=None, az=None,
                distLat=None, distLon=None, distAlt=None):
    """Calculate position of a distant point through one of several methods 
    Parameters
    ----------
    origLat : float
        geographic latitude of point of origin [degree]
    origLon : float
        geographic longitude of point of origin [degree]
    origAlt : float
        altitude of point of origin [km]
    dist : Optional[float]
        distance to point [km]
    el : Optional[float]
        elevation [degree]
    az : Optional[float]
        azimuth [degree]
    distLat : Optional[float]
        latitude [degree] of distant point
    distLon : Optional[float]
        longitude [degree] of distant point
    distAlt : Optional[float]
        altitide [km] of distant point
    Returns
    -------
    dictOut : (dict of floats)
        A dictionary containing the information about the origin and remote
        points, as well as their relative positions.  The keys are:
        origLat - origin latitude in degrees,
        origLon - origin longitude in degrees
        origAlt - origin altitude in km
        distLat - distant latitude in degrees
        distLon - distant longitude in degrees 
        distAlt - distant altitude in km
        az - azimuthal angle between origin and distant locations in degrees
        el - elevation angle between origin and distant locations in degrees
        dist - slant distance between origin and distant locaitons in km
        origRe - origin earth radius
        distRe - distant earth radius
    Notes
    -------
    Calculation methods
        - the coordinates and altitude of a distant point given a point of
          origin, distance, azimuth and elevation
        - the coordinates and distance of a distant point given a point of
          origin, altitude, azimuth and elevation
        - the distance, azimuth and elevation between a point of origin and a
          distant point
        - the distance, azimuth between a point of origin and a distant point
          and the altitude of said distant point given a point of origin,
          distant point and elevation angle.
    Input/output is in geodetic coordinates, distances are in km and angles in
    degrees.
    """
    
    # If all the input parameters (keywords) are set to 0, show a warning, and
    # default to fint distance/azimuth/elevation
    if dist is None and el is None and az is None:
        assert None not in [distLat, distLon, distAlt], \
            logging.error('Not enough keywords.')

        # Convert point of origin from geodetic to geocentric
        (gcLat, gcLon, origRe) = geodToGeoc(origLat, origLon)
        # Convert distant point from geodetic to geocentric
        (gcDistLat, gcDistLon, distRe) = geodToGeoc(distLat, distLon)
        # convert distant point from geocentric to global cartesian
        (pX, pY, pZ) = gspToGcar(gcDistLat, gcDistLon, distRe + distAlt)
        # convert pointing direction from global cartesian to local cartesian
        (dX, dY, dZ) = gcarToLcar(pX, pY, pZ, gcLat, gcLon, origRe+origAlt)
        # convert pointing direction from local cartesian to local spherical
        (gaz, gel, rho) = lspToLcar(dX, dY, dZ, inverse=True)
        # convert pointing azimuth and elevation to geodetic
        (lat, lon, Re, az, el) = geodToGeocAzEl(gcLat, gcLon, gaz, gel,
                                                inverse=True)
        dist = numpy.sqrt(dX**2 + dY**2 + dZ**2)

    elif distLat is None and distLon is None and distAlt is None:
        assert None not in [dist, el, az], logging.error('Not enough keywords.')

        # convert pointing azimuth and elevation to geocentric
        (gcLat, gcLon, origRe, gaz, gel) = geodToGeocAzEl(origLat, origLon, az,
                                                          el)
        # convert pointing direction from local spherical to local cartesian
        (pX, pY, pZ) = lspToLcar(gaz, gel, dist)
        # convert pointing direction from local cartesian to global cartesian
        (dX, dY, dZ) = gcarToLcar(pX, pY, pZ, gcLat, gcLon, origRe + origAlt,
                                  inverse=True)
        # Convert distant point from global cartesian to geocentric
        (gcDistLat, gcDistLon, rho) = gspToGcar(dX, dY, dZ, inverse=True)
        # Convert distant point from geocentric to geodetic
        (distLat, distLon, Re) = geodToGeoc(gcDistLat, gcDistLon, inverse=True)
        distAlt = rho - Re
        distRe = Re

    elif dist is None and distAlt is None and az is None:
        assert None not in [distLat, distLon, el], \
            logging.error('Not enough keywords')

        # Convert point of origin from geodetic to geocentric
        (gcLat, gcLon, origRe) = geodToGeoc(origLat, origLon)
        Dref = origRe + origAlt
        
        # convert point of origin from geocentric to global cartesian
        (pX, pY, pZ) = gspToGcar(gcLat, gcLon, Dref)
        
        # Convert distant point from geodetic to geocentric
        (gcDistLat, gcDistLon, distRe) = geodToGeoc(distLat, distLon)
        
        # convert distant point from geocentric to global cartesian
        (pdX, pdY, pdZ) = gspToGcar(gcDistLat, gcDistLon, Dref)
        
        # convert pointing direction from global cartesian to local cartesian
        (dX, dY, dZ) = gcarToLcar(pdX, pdY, pdZ, gcLat, gcLon, Dref)
        
        # convert pointing direction from local cartesian to local spherical
        (gaz, gel, rho) = lspToLcar(dX, dY, dZ, inverse=True)
        
        # convert pointing azimuth and elevation to geodetic
        (lat, lon, Re, az, el) = geodToGeocAzEl(gcLat, gcLon, gaz, gel,
                                                inverse=True)
        
        # convert pointing azimuth and elevation to geocentric
        (gcLat, gcLon, origRe, gaz, gel) = geodToGeocAzEl(origLat, origLon, az,
                                                          el)
        # calculate altitude and distance
        theta = numpy.arccos((pdX * pX + pdY * pY + pdZ * pZ) / Dref**2)
        distAlt = Dref * (numpy.cos(numpy.radians(gel)) /
                          numpy.cos(theta + numpy.radians(gel)) - 1.0)
        distAlt -= distRe - origRe
        dist = Dref * numpy.sin(theta) / numpy.cos(theta + numpy.radians(gel))

    elif distLat is None and distLon is None and dist is None:
        assert None not in [distAlt, el, az], \
            logging.error('Not enough keywords')

        # convert pointing azimuth and elevation to geocentric
        (gcLat, gcLon, origRe, gaz, gel) = geodToGeocAzEl(origLat, origLon, az,
                                                          el)
        
        # Calculate angles
        alpha = numpy.arcsin((origRe + origAlt) * numpy.cos(numpy.radians(gel)) /
                          (origRe + distAlt))
        theta = numpy.pi / 2.0 - alpha - numpy.radians(gel)
        
        # calculate distance
        dist = numpy.sqrt((origRe + origAlt)**2 + (origRe + distAlt)**2 - 2.0 *
                       (origRe + distAlt) * (origRe + origAlt) * numpy.cos(theta))
        
        # convert pointing direction from local spherical to local cartesian
        (pX, pY, pZ) = lspToLcar(gaz, gel, dist)
        
        # convert pointing direction from local cartesian to global cartesian
        (dX, dY, dZ) = gcarToLcar(pX, pY, pZ, gcLat, gcLon, origRe+origAlt,
                                  inverse=True)
        
        # Convert distant point from global cartesian to geocentric
        (gcDistLat, gcDistLon, rho) = gspToGcar(dX, dY, dZ, inverse=True)
        
        # Convert distant point from geocentric to geodetic
        (distLat, distLon, distRe) = geodToGeoc(gcDistLat, gcDistLon,
                                                inverse=True)
        distAlt = rho - distRe
    else:
        return dict()
    
    # Fill output dictionary
    dictOut = {'origLat': origLat, 'origLon': origLon, 'origAlt': origAlt,
               'distLat': distLat, 'distLon': distLon, 'distAlt': distAlt,
               'az': az, 'el': el, 'dist': dist, 'origRe': origRe,
               'distRe': distRe}
    
    return dictOut


def greatCircleMove(origLat, origLon, dist, az, alt=0.0, Re=6371.0):
    """Calculates the coordinates of an end point along a great circle path 
    given the original coordinates, distance, azimuth, and altitude.
    Parameters
    ----------
    origLat : float
        latitude [degree]
    origLon : float
        longitude [degree]
    dist : float
        distance [km]
    az : float
        azimuth [deg]
    alt : Optional[float]
        altitude [km] added to default Re = 6378.1 km (default=0.0)
    Re : Optional[float]
        Earth radius (default=6371.0)
    Returns
    -------
    latitude : (float)
        latitude in degrees
    longitude: (float)
        longitude in degrees
    """
    Re_tot = (Re + alt) * 1.0e3
    dist = dist * 1.0e3
    lat1 = numpy.radians(origLat) 
    lon1 = numpy.radians(origLon)
    az = numpy.radians(az)
    
    lat2 = numpy.arcsin(numpy.sin(lat1) * numpy.cos(dist / Re_tot) +
                     numpy.cos(lat1) * numpy.sin(dist / Re_tot) * numpy.cos(az))
    lon2 = lon1 + numpy.arctan2(numpy.sin(az) * numpy.sin(dist / Re_tot) * numpy.cos(lat1),
                             numpy.cos(dist / Re_tot)
                             - numpy.sin(lat1) * numpy.sin(lat2))

    # Convert everything to numpy arrays to make selective processing easier.
    ret_lat = numpy.degrees(lat2)
    ret_lon = numpy.degrees(lon2)
    
    ret_lat = numpy.array(ret_lat)
    if ret_lat.shape == ():
        ret_lat.shape = (1,)

    ret_lon = numpy.array(ret_lon)
    if ret_lon.shape == ():
        ret_lon.shape = (1,)

    # Put all longitudes on -180 to 180 domain.
    ret_lon = ret_lon % 360.0

    tf = ret_lon > 180.0
    ret_lon[tf] = ret_lon[tf] - 360.0

    return (ret_lat, ret_lon)


def greatCircleAzm(lat1, lon1, lat2, lon2):
    """Calculates the azimuth from the coordinates of a start point to and end
    point along a great circle path.
    Parameters
    ----------
    lat1 : float
        latitude [deg]
    lon1 : float
        longitude [deg]
    lat2 : float
        latitude [deg]
    lon2 : float
        longitude [deg]
    Returns
    -------
    azm : float
        azimuth [deg]
    """
    lat1 = numpy.radians(lat1)
    lon1 = numpy.radians(lon1)
    lat2 = numpy.radians(lat2)
    lon2 = numpy.radians(lon2)
    
    dlon = lon2 - lon1
    y = numpy.sin(dlon) * numpy.cos(lat2)
    x = numpy.cos(lat1) * numpy.sin(lat2) - numpy.sin(lat1) * numpy.cos(lat2) * numpy.cos(dlon)
    azm = numpy.degrees(numpy.arctan2(y,x))

    return azm


def greatCircleDist(lat1, lon1, lat2, lon2):
    """Calculates the distance in radians along a great circle path between two
    points.
    Parameters
    ----------
    lat1 : float
        latitude [deg]
    lon1 : float
        longitude [deg]
    lat2 : float
        latitude [deg]
    lon2 : float
        longitude [deg]
    Returns
    -------
    radDist : float
        distance [radians]
    """
    lat1 = numpy.radians(lat1)
    lon1 = numpy.radians(lon1)
    lat2 = numpy.radians(lat2)
    lon2 = numpy.radians(lon2)

    dlat = (lat2 - lat1) / 2.0
    dlon = (lon2 - lon1) / 2.0
    a = numpy.sin(dlat)**2 + numpy.cos(lat1) * numpy.cos(lat2) * numpy.sin(dlon)**2
    radDist = 2.0 * numpy.arctan2(numpy.sqrt(a), numpy.sqrt(1.0 - a))

    return radDist
