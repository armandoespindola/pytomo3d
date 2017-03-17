#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Utility functions for double difference

:copyright:
   Ridvan Orsvuran (orsvuran@geoazur.unice.fr), 2017
:license:
    GNU Lesser General Public License, version 3 (LGPLv3)
    (http://www.gnu.org/licenses/lgpl-3.0.en.html)
"""


from obspy.geodetics import locations2degrees


def get_dist_in_km(lat_i, lon_i, lat_j, lon_j):
    """Get distance between to points in kilometers

    :param lat_i: Latitude of first point
    :type lat_i: float
    :param lon_i: Longitude of first point
    :type lon_i: float
    :param lat_j: Latitude of first point
    :type lat_j: float
    :param lon_j: Longitude of first point
    :type lon_j: float
    :returns: distance in km
    :rtype: float

    """
    dist_deg = locations2degrees(lat_i, lon_i, lat_j, lon_j)
    return dist_deg*111


def deconstruct_winname(winname):
    """Deconstruct window name to important bits

    :param winname: window name string (e.g. net.sta.loc.comp:winid)
    :type winname: str
    :returns: (net.sta, net.sta.loc.comp, winid)
    :rtype: tuple

    """
    compname, win_id = winname.split(":")
    staname = ".".join(compname.split(".")[:2])
    return staname, compname, int(win_id)


def get_stanames_of_pair(pair):
    """Extract station names from pair object

    :param pair: pair dict which contains window_id_i and window_id_j keys.
    :type pair: dict
    :returns: station names of the pair
    :rtype: tuple

    """
    sta_i, _, _ = deconstruct_winname(pair["window_id_i"])
    sta_j, _, _ = deconstruct_winname(pair["window_id_j"])
    return sta_i, sta_j
