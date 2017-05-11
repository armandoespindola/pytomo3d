# -*- coding: utf-8 -*-
"""Methods for double difference adjoint source calculation

:copyright:
   Ridvan Orsvuran (orsvuran@geoazur.unice.fr), 2017
:license:
    GNU Lesser General Public License, version 3 (LGPLv3)
    (http://www.gnu.org/licenses/lgpl-3.0.en.html)
"""
from .pairing import get_stanames_of_pair
from .windows import get_adj_window
from pyadjoint import AdjointSource, calculate_adjoint_source_DD

import numpy as np
import obspy


def calculate_adjoint_pair(pair, adj_src_type, config, windows, obsd, synt):
    """Calculate adjoint sources for pair

    :param pair: pair object
    :type pair: dict
    :param adj_src_type: adjoint source type
    :type adj_src_type: str
    :param config: adjoint source config
    :type config: pyadjoint.Config
    :param windows: windows
    :type windows: dict
    :param obsd: obsd data dict
    :type obsd: dict
    :param synt: synt data dict
    :type synt: dict
    :returns: named adjoint sources
    :rtype: dict

    """
    win_id_i = pair["window_id_i"]
    win_id_j = pair["window_id_j"]

    win_i = get_adj_window(windows, win_id_i)
    win_j = get_adj_window(windows, win_id_j)

    sta_i, sta_j = get_stanames_of_pair(pair)

    adj_i, adj_j = calculate_adjoint_source_DD(adj_src_type,
                                               obsd[sta_i], synt[sta_i],
                                               obsd[sta_j], synt[sta_j],
                                               config,
                                               win_i, win_j)

    return {(win_id_i, win_id_j): (adj_i, adj_j)}


def calculate_measure_pair(pair, adj_src_type, config, windows, obsd, synt):
    """Calculate adjoint sources for pair

    :param pair: pair object
    :type pair: dict
    :param adj_src_type: adjoint source type
    :type adj_src_type: str
    :param config: adjoint source config
    :type config: pyadjoint.Config
    :param windows: windows
    :type windows: dict
    :param obsd: obsd data dict
    :type obsd: dict
    :param synt: synt data dict
    :type synt: dict
    :returns: named adjoint sources
    :rtype: dict

    """
    win_id_i = pair["window_id_i"]
    win_id_j = pair["window_id_j"]

    win_i = get_adj_window(windows, win_id_i)
    win_j = get_adj_window(windows, win_id_j)

    sta_i, sta_j = get_stanames_of_pair(pair)

    meas_i, meas_j = calculate_adjoint_source_DD(adj_src_type,
                                                 obsd[sta_i], synt[sta_i],
                                                 obsd[sta_j], synt[sta_j],
                                                 config,
                                                 win_i, win_j,
                                                 adjoint_src=False)

    meas_i = meas_i.measurement[0]
    meas_j = meas_j.measurement[0]

    meas_i.pop("ddt_w", None)
    meas_j.pop("ddt_w", None)

    return {(win_id_i, win_id_j): (meas_i, meas_j)}


def add_adjoint_sources(a, b):
    """Add two adjoint sources

    :param a: first adjoint aource
    :type a: pyadjoint.AdjointSource
    :param b: second adjoint_source
    :type b: pyadjoint.AdjointSource
    :returns: summed adjoint source
    :rtype: pyadjoint.AdjointSource

    """
    c = AdjointSource(a.adj_src_type,
                      a.misfit + b.misfit,
                      a.dt,
                      a.min_period,
                      a.max_period,
                      a.component,
                      a.measurement,
                      a.adjoint_source+b.adjoint_source,
                      a.network,
                      a.station,
                      a.location,
                      a.starttime)

    return c


def multiply_adjoint_source(c, a):
    """Add two adjoint sources

    :param a: adjoint aource
    :type a: pyadjoint.AdjointSource
    :param b: coefficient
    :type b: float
    :returns: multiplied adjoint source
    :rtype: pyadjoint.AdjointSource

    """
    c = AdjointSource(a.adj_src_type,
                      c*a.misfit,
                      a.dt,
                      a.min_period,
                      a.max_period,
                      a.component,
                      a.measurement,
                      c*a.adjoint_source,
                      a.network,
                      a.station,
                      a.location,
                      a.starttime)

    return c


def asdf_adj_to_adjoint(a):
    """Convert adjoint data read from asdf to pyadjoint adjoint source

    :param a: asdf auxiliary data
    :returns: adjoint source
    :rtype: pyadjoint.AdjointSource

    """
    try:
        net, sta = a.parameters["station_id"].split(".")
    except:
        net, sta = a.parameters["station_id"].split("_")

    c = AdjointSource(a.parameters["adjoint_source_type"],
                      a.parameters["misfit"],
                      a.parameters["dt"],
                      a.parameters["min_period"],
                      a.parameters["max_period"],
                      a.parameters["component"],
                      [],
                      np.array(a.data),
                      net,
                      sta,
                      a.parameters["location"],
                      obspy.UTCDateTime(a.parameters["starttime"]))

    loc = {"local_depth": a.parameters["depth_in_m"],
           "latitude": a.parameters["latitude"],
           "longitude": a.parameters["longitude"],
           "elevation": a.parameters["elevation_in_m"]}

    return c, loc
