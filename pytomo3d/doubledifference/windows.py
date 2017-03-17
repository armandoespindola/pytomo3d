# -*- coding: utf-8 -*-
"""Helper methods for double difference window operations

:copyright:
   Ridvan Orsvuran (orsvuran@geoazur.unice.fr), 2017
:license:
    GNU Lesser General Public License, version 3 (LGPLv3)
    (http://www.gnu.org/licenses/lgpl-3.0.en.html)
"""

from obspy import UTCDateTime


def component_based_windows_data(windows_data):
    """Convert station based windows data to component based one

    Converts
        "net.sta": {"net.sta.loc.comp": [{}, {}], ...}, ...
    into
        "comp": {"net.sta.loc.comp:0":{}, ...}, ...

    :param windows_data: station based windows data
    :type windows_data: dict
    :returns: component based windows data
    :rtype: dict

    """
    all_windows = {}
    for station, comps in windows_data.iteritems():
        for comp, windows in comps.iteritems():
            c = comp[-1]
            if c not in all_windows:
                all_windows[c] = {}
            for i, window in enumerate(windows):
                name = ":".join([comp, str(i)])
                all_windows[c][name] = window

    return all_windows


def convert_to_sta_based_windows(comp_based_windows):
    """Convert component based windows data to station based one

    :param comp_based_windows: component based windows
    :type comp_based_windows: dict
    :returns: station based windows data
    :rtype: dict

    """
    windows_data = {}
    for comp, comp_windows in comp_based_windows.iteritems():
        for wname, wdata in comp_windows.iteritems():
            compname = wname.split(":")[0]
            staname = ".".join(compname.split(".")[:2])
            # if there is a data for this station
            if staname in windows_data:
                # If there is already a window for this component
                if compname in windows_data[staname]:
                    windows_data[staname][compname].append(wdata)
                else:
                    windows_data[staname][compname] = [wdata]
            else:
                windows_data[staname] = {
                    compname: [wdata]
                }
    return windows_data


def get_windowed_trace(trace, windows_data, winname):
    """Window the trace

    :param trace: trace to window
    :type trace: obspy.Trace
    :param windows_data: windows dict with window names as keys
    :type windows_data: dict
    :param winname: window name (net.sta.loc.comp:winid)
    :type winname: str
    :returns: windowed trace data
    :rtype: numpy.array

    """
    window = windows_data[winname]
    windowed_trace = trace.slice(UTCDateTime(window["absolute_starttime"]),
                                 UTCDateTime(window["absolute_endtime"]))
    windowed_trace.taper(max_percentage=0.05)
    return windowed_trace.data


def get_adj_window(windows_data, winname):
    """Get window for using in pyadjoint

    :param windows_data: windows dict with window names as keys
    :type windows_data: dict
    :param winname: window name (net.sta.loc.comp:winid)
    :type winname: str
    :returns: windows list
    :rtype: list

    """
    window = windows_data[winname]
    return [[window["relative_starttime"],
             window["relative_endtime"]]]


def filter_paired_windows(windows, pairs):
    """Filter windows using pairs data

    :param windows: component based windows data
    :param pairs: pairs data
    :returns: paired_windows, single_windows
    :rtype: tuple

    """
    paired_windows = {}
    single_windows = {}
    for comp, comp_windows in windows.iteritems():
        paired_windows[comp] = {}
        single_windows[comp] = {}
        for wname, wdata in comp_windows.iteritems():
            is_paired = False
            for pair in pairs[comp]:
                if pair["window_id_i"] == wname or \
                   pair["window_id_j"] == wname:
                    is_paired = True
                    break
            if is_paired:
                paired_windows[comp][wname] = wdata
            else:
                single_windows[comp][wname] = wdata
    return paired_windows, single_windows
