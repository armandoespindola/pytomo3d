# -*- coding: utf-8 -*-
"""Methods for double difference pairing

:copyright:
   Ridvan Orsvuran (orsvuran@geoazur.unice.fr), 2017
:license:
    GNU Lesser General Public License, version 3 (LGPLv3)
    (http://www.gnu.org/licenses/lgpl-3.0.en.html)
"""

from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals


from .utils import get_dist_in_km, get_stanames_of_pair, deconstruct_winname
from .windows import get_windowed_trace

from itertools import combinations

import numpy as np
from collections import Counter


def create_all_pairs(windows):
    """Create all possible pairs from windows

    :param windows: component based windows data
    :type windows: dict
    :returns: pairs
    :rtype: dict

    """
    pairs = {}
    for comp, comp_windows in windows.iteritems():
        pairs[comp] = [{"window_id_i": i,
                        "window_id_j": j}
                       for i, j in combinations(comp_windows.keys(), 2)]
    return pairs


def _reduce_pairs(pairing_func, pairing_key, old_pairs=None):
    # General function for pairing
    # It reduces pairs to ones which satisfies pairing_func
    pairs = {}
    for comp, comp_pairs in old_pairs.iteritems():
        pairs[comp] = []
        for pair in comp_pairs:
            sta_i, sta_j = get_stanames_of_pair(pair)
            # Do not pair the station with itself
            if sta_i != sta_j:
                pair_value, is_paired = pairing_func(pair["window_id_i"],
                                                     pair["window_id_j"])
                if is_paired:
                    new_pair = pair.copy()
                    new_pair[pairing_key] = pair_value
                    pairs[comp].append(new_pair)
    return pairs


def close_pairs(locations, threshold, old_pairs):
    """Pair the data using geographical closeness

    :param locations: locations of stations
    :type locations: dict
    :param threshold: threshold distance in km
    :type threshold: float
    :param old_pairs: previous pairs
    :type old_pairs: dict
    :returns: pairs
    :rtype: dict

    """

    def get_comp_name(wname):
        return wname.split(":")[0].replace("BHR", "BHZ").replace("BHT", "BHZ")

    def pair_by_distance(pair_i, pair_j):
        comp_i = get_comp_name(pair_i)
        comp_j = get_comp_name(pair_j)
        distance = get_dist_in_km(locations[comp_i]["latitude"],
                                  locations[comp_i]["longitude"],
                                  locations[comp_j]["latitude"],
                                  locations[comp_j]["longitude"])
        is_paired = distance < threshold
        return distance, is_paired

    return _reduce_pairs(pair_by_distance, "distance", old_pairs)


def similar_pairs(data, threshold, old_pairs):
    """Pair the data using cc similarity

    :param data: windowed data
    :type data: station
    :param threshold: threshold cc value
    :type threshold: float
    :param old_pairs: previous pairs
    :type old_pairs: dict
    :returns: pairs
    :rtype: dict

    """

    def pair_by_similarity(pair_i, pair_j):
        data_i = data[pair_i]
        data_j = data[pair_j]
        corr = np.correlate(data_i, data_j, "full")
        ratio = max(corr)/np.sqrt(sum(data_i**2)*sum(data_j**2))
        # print(np.sqrt(sum(data_i**2)*sum(data_j**2)))
        is_paired = abs(ratio) > threshold
        return ratio, is_paired

    return _reduce_pairs(pair_by_similarity, "cc_similarity", old_pairs)


def phase_pairs(windows, phases, old_pairs):
    """Data
    """

    def pair_by_phase(pair_i, pair_j):
        _, comp_i, _ = deconstruct_winname(pair_i)
        comp = comp_i[-1]
        phases_i = phase_data[comp][pair_i]
        phases_j = phase_data[comp][pair_j]
        # Now look at the phases
        is_paired = False
        paired_phase = ""
        for phase in phases:
            if phase in phases_i and phase in phases_j:
                is_paired = True
                paired_phase = phase
                break
        return paired_phase, is_paired

    phase_data = {}
    for comp, comp_windows in windows.iteritems():
        phase_data[comp] = {}
        for winname, window in comp_windows.iteritems():
            phase_data[comp][winname] = [phase["phase_name"]
                                         for phase in window["phase_arrivals"]]

    return _reduce_pairs(pair_by_phase, "paired_phase", old_pairs)


def _get_windowed_traces(traces, windows, pairs):
    # Get the data for paired windows by windowing the traces
    windowed_traces = {}
    for comp, comp_pairs in pairs.iteritems():
        for winnames in comp_pairs:
            for winname in [winnames["window_id_i"], winnames["window_id_j"]]:
                if winname not in windowed_traces:
                    compname = winname.split(":")[0]
                    windowed_traces[winname] = get_windowed_trace(
                        traces[compname],
                        windows[compname[-1]],
                        winname)
    return windowed_traces


def find_weights(pairs):
    for comp, comp_pairs in pairs.iteritems():
        all_paired_windows = []
        for pair in comp_pairs:
            all_paired_windows.append(pair["window_id_i"])
            all_paired_windows.append(pair["window_id_j"])
        counts = Counter(all_paired_windows)
        for pair in comp_pairs:
            w_i = 1/(counts[pair["window_id_i"]] + 1)
            w_j = 1/(counts[pair["window_id_i"]] + 1)
            # w_i = 1/len(comp_pairs)
            # w_j = 1/len(comp_pairs)
            pair["weight_i"] = w_i
            pair["weight_j"] = w_j


def find_pairs(windows, traces=None, locations=None,
               closeness_on=False, closeness_threshold=250,
               similarity_on=False, similarity_threshold=0.9,
               phase_on=False, phases=[]):
    """Find pairs from windows

    :param windows: component based windows data
    :type windows: dict
    :param traces: traces as a dict
    :type traces: dict
    :param locations: location data as a dict
    :type locations: dict
    :param closeness_on: flag for pairing by distance
    :type closeness_on: bool
    :param closeness_threshold: maximum distance for pairing
    :type closeness_threshold: float
    :param similarity_on: flag for pairing by cc similarity
    :type similarity: bool
    :param similarity_threshold: minimum cc similarity value
    :type similarity_threshold: float
    :returns: pairs
    :rtype: dict

    """
    # Create all possible pairs first, and then eliminate pairs based
    # on conditions.
    pairs = create_all_pairs(windows)

    if closeness_on:
        pairs = close_pairs(locations, closeness_threshold, pairs)

    if phase_on:
        pairs = phase_pairs(windows, phases, pairs)

    if similarity_on:
        windowed_traces = _get_windowed_traces(traces, windows, pairs)
        pairs = similar_pairs(windowed_traces, similarity_threshold, pairs)

    find_weights(pairs)
    return pairs
