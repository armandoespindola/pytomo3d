# -*- coding: utf-8 -*-

from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals


from pytomo3d.doubledifference.pairing import (create_all_pairs,
                                               close_pairs,
                                               similar_pairs,
                                               find_pairs)
from pytomo3d.doubledifference.windows import component_based_windows_data

import numpy as np

import pytest


@pytest.fixture
def windows():
    windows_data = {
        "AA.AAA": {
            "AA.AAA.00.BHZ": [
                {"absolute_starttime": "1970-01-01T00:00:05.000000Z",
                 "absolute_endtime": "1970-01-01T00:00:15.000000Z"},
                {"absolute_starttime": "1970-01-01T00:00:55.000000Z",
                 "absolute_endtime": "1970-01-01T00:01:05.000000Z"}
            ],
            "AA.AAA.00.BHE": [],
            "AA.AAA.00.BHN": [
                {"absolute_starttime": "1970-01-01T00:00:05.000000Z",
                 "absolute_endtime": "1970-01-01T00:00:15.000000Z"}
            ]
        },
        "AA.BBB": {
            "AA.BBB.00.BHZ": [
                {"absolute_starttime": "1970-01-01T00:00:05.000000Z",
                 "absolute_endtime": "1970-01-01T00:00:15.000000Z"}
            ]
        },
        "AA.CCC": {
            "AA.CCC.00.BHZ": [
                {"absolute_starttime": "1970-01-01T00:00:55.000000Z",
                 "absolute_endtime": "1970-01-01T00:01:05.000000Z"}
            ]
        }
    }
    return component_based_windows_data(windows_data)


@pytest.fixture
def locations():
    locs = {
        "AA.AAA.00.BHZ": {
            "latitude": 41.00892,
            "longitude": 28.97644
        },
        "AA.AAA.00.BHN": {
            "latitude": 41.00892,
            "longitude": 28.97644
        },
        "AA.BBB.00.BHZ": {
            "latitude": 40.19041,
            "longitude": 29.06158
        },
        "AA.CCC.00.BHZ": {
            "latitude": 40.07179,
            "longitude": 29.50928
        }
    }
    return locs


@pytest.fixture
def windowed_data():
    sample_data = {}
    times = np.linspace(-50, 50, 100)
    data = [
        np.exp(-times),
        np.diff(np.exp(-times)),
        np.exp(-(times-10)),
        np.exp(-(times-30))
    ]
    all_windows = windows()
    for window, window_data in zip(all_windows["Z"].keys(), data):
        sample_data[window] = window_data
    return sample_data


@pytest.fixture
def trace_data():
    from obspy import Trace
    from obspy import UTCDateTime

    starttime = UTCDateTime(0)
    dt = 0.1
    stats = {"sampling_rate": 1/dt,
             "starttime": starttime}
    times = np.arange(0, 100, dt)
    traces = {}

    data = np.exp(-(times-10)**2)
    data[:-1] += np.diff(np.exp(-(times-60)**2))
    traces["AA.AAA.00.BHZ"] = Trace(data, stats)

    data = np.exp(-(times-12)**2)
    traces["AA.BBB.00.BHZ"] = Trace(data, stats)

    data = np.zeros(len(times))
    data[:-1] += np.diff(np.exp(-(times-63)**2))
    traces["AA.CCC.00.BHZ"] = Trace(data, stats)

    return traces


def test_create_all_pairs(windows):
    all_pairs = create_all_pairs(windows)
    assert len(all_pairs["Z"]) == 6
    assert len(all_pairs["E"]) == 0
    assert len(all_pairs["N"]) == 0


def test_close_pairs(windows, locations):
    all_pairs = create_all_pairs(windows)
    pairs = close_pairs(locations, threshold=100,
                        old_pairs=all_pairs)
    assert len(pairs["Z"]) == 3


def test_similar_pairs(windows, windowed_data):
    all_pairs = create_all_pairs(windows)
    pairs = similar_pairs(windowed_data,
                          threshold=0.90,
                          old_pairs=all_pairs)
    assert len(pairs["Z"]) == 3


def test_find_pairs_both_restriction_on(trace_data, windows, locations):
    pairs = find_pairs(windows, trace_data, locations,
                       closeness_on=True, closeness_threshold=100,
                       similarity_on=True, similarity_threshold=0.9)
    # only pair is AAA:0 and BBB:0
    assert len(pairs["Z"]) == 1
    assert len(pairs["E"]) == 0
    assert len(pairs["N"]) == 0


def test_find_pairs_both_restriction_off(trace_data, windows, locations):
    pairs = find_pairs(windows, trace_data, locations,
                       closeness_on=False, closeness_threshold=100,
                       similarity_on=False, similarity_threshold=0.9)
    # all possible pairs (4C2)
    assert len(pairs["Z"]) == 6
    assert len(pairs["E"]) == 0
    assert len(pairs["N"]) == 0


def test_find_pairs_only_closeness_on(trace_data, windows, locations):
    pairs = find_pairs(windows, trace_data, locations,
                       closeness_on=True, closeness_threshold=100,
                       similarity_on=False, similarity_threshold=0.9)
    # pairs are AAA:0, BBB:0; AAA:1, BBB:0; BBB:0, CCC:0
    assert len(pairs["Z"]) == 3
    assert len(pairs["E"]) == 0
    assert len(pairs["N"]) == 0


def test_find_pairs_only_similarity_on(trace_data, windows, locations):
    pairs = find_pairs(windows, trace_data, locations,
                       closeness_on=False, closeness_threshold=100,
                       similarity_on=True, similarity_threshold=0.9)
    # pairs are AAA:0, BBB:0; AAA:1; CCC:0
    assert len(pairs["Z"]) == 2
    assert len(pairs["E"]) == 0
    assert len(pairs["N"]) == 0
