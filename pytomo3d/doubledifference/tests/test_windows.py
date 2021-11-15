# -*- coding: utf-8 -*-

from pytomo3d.doubledifference.windows import (component_based_windows_data,
                                               filter_paired_windows,
                                               convert_to_sta_based_windows,
                                               get_adj_window)


import pytest


@pytest.fixture
def windows_data():
    return {
        "AA.AAA": {
            "AA.AAA.00.BHZ": [
                {"relative_starttime": 0, "relative_endtime": 100},
                {"relative_starttime": 100, "relative_endtime": 200}
            ],
            "AA.AAA.00.BHE": [],
            "AA.AAA.00.BHN": [
                {"relative_starttime": 0, "relative_endtime": 100}
            ]
        },
        "AA.BBB": {
            "AA.BBB.00.BHZ": [
                {"relative_starttime": 0, "relative_endtime": 100}
            ]
        }
    }


@pytest.fixture
def comp_win_data(windows_data):
    return component_based_windows_data(windows_data)


def test_comp_based_windows_data(windows_data):
    comp_win_data = component_based_windows_data(windows_data)
    assert len(comp_win_data) == 3

    for comp in "ZEN":
        assert comp in comp_win_data

    z_keys = ["AA.AAA.00.BHZ:0", "AA.AAA.00.BHZ:1", "AA.BBB.00.BHZ:0"]
    n_keys = ["AA.AAA.00.BHN:0"]

    for z_key in z_keys:
        assert z_key in comp_win_data["Z"]
    for n_key in n_keys:
        assert n_key in comp_win_data["N"]


def test_filter_paired_windows(comp_win_data):
    pairs = {"Z": [{"window_id_i": "AA.AAA.00.BHZ:0",
                    "window_id_j": "AA.BBB.00.BHZ:0"}],
             "E": [],
             "N": []}
    paired, single = filter_paired_windows(comp_win_data,
                                           pairs)
    assert len(paired["Z"]) == 2
    assert len(paired["E"]) == 0
    assert len(paired["N"]) == 0

    assert len(single["Z"]) == 1
    assert len(single["E"]) == 0
    assert len(single["N"]) == 1


def test_get_adj_window(comp_win_data):
    window = get_adj_window(comp_win_data["Z"], "AA.AAA.00.BHZ:1")
    assert len(window) == 1
    assert window[0][0] == 100
    assert window[0][1] == 200


def test_convert_to_sta_based_windows(windows_data):
    comp_win_data = component_based_windows_data(windows_data)
    converted = convert_to_sta_based_windows(comp_win_data)

    # Test equality of the whole structure
    for sta1, sta2 in zip(sorted(windows_data.keys()),
                          sorted(converted.keys())):
        assert sta1 == sta2
