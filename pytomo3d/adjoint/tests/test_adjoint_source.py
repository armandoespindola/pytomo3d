import os
import inspect
import json
from obspy import read, Stream
from pyflex.window import Window
import pytomo3d.adjoint.adjoint_source as adj
import pytomo3d.adjoint.io as adj_io
import pytest
import matplotlib.pyplot as plt
# import pyadjoint.adjoint_source


def _upper_level(path, nlevel=4):
    """
    Go the nlevel dir up
    """
    for i in range(nlevel):
        path = os.path.dirname(path)
    return path


# Most generic way to get the data folder path.
TESTBASE_DIR = _upper_level(os.path.abspath(
    inspect.getfile(inspect.currentframe())), 4)
DATA_DIR = os.path.join(TESTBASE_DIR, "tests", "data")

obsfile = os.path.join(DATA_DIR, "proc", "IU.KBL.obs.proc.mseed")
synfile = os.path.join(DATA_DIR, "proc", "IU.KBL.syn.proc.mseed")
winfile = os.path.join(DATA_DIR, "window", "IU.KBL..BHR.window.json")


@pytest.fixture
def config_mt():
    config_file = os.path.join(DATA_DIR, "adjoint",
                               "multitaper.adjoint.config.yaml")
    return adj_io.load_adjoint_config_yaml(config_file)


@pytest.fixture
def obs_tr():
    return read(obsfile).select(channel="*R")[0]


@pytest.fixture
def syn_tr():
    return read(synfile).select(channel="*R")[0]


@pytest.fixture
def win_tr():
    with open(winfile) as fh:
        wins_json = json.load(fh)
    windows = []
    for _win in wins_json:
        windows.append(Window._load_from_json_content(_win))

    return windows


def test_calculate_adjsrc_on_trace_raises_if_obs_is_not_trace(
        syn_tr, win_tr, config_mt):
    obs_tr = []
    with pytest.raises(ValueError):
        adj.calculate_adjsrc_on_trace(obs_tr, syn_tr, win_tr, config_mt,
                                      adj_src_type="multitaper_misfit")


def test_calculate_adrjrc_on_trace_raises_if_syn_is_not_trace(
        obs_tr, win_tr, config_mt):
    syn_tr = []
    with pytest.raises(ValueError):
        adj.calculate_adjsrc_on_trace(obs_tr, syn_tr, win_tr, config_mt,
                                      adj_src_type="multitaper_misfit")


def test_calculate_adjsrc_on_trace_raises_if_config_is_not_config(
        obs_tr, syn_tr, win_tr):
    config_mt = []
    with pytest.raises(ValueError):
        adj.calculate_adjsrc_on_trace(obs_tr, syn_tr, win_tr, config_mt,
                                      adj_src_type="multitaper_misfit")


def test_calculate_adjsrc_on_trace_raises_bad_windows_shape(
        obs_tr, syn_tr, config_mt):
    win_tr = []
    with pytest.raises(ValueError):
        adj.calculate_adjsrc_on_trace(obs_tr, syn_tr, win_tr, config_mt,
                                      adj_src_type="multitaper_misfit")


def test_calculate_adjsrc_on_trace_figure_mode_none_figure_dir(
        obs_tr, syn_tr, win_tr, config_mt):
    plt.switch_backend('agg')
    adjsrc = adj.calculate_adjsrc_on_trace(
        obs_tr, syn_tr, win_tr, config_mt, adj_src_type="multitaper_misfit",
        figure_mode=True)
    assert adjsrc


# def test_calculate_adjsrc_on_trace_waveform_misfit_produces_adjsrc():
#    obs, syn, win_time = setup_calculate_adjsrc_on_trace_args()
#    config = load_config_waveform()

#    adjsrc = adj.calculate_adjsrc_on_trace(
#        obs, syn, win_time, config, adj_src_type="waveform_misfit",
#        adjoint_src_flag=True, figure_mode=False)
#    assert adjsrc


def test_calculate_adjsrc_on_trace_multitaper_misfit_produces_adjsrc(
        obs_tr, syn_tr, win_tr, config_mt):
    adjsrc = adj.calculate_adjsrc_on_trace(
        obs_tr, syn_tr, win_tr, config_mt, adj_src_type="multitaper_misfit",
        adjoint_src_flag=True, figure_mode=False)
    assert adjsrc


# def test_calculate_adjsrc_on_trace_traveltime_misfit_produces_adjsrc():
#    obs, syn, win_time = setup_calculate_adjsrc_on_trace_args()
#    config = load_config_traveltime()
#
#    adjsrc = adj.calculate_adjsrc_on_trace(
#        obs, syn, win_time, config, adj_src_type="cc_traveltime_misfit",
#        adjoint_src_flag=True, figure_mode=False)
#    assert adjsrc


@pytest.fixture
def obs_st():
    return Stream(traces=[read(obsfile).select(channel="*R")[0]])


@pytest.fixture
def syn_st():
    return Stream(traces=[read(synfile).select(channel="*R")[0]])


@pytest.fixture
def win_st(obs_st):
    with open(winfile) as fh:
        wins_json = json.load(fh)

    return {obs_st[0].id: wins_json}


def test_calculate_adjsrc_on_stream_raises_if_obs_is_not_stream(
        syn_st, win_st, config_mt):
    obs_st = []
    with pytest.raises(ValueError):
        adj.calculate_adjsrc_on_stream(obs_st, syn_st, win_st, config_mt,
                                       adj_src_type="multitaper_misfit")


def test_calculate_adjsrc_on_stream_raises_if_syn_is_not_stream(
        obs_st, win_st, config_mt):
    syn_st = []
    with pytest.raises(ValueError):
        adj.calculate_adjsrc_on_stream(obs_st, syn_st, win_st, config_mt,
                                       adj_src_type="multitaper_misfit")


def test_calculate_adjsrc_on_stream_raises_if_config_is_not_config(
        obs_st, syn_st, win_st):
    config_mt = []
    with pytest.raises(ValueError):
        adj.calculate_adjsrc_on_stream(obs_st, syn_st, win_st, config_mt,
                                       adj_src_type="multitaper_misfit")


def test_calculate_adjsrc_on_stream_raises_if_windows_is_empty(
        obs_st, syn_st, config_mt):
    windows = None
    ret = adj.calculate_adjsrc_on_stream(obs_st, syn_st, windows, config_mt,
                                         adj_src_type="multitaper_misfit")
    assert ret is None
    windows = {}
    ret = adj.calculate_adjsrc_on_stream(obs_st, syn_st, windows, config_mt,
                                         adj_src_type="multitaper_misfit")
    assert ret is None


# def test_calculate_adjsrc_on_stream_multitaper_misfit_produces_adjsrc():
#    obs, syn, windows = setup_calculate_adjsrc_on_stream_args()
#    config = load_config_traveltime()
#
#    adjsrc = adj.calculate_adjsrc_on_stream(
#        obs, syn, windows, config, adj_src_type="multitaper_misfit",
#        adjoint_src_flag=True, figure_mode=False)
#    assert adjsrc


# def test_calculate_adjsrc_on_stream_waveform_misfit_produces_adjsrc():
#    obs, syn, windows = setup_calculate_adjsrc_on_stream_args()
#    config = load_config_traveltime()
#
#    adjsrc = adj.calculate_adjsrc_on_stream(
#        obs, syn, windows, config, adj_src_type="waveform_misfit",
#        adjoint_src_flag=True, figure_mode=False)
#    assert adjsrc


# def test_calculate_adjsrc_on_stream_traveltime_misfit_produces_adjsrc():
#    obs, syn, windows = setup_calculate_adjsrc_on_stream_args()
#    config = load_config_traveltime()
#
#    adjsrc = adj.calculate_adjsrc_on_stream(
#        obs, syn, windows, config, adj_src_type="cc_traveltime_misfit",
#        adjoint_src_flag=True, figure_mode=False)
#    assert adjsrc


def test_measure_adjoint_on_stream():
    pass
