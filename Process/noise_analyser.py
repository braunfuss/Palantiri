
from pyrocko.client import catalog

import logging
import numpy as num
from pyrocko.guts import Int, Bool, Float, String
from pyrocko.gf.meta import OutOfBounds
from grond.dataset import NotFound

logger = logging.getLogger('grond.analysers.NoiseAnalyser')


guts_prefix = 'grond'


def get_phase_arrival_time(engine, source, station, wavename, store_id):
    """
    Get arrival time from Greens Function store for respective
    :class:`pyrocko.gf.seismosizer.Target`,
    :class:`pyrocko.gf.meta.Location` pair.

    Parameters
    ----------
    engine : :class:`pyrocko.gf.seismosizer.LocalEngine`
    source : :class:`pyrocko.gf.meta.Location`
        can be therefore :class:`pyrocko.gf.seismosizer.Source` or
        :class:`pyrocko.model.Event`
    target : :class:`pyrocko.gf.seismosizer.Target`
    wavename : string
        of the tabulated phase_def that determines the phase arrival

    Returns
    -------
    scalar, float of the arrival time of the wave
    """
    store = engine.get_store(store_id)
    dist = station.distance_to(source)
    depth = source.depth
    return store.t(wavename, (depth, dist)) + source.time


def seismic_noise_variance(traces, engine, event, stations,
                           nwindows, pre_event_noise_duration,
                           check_events, phase_def, store_id):
    '''
    Calculate variance of noise (half an hour) before P-Phase onset, and check
    for other events interfering

    Parameters
    ----------
    data_traces : list
        of :class:`pyrocko.trace.Trace` containing observed data
    engine : :class:`pyrocko.gf.seismosizer.LocalEngine`
        processing object for synthetics calculation
    event : :class:`pyrocko.meta.Event`
        reference event from catalog
    targets : list
        of :class:`pyrocko.gf.seismosizer.Targets`
    nwindows : integer
        if given, windowed analysis of noise, else
        variance is calculated on the entire pre-event
        noise
    phase_def : :class:'pyrocko.gf.Timing'
    arrivals : list
        of :class'pyrocko.gf.Timing' arrivals of waveforms
        at station

    Returns
    -------
    :class:`numpy.ndarray`
    '''

    var_ds = []
    global_cmt_catalog = catalog.GlobalCMT()
    var_ds = []
    ev_ws = []
    for tr, station in zip(traces, stations):
        stat_w = 1.

        if tr is None:
            var_ds.append(0.)
            ev_ws.append(0.)
        else:

            arrival_time = get_phase_arrival_time(
                engine=engine, source=event,
                station=station, wavename=phase_def,
                store_id=store_id)
            if check_events:
                events = global_cmt_catalog.get_events(
                    time_range=(
                        arrival_time-pre_event_noise_duration-50.*60.,
                        arrival_time),
                    magmin=5.,)
                for ev in events:
                    try:
                        arrival_time_pre = get_phase_arrival_time(
                            engine=engine,
                            source=ev,
                            station=station,
                            wavename=phase_def,
                            store_id=store_id)

                        if arrival_time_pre > arrival_time \
                                - pre_event_noise_duration \
                                and arrival_time_pre < arrival_time:

                            stat_w = 0.
                            logger.info(
                                'Noise analyser found event %s phase onset of '
                                '%s for %s' % (
                                    ev.name, phase_def, station))

                        if arrival_time_pre > arrival_time-30.*60.\
                                and arrival_time_pre < arrival_time - \
                                pre_event_noise_duration:
                            stat_w *= 0.5
                            logger.info(
                                'Noise analyser found event %s possibly '
                                'contaminating the noise' % ev.name)

                            # this should be magnitude dependent
                    except Exception:
                        pass
            ev_ws.append(stat_w)

            if nwindows == 1:
                vtrace_var = num.nanvar(tr.ydata)
                var_ds.append(vtrace_var)
            else:
                win = arrival_time - (arrival_time -
                                      pre_event_noise_duration)
                win_len = win/nwindows
                v_traces_w = []
                for i in range(0, nwindows):
                    vtrace_w = tr.chop(
                        tmin=win+win_len*i,
                        tmax=arrival_time+win_len*i+1,
                        inplace=False)
                    v_traces_w.append(vtrace_w.ydata)
                v_traces_w = num.nanmean(v_traces_w)
                var_ds.append(v_traces_w)
    var_ds = num.array(var_ds, dtype=num.float)
    ev_ws = num.array(ev_ws, dtype=num.float)
    return var_ds, ev_ws




def analyse(traces, engine, event, stations,
            pre_event_noise_duration, store_id, nwindows=1,
            check_events='True', phase_def='P'):

    if pre_event_noise_duration == 0:
        return

    var_ds, ev_ws = seismic_noise_variance(
        traces, engine, event, stations,
        nwindows, pre_event_noise_duration,
        check_events, phase_def, store_id)
    norm_noise = num.nanmedian(var_ds)
    if norm_noise == 0:
            logger.info('Noise Analyser returned a weight of 0 for \
                        all stations')

    weights = norm_noise/var_ds
    if check_events:
        weights = weights*ev_ws
    return weights
