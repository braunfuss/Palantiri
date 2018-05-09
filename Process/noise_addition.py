from pyrocko.client import catalog
import logging
import numpy as num
from pyrocko.guts import Int, Bool, Float, String
from pyrocko.gf.meta import OutOfBounds
from grond.dataset import NotFound
from pyrocko import io, trace, pile

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


def add_noise(traces, engine, event, stations, store,
              phase_def='P'):
    '''
    Calculate variance of noise (half an hour) before P-Phase onset, and check
    for other events interfering

    Parameters
    ----------
    data_traces : list
        of :class:`pyrocko.trace.Trace` containing observed data
    event : :class:`pyrocko.meta.Event`
        reference event from catalog
    phase_def : :class:'pyrocko.gf.Timing'

    Returns
    -------
    :class:`numpy.ndarray`
    '''
    noised_traces = []
    for tr, station in zip(traces, stations):

        if tr is None:
            pass
        else:
            arrival_time = get_phase_arrival_time(
                engine=engine, source=event,
                station=station, wavename=phase_def,
                store_id=store_id)
            extracted = tr.chop(tr.tmin, arrival_time-10,
                                inplace=False)

            mean = num.mean(extracted.ydata)
            var = num.var(extracted.ydata)
            noise_data = num.random.normal(loc=mean,
            scale=var, size=num.shape(tr.ydata))
            tr.add(noise_data)
            noised_traces.append(tr)
    return noised_traces
