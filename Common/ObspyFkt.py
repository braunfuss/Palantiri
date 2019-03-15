import obspy
import obspy.core
from obspy.geodetics import locations2degrees
import Basic


def loc2degrees(a, b):

    if type(a) is dict:
        a1 = Basic.dictToLocation(a)
    else:
        a1 = a

    if type(b) is dict:
        b1 = Basic.dictToLocation(b)
    else:
        b1 = b

    delta = locations2degrees(float(a1.lat), float(a1.lon), float(b1.lat),
                              float(b1.lon))
    return delta


def obs_TravelTimes(delta1, depth1):

    model = obspy.taup.TauPyModel(model='ak135')
    return model.get_travel_times(distance_in_degree=delta1,
                                  source_depth_in_km=float(depth1))
