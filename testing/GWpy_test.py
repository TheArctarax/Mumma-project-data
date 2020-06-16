import gwpy
from gwosc.datasets import event_gps
gps = event_gps('GW150914')
print(gps)
segment = (int(gps)-5, int(gps)+5)
print(segment)
from gwpy.timeseries import TimeSeries
ldata = TimeSeries.fetch_open_data('L1', *segment, verbose=True)
print(ldata)
