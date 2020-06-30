from gwosc.datasets import event_gps
from gwpy.timeseries import TimeSeries


gps = event_gps('GW170817')


# fetch data and compute ASDs for each involved detector around GW170817's merger
ldata = TimeSeries.fetch_open_data('L1', int(gps)-512, int(gps)+512)
lasd=ldata.asd(fftlength=4, method='median')

hdata = TimeSeries.fetch_open_data('H1', int(gps)-512, int(gps)+512)
hasd=hdata.asd(fftlength=4, method='median')

vdata = TimeSeries.fetch_open_data('V1', int(gps)-512, int(gps)+512)
vasd=vdata.asd(fftlength=4, method='median')


# plot ASD's on a single set of axes
plot = lasd.plot()
ax = plot.gca()
ax.set_xlim(10, 1400)
ax.set_ylim(5e-24, 1e-20)
plot.show(warn=False)

