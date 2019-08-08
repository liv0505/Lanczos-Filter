"""
Applying a bandpass filter
(the difference of two lowpass lanczos filters)
to a time-series.
==================================

This example demonstrates low pass filtering a time-series by applying a
weighted running mean over the time dimension.

The time-series used here is the EAR5 Reanalysis hourly 850hpa vorticity,
which is first averaged to daily data, and then filtered using two different
Lanczos filters, one to filter out time-scales of less than 3 years and one
to filter out time-scales of less than 10 years.

References
----------

    Duchon C. E. (1979) Lanczos Filtering in One and Two Dimensions.
    Journal of Applied Meteorology, Vol 18, pp 1016-1022.

"""

import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, savefig
import matplotlib.colors
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER



def low_pass_weights(window, cutoff):
    """Calculate weights for a low pass Lanczos filter.

    Args:

    window: int
        The length of the filter window.

    cutoff: float
        The cutoff frequency in inverse time steps.

    """
    order = ((window - 1) // 2 ) + 1
    nwts = 2 * order + 1
    w = np.zeros([nwts])
    n = nwts // 2
    w[n] = 2 * cutoff
    k = np.arange(1., n)
    sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
    firstfactor = np.sin(2. * np.pi * cutoff * k) / (np.pi * k)
    w[n-1:0:-1] = firstfactor * sigma
    w[n+1:-1] = firstfactor * sigma
    return w[1:-1]




def main():
    # window length for filters
    window = 50

    # construct 3 days and 10 days low pass filters
    hfw = low_pass_weights(window, 1. / 3.)
    lfw = low_pass_weights(window, 1. / 10.)
    weight_high = xr.DataArray(hfw, dims = ['window'])
    weight_low = xr.DataArray(lfw, dims = ['window'])

    # Load the hourly 850hPa vorticity
    fname = '/Users/serene.meng/EAR5_Test/tcseed/vort_850_06_to_11_2016'
    ds = xr.open_dataset(fname)

    '''
    Results are presented for Jul–Oct although data in June and November
    are also needed because the time filter (see section 2b) requires
    extra data at the beginning and end of each year’s time series.
    '''
    hourly_data = ds.sel(time=slice("2016-06-07T00:00:00","2016-11-24T23:00:00"))['vo']

    # Since the EAR5 data are hourly, we first average them into daily data.
    daily_data = hourly_data.resample(time='24H').mean(dim = 'time')

    # apply the filters using the rolling method with the weights
    lowpass_hf = daily_data.rolling(time = len(hfw), center = True).construct('window').dot(weight_high)
    lowpass_lf = daily_data.rolling(time = len(lfw), center = True).construct('window').dot(weight_low)

    # the bandpass is the difference of two lowpass filters.
    bandpass = lowpass_hf - lowpass_lf

    '''
    Okay, here I calculate the variance of 3-10 day bandpass filtered
    850hPa vorticity (Jul to Oct), which could be defined as
    TC seed index:  "pre-TC synoptic-scale disturbances"

    References
    ----------
    Li, T., Kwon, M., Zhao, M., Kug, J.S., Luo, J.J. and Yu, W., 2010.
    Global warming shifts Pacific tropical cyclone location.
    Geophysical Research Letters, 37(21).

    '''
    dvar = xr.DataArray.var(bandpass, dim = 'time', skipna = True)
    print(dvar)



    # plot the TC Seed Index
    '''
    for i in list(range(int(window/2)-1,len(bandpass.time)-int(window/2)+1)):
    #for i in list(range(int(window/2)-1,int(window/2)+1)):

        norm = matplotlib.colors.Normalize(vmin = 0, vmax = 0.00005)
        fig = plt.figure(figsize = (12,4.9), dpi = 300)
        ax = plt.axes(projection = ccrs.PlateCarree(central_longitude = 180))

        contour = bandpass.isel(time = i).plot(ax = ax, transform = ccrs.PlateCarree(), cmap = 'Blues',
                            norm = norm,
                            linewidth = 0, antialiased = False)
        ax.set_global()
        gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels=True, alpha=0.5, linewidth = 0.1, color = 'gray')
        ax.coastlines(linewidth = 0.3)
        #ax.set_title('The variance of 3-10 day bandpass filtered 850hPa vorticity')
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xlocator = mticker.FixedLocator([-180,-135,-90, -45, 0, 45, 90, 135])
        gl.ylocator = mticker.FixedLocator([-90,-60,-30, 0, 30, 60, 90])
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 10, 'color': 'gray'}
        gl.ylabel_style = {'size': 10,'color': 'gray'}

        plt.show()
        fig.savefig('/Users/serene.meng/EAR5_Test/tcseed/everyday_2017/' + str(i) + '.png', figsize = (12,4.9), dpi = 300)

    '''

    # Plot the daily bandpass vorticity
    norm = matplotlib.colors.Normalize(vmin = 0, vmax = 2e-9)


    fig = plt.figure(figsize = (12,4.9), dpi = 300)

    ax = plt.axes(projection = ccrs.PlateCarree(central_longitude = 180))

    contour = dvar.plot(ax = ax, transform = ccrs.PlateCarree(),
                        cmap = 'Blues',
                        norm = norm,
                        linewidth = 0, antialiased = False)

    ax.set_global()
    gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels=True, alpha=0.5, linewidth = 0.1, color = 'gray')
    ax.coastlines(linewidth = 0.3)
    #ax.set_title('The variance of 3-10 day bandpass filtered 850hPa vorticity')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlocator = mticker.FixedLocator([-180,-135,-90, -45, 0, 45, 90, 135])
    gl.ylocator = mticker.FixedLocator([-90,-60,-30, 0, 30, 60, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'gray'}
    gl.ylabel_style = {'size': 10,'color': 'gray'}

    plt.show()
    fig.savefig('/Users/serene.meng/EAR5_Test/tcseed/3-10day/2016JASO_3-10d_wl50.png', figsize = (12,4.9), dpi = 300)

if __name__ == '__main__':
    main()
