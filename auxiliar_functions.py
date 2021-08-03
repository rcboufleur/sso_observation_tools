from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.visualization as vis
from astropy.wcs import WCS
from astroquery.jplhorizons import Horizons
from astroquery.skyview import SkyView
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.patches import Rectangle

def ephem(start_time, stop_time, step, utc=0, obj_id='301', loc='874'):
    '''Retrieves ephemeris for a particular body at a particular location
    '''
    start, stop =  Time(start_time, format='isot'), Time(stop_time, format='isot')
    startjd, stopjd = Time(start.to_value('jd')-utc/24, format='jd'), Time(stop.to_value('jd')-utc/24, format='jd')
    start_time, stop_time = startjd.to_value('isot'), stopjd.to_value('isot')
    
    try:
        if (obj_id.upper() == 'MOON') or (obj_id.upper() == 'LUNA'):
            obj = Horizons(id='301', id_type='majorbody', location=loc, epochs={'start':start_time, 'stop':stop_time, 'step':step})
            return obj.ephemerides(quantities='1,4,9,24,25')
        else:
            obj = Horizons(id=obj_id, id_type='smallbody', location=loc, epochs={'start':start_time, 'stop':stop_time, 'step':step})
            return obj.ephemerides(quantities='1,4,9,24,25')
    except:
        obj = Horizons(id=obj_id, id_type='majorbody', location=loc, epochs={'start':start_time, 'stop':stop_time, 'step':step})
        return obj.ephemerides(quantities='1,4,9,24,25')

def reorder_table(astropy_table):
    ''' Reorders the main table deleting some columns and adding ra and dec as hmsdms
    '''
    c = SkyCoord(ra=astropy_table['RA'].data.data*u.degree, dec=astropy_table['DEC'].data.data*u.degree, frame='icrs')
    tmp = c.to_string('hmsdms')
    tmp = [s.replace('h', ' ') for s in tmp]
    tmp = [s.replace('m', ' ') for s in tmp]
    tmp = [s.replace('s', ' ') for s in tmp]
    tmp = [s.replace('d', ' ') for s in tmp]
    
    astropy_table.remove_column('targetname')
    astropy_table.remove_column('datetime_jd')
    try:
        astropy_table.remove_column('H')
        astropy_table.remove_column('G')
    except:
        pass
    astropy_table.remove_column('surfbright')
    astropy_table.remove_column('solar_presence')    
    astropy_table.remove_column('RA')
    astropy_table.remove_column('DEC')
    astropy_table.remove_column('AZ')
    astropy_table.add_column(tmp, name='RA_DEC', index=1)
    astropy_table['RA_DEC'].unit = 'hmsdms'
    return astropy_table

def plot_image_and_location(astropy_table, survey='DSS', pixels='512', radius=6.5*u.arcmin, body_path=False, save=''):
    ''' Downloads an imagette from data provided from the astroquery jpl horizons search
        using the SkyView class from astroquery
    '''
    
    import warnings
    pos_idx = np.argmax(astropy_table['EL'].data.data)
    position = str(astropy_table['RA'].data.data[pos_idx])+' '+str(astropy_table['DEC'].data.data[pos_idx])
    coordinates='icrs'
    warnings.filterwarnings('ignore')
    hdu = SkyView.get_images(position=position, radius=radius, coordinates=coordinates, survey=survey, pixels=pixels)[0][0]
    warnings.filterwarnings('default')
    # Tell matplotlib how to plot WCS axes
    fig = plt.figure(figsize=(9,9))
    wcs = WCS(hdu.header)
    ax = plt.gca(projection=wcs)

    # Plot the image
    interval = vis.ZScaleInterval(nsamples=1000, contrast=0.25, max_reject=0.5, min_npixels=5, krej=2.5, max_iterations=5)
    vmin,vmax = interval.get_limits(hdu.data)
    norm = vis.ImageNormalize(vmin=vmin, vmax=vmax, stretch=vis.SinhStretch())
    
    ax.imshow(hdu.data, cmap='gray', norm=norm)

    # plot bodu position at highest elevation on the image
    ax.scatter(astropy_table['RA'].data.data[pos_idx], astropy_table['DEC'].data.data[pos_idx], transform=ax.get_transform('icrs'), 
               s=1500, edgecolor='chartreuse', facecolor='none', linestyle='--')

    # plot body path over image
    if body_path:
        ax.plot(astropy_table['RA'].data.data, astropy_table['DEC'].data.data, 'r--', transform=ax.get_transform('icrs'))

    ax.coords.grid(True, color='white', ls='dotted')
    ax.coords[0].set_axislabel('Right Ascension')
    ax.coords[1].set_axislabel('Declination')

    if (len(save) > 0):
        plt.savefig(save, dpi=200)
    plt.show()
    return None

def plot_elevation_curves(prioritary='', secondary='', start_time='', stop_time='', step='', utc=0, upper_lim=95, lower_lim=30, sun=True, moon=True, save=''):
    '''Plot elevation curves for a given set of bodies
    '''
    # plot elevation curves of the selected bodies
    fig, ax = plt.subplots(1, 1, constrained_layout=True, figsize=(12,8))

    for obj in prioritary:
        eph = ephem(start_time, stop_time, step, utc=utc, obj_id=obj)
        tmp =  Time(eph['datetime_jd'].data.data+(utc/24), format='jd')
        tmp = tmp.to_value('iso')
        time = [datetime.fromisoformat(t) for t in tmp]
        ax.plot(time, eph['EL'].data.data,'-')
        txt_idx = np.argmax(eph['EL'].data.data)
        ax.text(time[txt_idx], eph['EL'].data.data[txt_idx]+1, eph['targetname'].data.data[0], horizontalalignment='center', verticalalignment='bottom', rotation=90, fontsize=10)

    for obj in secondary:
        eph = ephem(start_time, stop_time, step, utc=utc, obj_id=obj)
        tmp =  Time(eph['datetime_jd'].data.data+(utc/24), format='jd')
        tmp = tmp.to_value('iso')
        time = [datetime.fromisoformat(t) for t in tmp]
        ax.plot(time, eph['EL'].data.data,'--', color='grey')
        txt_idx = np.argmax(eph['EL'].data.data)
        ax.text(time[txt_idx], eph['EL'].data.data[txt_idx]+1, eph['targetname'].data.data[0], horizontalalignment='center', verticalalignment='bottom', rotation=90, fontsize=10, alpha=0.5)
    
    #draw sun/astronomical twilight zone
    if sun == True:
        suneph = ephem(start_time, stop_time, step, utc=utc, obj_id='Sun')
        tmp =  Time(suneph['datetime_jd'].data.data+(utc/24), format='jd')
        tmp = tmp.to_value('iso')
        time = [datetime.fromisoformat(t) for t in tmp]
        idx = np.where(suneph['EL'].data.data < -18)
        idx2 = np.where(suneph['EL'].data.data < 0)
        ax.plot(time,suneph['EL'].data.data)
        ax.add_patch(Rectangle(( mdates.date2num(time[0]), lower_lim), mdates.date2num(time[idx[0][0]])-mdates.date2num(time[0]), upper_lim-lower_lim, 
                                facecolor = 'grey', alpha=0.5, label='Sun/Twilight')) #add rectangle to plot
        ax.add_patch(Rectangle(( mdates.date2num(time[idx[0][-1]]), lower_lim), mdates.date2num(time[-1])-mdates.date2num(time[idx[0][-1]]), upper_lim-lower_lim,
                                facecolor = 'grey', alpha=0.5)) #add rectangle to plot
        ax.plot([mdates.date2num(time[idx2[0][0]]),mdates.date2num(time[idx2[0][0]])],[lower_lim,upper_lim], 'r--', label='Sunset/Sunrise')
        ax.plot([mdates.date2num(time[idx2[0][-1]]),mdates.date2num(time[idx2[0][-1]])],[lower_lim,upper_lim], 'r--')

        
    #draw moon path twilight zone
    if moon == True:
        mooneph = ephem(start_time, stop_time, step, utc=utc, obj_id='moon')
        tmp =  Time(mooneph['datetime_jd'].data.data+(utc/24), format='jd')
        tmp = tmp.to_value('iso')
        time = [datetime.fromisoformat(t) for t in tmp]
        idx = np.where(mooneph['EL'].data.data > 0)
        ax.plot(time,mooneph['EL'].data.data, 'k-.', label='Moon')
        ax.plot([mdates.date2num(time[idx[0][0]]),mdates.date2num(time[idx[0][0]])],[lower_lim,upper_lim], 'k--', label='Moonset/Moonrise')
        
    ax.plot(time[0],[0], '--', color='grey', label='Non-Prioritary')
    locator = mdates.AutoDateLocator(minticks=6, maxticks=15)
    formatter = mdates.ConciseDateFormatter(locator)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)
    ax.set_ylim(lower_lim,upper_lim)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_xlabel('Time (UTC'+str(utc)+')', fontsize=14)
    ax.set_ylabel('Elevation (degrees)', fontsize=14)
    ax.margins(x=0.0) 
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),
          fancybox=True, shadow=True, ncol=5)
    if (len(save) > 0):
        plt.savefig(save, dpi=200)
    plt.show()