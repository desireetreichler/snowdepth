# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 14:36:36 2022

@author: desireet

This module contains functions to work with ICESat-2 data, elevation raster 
data (DEMs) and geodataframes. 
"""

# Developed and built on starting from topolib toolbox 
# - some functions are adjusted as the original produced errors
# - most functions are completely my own.


import h5py
import numpy as np
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point, Polygon
import rasterio
from rasterio import features
from pyproj import Proj   # assumes we only have one datum (WGMS84) 
from math import floor
import datetime
import time

#import rasterstats as rs

# DT added
import matplotlib.pyplot as plt


# TODO sort functions alphabetically and add docstring


# ---------------- convert ATL data to geodataframes -------------------------
#this function is from Ben
def ATL06_to_dict(filename, dataset_dict):
    """
        Read selected datasets from an ATL06 file

        Input arguments:
            filename: ATl06 file to read
            dataset_dict: A dictinary describing the fields to be read
                    keys give the group names to be read,
                    entries are lists of datasets within the groups
        Output argument:
            D6: dictionary containing ATL06 data.  Each dataset in
                dataset_dict has its own entry in D6.  Each dataset
                in D6 contains a list of numpy arrays containing the
                data
    """

    D6=[]
    pairs=[1, 2, 3]
    beams=['l','r']
    # open the HDF5 file
    with h5py.File(filename,'r') as h5f: # added mode to 'r', to suppress depreciation warning

        # loop over beam pairs
        for pair in pairs:
            # loop over beams
            for beam_ind, beam in enumerate(beams):
                # check if a beam exists, if not, skip it
                if '/gt%d%s/land_ice_segments' % (pair, beam) not in h5f:
                    continue
                # loop over the groups in the dataset dictionary
                temp={}
                for group in dataset_dict.keys():
                    for dataset in dataset_dict[group]:
                        DS='/gt%d%s/%s/%s' % (pair, beam, group, dataset)
                        #print(DS) # debugging
                        # since a dataset may not exist in a file, we're going to try to read it, and if it doesn't work, we'll move on to the next:
                        try:
                            temp[dataset]=np.array(h5f[DS])
                            # some parameters have a _FillValue attribute.  If it exists, use it to identify bad values, and set them to np.NaN
                            if '_FillValue' in h5f[DS].attrs:
                                fill_value=h5f[DS].attrs['_FillValue']
                                bad = temp[dataset]==fill_value
                                try:
                                    temp[dataset]=np.float64(temp[dataset])
                                    temp[dataset][bad]=np.NaN
                                except TypeError as e: # added this exception, as I got some type errors - temp[dataset] was 0.0
                                    pass
                            # a few attributes have 5 columns, corresponding to land, ocean, sea ice, land ice, inland water
                            try:
                                if len(temp[dataset][0]==5):
                                    for surftypenr, surftype in enumerate(['ocean','seaice','landice','inlandwater']):
                                        temp[dataset+'_'+surftype]=temp[dataset][:,surftypenr+1]
                                    temp[dataset]=temp[dataset][:,0] # default = land, first column
                            except TypeError as e: # added this exception, as I got some type errors - temp[dataset] was 0.0
                                pass    
                            except IndexError as e: # added this exception, as I got some *** IndexError: invalid index to scalar variable.
                                pass
                        except KeyError as e:
                            pass
                if len(temp) > 0:
                    # it's sometimes convenient to have the beam and the pair as part of the output data structure: This is how we put them there.
                    temp['pair']=np.zeros_like(temp['h_li'])+pair
                    temp['beam']=np.zeros_like(temp['h_li'])+beam_ind
                    #temp['filename']=filename
                    D6.append(temp)
    return D6



def ATL06_2_gdf(ATL06_fn,dataset_dict):
    """
    function to convert ATL06 hdf5 to geopandas dataframe, containing columns as passed in dataset dict
    Used Ben's ATL06_to_dict function
    """
    if ('latitude' in dataset_dict['land_ice_segments']) != True:
        dataset_dict['land_ice_segments'].append('latitude')
    if ('longitude' in dataset_dict['land_ice_segments']) != True:
        dataset_dict['land_ice_segments'].append('longitude')
    #use Ben's Scripts to convert to dict
    data_dict = ATL06_to_dict(ATL06_fn,dataset_dict)
    #print('converted to dict') # debugging
    #this will give us 6 tracks
    if len(data_dict) > 0: # added this to account for empty data_dicts, which occured
        i = 0
        for track in data_dict:
            #1 track
            #convert to datafrmae
            df = pd.DataFrame(track)
            try:
                df['p_b'] = str(track['pair'][0])+'_'+str(track['beam'][0])
            except:  # added, to account for errors - maybe where there is onluy one data point (?)
                df['p_b'] = str(track['pair'])+'_'+str(track['beam'])
    
            df['geometry'] = df.apply(point_covert,axis=1)
            if i==0:
                df_final = df.copy()
            else:
                df_final = df_final.append(df)
            i = i+1
        gdf_final = gpd.GeoDataFrame(df_final,geometry='geometry',crs='epsg:4326')  # changed from +init-version to avoid the upcoming warning
    else: 
        gdf_final = gpd.GeoDataFrame()
    return gdf_final

def get_ndv(ds):
    no_data = ds.nodatavals[0]
    if no_data == None:
        #this means no data is not set in tif tag, nead to cheat it from raster
        ndv = ds.read(1)[0,0]
    else:
        ndv = no_data
    return ndv


# ---------------------ATL03/ATL08 tools by Desiree -------------------------------------

# adapted for ATL08 by Desiree
def ATL08_to_dict(filename, dataset_dict):
    """
        Read selected datasets from an ATL06 file
        Input arguments:
            filename: ATl08 file to read
            dataset_dict: A dictinary describing the fields to be read
                    keys give the group names to be read,
                    entries are lists of datasets within the groups
        Output argument:
            D6: dictionary containing ATL08 data.  Each dataset in
                dataset_dict has its own entry in D6.  Each dataset
                in D6 contains a list of numpy arrays containing the
                data
    """

    D6=[]
    pairs=[1, 2, 3]
    beams=['l','r']
    # open the HDF5 file
    with h5py.File(filename,'r') as h5f: # added mode to 'r', to suppress depreciation warning
        # loop over beam pairs
        for pair in pairs:
            # loop over beams
            for beam_ind, beam in enumerate(beams):
                # check if a beam exists, if not, skip it
                if '/gt%d%s/land_segments' % (pair, beam) not in h5f:
                    continue
                # loop over the groups in the dataset dictionary
                temp={}
                for group in dataset_dict.keys():
                    for dataset in dataset_dict[group]:
                        DS='/gt%d%s/%s/%s' % (pair, beam, group, dataset)
                        # since a dataset may not exist in a file, we're going to try to read it, and if it doesn't work, we'll move on to the next:
                        try:
                            temp[dataset]=np.array(h5f[DS])
                            # some parameters have a _FillValue attribute.  If it exists, use it to identify bad values, and set them to np.NaN
                            if '_FillValue' in h5f[DS].attrs:
                                fill_value=h5f[DS].attrs['_FillValue']
                                bad = temp[dataset]==fill_value
                                try:
                                    temp[dataset]=np.float64(temp[dataset])
                                    temp[dataset][bad]=np.NaN
                                except TypeError as e: # added this exception, as I got some type errors - temp[dataset] was 0.0
                                    pass
                            # a few attributes have 5 columns, corresponding to either 30m segments or land, ocean, sea ice, land ice, inland water
                            try:
                                if len(temp[dataset][0]==5):
                                    if dataset in ['h_te_best_fit_20m']:
                                        for segnr, segstr in enumerate(['1','2','3','4','5']):
                                            temp[dataset+'_'+segstr]=temp[dataset][:,segnr]
                                        del temp[dataset]
                                    else: # default = land, first column
                                        for surftypenr, surftype in enumerate(['ocean','seaice','landice','inlandwater']):
                                            temp[dataset+'_'+surftype]=temp[dataset][:,surftypenr+1]
                                        temp[dataset]=temp[dataset][:,0] 
                            except TypeError as e: # added this exception, as I got some type errors - temp[dataset] was 0.0
                                pass    
                            except IndexError as e: # added this exception, as I got some *** IndexError: invalid index to scalar variable.
                                pass
                        except KeyError as e:
                            pass
                if len(temp) > 0:
                    # it's sometimes convenient to have the beam and the pair as part of the output data structure: This is how we put them there.
                    #a = np.zeros_like(temp['h_te_best_fit'])
                    #print(a)
                    temp['pair']=np.zeros_like(temp['h_te_best_fit'])+pair
                    temp['beam']=np.zeros_like(temp['h_te_best_fit'])+beam_ind
                    # RGT and cycle are also convenient. They are in the filename.
                    fs=filename.split('_')  
                    temp['RGT']=np.zeros_like(temp['h_te_best_fit'])+int(fs[-3][0:4])
                    temp['cycle']=np.zeros_like(temp['h_te_best_fit'])+int(fs[-3][4:6]) 
                    #temp['filename']=filename
                    D6.append(temp)
    return D6



def ATL08_2_gdf(ATL06_fn,dataset_dict):
    """
    function to convert ATL06 hdf5 to geopandas dataframe, containing columns as passed in dataset dict
    Used Ben's ATL06_to_dict function
    """
    if ('latitude' in dataset_dict['land_segments']) != True:
        dataset_dict['land_segments'].append('latitude')
    if ('longitude' in dataset_dict['land_segments']) != True:
        dataset_dict['land_segments'].append('longitude')
    #use Ben's Scripts to convert to dict
    
    data_dict = ATL08_to_dict(ATL06_fn,dataset_dict)
    #this will give us 6 tracks
    i = 0
    for track in data_dict:
        #1 track
        #convert to datafrmae
        df = pd.DataFrame(track)
        try:
            df['p_b'] = str(track['pair'][0])+'_'+str(track['beam'][0])
        except:  # added, to account for errors - maybe where there is onluy one data point (?)
            df['p_b'] = str(track['pair'])+'_'+str(track['beam'])

        df['geometry'] = df.apply(point_covert,axis=1)
        if i==0:
            df_final = df.copy()
        else:
            df_final = df_final.append(df)
        i = i+1
    gdf_final = gpd.GeoDataFrame(df_final,geometry='geometry',crs='epsg:4326')  # changed from +init-version to avoid the upcoming warning
    return gdf_final


  
# adapted for ATL03 by Desiree
def ATL03_to_dict(filename, dataset_dict=False, utmzone=False):
    """
        Read selected datasets from an ATL03 file
        Input arguments:
            filename: ATl03 file to read
            dataset_dict: A dictinary describing the fields to be read
                    keys give the group names to be read,
                    entries are lists of datasets within the groups. If not set, 
                    a standard minimum dict of lat, lon, height, time, quality flag is read
            utmzone: optional, if set, lat/lon are converted to the provided utmzone string (e.g.: 33N)
        Output argument:
            D3: dictionary containing ATL03 data.  Each dataset in
                dataset_dict has its own entry in D3.  Each dataset
                in D3 contains a list of numpy arrays containing the
                data
    """
#filename = ATL03_file[filenr]
#h5f=h5py.File(filename,'r') 
    if dataset_dict==False:
        dataset_dict = {'heights':['delta_time','lon_ph','lat_ph','h_ph','signal_conf_ph']}
     

    D3=[]
    pairs=[1, 2, 3]
    beams=['l','r']
    # open the HDF5 file
    with h5py.File(filename,'r') as h5f: # added mode to 'r', to suppress depreciation warning
        # loop over beam pairs
        for pair in pairs:
            # loop over beams
            for beam_ind, beam in enumerate(beams):
                # check if a beam exists, if not, skip it
                if '/gt%d%s/heights' % (pair, beam) not in h5f:
                    continue
                # loop over the groups in the dataset dictionary
                temp={}
                for group in dataset_dict.keys():
                    if group in ['geolocation']: # check whether this data is available for each photon or for each 20m segment only.
                        segmentrate = 1
                    else: segmentrate = 0
                    for dataset in dataset_dict[group]:
                        DS='/gt%d%s/%s/%s' % (pair, beam, group, dataset)
                        # since a dataset may not exist in a file, we're going to try to read it, and if it doesn't work, we'll move on to the next:
                        try:
                            temp[dataset]=np.array(h5f[DS])
                            # some parameters have a _FillValue attribute.  If it exists, use it to identify bad values, and set them to np.NaN
                            if '_FillValue' in h5f[DS].attrs:
                                fill_value=h5f[DS].attrs['_FillValue']
                                bad = temp[dataset]==fill_value
                                try:
                                    temp[dataset]=np.float64(temp[dataset])
                                    temp[dataset][bad]=np.NaN
                                except TypeError as e: # added this exception, as I got some type errors - temp[dataset] was 0.0
                                    pass
                            # a few attributes have 5 columns, corresponding to land, ocean, sea ice, land ice, inland water
                            try:
                                if len(temp[dataset][0]==5):
                                    for surftypenr, surftype in enumerate(['ocean','seaice','landice','inlandwater']):
                                        temp[dataset+'_'+surftype]=temp[dataset][:,surftypenr+1]
                                    temp[dataset]=temp[dataset][:,0] # default = land, first column
                            except TypeError as e: # added this exception, as I got some type errors - temp[dataset] was 0.0
                                pass  
                            # some data is only available at the 20 m segment rate (or even less). Fix this by duplicating
                            if segmentrate:
                                DS='/gt%d%s/%s/%s' % (pair, beam, 'geolocation', 'segment_ph_cnt')
                                segment_ph_cnt=np.array(h5f[DS])
                                temp2 = np.concatenate([np.repeat(y,x) for x,y in zip(segment_ph_cnt,temp[dataset])])
                                temp[dataset]=temp2 # overwrite
                        except KeyError as e:
                            pass
                if len(temp) > 0:
                    # it's sometimes convenient to have the beam and the pair as part of the output data structure: This is how we put them there.
                    #a = np.zeros_like(temp['h_te_best_fit'])
                    #print(a)
                    temp['pair']=np.zeros_like(temp['h_ph'])+pair
                    temp['beam']=np.zeros_like(temp['h_ph'])+beam_ind
                    # add x_atc, y_atc
                    temp['x_atc'], temp['y_atc']=get_ATL03_x_atc(h5f, pair, beam,temp)
                    if utmzone: # default: false
                        # add utm XXX x/y
                        temp['utmx'], temp['utmy']=ATL03_to_utm(temp,utmzone)
                    #print(temp)
                    #temp['filename']=filename
                    D3.append(temp)
    return D3

def get_ATL03_x_atc(h5f, pair, beam, temp):
    # calculate the along-track and across-track coordinates for ATL03 photons
        #-- data and attributes for beam gtx
    x_atc=np.zeros_like(temp['h_ph'])+np.NaN
    y_atc=np.zeros_like(temp['h_ph'])+np.NaN
                #-- ATL03 Segment ID
    DS='/gt%d%s/%s/%s' % (pair, beam, 'geolocation','segment_id')
    Segment_ID = np.array(h5f[DS])
            #Segment_ID[gtx] = val['geolocation']['segment_id']
    n_seg = len(Segment_ID)#[gtx])
            #-- first photon in the segment (convert to 0-based indexing)
    DS='/gt%d%s/%s/%s' % (pair, beam, 'geolocation','ph_index_beg')
    Segment_Index_begin = np.array(h5f[DS])-1
            #Segment_Index_begin[gtx] = val['geolocation']['ph_index_beg'] - 1
            #-- number of photon events in the segment
    DS='/gt%d%s/%s/%s' % (pair, beam, 'geolocation','segment_ph_cnt')
    Segment_PE_count = np.array(h5f[DS])
            #Segment_PE_count[gtx] = val['geolocation']['segment_ph_cnt']
            #-- along-track distance for each ATL03 segment
    DS='/gt%d%s/%s/%s' % (pair, beam, 'geolocation','segment_dist_x')
    Segment_Distance = np.array(h5f[DS])
            #Segment_Distance[gtx] = val['geolocation']['segment_dist_x']
            #-- along-track length for each ATL03 segment
    DS='/gt%d%s/%s/%s' % (pair, beam, 'geolocation','segment_length')
    Segment_Length = np.array(h5f[DS])
            #Segment_Length[gtx] = val['geolocation']['segment_length']
            #-- Transmit time of the reference photon
    DS='/gt%d%s/%s/%s' % (pair, beam, 'geolocation','delta_time')
    delta_time = np.array(h5f[DS])
            #delta_time = val['geolocation']['delta_time']
    # -- distance between photons
    DS='/gt%d%s/%s/%s' % (pair, beam, 'heights','dist_ph_along')
    dist_ph_along = np.array(h5f[DS])
    DS='/gt%d%s/%s/%s' % (pair, beam, 'heights','dist_ph_across')
    dist_ph_across = np.array(h5f[DS])    
            #-- iterate over ATL03 segments to calculate 40m means
            #-- in ATL03 1-based indexing: invalid == 0
            #-- here in 0-based indexing: invalid == -1   
    segment_indices, = np.nonzero((Segment_Index_begin[:-1] >= 0) &
                (Segment_Index_begin[1:] >= 0))
    for j in segment_indices:
                #-- index for segment j
        idx = Segment_Index_begin[j]
        #-- number of photons in segment (use 2 ATL03 segments)
        c1 = np.copy(Segment_PE_count[j])
        c2 = np.copy(Segment_PE_count[j+1])
        cnt = c1 + c2
                    #-- time of each Photon event (PE)
        segment_delta_times = temp['delta_time'][idx:idx+cnt]
                    #-- Photon event lat/lon and elevation (WGS84)
        segment_heights = temp['h_ph'][idx:idx+cnt]
        segment_lats = temp['lat_ph'][idx:idx+cnt]
        segment_lons = temp['lon_ph'][idx:idx+cnt]
                    #-- Along-track and Across-track distances
        distance_along_X = np.copy(dist_ph_along[idx:idx+cnt])
        distance_along_X[:c1] += Segment_Distance[j]
        distance_along_X[c1:] += Segment_Distance[j+1]
        distance_along_Y = np.copy(dist_ph_across[idx:idx+cnt])
        x_atc[idx:idx+cnt]=distance_along_X
        y_atc[idx:idx+cnt]=distance_along_Y
        
    return x_atc, y_atc
            
            

        
def ATL03_to_utm(temp, utmzone):
    if utmzone[-1]=='N':
        northsouth='north'
    else:
        northsouth='south'
    #utmx=np.zeros_like(temp['h_ph'])+np.NaN
    #utmy=np.zeros_like(temp['h_ph'])+np.NaN
    # convert to x/y utm coordinates 
    myProj = Proj("+proj=utm +zone=%s, +%s +ellps=WGS84 +datum=WGS84 +units=m +no_defs" % (utmzone,northsouth)) ## assumes we only have one datum (WGMS84, otherwise use transform) and want UTMXXX
    utmx, utmy = myProj(temp['lon_ph'],temp['lat_ph']) # to convert back, add argument: , inverse=True
    
    return utmx, utmy


def ATL03_2_gdf(ATL03_fn,dataset_dict=False,aoicoords=False, filterbackground=False, utmzone=False, v=False):
    """
    function to convert ATL03 hdf5 to geopandas dataframe, containing columns as passed in dataset dict. Based on Ben's functions.

    additionally: aoicoords -> [xmin,ymin,xmax,ymax] if set, the data will be clipped to that bounding box.
                  filterbackground -> default False, otherwise set to remove values (e.g. [-2, -1,0, 1], only excluded if true for BOTH land (standard) signal_conf_ph and the one for land ice)
                  v: if set to True, the function prints the current file - keeping track of progress in a loop.
                  
                  signal_conf_ph meaning/values available: Confidence level associated with each photon event selected as signal. 
                  0=noise. 1=added to allow for buffer but algorithm classifies as background; 2=low; 3=med; 4=high). This parameter 
                  is a 5xN array where N is the number of photons in the granule, and the 5 rows indicate signal finding for each
                  surface type (in order: land, ocean, sea ice, land ice and inland water). Events not associated with a specific surface
                  type have a confidence level of ­1. Events evaluated as TEP returns have a confidence level of ­2. flag_values: ­2, ­1, 0, 1, 2, 3, 4 
                  flag_meanings : possible_tep not_considered noise buffer low medium high
    """
    if dataset_dict==False:
        dataset_dict = {'heights':['delta_time','lon_ph','lat_ph','h_ph','signal_conf_ph']}
    clip=False
    if type(aoicoords)!=type(False): #aoicoords.any():
       xmin,ymin,xmax,ymax=aoicoords
       clip=True

    if ('lat_ph' in dataset_dict['heights']) != True:
        dataset_dict['heights'].append('lat_ph')
    if ('lon_ph' in dataset_dict['heights']) != True:
        dataset_dict['heights'].append('lon_ph')
    #use the above function to convert to dict
    
    # verbose: keep track of current file
    if v: print('converting '+ ATL03_fn +' to gdf...')
    
    data_dict = ATL03_to_dict(ATL03_fn,dataset_dict, utmzone)
    if len(data_dict) > 0: # added this to account for empty data_dicts, which occured
        #this will give us 6 tracks
        i = 0
        for track in data_dict:
            #1 track
            #convert to datafrmae
            df = pd.DataFrame(track)
            # filter by aoi (owerwrite df)
            if clip:
                df=df.loc[(df['lon_ph']>xmin) & (df['lon_ph']<xmax) & (df['lat_ph']>ymin) & (df['lat_ph']<ymax)]
            # filter photons classified as background (0,1, for land / ice), owerwrite df - should this also include filtering not assigned (-1)? TBD
            if filterbackground:
                df=df[ (~df['signal_conf_ph'].isin(filterbackground)) | (~df['signal_conf_ph_landice'].isin(filterbackground)) ]
                # old filter method, discarded
                #df = df.loc[((df['signal_conf_ph'] >1 ) | (df['signal_conf_ph'] == -1 )) &  ((df['signal_conf_ph_landice'] >1)| (df['signal_conf_ph_landice'] == -1 ))]
            # add track/beam
            try:
                df['p_b'] = str(track['pair'][0])+'_'+str(track['beam'][0])
            except:  # added, to account for errors - maybe where there is onluy one data point (?)
                df['p_b'] = str(track['pair'])+'_'+str(track['beam'])
            # df['geometry'] = df.apply(point_covert,axis=1) # this did not work here??? but for ATL08 yes? strange.
            # df['geometry'] = [Point(xy) for xy in zip(df.lon_ph, df.lat_ph)] # this works, but is slow
            if i==0:
                df_final = df.copy()
            else:
                df_final = df_final.append(df)
            i = i+1
        try:
            gdf_final = gpd.GeoDataFrame(df_final,geometry=gpd.points_from_xy(x=df_final.lon_ph, y=df_final.lat_ph),crs='epsg:4326')  # changed from +init-version to avoid the upcoming warning
        except: # added this exception, as I got some errors that the df_final was not defined when there was a problem reading out height data for the granule - skip
            gdf_final = gpd.GeoDataFrame() # empty frame   
    else:
        gdf_final = gpd.GeoDataFrame() # empty frame  
    return gdf_final 


# wrapper to handle lists of granules, and add some custom parameters

def hdf2gdf(ATL_list, ATL_version=8, dataset_dict={}):
    """ 
    Load data from granules and add to a single gdf. Add ID number and date
    to data rows.
    Usage: hdf2gdf(ATL_list, ATL_version=8, dataset_dict={})
    ATL_list:     a list with granule filenames to load
    ATL_version:  8 or 3, default: 8
    dataset_dict: dictionary for which datasets to add to the dataframe. 
                  Default available for ATL08 and ALT03 (see function for details).
    Returns:      gdf: gdf with all points (rows) and parameters in data_dict (columns)
                  dt: list of dates with data
    """
    
    # check whether dataset_dict is provided
    if len(dataset_dict)==0:
        if ATL_version==8:
            dataset_dict = {'land_segments': ['delta_time', 'longitude', 'latitude', 'dem_h', 'dem_flag', 'msw_flag', 'n_seg_ph',
                                  'night_flag', 'rgt', 'segment_id', 'segment_landcover', 'segment_snowcover',  # 'surf_type',
                                  'segment_watercover', 'sigma_h', 'sigma_topo', 'quality', 'terrain_flg'],
                'land_segments/terrain': ['h_te_best_fit', 'h_te_best_fit_20m','h_te_interp', 'h_te_max', 'h_te_mean', 'h_te_median', 'h_te_min',
                                          'h_te_mode',
                                          'h_te_skew', 'h_te_std', 'h_te_uncertainty', 'n_te_photons', 'terrain_slope']}
        elif ATL_version==3:
                dataset_dict = {'heights': ['delta_time', 'lon_ph', 'lat_ph', 'h_ph', 'signal_conf_ph'],
                                'geolocation': ['segment_id', 'segment_ph_cnt']}
                
    # load data into a list of gdfs          
    if ATL_version==8:     
           gdf_list = [(ATL08_2_gdf(x, dataset_dict)) for x in ATL_list]  # used my own, as the toolbox one caused some errors - maybe for a dataset that had just one point?
    elif ATL_version==3:    
           gdf_list = [(ATL03_2_gdf(x, dataset_dict)) for x in ATL_list]  # used my own, as the toolbox one caused some errors - maybe for a dataset that had just one point?
           
    # concatenate to a single geodataframe 
    gdf = concat_gdf(gdf_list)
           
    # add ID nr, useful to merge different subsets
    gdf['ID'] = np.arange(1, len(gdf) + 1)
    # add date
    deltatime = np.floor(gdf['delta_time'].to_numpy() /3600/24)
    currentdates, datecounts=np.unique(deltatime,return_counts=True)
    dt =[datetime.date(2018,1,1)+datetime.timedelta(d) for d in currentdates] 
    date =[datetime.date(2018,1,1)+datetime.timedelta(d) for d in deltatime] 
    gdf['date']=[10000*x.year + 100*x.month +x.day for x in date]
    gdf['month']=[x.month for x in date]
    
    return gdf, dt

# ---- other tools, copied from topolib toolbox, used in some of the wrappers below/above -----------------------------

def dem2polygon(dem_file_name):
    """
        Take DEM and return polygon geodataframe matching the extent and coordinate system of the input DEM.
        
        Input parameters:
        dem_file_name: Absolute path to a DEM file
        
        Output parameters:
        dem_polygon: A polygon geodataframe matching the extent and coordinate system of the input DEM.

    """
    
    # read in dem using rasterio
    dem = rasterio.open(dem_file_name)
    
    # extact total bounds of dem
    bbox = dem.bounds
    
    # convert to corner points
    p1 = Point(bbox[0], bbox[3])
    p2 = Point(bbox[2], bbox[3])
    p3 = Point(bbox[2], bbox[1])
    p4 = Point(bbox[0], bbox[1])
    
    # extract corner coordinates
    np1 = (p1.coords.xy[0][0], p1.coords.xy[1][0])
    np2 = (p2.coords.xy[0][0], p2.coords.xy[1][0])
    np3 = (p3.coords.xy[0][0], p3.coords.xy[1][0])
    np4 = (p4.coords.xy[0][0], p4.coords.xy[1][0])
    
    # convert to polygon
    bb_polygon = Polygon([np1, np2, np3, np4])

    # create geodataframe
    dem_polygon_gdf = gpd.GeoDataFrame(gpd.GeoSeries(bb_polygon), columns=['geometry'])
    
    dem_polygon_gdf.crs = dem.crs
    
    return dem_polygon_gdf

def point_covert(row):
    geom = Point(row['longitude'],row['latitude'])
    return geom



def sample_near_nbor(ds,geom):
    """
    sample values from raster at the given ICESat-2 points
    using nearest neighbour algorithm
    Inputs are a rasterio dataset and a geodataframe of ice_sat2 points
    """
    # reproject the shapefile to raster projection
    x_min,y_min,x_max,y_max = ds.bounds
    geom = geom.to_crs(ds.crs)
    #filter geom outside bounds
    geom = geom.cx[x_min:x_max,y_min:y_max]
    X = geom.geometry.x.values
    Y = geom.geometry.y.values
    xy = np.vstack((X,Y)).T
    sampled_values = np.array(list(ds.sample(xy)))
    no_data = get_ndv(ds)
    sample = np.ma.fix_invalid(np.reshape(sampled_values,np.shape(sampled_values)[0]))
    sample = np.ma.masked_equal(sample,no_data)
    x_atc = np.ma.array(geom.x_atc.values,mask = sample.mask)
    return x_atc, sample

def buffer_sampler(ds,geom,buffer,val='median',ret_gdf=False):
    """
    sample values from raster at the given ICESat-2 points
    using a buffer distance, and return median/mean or a full gdf ( if return gdf=True)
    Inputs = rasterio dataset, Geodataframe containing points, buffer distance, output value = median/mean (default median)
    and output format list of x_atc,output_value arrays (default) or  full gdf
    """
    import rasterstats as rs
    ndv = get_ndv(ds)
    array = ds.read(1)
    gt = ds.transform
    stat = val
    geom = geom.to_crs(ds.crs)
    x_min,y_min,x_max,y_max = ds.bounds
    geom = geom.cx[x_min:x_max, y_min:y_max]
    geom['geometry'] = geom.geometry.buffer(buffer)
    json_stats = rs.zonal_stats(geom,array,affine=gt,geojson_out=True,stats=stat,nodata=ndv)
    gdf = gpd.GeoDataFrame.from_features(json_stats)
    if val =='median':
        gdf = gdf.rename(columns={'median':'med'})
        call = 'med'
    else:
        gdf = gdf.rename(columns={'mean':'avg'})
        call = 'avg'
    if ret_gdf:
        out_file = gdf
    else:
        out_file = [gdf.x_atc.values,gdf[call].values]
    return out_file


def concat_gdf(gdf_list):
    """
    concatanate geodataframes into 1 geodataframe
    Assumes all input geodataframes have same projection
    Inputs : list of geodataframes in same projection
    Output : 1 geodataframe containing everything having the same projection
    """
    #from https://stackoverflow.com/questions/48874113/concat-multiple-shapefiles-via-geopandas
    gdf = pd.concat([gdf for gdf in gdf_list]).pipe(gpd.GeoDataFrame)
    gdf.crs = (gdf_list[0].crs)
    if gdf.crs is None:
        gdf.crs='epsg:4326'  # sometimes the first geodataframe in a list may be empty, causing the result not to have a coordinate system.
    return gdf
    
    
### Various tools to work with geodataframes (Desiree) ###
#------------------------------------------------------------------------"
def overlappingBB(gdfA, gdfB):
    """ Check whether the total bounding boxes of two geodataframes overlap.
    Usage: overlappingBB(gdfA, gdfB)
    Returns True/False """
    xmin, ymin, xmax, ymax = gdfA.total_bounds
    xminB, yminB, xmaxB, ymaxB = gdfB.total_bounds
    result1 =  ((xmin < xmaxB) & (xmax > xminB) & (ymin < ymaxB) & (ymax > yminB) ) 
    result2 = ((xminB < xmax) & (xmaxB > xmin) & (yminB < ymax) & (ymaxB > ymin) )
    result = result1|result2
    return result

# ------------ Point Extract wrapper ----------------------------------

def ic2elevextract(gdf,dem_file,demname):
    """ Extract elevation data from a DEM at point locations from a geodataframe.
    The function is based on rasterio. Adds a new column to the gdf.
    Usage:  ic2elevextract(gdf,dem_file,demname)
    gdf:        point geodataframe
    dem_file:   location of a DEM file (will be opened with rasterio)
    demname:    name of the DEM, will be the column name for the extracted elevation"""
                           
    tSTART=time.time()
    t2=time.time()    
    
    # open the elevation file
    dem_ds = rasterio.open(dem_file)
    
    # clip to dem extent  
    xmin,ymin,xmax,ymax = dem_ds.bounds
    if (gdf.crs.srs in dem_ds.crs.data['init'])==False: # are they in the same projection?
        print('reproject to DEM crs...')    
        gdf_dem_extent = gdf.to_crs(dem_ds.crs).cx[xmin:xmax,ymin:ymax]   # cx: coordinate indexer
    else: 
        gdf_dem_extent = gdf.cx[xmin:xmax,ymin:ymax]   # cx: coordinate indexer
        
    # sample elevations 
    #print('sample elevations for '+demname+'...')
    gdf_dem_extent=elevextrRasterio(gdf_dem_extent,dem_ds,demname)
    # locate in original gdf
    gdf.loc[gdf['ID'].isin(gdf_dem_extent['ID']), demname]=gdf_dem_extent[demname]
    
    
    if 0: # debug plot
        fig,ax = plt.subplots(figsize=(10,10))
        clim = np.nanpercentile(gdf['h_EGM'].values,(10,90))
        rasterio.plot.show(dem_ds,ax=ax)
        ax.plot([xmin,xmax,xmax,xmin,xmin],[ymin,ymin,ymax,ymax,ymin]) # dem extent
        gdf.plot('h_EGM',ax=ax,s=3,legend=True,cmap='inferno',vmin=clim[0],vmax=clim[1],legend_kwds={'label':'Elevation [m]'})
            
    print('total processing time for this DEM %f' %(time.time()-tSTART))   
    return gdf


def elevextrRasterio(gdf_dem_extent, dem_ds, demname):
    # raster sampling: get nearest neighbour DEM elevation and add to geodataframe. Use all tracks and beams.
    print(f'sampling nearest neighbour with rasterio.sample method - for {len(gdf_dem_extent)} points...') # 5.6 sec
    t2=time.time()
    # THIS WAS IN LAT/LON, CAUSING FAILURE OF ELEVATION EXTRACTION
    coords = [(x,y) for x, y in zip(gdf_dem_extent.geometry.x, gdf_dem_extent.geometry.y)]
    # sample
    gdf_dem_extent[demname] =np.array([x[0] for x in dem_ds.sample(coords)],)
    # set strange values out of range to nan:
    gdf_dem_extent.loc[gdf_dem_extent[demname]<-200 ,demname] =np.nan
    gdf_dem_extent.loc[gdf_dem_extent[demname]>9000 ,demname] =np.nan
    # store also in the original dataframe
    #gdf.loc[gdf['ID'].isin(gdf_dem_extent['ID']), demname]=gdf_dem_extent[demname]
    print('elapsed time: %f' %(time.time()-t2))


def ic2dfelevextract(df,dem_file,demname, lon='longitude',lat='latitude', date=None):
    """ Extract elevation data from a DEM at point locations from a dataframe.
    The function is based on rasterio. Adds a new column to the df.
    CAUTION: the function does not check whether the projections agree! Use 
             ic2elevextract instead, if you need that. 
    Usage: df = ic2dfelevextract(df,dem_file,demname, lon='longitude',lat='latitude')
    df:        point dataframe
    dem_file:   location of a DEM file (will be opened with rasterio)
    demname:    name of the DEM, will be the column name for the extracted elevation
    lon:       optional, name of the column with lon/x coordinates. Default: 'longitude'
    lat:       optional, name of the column with lat/y coordinates. Default: 'latitude'
    date:      optional, one single date (integer) present in the 'date' column - 
                         elevextract done only for points of this date 
    """                           
    tSTART=time.time()
   
    # open the elevation file
    dem_ds = rasterio.open(dem_file)
    
    # ensure that the df has an ID column 
    if ~('ID' in df.columns):
        removeID = True
        df['ID'] = np.arange(1, len(df) + 1)
    else: 
        removeID = False

    # clip to dem extent  
    xmin,ymin,xmax,ymax = dem_ds.bounds
    if date != None: 
        df_dem_extent= df[(df[lon]>=xmin) & (df[lon]<=xmax) & (df[lat]>=ymin) & (df[lat]<=ymax) & (df['date']==date)]
    else:
        df_dem_extent= df[(df[lon]>=xmin) & (df[lon]<=xmax) & (df[lat]>=ymin) & (df[lat]<=ymax)]
        
    # sample elevations 
    #print('sample elevations for '+demname+'...')
    print(f'sampling {len(df_dem_extent)} points...')
    coords = [(x,y) for x, y in zip(df_dem_extent[lon], df_dem_extent[lat])]
    # sample
    df_dem_extent[demname] =np.array([x[0] for x in dem_ds.sample(coords)],)
    # set strange values out of range to nan:
    df_dem_extent.loc[df_dem_extent[demname]<-200 ,demname] =np.nan
    df_dem_extent.loc[df_dem_extent[demname]>9000 ,demname] =np.nan
    # locate in original gdf
    df.loc[df['ID'].isin(df_dem_extent['ID']), demname]=df_dem_extent[demname]
    
    if removeID:
        # never remove ID... seems not to work reliably
        print()
        #df=df.drop(columns={'ID'})
    
    if 0: # debug plot
        fig,ax = plt.subplots(figsize=(10,10))
        clim = np.nanpercentile(df[demname].values,(10,90))
        rasterio.plot.show(dem_ds,ax=ax)
        ax.plot([xmin,xmax,xmax,xmin,xmin],[ymin,ymin,ymax,ymax,ymin]) # dem extent
        plt.scatter(df.lon,df.lat,df.demname,ax=ax,s=3,legend=True,cmap='inferno',vmin=clim[0],vmax=clim[1],legend_kwds={'label':'Elevation [m]'})
            
    print('total processing time for this DEM %f' %(time.time()-tSTART))   
    return df



# def gdalslope(dem_file):
#     # use gdal to compute the slope of the input DEM file
#     #sp_file=os.path.splitext(dem_file)[0]+'_slope.tif'    
#     #! gdaldem slope $dem_file $sp_file
#     return sp_file

# def gdalaspect(dem_file):
#     # use gdal to compute the aspect of the input DEM file
#     as_file=os.path.splitext(dem_file)[0]+'_aspect.tif'    
#     #! gdaldem aspect $dem_file $as_file
#     return as_file

# def gdalhillshade(dem_file):
#     # use gdal to compute the hillshade of the input DEM file
#     hs_file=os.path.splitext(dem_file)[0]+'_hs.tif'    
#     #! gdaldem hillshade $dem_file $hs_file
#     return hs_file            

    
def elevextrNN(gdf_dem_extent, dem_ds, demname):
    # raster sampling from ICESat-2 hackathlon from ca 2019, takes ca. twice the time of elevextrRasterio
    print(f'sampling with nearest neighbour method - for {len(gdf_dem_extent)} points...') # 9 sec
    t2=time.time()
    del_time,elev = sample_near_nbor(dem_ds,gdf_dem_extent)
    gdf_dem_extent[demname] = elev # but somehow, the mask doesn't get transported into the geodataframe.
    # store also in the original dataframe
    #gdf.loc[gdf['ID'].isin(gdf_dem_extent['ID']), demname]=elev
    print('elapsed time: %f' %(time.time()-t2))
    return gdf_dem_extent



### --------- display geodataframes on raster ---------------------------

def gdf_on_raster(gdf,ds,ax,shp,hs_ds,cmap='inferno'): 
    """
    plot the ATL08 data points on top of your own raster, with RASTER coordinate system
    Parameters
    ----------
    gdf : ATL08 geodataframe
    ds : dem dataset
    ax : figure axes
    shp : shapefile
    hs_ds : hillshade dataset, optional
    cmap : colormap, optional. The default is 'inferno'.

    """
    gdf = gdf.to_crs(ds.crs)
    ## added shp functionality to function myself
    shp=shp.to_crs(ds.crs)
    xmin,ymin,xmax,ymax = ds.bounds #ds.bounds
    ndv = get_ndv(ds) # get no data values
    img = ds.read(1)
    img = np.ma.masked_less_equal(img,ndv) # mask them
    clim = np.nanpercentile(img.compressed(),(2,98)) # if not vompressed, all the masked values are still counted
    if hs_ds:
        hs = hs_ds.read(1)
        ndv = get_ndv(hs_ds)
        hs = np.ma.masked_less_equal(hs,ndv)
        ax.imshow(hs,cmap='gray',extent=[xmin,xmax,ymin,ymax])
        im = ax.imshow(img,alpha=0.6,cmap=cmap,extent=[xmin,xmax,ymin,ymax])
        print('colorbar caxis:'), print(clim)
    else:
        im = ax.imshow(img,cmap=cmap,vmin=clim[0],vmax=clim[1],extent=[xmin,xmax,ymin,ymax])
    #gdf.plot('h_te_best_fit',ax=ax,s=1,cmap='inferno',vmin=clim[0],vmax=clim[1])# brfore:('p_b',ax=ax,s=1)   # p_b   - pair-beam
    ## added shp functionality to function myself
    shp.plot(color='k',alpha=0.3, ax=ax)
    shp.boundary.plot(edgecolor="black", ax=ax)
    ###
    gdf.plot('p_b',ax=ax,s=1)   # p_b   - pair-beam
    plt.colorbar(im,ax=ax,extend='both',label='Elevation (m)')
    
    
def gdfsubset_on_raster(gdf,ds,ax,shp,hs_ds=False,subset=0.01,cmap='inferno',vartoplot='p_b'): 
    """
    plot a subset (fraction 2subset) of  ATL03 data points on top of your own raster, with RASTER coordinate system
    Parameters
    ----------
    gdf : ATL03 geodataframe
    ds : dem dataset
    ax : figure axes
    shp : shapefile
    hs_ds : hillshade dataset, optional (default: false)
    subset: fraction of datapoints to plot, default: 1% (0.01)
    cmap : colormap, optional. The default is 'inferno'.
    vartoplot: p_b, specify something else if you want

    """
    gdf = gdf.to_crs(ds.crs)
    slicestep=floor(1/subset)
    gdf=gdf[0::slicestep]
    ## added shp functionality to function myself
    shp=shp.to_crs(ds.crs)
    xmin,ymin,xmax,ymax = ds.bounds #ds.bounds
    ndv = get_ndv(ds) # get no data values
    img = ds.read(1)
    img = np.ma.masked_less_equal(img,ndv) # mask them
    clim = np.nanpercentile(img.compressed(),(2,98)) # if not vompressed, all the masked values are still counted
    if hs_ds:
        hs = hs_ds.read(1)
        ndv = get_ndv(hs_ds)
        hs = np.ma.masked_less_equal(hs,ndv)
        ax.imshow(hs,cmap='gray',extent=[xmin,xmax,ymin,ymax])
        im = ax.imshow(img,alpha=0.6,cmap=cmap,extent=[xmin,xmax,ymin,ymax])
        print('colorbar caxis:'), print(clim)
    else:
        im = ax.imshow(img,cmap=cmap,vmin=clim[0],vmax=clim[1],extent=[xmin,xmax,ymin,ymax])
    #gdf.plot('h_te_best_fit',ax=ax,s=1,cmap='inferno',vmin=clim[0],vmax=clim[1])# brfore:('p_b',ax=ax,s=1)   # p_b   - pair-beam
    ## added shp functionality to function myself
    shp.plot(color='k',alpha=0.3, ax=ax)
    shp.boundary.plot(edgecolor="black", ax=ax)
    ###
    gdf.plot(vartoplot,ax=ax,s=1)   # p_b   - pair-beam
    plt.colorbar(im,ax=ax,extend='both',label='Elevation (m)')


def gdf_on_shp(gdf,ax,shp,cmap='inferno', colname='h_te_best_fit', legend_kwds={}): 
    """
    plot the ATL08 'COLUMN_NAME' gdf on top of the shp, with GDF coordinate system
        Parameters
    ----------
    gdf : ATL08 geodataframe
    ax : figure axes
    shp : shapefile
    cmap : colormap, optional. The default is 'inferno'.
    colname : name of column to plot, optional. The default is 'h_te_best_fit'
    legend_kwds: geopandas colorbar legend keywords as, e.g., legend_kwds={'label': "Population by Country", 'orientation': "horizontal"})

    """
    #gdf = gdf.to_crs(shp.crs)
    ## added shp functionality to function myself
    shp=shp.to_crs(gdf.crs)
    xmin,ymin,xmax,ymax = shp.total_bounds #ds.bounds
    #ndv = np.min((np.array(gdf['h_te_best_fit'])) # get no data values = nan
    clim = np.nanpercentile(gdf[colname],(2,98)) 
    ## added shp functionality to function myself
    shp.plot(color='k',alpha=0.3, ax=ax)
    shp.boundary.plot(edgecolor="black",linewidth=0.5, ax=ax)
    #gdf.plot('p_b',ax=ax,s=1)   # p_b   - pair-beam
    gdf.plot(colname,ax=ax,s=1,cmap='inferno',vmin=clim[0],vmax=clim[1],legend=True,legend_kwds=legend_kwds)# brfore:('p_b',ax=ax,s=1)   # p_b   - pair-beam
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    
    
# THIS IS CURRENTLY NOT WORKING!
def gdfsubset_on_shp(gdf,ax,shp,subset=0.01, cmap='inferno', colname='h_ph', legend_kwds={}): 
    """
    plot the ATL03 'COLUMN_NAME' gdf on top of the shp, with GDF coordinate system
        Parameters
    ----------
    gdf : ATL03 geodataframe
    ax : figure axes
    shp : shapefile
    cmap : colormap, optional. The default is 'inferno'.
    colname : name of column to plot, optional. The default is 'h_ph'
    legend_kwds: geopandas colorbar legend keywords as, e.g., legend_kwds={'label': "Population by Country", 'orientation': "horizontal"})

    
    gdf=gdf_chongtar
    subset=0.0001

    """
    #gdf = gdf.to_crs(shp.crs)
    slicestep=floor(1/subset)
    gdfs=gdf[0::slicestep]
    ## added shp functionality to function myself
    shp=shp.to_crs(gdf.crs)
    xmin,ymin,xmax,ymax = shp.total_bounds #ds.bounds
    #ndv = np.min((np.array(gdf['h_te_best_fit'])) # get no data values = nan
    clim = np.nanpercentile(gdf[colname],(2,98)) 
    ## added shp functionality to function myself
    shp.plot(color='k',alpha=0.3, ax=ax)
    shp.boundary.plot(edgecolor="black",linewidth=0.5, ax=ax)
    #gdf.plot('p_b',ax=ax,s=1)   # p_b   - pair-beam
    gdfs.plot(colname,ax=ax,s=1,cmap='inferno',vmin=clim[0],vmax=clim[1],legend=True,legend_kwds=legend_kwds)# brfore:('p_b',ax=ax,s=1)   # p_b   - pair-beam
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
