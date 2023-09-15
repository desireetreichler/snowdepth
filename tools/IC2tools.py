# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 10:19:20 2022

@author: desireet

This module contains functions to work with ICESat-2 data -
same as gda_lib_DT, but this package doesn't require rasterio.

"""


import h5py
import numpy as np
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point, Polygon, LineString
#import rasterio
#from rasterio import features
from pyproj import Proj   # assumes we only have one datum (WGMS84) 
#from math import floor
import datetime
import os
from glob import glob
#import time
import pdb

#import rasterstats as rs

# DT added
#import matplotlib.pyplot as plt


# ---------------------ATL03/ATL08 to gdf converters by Desiree -------------------------------------
# # Developed according to a ATL06 example from topolib toolbox 

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
                                except TypeError:# as e: # added this exception, as I got some type errors - temp[dataset] was 0.0
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
                            except TypeError:# as e: # added this exception, as I got some type errors - temp[dataset] was 0.0
                                pass    
                            except IndexError:# as e: # added this exception, as I got some *** IndexError: invalid index to scalar variable.
                                pass
                        except KeyError:# as e:
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


def point_convert(row):
    geom = Point(row['longitude'],row['latitude'])
    return geom



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
            df['pb'] = df['pair']*10+df['beam']
            #df['p_b'] = str(track['pair'][0])+'_'+str(track['beam'][0])
        except:  # added, to account for errors - maybe where there is onluy one data point (?)
            df['pb'] = track['pair']*10+track['beam']
            #df['p_b'] = str(track['pair'])+'_'+str(track['beam'])

        df['geometry'] = df.apply(point_convert,axis=1)
        if i==0:
            df_final = df.copy()
        else:
            df_final = pd.concat((df_final,df))
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
                #print('/gt%d%s/' % (pair, beam))
                # loop over the groups in the dataset dictionary
                temp={}
                for group in dataset_dict.keys():
                    if group in ['geolocation']: # check whether this data is available for each photon or for each 20m segment only.
                        segmentlevel = 1
                    elif group in ['heights']:
                        segmentlevel = 0
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
                                except TypeError:# as e: # added this exception, as I got some type errors - temp[dataset] was 0.0
                                    pass
                            # a few attributes have 5 columns, corresponding to land, ocean, sea ice, land ice, inland water
                            try:
                                if len(temp[dataset][0]==5):
                                    for surftypenr, surftype in enumerate(['ocean','seaice','landice','inlandwater']):
                                        temp[dataset+'_'+surftype]=temp[dataset][:,surftypenr+1]
                                    temp[dataset]=temp[dataset][:,0] # default = land, first column
                            except TypeError:# as e: # added this exception, as I got some type errors - temp[dataset] was 0.0
                                pass  
                            # some data is only available at the 20 m segment rate (or even less). Fix this by duplicating
                            if segmentlevel==1:
                                DS='/gt%d%s/%s/%s' % (pair, beam, 'geolocation', 'segment_ph_cnt')
                                segment_ph_cnt=np.array(h5f[DS])
                                #pdb.set_trace(); # for tracks with little/no actual data, 'segment_ph_cnt' values seem to sometimes be really high (500 times higher than normal), making this duplication fail
                                temp2 = np.concatenate([np.repeat(y,x) for x,y in zip(segment_ph_cnt,temp[dataset])])
                                temp[dataset]=temp2 # overwrite
                        except KeyError:# as e:
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



def get_ATL03_x_atc(h5f, pair, beam, temp):
    # calculate the along-track and across-track coordinates for ATL03 photons
        #-- data and attributes for beam gtx
    x_atc=np.zeros_like(temp['h_ph'])+np.NaN
    y_atc=np.zeros_like(temp['h_ph'])+np.NaN
                #-- ATL03 Segment ID
    DS='/gt%d%s/%s/%s' % (pair, beam, 'geolocation','segment_id')
    #Segment_ID = np.array(h5f[DS])
            #Segment_ID[gtx] = val['geolocation']['segment_id']
    #n_seg = len(Segment_ID)#[gtx])
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
    #Segment_Length = np.array(h5f[DS])
            #Segment_Length[gtx] = val['geolocation']['segment_length']
            #-- Transmit time of the reference photon
    DS='/gt%d%s/%s/%s' % (pair, beam, 'geolocation','delta_time')
    #delta_time = np.array(h5f[DS])
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
        #segment_delta_times = temp['delta_time'][idx:idx+cnt]
                    #-- Photon event lat/lon and elevation (WGS84)
        #segment_heights = temp['h_ph'][idx:idx+cnt]
        #segment_lats = temp['lat_ph'][idx:idx+cnt]
        #segment_lons = temp['lon_ph'][idx:idx+cnt]
                    #-- Along-track and Across-track distances
        distance_along_X = np.copy(dist_ph_along[idx:idx+cnt])
        distance_along_X[:c1] += Segment_Distance[j]
        distance_along_X[c1:] += Segment_Distance[j+1]
        distance_along_Y = np.copy(dist_ph_across[idx:idx+cnt])
        x_atc[idx:idx+cnt]=distance_along_X
        y_atc[idx:idx+cnt]=distance_along_Y
        
    return x_atc, y_atc
            
            
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
            # check that all data have the same length - sometimes, the multiplication with segment_ph_cnt fails. Set these to Nan.
            nrdatapts = len(track['delta_time'])
            for key in track.keys():
                if len(track[key]) != nrdatapts:
                    track[key]= np.empty_like(track['delta_time'])*np.nan
                    if v: 
                        print(f'dropped {key}: wrong nr of data points')
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
                df['pb'] = df['pair']*10+df['beam']
                #df['p_b'] = str(track['pair'][0])+'_'+str(track['beam'][0])
            except:  # added, to account for errors - maybe where there is onluy one data point (?)
                df['pb'] = track['pair']*10+track['beam']
                #df['p_b'] = str(track['pair'])+'_'+str(track['beam'])
            # df['geometry'] = df.apply(point_covert,axis=1) # this did not work here??? but for ATL08 yes? strange.
            # df['geometry'] = [Point(xy) for xy in zip(df.lon_ph, df.lat_ph)] # this works, but is slow
            if i==0:
                df_final = df.copy()
            else:
                df_final=pd.concat([df,df_final]) # possible solution? Not quite yet it seems - FutureWarning: The frame.append method is deprecated and will be removed from pandas in a future version. Use pandas.concat instead.
                #df_final = df_final.append(df)
            i = i+1
        try:
            gdf_final = gpd.GeoDataFrame(df_final,geometry=gpd.points_from_xy(x=df_final.lon_ph, y=df_final.lat_ph),crs='epsg:4326')  # changed from +init-version to avoid the upcoming warning
        except: # added this exception, as I got some errors that the df_final was not defined when there was a problem reading out height data for the granule - skip
            gdf_final = gpd.GeoDataFrame() # empty frame   
    else:
        gdf_final = gpd.GeoDataFrame() # empty frame  
    return gdf_final 


# wrapper to handle lists of granules, and add some custom parameters

def hdf2gdf(ATL_list, ATL_version=8, dataset_dict={}, aoicoords=False, filterbackground=False, utmzone=False, v=False):
    """ 
    Load data from granules and add to a single gdf. Add ID number and date
    to data rows.
    Usage: hdf2gdf(ATL_list, ATL_version=8, dataset_dict={})
    ATL_list:     a list with granule filenames to load
    ATL_version:  8 or 3, default: 8
    dataset_dict: dictionary for which datasets to add to the dataframe. 
                  Default available for ATL08 and ALT03 (see function for details).
        additionally, for ATL03: aoicoords -> [xmin,ymin,xmax,ymax] if set, the data will be clipped to that bounding box.
                      filterbackground -> default False, otherwise set to remove values (e.g. [-2, -1,0, 1], only excluded if true for BOTH land (standard) signal_conf_ph and the one for land ice)
                      v: if set to True, the function prints the current file - keeping track of progress in a loop.
                      
                      signal_conf_ph meaning/values available: Confidence level associated with each photon event selected as signal. 
                      0=noise. 1=added to allow for buffer but algorithm classifies as background; 2=low; 3=med; 4=high). This parameter 
                      is a 5xN array where N is the number of photons in the granule, and the 5 rows indicate signal finding for each
                      surface type (in order: land, ocean, sea ice, land ice and inland water). Events not associated with a specific surface
                      type have a confidence level of ­1. Events evaluated as TEP returns have a confidence level of ­2. flag_values: ­2, ­1, 0, 1, 2, 3, 4 
                      flag_meanings : possible_tep not_considered noise buffer low medium high

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
           gdf_list = [(ATL03_2_gdf(x, dataset_dict, aoicoords=aoicoords, filterbackground=filterbackground, 
                                    utmzone=utmzone, v=v)) for x in ATL_list]  # used my own, as the toolbox one caused some errors - maybe for a dataset that had just one point?
           
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

    
    
### Various tools to work with geodataframes / ICESat-2 data (Desiree) ###
#------------------------------------------------------------------------"

def concat_gdf(gdf_list):
    """
    concatanate geodataframes into 1 geodataframe
    Assumes all input geodataframes have same projection
    Inputs : list of geodataframes in same projection
    Output : 1 geodataframe containing everything having the same projection
    """
    #from https://stackoverflow.com/questions/48874113/concat-multiple-shapefiles-via-geopandas
    gdf = pd.concat([gdf for gdf in gdf_list]).pipe(gpd.GeoDataFrame)
    try:
        gdf.crs = (gdf_list[0].crs)
        if gdf.crs is None:
            gdf.crs='epsg:4326'  # sometimes the first geodataframe in a list may be empty, causing the result not to have a coordinate system.
    except: 
        print('Warning: no CRS assigned')
        pass
    return gdf
    
def makebboxgdf(min_x, min_y, max_x, max_y, crs='EPSG:4326'):
    """ Creates a gdf from bounding box coordinates.
    Source: https://github.com/ICESAT-2HackWeek/2022-snow-dem-large/commits/main/notebooks/create_bbox.ipynb
    
    Usage: makebboxgdf(min_x, min_y, max_x, max_y, crs='EPSG:4326')
        min_x etc:  corner coordinates of the bounding box.
        crs:    coordinate reference system the coordinates are in. 
    Output: polygon as a geodataframe
    """
    # create polygon coordinates from values
    coords = [(min_x, min_y), (max_x, min_y), (max_x, max_y), (min_x, max_y), (min_x, max_y)]
    # convert bbox to polygon
    poly_geom = Polygon(coords)
    # apply crs and create geopandas
    polygon = gpd.GeoDataFrame(index=[0], crs=crs, geometry=[poly_geom]) 
    
    return polygon

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

def overlapping_totalbounds(tbA, tbB):
    """ Check whether the total bounding boxes of two spatial entities overlap.
    Usage: overlapping_totalbounds(tbA, tbB) - where tbA and B are arrays with 
    total_bounds: [xmin, ymin, xmax, ymax].
    Returns True/False """
    xmin, ymin, xmax, ymax = tbA
    xminB, yminB, xmaxB, ymaxB = tbB
    result1 =  ((xmin < xmaxB) & (xmax > xminB) & (ymin < ymaxB) & (ymax > yminB) ) 
    result2 = ((xminB < xmax) & (xmaxB > xmin) & (yminB < ymax) & (ymaxB > ymin) )
    result = result1|result2
    return result

def dfxyclip(df,xminmax, yminmax, x='x',y='y'):
    """ dclip = dfxyclip(df,xminmax, yminmax, x='x',y='y') 
    clips a dataframe with coordinate columns (specify column names with 
    optional arguments x, y) to the extents defined by xminmax, yminmax 
    (two arrays of size two)."""
    dfclip= df[(df[x]>=xminmax[0]) & (df[x]<=xminmax[1]) & (df[y]>=yminmax[0]) & (df[y]<=yminmax[1])]
    return dfclip

    
    
def getIC2dates(gdf):
    """dates, dateint, currentdates, datecounts, datemonth= getIC2dates(gdf) 
    takes a geodataframe with a column delta_time, and returns:
        dates: acquisition datetime dates for each point
        dateint: integer date YYYYMMDD
        currentdates: list of datetime dates within the gdf
        datecounts: how many points of each dt
        datemonth: the acquisition months of each point"""
    
    deltatime = np.floor(gdf['delta_time'].to_numpy() /3600/24)
    # plt.hist(deltatime)
    dates =[datetime.date(2018,1,1)+datetime.timedelta(d) for d in deltatime] 
    dateint=[10000*x.year + 100*x.month +x.day for x in dates]
    currentdates, datecounts=np.unique(dateint,return_counts=True)
    #dt =[datetime.date(2018,1,1)+datetime.timedelta(d) for d in currentdates]
    datemonth=[x.month for x in dates]
    return dates, dateint, currentdates, datecounts, datemonth

def points2linestring(gf,datecol = 'dateint', beamcol = 'pb',other_cols=[]): 	
    """Function to convert ICESat-2 ATL_08 data (or ATL_08 quicklook data)
    to a gdf with lines indicating the strips with data. The tool assumes 
    (and assigns) the crs epsg:4326, i.e. lat/lon in WGS84.
    Usage: gf_lines = points2linestring(gf, datecol = 'dateint', beamcol = 'pb', other_cols=[])
    Parameters: 
        gf      - input point geodataframe
        datecol - column name indicating the overpass date. Default: 'dateint'
        beamcol - column name indicating the beam identifier (six per overpass). Default: 'pb'
        other_cols - list of other column names to keep, e.g. RGT. (The value of 
                     the first point for each line is retained). Default: None
        """
    appender=[]
    dates=gf[datecol].unique()
    beams=gf[beamcol].unique()

    for date in dates:
        for beam in beams:
            gf_sub = gf[gf[datecol]==date]
            gf_sub = gf_sub[gf_sub[beamcol]==beam]
            
            xylist = [xy for xy in zip(gf_sub.geometry.x,gf_sub.geometry.y)]
            if len(xylist)<2:
                continue
            geom=LineString(xylist)
            
            #create the geodataframe with this Multiline and the date and beam information
            gdf=gpd.GeoDataFrame(data={
            'date':[date],
            'beam':[beam],
            'geometry':[geom]})          
            # add other columns
            if len(other_cols)>0:
                for col in other_cols:
                    gdf[col]=gf_sub[col].iloc[0]
            
            appender.append(gdf)                   
    #concat the created geodataframes, and set the crs.
    gf_lines = pd.concat(appender) 
    gf_lines.crs='epsg:4326'    
    return gf_lines    

### Related to coregistration
# ----------------------------------------------------------------
def savecoreginfo(destpath, coregobj, fileprefix = 0):
    """ savecoreginfo(destpath, coregobj, fileprefix=0)
        saves the coregistration information of coregobj produced by xdem: 
        - store the coreg coefficients in a textfile in destination folder destpath, 
        - and also move produced plots stored in the working dir (if any) there. 
        - if fileprefix = 1, treat the last part of the destpath as a file prefix rather than a new subfolder.
    """
    if fileprefix:
        destpath, prefix = os.path.split(destpath)
    else:
        prefix = ''
    isdir = os.path.exists(destpath)
    if not isdir: os.mkdir(destpath)
    with open (destpath+'/'+prefix+'coreg_matrix.csv', 'w') as f: np.savetxt(f,coregobj.to_matrix(),delimiter=";")
    with open (destpath+'/'+prefix+'coreg_meta.txt', 'w') as f: f.write(str(coregobj._meta))
    currd = os.getcwd()
    fs = glob(currd+'/it*')
    for f in fs:
        fn = os.path.basename(f)
        try:
            os.rename(f,destpath+'/'+prefix+fn)
        except FileExistsError as e:
            print("replacing existing file {destpath+'/'+fn}")
            os.remove(destpath+'/'+prefix+fn)
            os.rename(f,destpath+'/'+prefix+fn)
