# -*- coding: utf-8 -*-
"""
ICESAT-2 HACKATHON - Data access tutorial (03)
--------------------------------------------------
https://github.com/ICESAT-2HackWeek/ICESat2_hackweek_tutorials

Created on Wed Nov 27 14:17:53 2019

@author: desireet

This script is an adaptation from the ICESat-2 hackathon tutorial of 2019.
"""

import numpy as np
import requests
import getpass
import socket
import json
import zipfile
import io
import math
import os
#import shutil
import pprint
import time
import geopandas as gpd
import matplotlib.pyplot as plt
import fiona
#import h5py
#import re
from glob import glob
#import pickle
#import re

# To read KML files with geopandas, we will need to enable KML support in fiona (disabled by default)
fiona.drvsupport.supported_drivers['LIBKML'] = 'rw'
from shapely.geometry import Polygon, mapping
from shapely.geometry.polygon import orient
from statistics import mean
#from requests.auth import HTTPBasicAuth
from xml.etree import ElementTree as ET


def getIC2data(uid,email,pswd, short_name,start_date, end_date, foldername, requesttype, aoilist, start_time='00:00:00', end_time='23:59:59'):
    """Function to download ICESat-2 ATLAS data. Usage:
        getIC2data(uid,email,short_name,pswd, start_date, end_date, foldername, 
                   requesttype, aoilist, start_time='00:00:00', end_time='23:59:59')
        Downloads the requested data according to the parameters and stores it 
        in the provided foldername. Requires an EarthData account. The password 
        will be prompted.
        
        uid:    EarthData user name
        email:  EarthData email address, data download links are sendt there
        pswd:   the prompt below stopped to work
        short_name:     ATL product, e.g. 'ATL08'
        start_date / end_date:     date in format '2018-01-01'
        foldername:     e.g. granules/ATL08_locationXX_date. Is created if it
                        doesn't exist.
        requesttype:    'subsetting' (spatial) or 'native' (entire granules)
        aoilist:    1.) either a four elements bounding box list in 'W,S,E,N' format
                        ['LL_lon','LL_lat','UR_lon','UR_lat']
                    2.) or a two-entry list of two lists of longitude / latitude 
                        coordinates defining a polygon (not tested)
                    3.) or a one-element list containing the path to a lat/lon
                        shapefile (not tested)
        start_time /end_time:   default set to '00:00:00'/'23:59:59'
        
        Example for Désirée:
            getIC2data('himmmelblau','desiree.treichler@geo.uio.no', 'XXXX', 'ATL08',
                       '2018-09-01','2022-12-31', 'I:/science/ICESat-2/granules/ATL08_testfolder',
                       'subsetting', ['70','39','71','40'])
        """
    """ define all input here - extracted from the script below """
    
    # go to working directory, this could possibly be skipped?
    wd=os.getcwd() 
    os.chdir(wd)
    
    
    # Enter your Earthdata Login user name
    # uid = 'himmmelblau'
    # # Enter your email address associated with your Earthdata Login account
    # email = 'desiree.treichler@geo.uio.no'
    # pswd = getpass.getpass('Earthdata Login password: ')  # hint: Cheezy
    
    # # define which data
    # short_name = 'ATL08'   # Input data set ID (e.g. ATL06) of interest here, also known as "short name".
    
    # # Input start date in yyyy-MM-dd format
    # start_date = '2018-01-01'
    # # Input start time in HH:mm:ss format
    # start_time = '00:00:00'
    # # Input end date in yyyy-MM-dd format
    # end_date = '2022-03-01'
    # # Input end time in HH:mm:ss format
    # end_time = '23:59:59'
    
    # download only new granules - functionality added 3/21 (set to 0 if all found should be downloaded)
    downloadonlynew = 1
    # compare with a list of granules supplied, and download only those (also) in the list
    comparewithlist = 0 # NOT WORKING PROPERLY, DON'T TURN ON
    
    # Creating an output folder if the folder does not already exist. All data are unzpipped and moved to that folder.
    # foldername = '../granules/ATL08_Tavildara_20220208'
    path = foldername #str(os.getcwd() + '/' + foldername)
    if not os.path.exists(path):
        os.mkdir(path)
    
    # the subset_request_params are also stored in that folder:  subset_request_params.txt
    
    # define whether you want the subset data or run a native request
    # requesttype = 'subsetting' # subsetting or native
    
    # spatial location. select AOI (1), coordinate pairs (2) or a polygon (3) and corresponding input
    # aoi value used for CMR params below
    
    if len(aoilist)==4:
        aoi = '1'
    elif len(aoilist)==2:
        aoi = '2'
    elif len(aoilist)==1:
        aoi = '3'
    
    #aoi = '1'    
    
    if aoi == '1':
         # Input lower left longitude in decimal degrees
        LL_lon = aoilist[0] #'71.3'
         # Input lower left latitude in decimal degrees
        LL_lat = aoilist[1]#'39.45'
         # Input upper right longitude in decimal degrees
        UR_lon = aoilist[2]#'71.9'
         # Input upper right latitude in decimal degrees
        UR_lat = aoilist[3]#'39.85'
        
    if aoi == '2':
            #create list of x (longitude) values in decimal degrees
        x = aoilist[0]#[]
        #create list of y (latitude) values in decimal degrees
        y = aoilist[0]#[]
        xylist = list(zip(x, y))
    
    if aoi == '3': #Use geopandas to read in polygon file. A shapefile or geojson or other vector-based spatial data format could be substituted here.
    # see also further notes in the code below for alternatives for how to make a polygon!    
    
        # kml_filepath = str(os.getcwd() + '/pine_island_glims/glims_polygons.kml')
        shp_filepath = aoilist[0]#'../../DEM/eksport_183900_20190923/dom1/metadata/dom1_Klippefil.shp'
        
        #Return a GeoDataFrame object
        gdf = gpd.read_file(shp_filepath) # same with kml example , just substitute path
        gdf.head()
        # Note: the above example is in norwegian UTM, thus: set projection and convert to lat/lon
        #gdf.crs={'proj':'utm +zone=33 +ellips=WGS84 +datum=WGS84 +units=m +no_defs'}
        gdf = gdf.to_crs('epsg:4326')
        
        # We need to get from the geopandas GeoDataFrame object to an input that is readable by CMR.
        #Integer position based indexing of GeoDataFrame object to get it into a shapeply geometry object.
        poly = gdf.iloc[0].geometry
        
        # Simplify polygon. The larger the tolerance value, the more simplified the polygon.
        if max(np.shape(poly.exterior.coords)) > 20:  # check the size, simplify if above max nr (wild guess) of point pairs
            poly = poly.simplify(0.05, preserve_topology=False)
        # Orient counter-clockwise
        poly = orient(poly, sign=1.0)
        print(poly)
        #Format dictionary to polygon coordinate pairs for CMR polygon filtering
        polygon = ','.join([str(c) for xy in zip(*poly.exterior.coords.xy) for c in xy])
    
    
    # if comparewithlist:   # NOT WORKING
    #             # supply the list
    #     with open('chongtar/ROItracks'+'.pkl', 'rb') as f:  # Python 3: open(..., 'wb')
    #         ROItracks = pickle.load(f)
    #     orderlist= ['ATL03'+re.split('_004', re.split('ATL08',f)[2])[0] for f in ROItracks]  # remove everything up to ATL08 and the version, and add ATL03 instead
    
        
    
    """
    ______________________________________________________________________________
    Here are the steps you will learn in this tutorial (and that the script below does):
    
       1) Set up NASA Earthdata Login authentication for direct data access.
       2) Obtain data set metadata.
       3) Input data set search criteria using spatial and temporal filters. 
       4) Explore different methods of specifying spatial criteria including lat/lon bounds,
       polygon coordinate pairs, and polygon input from a Shapefile or KML.  
       5) Search for matching files and receive information about file volume.
       6) Obtain information about subsetting and reformatting capabilities for a specfic data set.
       7) Configure data request by "chunking" request by file number. 
       8) Submit a request and monitor the status of the request.
       9) Compare subsetted and original, "native" data outputs.
    
    """
    
    
    """ 1) Login data - generate token
        We will generate a token needed in order to access data using your Earthdata Login credentials, and we will 
        apply that token to the following queries. If you do not already have an Earthdata Login account, go to 
        http://urs.earthdata.nasa.gov to register. Your password will be prompted for privacy.
    """
    
    # Earthdata Login credentials
    # input above
    
    # Request token from Common Metadata Repository using Earthdata credentials
    token_api_url = 'https://cmr.earthdata.nasa.gov/legacy-services/rest/tokens'
    hostname = socket.gethostname()
    ip = socket.gethostbyname(hostname)
    
    data = {
        'token': {
            'username': uid,
            'password': pswd,
            'client_id': 'NSIDC_client_id',
            'user_ip_address': ip
        }
    }
    headers={'Accept': 'application/json'}
    response = requests.post(token_api_url, json=data, headers=headers)
    token = json.loads(response.content)['token']['id']
    print('Token:', token)
    
    
    """ 2) Select a data set of interest and explore resources available through NSIDC.
    Let's begin discovering ICESat-2 data by first inputting the data set of interest.
    See the ICESat-2 Data Sets page for a list of all ICESat-2 data set titles and IDs. Below we will input data 
    set ID ATL08, which is the ID for the ATLAS/ICESat-2 L3A Land and Canopy Height data set. (Original in tutorial: ATL06)
    > See an overview over datasets here: https://nsidc.org/data/icesat-2/data-sets
    with links to
    > user guides here: https://nsidc.org/data/atl08?qt-data_set_tabs=3#qt-data_set_tabs (and similar)
    > data dictionaries here: https://nsidc.org/sites/nsidc.org/files/technical-references/ATL08-data-dictionary-v001.pdf
    """
    # Input data set ID (e.g. ATL06) of interest here, also known as "short name".
    # short_name = input above
    
    # Determine the number and size of granules available within a time range and location.
    
    # Let's explore information about our data set. The Common Metadata Repository is queried to explore this information.
    
    search_params = {
        'short_name': short_name
    }
    
    cmr_collections_url = 'https://cmr.earthdata.nasa.gov/search/collections.json'
    response = requests.get(cmr_collections_url, params=search_params)
    results = json.loads(response.content)
    print('--- metadata for query: ---')
    pprint.pprint(results)
    print('--- ------------------- ---')
    
    
    #  We'll start by determining the most recent version number of our data set. 
    # Find all instances of 'version_id' in metadata and print most recent version number
    
    versions = [i['version_id'] for i in results['feed']['entry']]
    print(versions)
    latest_version = max(versions)
    print('latest version of dataset: ', latest_version)           # 002
    
    
    """ 3./4.) We will now find out how many data granules (files) exist over an area and time of interest.
    """
    
    # Input temporal range - according to the CMR API documentation:
    # https://cmr.earthdata.nasa.gov/search/site/docs/search/api.html#g-temporal
    # date range provided above
    
    temporal = start_date + 'T' + start_time + 'Z' + ',' + end_date + 'T' + end_time + 'Z'
    print('input temporal range time string: ', temporal)
    
    # There are three different options for inputting an area of interest to be applied to our granule search:
    # 1) Bounding Box 
    # 2) Polygon coordinate pairs 
    # 3) Spatial file input, including Esri Shapefile or KML/KMZ. 
    
     # 1) commented example code for 1: # # Bounding Box spatial parameter in 'W,S,E,N' format
    if aoi == '1':
        bounding_box = LL_lon + ',' + LL_lat + ',' + UR_lon + ',' + UR_lat
         # aoi value used for CMR params below
        aoi = '1'
        print(bounding_box)
    
    # # 2) commented example for polygon coordinate pair spatial parameter
    if aoi == '2':
        # # Polygon points need to be provided in counter-clockwise order. The last point should match the first point to close the polygon. 
        # # Input polygon coordinates as comma separated values in longitude latitude order, i.e. lon1, lat1, lon2, lat2, lon3, lat3, and so on.
        polygon = ','.join(map(str, list(sum(xylist, ()))))
        print(polygon)
        # # aoi value used for CMR params below
        aoi = '2'
    
    # 3) Use geopandas to read in polygon file. A shapefile or geojson or other vector-based spatial data format could be substituted here.
    if aoi == '3':
    
        if 0: # plot
            # simple visualisation 
            # Load "Natural Earth” countries dataset, bundled with GeoPandas
            world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
            
            # Overlay glacier outline
            f, ax = plt.subplots(1, figsize=(12, 6))
            world.plot(ax=ax, facecolor='lightgray', edgecolor='gray')
            gdf.plot(ax=ax, cmap='Set2')
            ax.set_ylim([50, 70])
            ax.set_xlim([0,20]);
    
        # We need to get from the geopandas GeoDataFrame object to an input that is readable by CMR.
        # CMR API documentation: "Polygon points are provided in counter-clockwise order. The last point 
        # should match the first point to close the polygon. The values are listed comma separated in 
        # longitude latitude order, i.e. lon1, lat1, lon2, lat2, lon3, lat3, and so on.`
        # -> simplify and reorder the GeoDataFrame object using the shapely package and convert the 
        # object back to a dictionary to be applied to the CMR polygon parameter. Simplification is 
        # needed in order to pass a reasonable request length to CMR. You may need to modify the 
        # simplification tolerance depending on the number of points of your polygon.
        # this is done above.
        
        # # Alternative option for spatial file: Post file to OGR service for conversion to CMR polygon format 
        # # Spatial file input, including Esri Shapefile or KML/KMZ - commented, as not used here: 
        # 
        # # POST shapefile or KML polygon to OGR for geojson conversion
        # url = 'http://ogre.adc4gis.com/convert'
        # shapefile = kml_filepath
        # files = {'upload': open(shapefile, 'rb')}
        # r = requests.post(url, files=files)
        # results = json.loads(r.content)
        # # Results is a dictionary representing a feature collection. List coordinates from the Polygon feature:
        # polygon_list = list(results['features'][0]['geometry']['coordinates'][0])     
        # # Remove z value from polygon list
        # for i in range(len(polygon_list)):
        #     del polygon_list[i][2] 
        # # Create shapely Polygon object for simplification and counter-clockwise ordering for CMR filtering
        # poly = Polygon(tuple(polygon_list))
        
        
        # Simplify polygon. The larger the tolerance value, the more simplified the polygon.
        # this is also done above
    
        print(polygon)
        # # aoi value used for CMR params below
        aoi = '3'
    
    
    
    
    """ 5) We will now populate dictionaries to be applied to our search query below based on spatial 
        and temporal inputs. For additional search parameters, see the The Common Metadata Repository API documentation.
    """
    
    # Create CMR parameters used for granule search. Modify params depending on bounding_box or polygon input.
    if aoi == '1':
    # bounding box input:
        search_params = {
        'short_name': short_name,
        'version': latest_version,
        'temporal': temporal,
        'page_size': 100,
        'page_num': 1,
        'bounding_box': bounding_box
        }
    else:
    # If polygon input (either via coordinate pairs or shapefile/KML/KMZ):
        search_params = {
        'short_name': short_name,
        'version': latest_version,
        'temporal': temporal,
        'page_size': 100,
        'page_num': 1,
        'polygon': polygon,
        }
    print('CMR search parameters: ', search_params)
    
    # Query number of granules using our criteria (paging over results)
    granule_search_url = 'https://cmr.earthdata.nasa.gov/search/granules'
    granules = []
    while True:
        response = requests.get(granule_search_url, params=search_params, headers=headers)
        results = json.loads(response.content)
        if len(results['feed']['entry']) == 0:
            # Out of results, so break out of loop
            break
        # Collect results and increment page_num
        granules.extend(results['feed']['entry'])
        search_params['page_num'] += 1
    
    # Get number of granules over my area and time of interest
    print('nr of granules over my area and time of interest: ', len(granules))
    
    # We could view this in the NASA Earthdata Search web interface, which relies on the 
    # same metadata, although their simplified polygon may differ slightly. With the same 
    # search criteria applied, we can view the  granules of ATL08 over the AOI.
    
    # added in 3/2021: check whether we already have these granules.
    if downloadonlynew:
        granules_all = granules
        
        existinglist=glob(foldername +'/*')
        existinglist=[ r[(len(foldername)+11):] for r in existinglist]
              
        notthere = [item for item in granules_all if np.setdiff1d(item['producer_granule_id'], existinglist)]
        print(len(granules_all)-len(notthere), ' of desired granules already downloaded. Skipping these. ')
        granules=notthere
    
    # # added in 4/2021: check which available granules are in the orderlist (e.g., to order some ATL03 that we identified as useful)
    # if comparewithlist:
    #     granules_all = granules
    #     there = [item for item in granules_all if re.split('_'+latest_version,item['producer_granule_id'])[0] in orderlist]
    #     print(len(orderlist)-len(there), ' of desired granules not found in granule list. Skipping these.')
    #     granules=there
    
    print('Granules to download: ',len(granules))
    
    # added in 8/2012, as this granule selection was not transferred to the API request: 
    granlist= [granules[x]['producer_granule_id'] for x in np.arange(0,len(granules))]
    # convert to a string input as with the coverage list?
    # using list comprehension 
    granlist2 = ', '.join([str(elem) for elem in granlist]) 
    #granlist2=granlist[0]
     
    
    # Now query the average size of those granules, as well as the total volume
    granule_sizes = [float(granule['granule_size']) for granule in granules]
    # Average size of granules in MB
    print('avg granule size, in MB: ', mean(granule_sizes))
    # Total volume in MB
    print('total size, in MB: ',sum(granule_sizes))
    # Although subsetting, reformatting, or reprojecting can alter the size of the granules, 
    # this "native" granule size can still be used to guide us towards the best download method 
    # to pursue, which we will come back to later on in this tutorial.
    
    
    
    """ 6) Select the subsetting and reformatting services enabled for your data set of interest.
        The NSIDC DAAC supports customization services on many of our NASA Earthdata mission 
        collections. Reformatting and subsetting are available on all Level-2 and -3 ICESat-2 data 
        sets. Let's discover the specific service options supported for this data set and select 
        which of these services we want to request.
    """  
    # We will start by querying the service capability to gather and select customization options.
    capability_url = f'https://n5eil02u.ecs.nsidc.org/egi/capabilities/{short_name}.{latest_version}.xml'
    print('capability_url: ', capability_url)
    
    # Create session to store cookie and pass credentials to capabilities url
    session = requests.session()
    s = session.get(capability_url)
    response = session.get(s.url,auth=(uid,pswd))
    
    root = ET.fromstring(response.content)
    
    # From the service capability XML, we can collect lists with each service option to gather service information.
    # collect lists with each service option:
    subagent = [subset_agent.attrib for subset_agent in root.iter('SubsetAgent')]
    
    # variable subsetting - list of all variables
    variables = [SubsetVariable.attrib for SubsetVariable in root.iter('SubsetVariable')]  
    variables_raw = [variables[i]['value'] for i in range(len(variables))]
    variables_join = [''.join(('/',v)) if v.startswith('/') == False else v for v in variables_raw] 
    variable_vals = [v.replace(':', '/') for v in variables_join]
    
    # reformatting - list of available options
    formats = [Format.attrib for Format in root.iter('Format')]
    format_vals = [formats[i]['value'] for i in range(len(formats))]
    format_vals.remove('')
    
    # reprojection only applicable on ICESat-2 L3B products, yet to be available!
    
    if 0:  # this caused troubles from autumn 2020, thus disabled. There shouldn't be reformatting options for level 2 products anyway.
        # reformatting options that support reprojection
        normalproj = [Projections.attrib for Projections in root.iter('Projections')] # edit: removed s, before: root.iter('Projection')
        normalproj_vals = []
        normalproj_vals.append(normalproj[0]['normalProj'])
        format_proj = normalproj_vals[0].split(',')
        format_proj.remove('')
        format_proj.append('No reformatting')
        
        #reprojection options - none so far
        projections = [Projection.attrib for Projection in root.iter('Projection')]
        proj_vals = []
        for i in range(len(projections)):
            if (projections[i]['value']) != 'NO_CHANGE' :
                proj_vals.append(projections[i]['value'])
                
        # reformatting options that do not support reprojection
        no_proj = [i for i in format_vals if i not in format_proj]
    
    
    
    # Let's confirm that subset services exist for our data set by reviewing the subagent list. 
    # If the list contains service information, we know that services are available. If not, we 
    # need to set the agent API parameter to NO to indicate that our request will bypass the 
    # subsetter. This will quickly send back the data "natively" without any customization applied.
    print('subagent: ',subagent)
    if len(subagent) < 1 :
        agent = 'NO'
        
        # output: [{'id': 'ICESAT2', 'spatialSubsetting': 'true', 'spatialSubsettingShapefile': 'true', 'temporalSubsetting': 'true', 'type': 'both', 'maxGransSyncRequest': '100', 'maxGransAsyncRequest': '2000'}]
        
    # More information is contained in the subagent list, including the maximum number of granules 
    # that we can request per order depending on our configuration. We'll come back to these options below.
    
    
    # We'll begin populating the subsetting and reformatting parameters used for our NSIDC API 
    #    request. In addition to the CMR information we queried above, the NSIDC API accepts Key-Value-Pairs 
    #    (KVPs) for subsetting and reformatting services.
    
    # Let's start with spatial subsetting. Recall that there are three options to filter our search results by spatial constraint:
    # 1) Bounding Box: Corresponding to the CMR bounding_box KVP
    # 2) Polygon coordinate pairs: Corresponding to the CMR polygon KVP
    # 3) Spatial file input, including Esri Shapefile or KML/KMZ: We simplified the file input to also be read by the CMR polygon KVP
    
    # We see above that spatialSubsetting is true and spatialSubsettingShapefile is true. Therefore the same 
    # filtering options can be applied to our subset constraint, with unique KVPs for the subsetting service. 
    # Note that both the search bounding_box KVP and bbox subsetting KVP need to be included in order to search 
    # and subset against the same granules.
    # 1) Bounding Box: bbox subset KVP
    # 2) Polygon coordinate pairs: bounding_shape subset KVP in GeoJSON format.
    # 3) Spatial file input: The file can be read directly by the subsetter without 
    #    simplification. This file will be posted to the API endpoint, so we don't need to specify an additional subset KVP here.
        
    # -> Because we're pursuing option 3), we don't need to provide an additional subset parameter. 
    
    ## Commented code for options 1/2:
    ## Bounding box subsetting (bbox) in same format as bounding_box
    if aoi == '1':
        bbox = bounding_box
    
    ## Polygon coordinate pair subsetting in GeoJSON format. Or for simplicity, get polygon bounds to be used as bounding box input
    if aoi == '2':
        ## Create shapely Polygon object from x y list
        p = Polygon(tuple(xylist))
            ## Extract the point values that define the perimeter of the polygon
        bounds = p.bounds
        bbox = ','.join(map(str, list(bounds)))
    
    
    # Temporal subsetting is next, since we saw above that temporalSubsetting is true. 
    # We filtered data over XXX-XXX and we can also subset the data to those dates if desired.
    # The time KVP is used to subset temporally. This can be entered in the following formats:
    # time=yyyy-mm-dd,yyyy-mm-dd
    # time=yyy-mm-ddThh:MM:ss,yyy-mm-ddThh:MM:ss
    
    # Temporal subsetting KVP
    timevar = start_date + 'T' + start_time + ',' + end_date + 'T' + end_time
    print('temporal subsetting: ', timevar)
    
    # Next, let's explore the refrmatting and reprojection options available.
    print('exploting reformatting and reprojection options: format_vals, proj_vals')
    print(format_vals)
    # These options can be inputted into the API request exactly as printed in the list, with quotes 
    # removed, using the format= Key-Value-Pair. For example:
    # format=TABULAR_ASCII
    # We will be exploring the data in its native HDF5 format so we won't pursue this option in this tutorial.
    # print(proj_vals)   # commented as this is currently disabled above
    # none, currently - will be from L3B
    print('nr of variables, stored in coveragelist: ', len(variable_vals)) # nr of variables
    coveragelist=variable_vals
    #coveragelist = pprint.pprint(variable_vals) # print all
    
    # And we can enter a list of variables to subset separated by comma using the coverage key. All forward 
    # slashes need to be included to indicate HDF group hierarchy.
     
    # coverage = '/ancillary_data/atlas_sdp_gps_epoch,\
    # /gt1l/land_ice_segments/atl06_quality_summary,\
    # /gt1l/land_ice_segments/delta_time'
    
    # in the script playground/data_access_tutorial.py there is the full list for WHICH DATASET SHORT NAME?? I think it's 08
    #coveragelist = ['/gt1l/land_segments/sigma_h',
    # ...
    
    # convert to a string input 
    # using list comprehension 
    coverage = ',\ '.join([str(elem) for elem in coveragelist]) 
    print('by default, we select ALL these values. Adapt the script if you want it different...')
    
    # coverage = '/ancillary_data/atlas_sdp_gps_epoch,\
    # /gt1l/land_ice_segments/atl06_quality_summary,\
    # /gt1l/land_ice_segments/delta_time'
      
    # seems to be OK now
    # print(coverage)  
    
    
    
    """ 7) Request data from the NSIDC data access API.
        We will now set up our data download request. The data access and service API (labeled EGI below) 
        incorporates the CMR parameters that we explored above, plus customization service parameters as well 
        as a few configuration parameters.
    """
    # As described above, the API is structured as a URL with a base plus individual key-value-pairs (KVPs) 
    # separated by ‘&’. The base URL of the NSIDC API is: 
    # Set NSIDC data access base URL
    base_url = 'https://n5eil02u.ecs.nsidc.org/egi/request'
    
    # Let's go over the configuration parameters:
    #    request_mode:   is "synchronous" by default, meaning that the request relies on a direct, 
    #                    continous connection between you and the API endpoint. Outputs are directly 
    #                   downloaded, or "streamed" to your working directory. For this tutorial, we will 
    #                   set the request mode to asynchronous, which will allow concurrent requests to be 
    #                   queued and processed without the need for a continuous connection.
    #                   Use the streaming request_mode with caution: While it can be beneficial to stream 
    #                   outputs directly to your local directory, note that timeout errors can result 
    #                   depending on the size of the request, and your request will not be queued in the 
    #                   system if NSIDC is experiencing high request volume. For best performance, I recommend 
    #                   setting page_size=1 to download individual outputs, which will eliminate extra time
    #                    needed to zip outputs and will ensure faster processing times per request. An example 
    #                   streaming request loop is available at the bottom of the tutorial below.
    #    page_size:      
    #    page_num:  Recall that we queried the total number and volume of granules prior to applying customization 
    #               services. page_size and page_num can be used to adjust the number of granules per request up to 
    #               a limit of 2000 granules for asynchronous, and 100 granules for synchronous (streaming). For now, 
    #               let's select 10 granules to be processed in each zipped request. For ATL06, the granule size can 
    #               exceed 100 MB so we want to choose a granule count that provides us with a reasonable zipped download size.
    
    # Set number of granules requested per order, which we will initially set to 10.
    page_size = 10
    if short_name == 'ATL03': # I had up to 7 individual zip files for ATL03 data. Seems a better idea to cut down the page size there!
        page_size = 3
    #Determine number of pages basd on page_size and total granules. Loop requests by this value
    page_num = math.ceil(len(granules)/page_size)
    #Set request mode. 
    request_mode = 'async'
    #Create config dictionary
    config_params = {
        'request_mode': request_mode, 
        'page_size': page_size,  
        'token': token, 
        'email': email,   
    }
    # Determine how many individual orders we will request based on the number of granules requested
    print('given our page_size of ', page_size, ' , we will request ', page_num , ' individual orders.')
    print('--- ------------------- ---')
    
    # After all of these KVP inputs, what does our request look like? Here's a summary of all possible KVPs that
    # we explored, both for CMR searching and for the subsetter:
    
    #CMR search keys (search_params dictionary):¶
    #    short_name=
    #    version=
    #    temporal=
    #    bounding_box=
    #    polygon=
    #
    #Customization service keys:
    
    #    time=
    #    bbox=
    #    bounding_shape=
    #    format=
    #    projection=
    #    projection_parameters=
    #    Coverage=
    #
    #No customization (access only):
    
    #    agent=
    #    include_meta=
    #        Y by default. N for No metadata requested.
    #
    #Request configuration keys:
    
    #    request_mode=
    #    page_size=
    #    page_num=
    #    token=
    #    email=
    
    # If we were to create an API request based on our request parameters , here's what we end up with:
    #Print API base URL + request parameters
    try: 
        API_request = f'{base_url}?short_name={short_name}&version={latest_version}&temporal={temporal}&time={timevar}&bounding_box={bounding_box}&Coverage={coverage}&request_mode={request_mode}&page_size={page_size}&page_num={page_num}&token={token}&email={email}'
        print('API request')
        #print('API request: ', API_request)
    except: print('API request not possible, but subset/native should be OK... testing these')
    
    # Adding customization parameter dictionary 
    custom_params = {
        'time': timevar,
        'Coverage': coverage,
        'bbox': bbox} 
    
    
    # Creating final request parameter dictionary with search, config, and customization parameters.
    subset_request_params = {**search_params, **config_params, **custom_params}
    print('subset request params (print disabled)')
    #print('subset request params: ', subset_request_params)
    print('--- ------------------- ---')
    
    
    # We'll also request the same data but without any subsetting services applied. Let's create another "no processing" dictionary to specify agent=NO for no processing (native data request) and option to include metadata, along with a new request parameter dictionary with CMR search and configuration parameters.
    no_processing_params = {
        'agent' : 'NO',
        'include_meta' : 'Y'}
    # Creating additional final request parameter dictionary with search, config, and no processing parameters.
    native_request_params = {**search_params, **config_params, **no_processing_params}
    print('alternative - native request params: ', native_request_params)
    print('--- ------------------- ---')
    
    
    
    """ 8) Request Data
        Finally, we'll download the data directly to this notebook directory in a new 
        Outputs folder. The progress of each order will be reported.
    """    
    
    # foldername is defines above
    
    # First we'll submit our native request without subsetting services:
    
    # Request data service for each page number, and unzip outputs
    if requesttype == 'native':
        for i in range(page_num):
            page_val = i + 1
            print('Order: ', page_val)
            native_request_params.update( {'page_num': page_val} )
            
        # For all requests other than spatial file upload, use get function
            request = session.get(base_url, params=native_request_params)
            
            print('Request HTTP response: ', request.status_code)
        
        # Raise bad request: Loop will stop for bad response code.
            request.raise_for_status()
            print('Order request URL: ', request.url)
            esir_root = ET.fromstring(request.content)
            #print('Order request response XML content: ', request.content)  # this is potentially very long
        
        #Look up order ID
            orderlist = []   
            for order in esir_root.findall("./order/"):
                orderlist.append(order.text)
            orderID = orderlist[0]
            print('order ID: ', orderID)
        
        #Create status URL
            statusURL = base_url + '/' + orderID
            print('status URL: ', statusURL)
        
        #Find order status
            request_response = session.get(statusURL)    
            print('HTTP response from order response URL: ', request_response.status_code)
            
        # Raise bad request: Loop will stop for bad response code.
            request_response.raise_for_status()
            request_root = ET.fromstring(request_response.content)
            statuslist = []
            for status in request_root.findall("./requestStatus/"):
                statuslist.append(status.text)
            status = statuslist[0]
            print('Data request ', page_val, ' is submitting...')
            print('Initial request status is ', status)
        
        #Continue loop while request is still processing
            while status == 'pending' or status == 'processing': 
                print('Status is not complete. Trying again.')
                time.sleep(10)
                loop_response = session.get(statusURL)
        
        # Raise bad request: Loop will stop for bad response code.
                loop_response.raise_for_status()
                loop_root = ET.fromstring(loop_response.content)
        
        #find status
                statuslist = []
                for status in loop_root.findall("./requestStatus/"):
                    statuslist.append(status.text)
                status = statuslist[0]
                print('Retry request status is: ', status)
                if status == 'pending' or status == 'processing':
                    continue
        
        #Order can either complete, complete_with_errors, or fail:
        # Provide complete_with_errors error message:
            if status == 'complete_with_errors' or status == 'failed':
                messagelist = []
                for message in loop_root.findall("./processInfo/"):
                    messagelist.append(message.text)
                print('error messages:')
                pprint.pprint(messagelist)
        
        # Download zipped order if status is complete or complete_with_errors
            if status == 'complete' or status == 'complete_with_errors':
                downloadURL = 'https://n5eil02u.ecs.nsidc.org/esir/' + orderID + '.zip'
                print('Zip download URL: ', downloadURL)
                print('Beginning download of zipped output...')
                zip_response = session.get(downloadURL)
                # Raise bad request: Loop will stop for bad response code.
                zip_response.raise_for_status()
                with zipfile.ZipFile(io.BytesIO(zip_response.content)) as z:
                    z.extractall(path)
                print('Data request', page_val, 'is complete.')
            else: print('Request failed.')
        
    # Let's run our request loop again, this time with subsetting services applied. We will post the KML file directly to the API:   
    if requesttype == 'subsetting' and comparewithlist:
            # an individual request for each granule, and unzip outputs
            # Adding customization parameter dictionary 
            
        config_params = {
            'request_mode': request_mode, 
            'page_size': page_size,  
            'token': token, 
            'email': email} 
    
        for i in np.arange(0,len(granlist)):
            page_val = i + 1
            print('Order: ', page_val)
            
            custom_params = {
                'producer_granule_id': granlist[i]} 
    
    
            # Creating final request parameter dictionary with search, config, and customization parameters.
            subset_request_params = {**search_params, **config_params, **custom_params}
            
        # Post polygon to API endpoint for polygon subsetting to subset based on original, non-simplified KML file
            if aoi == '3':
                shape_post = {'shapefile': open(shp_filepath, 'rb')}
                request = session.post(base_url, params=subset_request_params, files=shape_post) 
            else:
            # # FOR ALL OTHER REQUESTS THAT DO NOT UTILIZED AN UPLOADED POLYGON FILE, USE A GET REQUEST INSTEAD OF POST:
                request = session.get(base_url, params=subset_request_params)
            
            print('Request HTTP response: ', request.status_code)
        
        # Raise bad request: Loop will stop for bad response code.
            request.raise_for_status()
            # print('Order request URL: ', request.url)   # this is very long
            esir_root = ET.fromstring(request.content)
            print('Order request response XML content - print suppressed (too wordy)')
            #print('Order request response XML content: ', request.content)
        
        # Look up order ID
            orderlist = []   
            for order in esir_root.findall("./order/"):
                orderlist.append(order.text)
            orderID = orderlist[0]
            print('order ID: ', orderID)
        
        # Create status URL
            statusURL = base_url + '/' + orderID
            print('status URL: ', statusURL)
        
        # Find order status
            request_response = session.get(statusURL)    
            print('HTTP response from order response URL: ', request_response.status_code)
            
        # Raise bad request: Loop will stop for bad response code.
            request_response.raise_for_status()
            request_root = ET.fromstring(request_response.content)
            statuslist = []
            for status in request_root.findall("./requestStatus/"):
                statuslist.append(status.text)
            status = statuslist[0]
            print('Data request ', page_val, ' is submitting...')
            print('Initial request status is ', status)
        
        # Continue to loop while request is still processing
            while status == 'pending' or status == 'processing': 
                print('Status is not complete. Trying again.')
                time.sleep(30)
                loop_response = session.get(statusURL)
        
        # Raise bad request: Loop will stop for bad response code.
                loop_response.raise_for_status()
                loop_root = ET.fromstring(loop_response.content)
        
        # Find status
                statuslist = []
                for status in loop_root.findall("./requestStatus/"):
                    statuslist.append(status.text)
                status = statuslist[0]
                print('Retry request status is: ', status)
                if status == 'pending' or status == 'processing':
                    continue
        
        # Order can either complete, complete_with_errors, or fail:
        # Provide complete_with_errors error message:
            if status == 'complete_with_errors' or status == 'failed':
                messagelist = []
                for message in loop_root.findall("./processInfo/"):
                    messagelist.append(message.text)
                print('error messages:')
                pprint.pprint(messagelist)
        
    # original: does not consider that there could be several zip files!    
    #    # Download zipped order if status is complete or complete_with_errors
    #        if status == 'complete' or status == 'complete_with_errors':
    #            downloadURL = 'https://n5eil02u.ecs.nsidc.org/esir/' + orderID + '.zip'
    #            print('Zip download URL: ', downloadURL)
    #            print('Beginning download of zipped output...')
    #            zip_response = session.get(downloadURL)
    #            # Raise bad request: Loop will stop for bad response code.
    #            zip_response.raise_for_status()
    #            with zipfile.ZipFile(io.BytesIO(zip_response.content)) as z:
    #                z.extractall(path)
    #            print('Data request', page_val, 'is complete.\n --------')
    #        else: print('Request failed.\n --------')
            
        # Download zipped order if status is complete or complete_with_errors
            if status == 'complete' or status == 'complete_with_errors':
                downloadURL = 'https://n5eil02u.ecs.nsidc.org/esir/' + orderID + '.zip'
                print('Zip download URL: ', downloadURL)
                print('Beginning download of zipped output...')
                zip_response = session.get(downloadURL)
                # Raise bad request: Loop will stop for bad response code.
                zip_response.raise_for_status()
                with zipfile.ZipFile(io.BytesIO(zip_response.content)) as z:
                    z.extractall(path)
                if short_name == 'ATL03':
                    # try whether there were further zip files!
                    for zipnr in np.arange(2,page_size):
                       downloadURL = 'https://n5eil02u.ecs.nsidc.org/esir/' + orderID + '.zip?' + str(zipnr)
                       try:
                           print('Try - beginning download of zip download URL: ', downloadURL)
                           zip_response = session.get(downloadURL)
                            # Raise bad request: Loop will stop for bad response code.
                           zip_response.raise_for_status()
                           with zipfile.ZipFile(io.BytesIO(zip_response.content)) as z:
                               z.extractall(path)
                       except: 
                           print('failed, assuming zip download URL does not exist.')
                           break # breaks out of the for loop, ignoring further pages
                print('Data request', page_val, 'is complete.\n --------')
            else: print('Request failed.\n --------')    
        
    print('----------')
    print('DONE. Please double-check whether all files are downloaded and unzipped.')    
    print('Downloading several zip files per page currently only attempted for ATL03!')    
    
    
    if requesttype == 'subsetting' and not comparewithlist:
            # Request data service for each page number, and unzip outputs
        
        for i in range(page_num):
            page_val = i + 1
            print('Order: ', page_val)
            subset_request_params.update( {'page_num': page_val} )
            
        # Post polygon to API endpoint for polygon subsetting to subset based on original, non-simplified KML file
            if aoi == '3':
                shape_post = {'shapefile': open(shp_filepath, 'rb')}
                request = session.post(base_url, params=subset_request_params, files=shape_post) 
            else:
            # # FOR ALL OTHER REQUESTS THAT DO NOT UTILIZED AN UPLOADED POLYGON FILE, USE A GET REQUEST INSTEAD OF POST:
                request = session.get(base_url, params=subset_request_params)
            
            print('Request HTTP response: ', request.status_code)
        
        # Raise bad request: Loop will stop for bad response code.
            request.raise_for_status()
            # print('Order request URL: ', request.url)   # this is very long
            esir_root = ET.fromstring(request.content)
            print('Order request response XML content - print suppressed (too wordy)')
            #print('Order request response XML content: ', request.content)
        
        # Look up order ID
            orderlist = []   
            for order in esir_root.findall("./order/"):
                orderlist.append(order.text)
            orderID = orderlist[0]
            print('order ID: ', orderID)
        
        # Create status URL
            statusURL = base_url + '/' + orderID
            print('status URL: ', statusURL)
        
        # Find order status
            request_response = session.get(statusURL)    
            print('HTTP response from order response URL: ', request_response.status_code)
            
        # Raise bad request: Loop will stop for bad response code.
            request_response.raise_for_status()
            request_root = ET.fromstring(request_response.content)
            statuslist = []
            for status in request_root.findall("./requestStatus/"):
                statuslist.append(status.text)
            status = statuslist[0]
            print('Data request ', page_val, ' is submitting...')
            print('Initial request status is ', status)
        
        # Continue to loop while request is still processing
            while status == 'pending' or status == 'processing': 
                print('Status is not complete. Trying again.')
                time.sleep(30)
                loop_response = session.get(statusURL)
        
        # Raise bad request: Loop will stop for bad response code.
                loop_response.raise_for_status()
                loop_root = ET.fromstring(loop_response.content)
        
        # Find status
                statuslist = []
                for status in loop_root.findall("./requestStatus/"):
                    statuslist.append(status.text)
                status = statuslist[0]
                print('Retry request status is: ', status)
                if status == 'pending' or status == 'processing':
                    continue
        
        # Order can either complete, complete_with_errors, or fail:
        # Provide complete_with_errors error message:
            if status == 'complete_with_errors' or status == 'failed':
                messagelist = []
                for message in loop_root.findall("./processInfo/"):
                    messagelist.append(message.text)
                print('error messages:')
                pprint.pprint(messagelist)
        
    # original: does not consider that there could be several zip files!    
    #    # Download zipped order if status is complete or complete_with_errors
    #        if status == 'complete' or status == 'complete_with_errors':
    #            downloadURL = 'https://n5eil02u.ecs.nsidc.org/esir/' + orderID + '.zip'
    #            print('Zip download URL: ', downloadURL)
    #            print('Beginning download of zipped output...')
    #            zip_response = session.get(downloadURL)
    #            # Raise bad request: Loop will stop for bad response code.
    #            zip_response.raise_for_status()
    #            with zipfile.ZipFile(io.BytesIO(zip_response.content)) as z:
    #                z.extractall(path)
    #            print('Data request', page_val, 'is complete.\n --------')
    #        else: print('Request failed.\n --------')
            
        # Download zipped order if status is complete or complete_with_errors
            if status == 'complete' or status == 'complete_with_errors':
                downloadURL = 'https://n5eil02u.ecs.nsidc.org/esir/' + orderID + '.zip'
                print('Zip download URL: ', downloadURL)
                print('Beginning download of zipped output...')
                zip_response = session.get(downloadURL)
                # Raise bad request: Loop will stop for bad response code.
                zip_response.raise_for_status()
                with zipfile.ZipFile(io.BytesIO(zip_response.content)) as z:
                    z.extractall(path)
                if short_name == 'ATL03':
                    # try whether there were further zip files!
                    for zipnr in np.arange(2,page_size):
                       downloadURL = 'https://n5eil02u.ecs.nsidc.org/esir/' + orderID + '.zip?' + str(zipnr)
                       try:
                           print('Try - beginning download of zip download URL: ', downloadURL)
                           zip_response = session.get(downloadURL)
                            # Raise bad request: Loop will stop for bad response code.
                           zip_response.raise_for_status()
                           with zipfile.ZipFile(io.BytesIO(zip_response.content)) as z:
                               z.extractall(path)
                       except: 
                           print('failed, assuming zip download URL does not exist.')
                           break # breaks out of the for loop, ignoring further pages
                print('Data request', page_val, 'is complete.\n --------')
            else: print('Request failed.\n --------')    
        
    print('----------')
    print('DONE. Please double-check whether all files are downloaded and unzipped.')    
    print('Downloading several zip files per page currently only attempted for ATL03!')    
    
        
        
    # now move all the files from their subfolders to the main folder, and delete all the folders
    print('--------------> moving the files from their subfolders to the main folder...')
    folderlist=glob(foldername +'/*')
    for f in folderlist:
        filelist=glob(f+'/*')
        for f2 in filelist:
            f3=os.path.basename(f2)
            os.rename(f2,foldername + '/' + f3)
            # checking whether the folder is empty or not
            if len(os.listdir(f)) == 0:
                # removing the file using the os.remove() method
                os.rmdir(f)
            else:
                # messaging saying folder not empty
                print("Folder is not empty")
            
    # write the subset_request_params to a text file in the folder
    with open(foldername+ '/subset_request_params.txt', 'w') as file:
         file.write(json.dumps(subset_request_params)) # use `json.loads` to do the reverse
        