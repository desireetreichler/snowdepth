import h5py
import pandas as pd
import geopandas as gpd
import numpy as np

from shapely.geometry import LineString, MultiLineString



#read ICESat-2 [Quicklook] data into a GeoDataFrame.
def read_files(files_list, files_fld, bounds):

    # Load the ICESat-2 data and append them to an appender
    appender=[]

    # Loop trough all the files.

    for is2_file in files_list:


            # open the file
            with h5py.File(files_fld + is2_file, 'r') as f:
                

                # load into the appender latitude, longitude, best fit and ancillary data as date and beam name
                try:
                    n= f['gt2l/land_segments/longitude'][:].shape[0]
                    appender.append(pd.DataFrame(data={
                                                'date':np.repeat(np.datetime64(f['ancillary_data/data_start_utc'][0][:10]),n),
                                                 'lat': f['gt2l/land_segments/latitude'][:],
                                                 'lon': f['gt2l/land_segments/longitude'][:],
                                                 'h_te_best_fit': f['gt2l/land_segments/terrain/h_te_best_fit'][:],
                                                 'beam': np.repeat('2l',n)})
                                   )
                except:
                    pass
                
                try:
                    n= f[('gt2r/land_segments/longitude')][:].shape[0]
                    appender.append(pd.DataFrame(data={
                                                'date':np.repeat(np.datetime64(f['ancillary_data/data_start_utc'][0][:10]),n),
                                                 'lat': f['gt2r/land_segments/latitude'][:],
                                                 'lon': f['gt2r/land_segments/longitude'][:],
                                                 'h_te_best_fit': f['gt2r/land_segments/terrain/h_te_best_fit'][:],
                                                 'beam': np.repeat('2r',n)})
                                   )
                except:
                    pass
                
                try:
                    
                    n= f['gt1r/land_segments/longitude'][:].shape[0]        
                    appender.append(pd.DataFrame(data={
                                'date':np.repeat(np.datetime64(f['ancillary_data/data_start_utc'][0][:10]),n),
                                 'lat': f['gt1r/land_segments/latitude'][:],
                                 'lon': f['gt1r/land_segments/longitude'][:],
                                 'h_te_best_fit': f['gt1r/land_segments/terrain/h_te_best_fit'][:],
                                 'beam':np.repeat('1r',n)}))
                except:
                    pass
                
                try:

                    n= f['gt1l/land_segments/longitude'][:].shape[0]
                    appender.append(pd.DataFrame(data={
                                        'date':np.repeat(np.datetime64(f['ancillary_data/data_start_utc'][0][:10]),n),
                                         'lat': f['gt1l/land_segments/latitude'][:],
                                         'lon': f['gt1l/land_segments/longitude'][:],
                                         'h_te_best_fit': f['gt1l/land_segments/terrain/h_te_best_fit'][:],
                                         'beam':np.repeat('1l',n)})
                                    )
                except:
                    pass
                
                try:
                    
                    n= f['gt3l/land_segments/longitude'][:].shape[0]        
                    appender.append(pd.DataFrame(data={
                                'date':np.repeat(np.datetime64(f['ancillary_data/data_start_utc'][0][:10]),n),
                                 'lat': f['gt3l/land_segments/latitude'][:],
                                 'lon': f['gt3l/land_segments/longitude'][:],
                                 'h_te_best_fit': f['gt3l/land_segments/terrain/h_te_best_fit'][:],
                                 'beam':np.repeat('3l',n)})
                                   )
                except:
                    pass
                                
                try:
                    n= f['gt3r/land_segments/longitude'][:].shape[0]        
                    appender.append(pd.DataFrame(data={
                                'date':np.repeat(np.datetime64(f['ancillary_data/data_start_utc'][0][:10]),n),
                                 'lat': f['gt3r/land_segments/latitude'][:],
                                 'lon': f['gt3r/land_segments/longitude'][:],
                                 'h_te_best_fit': f['gt3r/land_segments/terrain/h_te_best_fit'][:],
                                 'beam':np.repeat('3r',n)})
                                   )
                except:
                    pass
                                
    #concatenate the appender
    is2 = pd.concat(appender)
    
    #create a GeoDataFrame from the normal DataFrame
    gf = gpd.GeoDataFrame(is2, geometry=gpd.points_from_xy(is2.lon, is2.lat), crs='epsg:7661')
    
    minx, miny, maxx, maxy = bounds
	
    gf=gf[((gf.lat > miny) & (gf.lat < maxy))]
    gf=gf[((gf.lon > minx) & (gf.lon < maxx))]

    return gf



# GeoDataFrame Points --> GeoDataFrame with MultiLineString

def points2multilines(gf, step = 10, datecol = 'date', beamcol = 'beam'):
	
    """Function to convert ICESat-2 ATL_08 data (or ATL_08 quicklook data)
    to a gdf with lines indicating the strips with data. 
    Usage: gf_lines = points2multilines(gf, datecol = 'date', beamcol = 'beam')
    Parameters: 
        gf      - input point geodataframe
		step    - distance (in km) between the points in the resulting multiline
        datecol - column name indicating the overpass date. Default: 'date'
        beamcol - coulmn name indicating the beam identifier (six per overpass). Default: 'beam'
        """
	
    # define the variable where data will be appended
    appender=[]



    # iterate through the dates and beams we have
    dates=gf[datecol].unique()
    beams=gf[beamcol].unique()


    for date in dates:
        for beam in beams:

            # select the relevant data 
            gf_sub = gf[gf[datecol]==date]
            gf_sub = gf_sub[gf_sub[beamcol]==beam]

            if gf_sub.shape[0] > 2:
                # sort the data from north to south
                try:
                    gf_sub.sort_values(by='lat',axis=0, inplace=True)
                except:
                    gf_sub.sort_values(by='latitude',axis=0, inplace=True)
                    
                # and reset the index, so that we can iterate following the north-south order
                gf_sub.index=np.array(range(gf_sub.shape[0]))

                # start_index is the first point of the "line"
                start_index= gf_sub.index[0]
                
                # define an empty variable were we can append lines
                multiline_i=[]

                # iterate trough the points
                for i in gf_sub.index[:-1]:
                    
                    dist=gf_sub.loc[i].geometry.distance(gf_sub.loc[i+1].geometry)
                    
                    dist_from_last=gf_sub.loc[i].geometry.distance(gf_sub.loc[start_index].geometry)

                    # check if the following point is farer than about 100 m 

                    if dist > 0.001 :

                        #if so, the line should end at THIS point
                        geom=LineString([gf_sub.loc[start_index].geometry,gf_sub.loc[i].geometry])

                        # the first point of the next line should be the following point
                        start_index=i+1

                        #add this line to the line list
                        multiline_i.append(geom)
                        
                    # if this point is further than the "step" distance from the last point in the multiline:
                    else:
                        if dist_from_last > (step /111):
                        
                            #if so, the line should end at THIS point
                            geom=LineString([gf_sub.loc[start_index].geometry,gf_sub.loc[i].geometry])
                            
                            # the first point of the next line should be the same point (no gap wanted)
                            start_index=i

                            #add this line to the line list
                            multiline_i.append(geom)



                
                dist=gf_sub.loc[i].geometry.distance(gf_sub.loc[i-1].geometry)
                
                if dist<0.001 and start_index != gf_sub.index[-1]:

                    geom=LineString([gf_sub.loc[start_index].geometry,gf_sub.loc[i-1].geometry])
                    multiline_i.append(geom)
                
                #create a MultiLine with the list we have
                geom=MultiLineString(multiline_i)        

                #create the geodataframe with this Multiline and the date and beam information
                gdf=gpd.GeoDataFrame(data={
                'date':[date],
                'beam':[beam],
                'geometry':[geom]})
                
                appender.append(gdf)
                   
    #concat the created geodataframes, and set the crs.
    gf_lines = pd.concat(appender) 
    gf_lines.crs='epsg:4326'    

    return gf_lines
