# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 10:24:21 2022

@author: desireet

This script contains functions to handle and predict actual ground tracks from RGT tracks. 
"""
import os
#import sys
from glob import glob
import re
import pickle
import geopandas as gpd
import pandas as pd
import fiona
import time
import numpy as np
import datetime
from pyproj import CRS
from shapely.geometry import Point, LineString, Polygon
from geopandas import GeoSeries
import matplotlib.pyplot as plt
from zipfile import ZipFile
try:
    import contextily as ctx
except: pass

try:  
    gpd.io.file.fiona.drvsupport.supported_drivers['KML']='rw' # native gpd has no functionality to read kml
    gpd.io.file.fiona.drvsupport.supported_drivers['KMZ']='rw' # native gpd has no functionality to read kml
except:
    from fiona.drvsupport import supported_drivers # fix from https://stackoverflow.com/questions/72960340/attributeerror-nonetype-object-has-no-attribute-drvsupport-when-using-fiona
    supported_drivers['LIBKML'] = 'rw'
try:
    from fiona.drvsupport import supported_drivers # fix from https://stackoverflow.com/questions/72960340/attributeerror-nonetype-object-has-no-attribute-drvsupport-when-using-fiona
    supported_drivers['LIBKML'] = 'rw'
except: 
    gpd.io.file.fiona.drvsupport.supported_drivers['KML']='rw' # native gpd has no functionality to read kml


# import custom modules
# if os.name=='posix':
#     mytoolspath='/uio/lagringshotell/geofag/icemass/icemass-users/desireet/science/ICESat-2/scripts/mytools/'
# if os.name=='nt':
#     mytoolspath = 'N:/science/ICESat-2/scripts/mytools/'
# sys.path.append(mytoolspath)
    



def matchRGTdata(clippedRGT, gdf_08, wd,cycle,ll=0, p=0 ):
    """ USAGE: GT_ortho, GT_xshift, rowswithcontent = matchRGTdata(clippedRGT, gdf_08, wd,cycle,ll=0, p=0)
    
        Matches RGT (clipped to AOI) and ATL08 data (from the same area)
        to check how far away the GT were from the RGT. 
        wd: In this directory, a folder called RGT will be made. The function 
        always plots at least one figure per RGT with the derived/predicted GT.
        c: cycle number
        ll: if set to 1, the analysis is done both in UTM (standard) and lat/lon.
        p=1: verbose output (some stats) and some more plotting.
        
        OUTPUT: GT_xxx: gdfs with the GT distances (and locations) for all six beams, 
        for mean / min / max distances each, for two methods: orthogonally to 
        the RGT (GT_ortho) and a parrallel shift in x (GT_xshift)
        rowswithcontent: a list of the row indices of the RGT where we found corresponding ATL08 data.
    """
    # RGT dates & numbers
    RGTdates = clippedRGT['track']
    RGTdates = [datetime.datetime.strptime(x[-15:-4], '%d-%b-%Y').date() for x in RGTdates]
    RGTdeltatime = [d-datetime.date(2018,1,1) for d in RGTdates] 
    RGTdeltatime = [d.days for d in RGTdeltatime] 
    RGTnrs = clippedRGT['RGT']
    # atl08 dates & RGT numbers
    deltatime08 = np.floor(gdf_08['delta_time'].to_numpy() /3600/24)
    deltatimes, datecounts=np.unique(deltatime08[:],return_counts=True)
    try:
        RGT08 = gdf_08['RGT'] 
        noRGTnr = False
    except: # this only works if the atl08 gdf has the RGT nr - not the case for the oldest downloaded tracks (older version of the scripts...)
        noRGTnr = True
            
    #dt =[datetime.date(2018,1,1)+datetime.timedelta(d) for d in deltatimes] 
    #date08 =[datetime.date(2018,1,1)+datetime.timedelta(d) for d in deltatime08] 
    
    # what is the current UTM zone?
    temp = clippedRGT.to_crs('EPSG:4326')
    crs = getUTMcrs(temp)
    
    rowswithcontent=[]
    gdf_ortho_list = []
    gdf_xshift_list = []
    for d in np.arange(0,len(RGTdeltatime)):
        #d = 0
        currdt=RGTdeltatime[d]
        currRGTnr = RGTnrs[d]
        incurrentdate = (deltatime08 == currdt)
        print(RGTdates[d], 'data points: ',sum(incurrentdate))
        #print(RGTdeltatime[d])
        if sum(incurrentdate)>0:
            # also the RGT number needs to match, as we could have two tracks from the same day
            if noRGTnr==False:
                iscurrentRGT= (RGT08==currRGTnr)
                incurrentdate = iscurrentRGT&incurrentdate                       
            gdf_curr = gdf_08[incurrentdate]
            # check that RGT and gdf are in same coordinate system
            gdfcrs=gdf_curr.crs.to_authority()
            rgtcrs=clippedRGT.crs.to_authority()
            if (gdfcrs[1]==rgtcrs[1]) == False:
                print('crs between RGT and data points do not agree! Converting to ',crs)
                gdf_curr=gdf_curr.to_crs(crs)
                clippedRGT = clippedRGT.to_crs(crs)
                
            xmin,ymin,xmax,ymax = gdf_curr.total_bounds #ds.bounds 
            #currRGTxy= clippedRGT.loc[[d],'geometry']
            currRGT= clippedRGT.loc[[d]]
            currRGT=currRGT.reset_index(drop=True)
            if (ymin<currRGT.total_bounds[3]) &  (ymax>currRGT.total_bounds[1]) & (xmax>currRGT.total_bounds[0]) & (xmin<currRGT.total_bounds[2]):
                print('match!')
                rowswithcontent.append(d)
                if p:
                    if (wd+'\\RGT' in glob(wd+'/*')) == False:
                        os.mkdir(wd+'/RGT') # make a RGT dir in the local directory
                    fig,ax = plt.subplots(figsize=(5,15))
                    #clim = np.nanpercentile(gdf_08['h_te_best_fit'].values,(2,98))
                    #fodar.to_crs(utm42n).boundary.plot(edgecolor="blue",linewidth=0.3, ax=ax)
                    #hefshp.boundary.plot(edgecolor="black",linewidth=0.5, ax=ax)
                    currRGT.plot(ax=ax, linewidth=3)
                    gdf_curr.plot('h_te_best_fit',ax=ax,s=0.1,legend=True,cmap='inferno',legend_kwds={'label':'Elevation [m]'})
                    try:
                        ctx.add_basemap(ax=ax, crs=gdf_curr.crs) # contextily map using our desired crs
                    except: 
                        pass
                    fig.suptitle(d)
                    ax.set_xlim((xmin, xmax))
                    ax.set_ylim((ymin, ymax))
                    fig.suptitle(clippedRGT.loc[d,'track'])#+' in FODAR areas (UTM42N)')
                    plt.savefig(wd+'/RGT/map_'+clippedRGT.loc[d,'track'][:-4]+'.png')
                
                # look at one RGT line point. How far away are the corresponding ATL points? both in x, lon and orthogonally?
                if ll: # latlon
                    GT_ortholl, GT_xshiftll = GTdiststats(currRGT.to_crs('EPSG:4326'), gdf_curr.to_crs('EPSG:4326'), cycle, wd+'/RGT/ll_', p)                        
                if 1: # UTM
                    GT_ortho, GT_xshift = GTdiststats(currRGT, gdf_curr, cycle, wd+'/RGT/', p)
                # store, unless the above returned none (may happen for cycle one)
                if (GT_ortho is None)==False:
                    gdf_ortho_list.append(GT_ortho)
                    gdf_xshift_list.append(GT_xshift)
    if len(gdf_ortho_list) >0:
        GT_ortho = concat_gdf(gdf_ortho_list)
        GT_xshift = concat_gdf(gdf_xshift_list)
    else:
        print('no data in cycle ',)
        return None, None, rowswithcontent
    return GT_ortho, GT_xshift, rowswithcontent

                
def GTdiststats(currRGT, gdf_curr, cycle, figprefix, p):
    """USAGE: GT_ortho, GT_xshift = GTdiststats(currRGT, gdf_curr, cycle, figprefix, p)
        Checks how far some actual data are away from the provided RGT.    
        currRGT and gdf_curr are a RGT line and the corresponding (acquired)
        ATL08 data in the SAME crs. figprefix is path where output figure with 
        the derived/predicted RGT is saved (for example the wd). p=1: verbose, 
        and some more plotting.
        
        OUTPUT: gdfs with the GT distances (and locations) for all six beams, 
        for mean / min / max distances each, for two methods: orthogonally to 
        the RGT (GT_ortho) and a parrallel shift in x (GT_xshift)"""
    RGTcoords = list(currRGT.geometry[0].coords)
    # make sure the RGT column contains numbers
    try:    
        currRGTnr = int(currRGT.loc[[0],'RGT'])
    except: 
        currRGTnr = currRGT.loc[[0],'RGT']
        try:
            currRGTnr = int(currRGTnr[0][4:])
        except: 
            currRGTnr = int(currRGTnr[0])
    # are we in UTM or lat lon?
    currcrs=currRGT.crs.to_authority()
    if currcrs[1] in ['4326']:
        dythreshold = 0.0006
    else: 
        dythreshold = 60
    # first, some statistics
    dxs = [] # in UTM x
    ds = []  # orthogonally
    
    # loop through all points in the RGT
    for r in np.arange(0,len(RGTcoords)):
        # take the second one, for example
        #r=1
        x = RGTcoords[r][0]
        y = RGTcoords[r][1]
        if r<(len(RGTcoords)-1):
            x1 = RGTcoords[r+1][0]
            y1 = RGTcoords[r+1][1]
            drx = x1-x
            dry = y1-y
            alpha = np.arctan(dry/drx)
            ascdesc=(y1-y)/abs(y1-y) # -1: desc, 1: asc
        
        # distance to left and right, in UTM
        for pair in [1,2,3]:              
            #if p: 
                #fig,ax = plt.subplots(figsize=(5,15))
            for b in [0,1]:
                gdf_pb=gdf_curr[(gdf_curr['pair']==pair) & (gdf_curr['beam']==b)]
                gdf_pb=gdf_pb.reset_index()
                absdy=list(abs(gdf_pb.geometry.y-y))
                #dy=list(gdf_pb.geometry.y-y)
                try:
                    dymin=min(absdy)
                except: # the list may be empty
                    continue
                #print(dymin)
                if dymin > dythreshold: # in lat: ca 0.0006, in m ca 60m
                    #print('dy too big, continue')
                    continue
                ind = absdy.index(dymin)
                # point with "same" x
                xmin = gdf_pb.geometry.x[ind]
                dx=xmin-x
                dxs.append([dx, pair, b, x+dx,y,x, y,currRGTnr])
                if 0: # debugging plot
                    #print('distance to the side, in m: ', round(dx))
                    fig,ax = plt.subplots(figsize=(15,5))
                    ax.plot(x,y,'o-',c='blue',linewidth=1, label='RGT')
                    ax.plot(x+dx,y,'*-',c='red',linewidth=1,label='predicted')
                    ax.plot(gdf_pb.geometry.x[ind-1:ind+2],gdf_pb.geometry.y[ind-1:ind+2], 'x',c='orange',markersize=10, label='ATL08')
                    ax.plot(gdf_pb.geometry.x[ind],gdf_pb.geometry.y[ind], '+',c='green',markersize=10, label='ATL08 match')
                    ax.legend()
                # orthogonal point
                if r<(len(RGTcoords)):
                    # find the two neighbouring points and interpolate Rx
                    try:
                        if absdy[ind-1] > absdy[ind+1]:
                            # if larger, take ind+1
                            inda = ind
                            indb = ind+1
                        else:
                            # if smaller, take ind-1 
                            inda = ind-1
                            indb = ind
                            if inda <0: raise Exception 
                    except:
                        #print('outside data range')
                        continue
                    A1y=gdf_pb.geometry.y[inda]
                    A2y=gdf_pb.geometry.y[indb]
                    A1x=gdf_pb.geometry.x[inda]
                    A2x=gdf_pb.geometry.x[indb]
                    # we want to find A0, which has coordinates A0x,Ry
                    dax0=(A2x-A1x)/(A2y-A1y)*(y-A1y)
                    A0x=A1x + dax0
                    dax = A0x - x  # should-be dist to Rx,y in x UTM
                    # find the orthogonal distance s
                    s = np.sin(alpha)*dax
                    # coordinates of the orthogonal point
                    Asx = np.sin(alpha)*s + x
                    Asy = -np.cos(alpha)*s + y
                    # asc/desc both work with the same code 
                    """if ascdesc ==-1: #descending case
                        s = np.sin(alpha)*dax
                        # coordinates of the orthogonal point
                        Asx = np.sin(alpha)*s + x
                        Asy = -np.cos(alpha)*s + y
                    else: # ascending
                        s = np.sin(alpha)*dax
                        # coordinates of the orthogonal point
                        Asx = np.sin(alpha)*s + x
                        Asy = -np.cos(alpha)*s + y"""
                    # store the point    
                    ds.append([s, pair, b, Asx, Asy, x, y, currRGTnr, ascdesc])
                    if 0: # debuggiing plot
                        print('distance to the side: ', dx)
                        fig,ax = plt.subplots(figsize=(15,5))
                        #ax.plot(RGTc[:,0],RGTc[:,1],'o-',c='blue',linewidth=1, label='RGT')
                        ax.plot(x,y,'o-',c='blue',linewidth=1,label='RGT pt')
                        ax.plot(gdf_pb.geometry.x[ind-5:ind+5],gdf_pb.geometry.y[ind-5:ind+5], 'x',c='orange',markersize=10, label='ATL08')
                        #ax.plot(A1x,A1y, '+',c='green',markersize=10, label='A1')
                        #ax.plot(A2x,A2y, '+',c='darkgreen',markersize=10, label='A2')
                        ax.plot([A1x,A2x],[A1y,A2y], '-+',c='darkgreen',markersize=10, label='A1-A2')
                        ax.plot(A0x,y, '*',c='green',markersize=10, label='A0 at Rx')
                        ax.plot(Asx,Asy, '*',c='red',markersize=10, label='predicted AS, orthogonal to RGT')
                        ax.legend()                    
                        
    # look at the resulting dx, in x (UTM) and store them
    dxs=np.array(dxs)
    ds=np.array(ds)
    if len(ds)==0: # this may happen especially for cycle 1
        print('no ATL08 data to match this RGT')
        return None, None
    
    for case in [1,2]:
        s = []
        smin = []
        smax= []
        if case ==1: # ortho
            dd=ds
            if p:
                print('pair/beam/ s orthogonal dist to RGT, mean, min, max')
        if case ==2: # x-shift
            dd=dxs
            if p:
                print('pair/beam/ x shift dist to RGT, mean, min, max')
        for pair in [1,2,3]:
            for b in [0,1]:
                I=(dd[:,1]==pair) & (dd[:,2]==b)
                if sum(I) ==0:
                    ss = np.nan
                    ssmin = np.nan
                    ssmax = np.nan
                else:
                    ss= np.nanmean(dd[I,0])
                    ssmin = np.nanmin(dd[I,0])
                    ssmax = np.nanmax(dd[I,0])
                if p:
                    print(pair, b, ss, ssmin, ssmax)
                # store the s orthogonal distance values of all RGTs    
                s.append(ss)
                smin.append(ssmin)
                smax.append(ssmax)  
        if case ==1:
            GT_ortho = predictGTortho(currRGT, s, smin, smax)
        if case ==2:
            GT_xshift = predictGTxshift(currRGT,s, smin, smax)
    
    # add the cycle to the RGT
    GT_ortho['cycle']=cycle
    GT_xshift['cycle']=cycle
    
    # plot the result
    fig,ax = plt.subplots(figsize=(10,25))
    ax.plot(gdf_curr.geometry.x,gdf_curr.geometry.y, 'x',c='orange',markersize=10, label='ATL08')
    # geopandas plots lines without markers
    GT_xshift.plot(color='darkred',marker = 'o',ax=ax,markersize=10, label='pred. x-dist')
    GT_ortho.plot(color='red',marker = '*',ax=ax,markersize=10, label='pred. s ortho')
    # add the points
    RGTline=np.array(RGTcoords[:])
    for pair in [1,2,3]:
        for b in [0,1]:
            I=(dxs[:,1]==pair) & (dxs[:,2]==b)                        
            ax.plot(dxs[I,3],dxs[I,4], 'o',c='darkred',markersize=10)#, label='pred. x-dist')
            I=(ds[:,1]==pair) & (ds[:,2]==b)                        
            ax.plot(ds[I,3],ds[I,4], '*',c='red',markersize=10)#, label='pred. s ortho')
    ax.plot(RGTline[:,0], RGTline[:,1],'-',c='blue',label='RGT')
    ax.plot(ds[I,5],ds[I,6], '.-',c='blue',markersize=10, label='RGT pts used')
    """ax.plot(np.array(Asx),np.array(Asy), '.',c='green',markersize=10, label='GT ortho pred')
    ax.plot(Asxmin,Asymin, '.',c='cyan',markersize=10, label='GT min')
    ax.plot(Asxmax,Asymax, '.',c='cyan',markersize=10, label='GT max')"""
    #ctx.add_basemap(ax=ax, crs=utm42n) # contextily map using our desired crs
    #fodar.to_crs(utm42n).boundary.plot(edgecolor="blue",linewidth=0.3, ax=ax)
    ax.legend() 
    fig.suptitle('RGT '+str(currRGTnr)+', cycle '+str(cycle))
    plt.savefig(figprefix+currRGT.loc[0,'track'][:-4]+'_predicted.png')
    
    return GT_ortho, GT_xshift
   


def predictGTortho(currRGT, s, smin=-1, smax=-1):
    """ USAGE: GT_gdf = predictGTortho(currRGT, s, smin = -1, smax = -1)
        Predicts GT from a given RGT  and the distances s (m or deg)
        ORTHOGONALLY to the RGT. If smin/smax are -1, they are not computed. If
        provided, they must have length 6, otherwise they are replaced 
        by standard values 100m / 0.001 deg in both directions.
        OUTPUT: gdf with estimated mean, min and max GT (3x6 lines)"""
    RGTcoords = list(currRGT.geometry[0].coords)
    if len(RGTcoords)<2:
        raise Exception('RGT has to contain at least 2 points (ca 7km length) to compute its angle.')
    RGTcoords=np.array(RGTcoords)
    RGTx=RGTcoords[:,0]
    RGTy=RGTcoords[:,1]
    try:    
        currRGTnr = int(currRGT.loc[[0],'RGT'])
    except: 
        currRGTnr = currRGT.loc[[0],'RGT']
        try:
            currRGTnr = int(currRGTnr[0][4:])
        except: 
            currRGTnr = int(currRGTnr[0])
    try:
        currdate = currRGT.loc[[0],'date'][0]
    except:
        currdate = np.nan
    try:
        currtrack = currRGT.loc[[0],'track'][0]
    except: currtrack = ''
    try:
        currcycle = currRGT.loc[[0],'cycle'][0]
    except: 
        currcycle = np.nan
            
    if len(s)<6:
        s=[3230, 3140, 45, -45, -3140, -3230]
    s= np.array(s)   
    sdict = {'mean':s}
    
    # are we in latlon?
    currcrs=currRGT.crs.to_authority()
    if currcrs[1] in ['4326']:
        sd = 0.001
    else: 
        sd = 100
        
    if (smin ==-1)==False:
        if (len(smin)<6):
            smin = s-sd
        sdict['min']=smin
    if (smax ==-1)==False:
        if (len(smax)<6):
            smax = s+sd
        sdict['max']=smax
        
    ascdesc=(RGTy[1]-RGTy[0])/abs(RGTy[1]-RGTy[0]) # -1: desc, 1: asc
    RGTdx=RGTx[1:]-RGTx[:-1]
    RGTdy=RGTy[1:]-RGTy[:-1]
    alpha = np.mean(np.arctan(RGTdy/RGTdx)) # varied at 3rd place after dec pt. e.g 1.5058, 1.5050
    
    # make a list of dictionaries with GT lines
    linelist = []
    #sdict = {'mean':s, 'min':smin, 'max':smax}
    # coordinates of the orthogonal point
    for k in np.arange(0,6):
        for tp, ss in sdict.items() :
            # s
            j=ss[k]
            Axj = np.sin(alpha)*j + RGTx
            Ayj = -np.cos(alpha)*j + RGTy
            points = GeoSeries(map(Point, zip(Axj, Ayj)))
            try:
                line = LineString(points.tolist())
                linelist.append({'pb':k, 's':j, 'RGT':currRGTnr, 'date':currdate,'cycle':currcycle,'track':currtrack,'stype':tp,'ascending':ascdesc,'predicted':'ortho', 'geometry':line })
            except:
                pass
    GT_gdf = gpd.GeoDataFrame(linelist)
    #GT_gdf.plot()
    # to return individual lines instead:
    #Asx = [np.sin(alpha)*j + RGTx for j in s]# with alpha variation: decimeter changes
    #Asy = [-np.cos(alpha)*j + RGTy for j in s]## with alpha variation: ca 2m differences
    return(GT_gdf)

def predictGTxshift(currRGT, s, smin=-1, smax=-1):
    """ USAGE: GT_gdf = predictGTxshift(currRGT, s, smin = -1, smax = -1)
        Predicts GT from a given RGT and the distances s - shifted 
        in x (or lon) direction. If smin/smax are -1, they are not computed. If
        provided, they must have length 6, otherwise they are replaced 
        by standard values 100m / 0.001 deg in both directions.
        OUTPUT: gdf with estimated mean, min and max GT (3x6 lines)"""
    RGTcoords = list(currRGT.geometry[0].coords)
    if len(RGTcoords)<2:
        raise Exception('RGT has to contain at least 2 points (ca 7km length) to compute its angle.')
    RGTcoords=np.array(RGTcoords)
    RGTx=RGTcoords[:,0]
    RGTy=RGTcoords[:,1]
    try:    
        currRGTnr = int(currRGT.loc[[0],'RGT'])
    except: 
        currRGTnr = currRGT.loc[[0],'RGT']
        try:
            currRGTnr = int(currRGTnr[0][4:])
        except: 
            currRGTnr = int(currRGTnr[0])
    try:
        currdate = currRGT.loc[[0],'date']
    except:
        currdate = np.nan
    try:
        currtrack = currRGT.loc[[0],'track']
    except: currtrack = ''
    try:
        currcycle = currRGT.loc[[0],'cycle']
    except: 
        currcycle = np.nan
        
    if len(s)<6:
        raise Exception(('distance s from RGT has the wrong size'))
    s= np.array(s)    
    sdict = {'mean':s}
    
    # are we in latlon?
    currcrs=currRGT.crs.to_authority()
    if currcrs[1] in ['4326']:
        sd = 0.001
    else: 
        sd = 100
    if (smin ==-1)==False:
        if (len(smin)<6):
            smin = s-sd
        sdict['min']=smin
    if (smax ==-1)==False:
        if (len(smax)<6):
            smax = s+sd
        sdict['max']=smax
    ascdesc=(RGTy[1]-RGTy[0])/abs(RGTy[1]-RGTy[0]) # -1: desc, 1: asc
    
    # make a list of dictionaries with GT lines
    linelist = []
    #sdict = {'mean':s, 'min':smin, 'max':smax}
    # coordinates of the tracks shifted parallel in x direction
    for k in np.arange(0,6):
        for tp, ss in sdict.items() :
            # s
            j=ss[k]
            Axj = j + RGTx
            Ayj = RGTy
            points = GeoSeries(map(Point, zip(Axj, Ayj)))
            line = LineString(points.tolist())
            linelist.append({'pb':k, 's':j, 'RGT':currRGTnr, 'date':currdate,'cycle':currcycle,'track':currtrack,'stype':tp,'ascending':ascdesc,'predicted':'xshift', 'geometry':line })
    GT_gdf = gpd.GeoDataFrame(linelist)
    
    return(GT_gdf)


def cleanRGT(clippedRGT):
    """ USAGE: clippedRGT=cleanRGT(clippedRGT)
        explodes multilines, remove horizontal ones, resets index, adds some info etc."""
    clippedRGT = clippedRGT.reset_index(drop=True)
    clippedRGT=clippedRGT.explode() #  there could be several multilinestrings
    # check whether we have RGT as a number
    try:
        if type(clippedRGT['RGT'][0][0])!=int:
            raise Exception 
    except: 
        # add RGT number
        try: # sometimes already renamed
            RGTnr = clippedRGT['Name']
        except:
            try:
                RGTnr = clippedRGT['RGT']
                try:
                    RGTnr=RGTnr.str.split()
                    RGTnr = [int(x[1]) for x in RGTnr]
                except:
                    pass
            except:
                RGTnr = clippedRGT['id']
    
        clippedRGT['RGT']= RGTnr
        
    try:
        clippedRGT.drop(['Name'], axis=1)
    except: 
        pass
    try:
        clippedRGT.drop(['id'], axis=1)
    except: 
        pass
    clippedRGT['linenr']=np.arange(len(clippedRGT))
    clippedRGT=clippedRGT.set_index('linenr')
    # clippedRGT=clippedRGT.set_index('RGT')
    #clippedRGT['RGT']= RGTnr
    
    # convert to cartesian 
    crs = getUTMcrs(clippedRGT)
    clippedRGT=clippedRGT.to_crs(crs)
    # remove the horizontal lines
    clippedRGT = dropfakeRGTs(clippedRGT,v=0)
    clippedRGT = clippedRGT.reset_index(drop=True)
    # if p:
    #     fig,ax = plt.subplots(figsize=(15,5))
    #     fodar.to_crs(currUTM).boundary.plot(edgecolor="blue",linewidth=0.3, ax=ax)
    #     #fodar.boundary.plot(edgecolor="blue",linewidth=0.3, ax=ax)
    #     clippedRGT.to_crs(currUTM).plot(ax=ax, linewidth=3)
    #     ctx.add_basemap(ax=ax, crs=currUTM) # contextily map using our desired crs
    return clippedRGT


def dropfakeRGTs(clippedRGT,v=1, polar = False):
    """ USAGE: clippedRGT = dropfakeRGTs(clippedRGT,v=1)
    drops RGTs with an ascending/descanding angle < 2 deg (< 0.01 deg in polar case)- these are artifacts 
    (where the RGT crossed 0 lon)
    v=1: verbose
    polar = False"""
    anglethres = 2
    if polar:
        anglethres = 0.01
    for index,row in clippedRGT.iterrows():
        #print(index,row)
        # linestrings:
        try:             
            coords = list(row.geometry.coords)
            angle = (coords[0][1]-coords[1][1])/(coords[0][0]-coords[1][0])
            # multilinestrings: try to avoid (by exploding the gdf before)
            # average all angles - assume there is no occasion with both valild and invalid (crossing horizontally) RGTs
        except:    
            angle=list()
            for g in range(len(row.geometry.geoms)):
                coords = list(row.geometry.geoms[g].coords)
                newangle=(coords[0][1]-coords[1][1])/(coords[0][0]-coords[1][0])
                angle.append(newangle)
            angle=np.mean(angle)
        if v:
            print(angle)
        
        if abs(angle)<anglethres:
            clippedRGT.drop(index, inplace=True)
            if v:
                print('-deleted')
        # add new index 
    clippedRGT = clippedRGT.reset_index(drop=True)
    return(clippedRGT)




def clipRGTcycle(RGTfolder, shp, cycle=-1):
    """ USAGE: clippedRGT = clipRGTcycle(RGTfolder, shp, cycle=-1)
    this function takes the folder in the RGTfolder with the highest cycle number (default),
    and clips all containing kml files with the provided shapefile (a gdf). 
    Alternatively, state another cycle.
    Takes about 10 minutes to run through the folder!!! Not very efficient...
    Returns: a geodataframe with clipped track lines, and a list with track names."""
    RGTcycles = sorted(glob(RGTfolder+'/IS2*'))
    # sort to find the highest cycle number
    howmany=len(RGTcycles)-9
    RGTcycles = RGTcycles[howmany:]+RGTcycles[:howmany]
    #print(RGTcycles)
    if type(cycle)==int:
        if cycle<=0:
            c=len(RGTcycles)-1
        else: 
            c=cycle-1
    else:   # if several are provided
        c=cycle-1
    selectedRGTs = [RGTcycles[c]]
    
    # convert shp to UTM and add a 5km buffer, reconvert to lat/lon
    crs = getUTMcrs(shp)
    buffer = shp.to_crs(crs).buffer(5000)
    bufferll = buffer.to_crs('EPSG:4326')
        
    clippedRGT=gpd.GeoDataFrame
    tracklist = []
    first = True
    t=time.time()
    for c in selectedRGTs: # take the last cycle 
        #c = selectedRGTs[0]
        RGTtracks = glob(c+'/*kml')
        #track = RGTtracks[0]
        
        for track in RGTtracks: 
            # track = RGTtracks[0]
            # I found out that geopandas only reads the first layer in the kml. To know the name of the layer we want, first, list the entire content
            #f=fiona.listlayers(RGTfolder+c+'/'+track)
            #IS2 = [s for s in f if 'IS2' in s]# the line is called IS2... find with pattern matching, or f[-1]), should always be the last (?))
            #kml = gpd.read_file(RGTfolder+c+'/'+track,driver='KML',layer=IS2[0]) # this works
            # try to speed up the above by guessing the layer name - it seems identical to the kml name:
            basename, trackname = os.path.split(track)
            try: kml = gpd.read_file(track,driver='KML',layer=trackname)   
                # success! this returns a linestring
            except: # the first few folders/kmls have a different structure, catch for these:
                f=fiona.listlayers(track)
                IS2 = [s for s in f if 'IS2' in s]
                kml = gpd.read_file(track,driver='KML',layer=IS2[0]) 
            currtrack = gpd.clip(kml, bufferll)
            if currtrack.empty==False:
                print(track)
                currtrack['track']= trackname
                tracklist.append(track)
                try: 
                    clippedRGT=clippedRGT.append(currtrack)
                except: # appending to the initialised variable somehow does not work, replace the first time a track is found (an dan error thrown)
                    print('could not append track')    
                    print(trackname)
                    if first==True: 
                        clippedRGT=currtrack
                        first = False
                    else: break
    print('elapsed time: %f' %(time.time()-t))   
    return clippedRGT


def clipRGTgdf(RGTfolder, shp, cycle=-1):
    """ USAGE: clippedRGT = clipRGTgdf(RGTfolder, shp, cycle=-1)
    this function takes the RGT gdf with the highest cycle number (default),
    and clips all containing RGTs with the provided shapefile (a gdf). 
    Alternatively, state another cycle.
    Takes about 10 minutes to run through the folder!!! Not very efficient...
    Returns: a geodataframe with clipped track lines, and a list with track names."""
    t=time.time()
    RGTcycles = sorted(glob(RGTfolder+'/*.json'))
    # sort to find the highest cycle number
    howmany=len(RGTcycles)-8
    RGTcycles = [RGTcycles[0]]+RGTcycles[howmany:]+RGTcycles[1:howmany]
    #print(RGTcycles)
    if type(cycle)==int:
        if cycle<=0:
            c=len(RGTcycles)-1
        else: 
            c=cycle-1
    else:   # if several are provided
        c=cycle-1
    
    selectedRGTs = RGTcycles[c]
    basen, filen = os.path.split(selectedRGTs)
    try: # picle is faster but fails in some cases
        with open(basen+'/'+filen[:-4]+'pkl', 'r') as f: RGT_gdf = pickle.load(f)
    except: # read the geojson object instead
        with open(basen+'/'+filen, 'r') as f: RGT_gdf = gpd.read_file(f)
    #    with open('N:/science/ICESat-2/tracks/RGT/RGTc1.pkl', 'rb') as f: RGT_gdf=pickle.load(f)
    
    # check whether the dataframe has the correct cycle numbers
    if RGT_gdf['cycle'][0] != c+1:
        print('The requested cycle number does not match the gdf content.')
        return
    
    # convert shp to UTM and add a 5km buffer, reconvert to lat/lon
    crs = getUTMcrs(shp)
    buffer = shp.to_crs(crs).buffer(5000)
    bufferll = buffer.to_crs('EPSG:4326')
    RGT_gdf=RGT_gdf.to_crs(bufferll.crs)
        
    clippedRGT = gpd.clip(RGT_gdf,bufferll)
        
    print('elapsed time: %f' %(time.time()-t))   
    return clippedRGT



def loadRGTtimestamps(clippedRGT, RGTfolder, cycle,shp, v=1):
    """ USAGE: loadRGTtimestamps(clippedRGT, RGTfolder, cycle,shp, v=0)
    adds a timestamp string to each RGT of the specified cycle number. The function
    takes the average/middle timestamps from all timestamp points within the 
    input shapefile shp, a gdf with a square bounding box.
    The input RGT gdf (previously clipped with shp) needs to have an integer 
    column named "RGT" with the RGT numbers (run e.g. cleanRGT to achieve this).
    Note that each RGT can only be once in the RGT_gdf!! (i.e. no strange clipping shapes)
    v=1: verbose, prints all loaded tracks etc.
    Very slow/inefficient function, as it loads all timestamps from the kml..."""
    # find the RGT nrs of all kml files of the cycle in question
    #cleanedRGT=cleanRGT(clippedRGT)
    ti=time.time()
    RGTcycles = sorted(glob(RGTfolder+'/IS2*'))
    cyclesplit = [re.split('cycle', x, maxsplit=0, flags=0) for x in RGTcycles ]
    cyclenrs = [re.split('_date', x[1], maxsplit=0, flags=0) for x in cyclesplit ]
    cyclenrs = [int(x[0] ) for x in cyclenrs ]
    cyclematches = [x ==cycle for x in cyclenrs]
    RGTcycle = [i for (i,v) in zip(RGTcycles,  cyclematches) if v] 
            
    # find the kml tracks
    RGTtracks = glob(RGTcycle[0]+'/*')
    RGTsplit = [re.split('RGT_', x, maxsplit=0, flags=0) for x in RGTtracks ]
    RGTnrs = [int(x[1][:4]) for x in RGTsplit ]
    # match
    RGTmatches = [(x in clippedRGT.RGT.values) for x in RGTnrs]
    RGTtrackstoload = [i for (i,v) in zip(RGTtracks,  RGTmatches) if v] 
    
    #timestamp_list=pd.DataFrame()
    
    # convert shp to UTM and add a 10km buffer, reconvert to lat/lon
    crs = getUTMcrs(shp)
    buffer = shp.to_crs(crs).buffer(510000)
    bufferll = buffer.to_crs('EPSG:4326')
            
    for k in np.arange(len(RGTtrackstoload)): 
        track = RGTtrackstoload[k]
        timestamps_gdf=gpd.GeoDataFrame()
        first = True
        t=time.time()       
        basename, trackname = os.path.split(track)
        if v:
            print('Finding timestamp of track '+trackname)
        f=fiona.listlayers(track)
        for t in f:
            if t[:3]=='RGT':
                kml = gpd.read_file(track,driver='KML',layer=t) 

                # append
                try: 
                    timestamps_gdf=pd.concat([timestamps_gdf,kml])
                except: # appending to the initialised variable somehow does not work, replace the first time a track is found (an dan error thrown)
                    if v: 
                        print('could not append timestamp')    
                        print(trackname)
                    if first==True: 
                        timestamps_gdf=kml
                        first = False
                    else: break
        
        currrgtvalue = int(t.split(' ')[1])
        # clip
        timestamps_gdf=timestamps_gdf.to_crs(bufferll.crs)
        timestamps_gdf = timestamps_gdf.reset_index(drop=True)
        clippedtimestamps= gpd.clip(timestamps_gdf,bufferll)
        
        # find timestamp
        timesplit = [re.split(' DOY', x, maxsplit=0, flags=0) for x in clippedtimestamps.Description.values ]
        timesplit = [x[0][-20:] for x in timesplit]
        # take the middle value
        timestamp = timesplit[int(np.floor(len(timesplit)/2))]
        
        #add to the clipped gdf with same RGT:
        rgtI = clippedRGT.RGT==currrgtvalue
        clippedRGT.loc[rgtI,'timestamp']=timestamp

    if v:
        print('elapsed time: %f' %(time.time()-ti))   
    return clippedRGT







def clipGTgdf(GTfolder, shp):
    """ USAGE: clippedRGT = clipRGTgdf(GTfolder, shp)
    this function takes the GT gdf,
    and clips all containing GTs with the provided shapefile (a gdf). 
    Returns: a geodataframe with clipped track lines, and a list with track names."""
    t=time.time()
    selectedGTs = sorted(glob(GTfolder+'/*.json'))
        
    # check whether the there are several json files (there should be only one)
    if len(selectedGTs) != 1:
        print('There should be (only) one .json file in this folder.')
        return
    
    basen, filen = os.path.split(selectedGTs[0])
    try: # picle is faster but fails in some cases
        with open(basen+'/'+filen[:-4]+'pkl', 'r') as f: GT_gdf = pickle.load(f)
    except: # read the geojson object instead
        with open(basen+'/'+filen, 'r') as f: GT_gdf = gpd.read_file(f)
    
    # convert shp to UTM and add a 5km buffer, reconvert to lat/lon
    crs = getUTMcrs(shp)
    buffer = shp.to_crs(crs).buffer(5000)
    bufferll = buffer.to_crs('EPSG:4326')
    GT_gdf=GT_gdf.to_crs(bufferll.crs)
        
    clippedGT = gpd.clip(GT_gdf,bufferll)
    
    # remove the horizontal GTs (artefacts of last to first kml point)
    
    print('elapsed time: %f' %(time.time()-t))   
    return clippedGT



def getUTMcrs(geodataset):
    """ USAGE: crs = getUTMcrs(geodataset) - crs in correct UTM zone for data in lat/lon"""
    utm = np.round((183 + np.mean([geodataset.total_bounds[0],geodataset.total_bounds[2]]))/6)
    if utm > 60:
        temp = geodataset.to_crs('EPSG:4326')
        utm = np.round((183 + np.mean([temp.total_bounds[0],temp.total_bounds[2]]))/6)
        #raise Exception(('UTM zone cannot be computed. Already in UTM..?'))
    if geodataset.total_bounds[1]<0: # south
       crs = CRS.from_string('+proj=utm +zone='+str(int(utm))+'+south')
    else:
       crs = CRS.from_string('+proj=utm +zone='+str(int(utm))) 
    return crs

def getUTMcrslatlon(lat,lon):
    """ USAGE: crs = getUTMcrslatlon(lat,lon) - crs in correct UTM zone for data in lat/lon"""
    utm = np.round((183 + lon)/6)
    if utm > 60:
        print('Coordinates need to be in lat/lon.')
        raise Exception(('UTM zone cannot be computed. Already in UTM..?'))
    if lat<0: # south
       crs = CRS.from_string('+proj=utm +zone='+str(int(utm))+'+south')
    else:
       crs = CRS.from_string('+proj=utm +zone='+str(int(utm))) 
    return crs

def sstats(GT_all, slist=['mean','min','max']):
    #GT_all=GT_ortho_all
    for ss in slist:
        try:
            GT_type = GT_all[GT_all['stype'] == ss]
        except:
            GT_type = GT_all[GT_all['type'] == ss]

        if len(GT_type)<1:
                continue
        print('-------'+ss+' - median, mean, min, max:')
        for ids in [0,1,2,3,4,5]:
            try:
                GT=GT_type[GT_type['pb']==ids]
            except:
                GT=GT_type[GT_type['id']==ids]

            if len(GT)<1:
                continue
            try:
                med = round(np.nanmedian(GT.s))
                m = round(np.nanmean(GT.s))
                mi = round(np.nanmin(GT.s))
                ma = round(np.nanmax(GT.s))
                print(ids, med,m, mi, ma)
            except:
                print(ids, '-','-','-','-')
             
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


""" convert rgts to gdf """

# mid-latitudes:
def convertRGTtogdf(cycles, RGTfolder='lagringshotell.uio.no/geofag/projects/snowdepth/ICESat-2/tracks/RGT/'):
    """ USAGE: convertRGTtogdf(cycles, RGTfolder='lagringshotell.uio.no/geofag/projects/snowdepth/ICESat-2/tracks/RGT/')
        Reads kml files and converts all tracks of one cycle to one geodataframe.
        The geodataframe is saved as geojson and pkl files. 
    Input:    
    cycles:     list of [cycle(s)] you want to convert from kml to gdf
    RGTfolder:  base folder where the RGT kmls are stored, 
                default 'N:/science/ICESat-2/tracks/RGT/' has to be adapted!
    Output:     printed file path -> to read: 
                with open(RGTfolder+'RGTc'+str(cycle)+'.pkl', 'rb') as f: RGT_gdf=pickle.load(f)"""
    
    if RGTfolder[-1]!='/':
        RGTfolder = RGTfolder + '/'
    
    # find how many cycle folders we have
    #RGTfolder = 'N:/science/ICESat-2/tracks/RGT/'
    RGTcycles = glob(RGTfolder+'/IS2*')
    # sort
    howmany=len(RGTcycles)-9
    RGTcycles = RGTcycles[howmany:]+RGTcycles[:howmany]
    
    # the loop below fails if we only have one int, put into a list:
    if type(cycles)==int:
        cycles=[cycles]
        
    # loop through cycles
    for cycle in cycles:#np.arange(0,len(RGTcycles)):
        print('convert cycle '+str(cycle))
        #c = RGTcycles[c]
        c = RGTcycles[cycle-1]
        # list the tracks
        RGTtracks = glob(c+'/*kml')
        RGT_gdf = loadRGTs(RGTtracks,v=1)
        #tracklist=RGT_gdf[1]
        #RGT_gdf=RGT_gdf[0]
        # clean up a bit
        
        RGT_gdf['cycle']=cycle
        RGT_gdf = RGT_gdf.reset_index(drop=True)
        
        # save
        with open(RGTfolder+'RGTc'+str(cycle)+'.json', 'w') as f: f.write(RGT_gdf.to_json())
        with open(RGTfolder+'RGTc'+str(cycle)+'.pkl', 'wb') as f: pickle.dump(RGT_gdf,f) 
        print('saved: '+RGTfolder+'RGTc'+str(cycle)+'.pkl/json')
    # to read them: 
    #for cycle in np.arange(0,len(RGTcycles)):
    #    with open(RGTfolder+'RGTc'+str(cycle)+'.pkl', 'rb') as f: RGT_gdf=pickle.load(f)
    
    


def convertkmlstogdf(kmlfolder): # polar
    """ USAGE: convertkmlstogdf(kmlfolder)
        Reads kml files in the folder and converts all tracks to one geodataframe.
        The geodataframe is saved as geojson and pkl files. 
    Input:    
    kmlfolder:  base folder where the kmls are stored, 
    Output:     printed file path -> to read: 
                read: with open(kmlfolder+GTname+'.pkl', 'rb') as f: GT_gdf=pickle.load(f) """
    
    if kmlfolder[-1]!='/':
        kmlfolder = kmlfolder + '/'
    
    # find how many cycle folders we have
    #RGTfolder = 'N:/science/ICESat-2/tracks/RGT/'
    GTs = glob(kmlfolder+'*kmz')
    
    # initialise geodataframe   
    GT_gdf=gpd.GeoDataFrame()
    #tracklist = []
    first = True
    t=time.time()    
    
    # loop through GTs
    for GT in GTs:#np.arange(0,len(RGTcycles)):
        print('convert GT '+str(GT))
        # unzip kmz
        kmz = ZipFile(GT,'r')
        kmls = kmz.namelist()
        # there should only be one
        kml = kmz.open(kmls[0])
        # list the tracks
        curr_gdf = gpd.read_file(kml,driver='KML')  
            # clean up a bit
        curr_gdf['RGT']=curr_gdf['Name']
        curr_gdf['GT']=GT[-8:-4]
        curr_gdf = curr_gdf.reset_index(drop=True)
        
        # append
        try: 
            GT_gdf=GT_gdf.append(curr_gdf)
        except: # appending to the initialised variable somehow does not work, replace the first time a track is found (an dan error thrown)
            print('could not append track')    
            print(GT)
            if first==True: 
                GT_gdf=curr_gdf
                first = False
            else: break
        print('elapsed time: %f' %(time.time()-t))   
        
    # save
    basename, GTname = os.path.split(kmlfolder[:-1])
    with open(kmlfolder+GTname+'.json', 'w') as f: f.write(GT_gdf.to_json())
    with open(kmlfolder+GTname+'.pkl', 'wb') as f: pickle.dump(GT_gdf,f) 
    print('saved: '+kmlfolder+GTname+'.pkl/json')
    # to read them: 
    #    with open(kmlfolder+GTname+'.pkl', 'rb') as f: GT_gdf=pickle.load(f)



def loadRGTs(RGTtracks,v=0):
    """ USAGE: RGT_gdf= loadRGTs(RGItracks,v=1)
    stores all kml files in RGItracks (list of paths) in a geodataframe.
    v=1: verbose, prints all loaded tracks etc.
    Takes about 10 minutes to run through one RGT folder!!! Not very efficient...
    Returns: a geodataframe with RGT track lines."""
    #print(RGTcycles)
    # if type(cycle)==int:
    #     if cycle<=0:
    #         c=len(RGTcycles)-1
    #     else: 
    #         c=cycle-1
    # else:   # if several are provided
    #     c=cycle-1
    # selectedRGTs = [RGTcycles[c]]
    
    # convert shp to UTM and add a 5km buffer, reconvert to lat/lon
    # crs = getUTMcrs(shp)
    # buffer = shp.to_crs(crs).buffer(5000)
    # bufferll = buffer.to_crs('EPSG:4326')
        
    RGT_gdf=gpd.GeoDataFrame()
    #tracklist = []
    first = True
    t=time.time()
    #for c in selectedRGTs: # take the last cycle 
        #c = selectedRGTs[0]
        #RGTtracks = glob(RGTfolder+c+'/*kml')
        #track = RGTtracks[0]
        
    for track in RGTtracks: 
        # track = RGTtracks[0]
        # I found out that geopandas only reads the first layer in the kml. To know the name of the layer we want, first, list the entire content
        #f=fiona.listlayers(RGTfolder+c+'/'+track)
        #IS2 = [s for s in f if 'IS2' in s]# the line is called IS2... find with pattern matching, or f[-1]), should always be the last (?))
        #kml = gpd.read_file(RGTfolder+c+'/'+track,driver='KML',layer=IS2[0]) # this works
        # try to speed up the above by guessing the layer name - it seems identical to the kml name:
        basename, trackname = os.path.split(track)
        try: kml = gpd.read_file(track,driver='KML',layer=trackname)   
            # success! this returns a linestring
        except: # the first few folders/kmls have a different structure, catch for these:
            f=fiona.listlayers(track)
            IS2 = [s for s in f if 'IS2' in s]
            kml = gpd.read_file(track,driver='KML',layer=IS2[0]) 
            
        #currtrack = gpd.clip(kml, bufferll)
        if kml.empty==False:
            if v:
                print(track)
            kml['track']= trackname
            kml=kml.rename(columns={'Name':'RGT'})
            kml=kml.drop(['Description'],axis=1)
            #geojson can't store datetime
            #currdate = datetime.datetime.strptime(trackname[-15:-4],'%d-%b-%Y')
            currdate = trackname[-15:-4]
            kml['date']=currdate
            #tracklist.append(track)
            try: 
                RGT_gdf = pd.concat([RGT_gdf, kml])
                #RGT_gdf=RGT_gdf.append(kml) # this caused a future warning
            except: # appending to the initialised variable somehow does not work, replace the first time a track is found (an dan error thrown)
                if v: 
                    print('could not append track')    
                    print(trackname)
                if first==True: 
                    RGT_gdf=kml
                    first = False
                else: break
    if v:
        print('elapsed time: %f' %(time.time()-t))   
    return RGT_gdf#   , tracklist







    