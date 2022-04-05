import numpy as np
import h5py
from datetime import datetime
from pyproj import Proj   # assumes we only have one datum (WGMS84) and want UTM33N


def ATL08_to_dict(filename, dataset_dict, ancillary_dict, *args):
    """
        Read selected datasets from an ATL08 file

        Input arguments:
            filename: ATl08 file to read
            dataset_dict: For all six profiles, a dictinary describing the fields to be read
                    keys give the group names to be read, 
                    entries are lists of datasets within the groups
            ancillary_dict: A dictionary describing the ancillary data to be read (i.e. for the entire granule)
            'v':    verbose mode
        Output argument:
            D6: dictionary containing ATL06 data.  Each dataset in 
                dataset_dict has its own entry in D6.  Each dataset 
                in D6 contains a list of numpy arrays containing the 
                data, plus the data from the ancillary_dict
                
        Adapted from corresponding ATL06 hackathon week script..        
    """
    otherargs= args
    testword="v"
    if testword in otherargs:
        print('verbose mode on')
        verbose = 1
    else: verbose = 0
    
    D6=[]
    pairs=[1, 2, 3]
    beams=['l','r']
    # open the HDF5 file
    with h5py.File(filename,'r') as h5f:   # open hp5 file in read mode
        
        # debug: look at the file
        # inspect HDF5 from python
        # print(h5f.keys())
    # view full structure
#        def print_attrs(name, obj):
#            print(name)
#            for key,val in obj.attrs.items():
#                print("    %s: %s" % (key, val))
#        h5f.visititems(print_attrs)
        
        
        
        # loop over the groups in the ancillary dictionary
        #print('ANCILLARY DATA')
        anc={}
        for group in ancillary_dict.keys():
            if verbose: print('reading /%s' % (group))
            for dataset in ancillary_dict[group]:
                DS='/%s/%s' % (group,dataset)
                #print(DS)
                anc = read_group(anc, h5f, DS, dataset)
        #print(anc)
        # loop over beam pairs
        for pair in pairs:
            # loop over beams
            for beam_ind, beam in enumerate(beams):
                # check if a beam exists, if not, skip it
                if '/gt%d%s/land_segments' % (pair, beam) not in h5f:
                    if verbose: print('/gt%d%s/land_segments missing' % (pair, beam))
                    continue
                
                if verbose: print('reading /gt%d%s/land_segments' % (pair, beam))
                # loop over the groups in the dataset dictionary
                temp={}
                for group in dataset_dict.keys():
                    for dataset in dataset_dict[group]:
                        DS='/gt%d%s/%s/%s' % (pair, beam, group, dataset)
                        temp = read_group(temp, h5f, DS, dataset)
                if len(temp) > 0:
                    # it's sometimes convenient to have the beam and the pair as part of the output data structure: This is how we put them there.
                    temp['pair']=np.zeros_like(temp['latitude'])+pair
                    temp['beam']=np.zeros_like(temp['latitude'])+beam_ind
                    temp['filename']=filename
                    for anci in anc.keys():
                        #print(anci)
                        temp[anci]=anc[anci]
                    #convert lat/lon to UTM 33N
                    utmx, utmy = to_utm(temp)
                    temp['utmx']=utmx
                    temp['utmy']=utmy

                    D6.append(temp)
    return D6

#    rdata = np.empty((DIM0,), dtype=npdtype)
#    dset.read(h5py.h5s.ALL, h5py.h5s.ALL, rdata)
#
#    # Output the data to the screen.
#    print("%s:" % DATASET)
#    print(rdata)


def read_group(temp, h5f, DS, dataset):
     # since a dataset may not exist in a file, we're going to try to read it, and if it doesn't work, we'll move on to the next:                        
    try:
        # check for non-numeric content
        try: STR = [x.decode() for x in h5f[DS]]
        except: STR= [1]; 
        # handle non-numeric content
        if  isinstance(STR[0], (str)):
            #print('found a string')
            #print(STR[0])
            # utc = datetime.utcnow()
            DS = datetime.strptime(STR[0],'%Y-%m-%dT%H:%M:%S.%fZ') # e.g. '2011-01-21 02:37:21'
            temp[dataset]=DS
        else: 
            #print('data found')
            #print(h5f[DS])
            temp[dataset]=np.array(h5f[DS])
            # some parameters have a _FillValue attribute.  If it exists, use it to identify bad values, and set them to np.NaN
            if '_FillValue' in h5f[DS].attrs:
                #print('fillvalue')
                fill_value=h5f[DS].attrs['_FillValue']
                #print('fillvalue found')
                #try to fill in fill value
                temp[dataset][temp[dataset]==fill_value]=np.NaN
    except KeyError as e:
        print('key error:')
        print(DS)
        #print(fill_value)
        raise KeyError(e) #was: pass
    except ValueError as e:
        #print(str(e))
        if str(e) == 'cannot convert float NaN to integer':
            pass # seems the cover masks are integer. Then we can't replace the fill value with a float NaN. Thus, leave as is (original fill value)
        else: 
            print('value error:')
            print(DS)
            print(fill_value)
            raise ValueError(e)   
    return temp

def to_utm(temp):
    # convert to x/y utm coordinates 
    lon=temp['longitude']
    lat =temp['latitude']
    
    myProj = Proj("+proj=utm +zone=33N, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs") ## assumes we only have one datum (WGMS84, otherwise use transform) and want UTM33N
    utmx, utmy = myProj(lon, lat) # to convert back, add argument: , inverse=True
    
    return utmx, utmy
