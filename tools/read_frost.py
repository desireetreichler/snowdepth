import requests, os, glob
import pandas as pd


def read_frost(sources, elements, start, end):
    """ Reads data from Met.no's frost.met.no.
    Based on example from https://frost.met.no/python_example.html.

    Parameters:
    sources: list of station name (list)
    elements: list of element (list)
    start: starttime in formay yyy-mm-dd (str)
    end: endtime in formay yyy-mm-dd (str)

    Returns:
    Dataframe with requested data
    """
    client_id = os.getenv('FROST_API_CLIENTID')

    sep = ','
    time_sep = '/'
    sources = sep.join(sources)
    elements = sep.join(elements)
    referencetime = time_sep.join([start,end])

    endpoint = 'https://frost.met.no/observations/v0.jsonld'
    parameters = {'sources': sources,
                'elements': elements,
                'referencetime': referencetime}

    # Issue an HTTP GET request
    r = requests.get(endpoint, parameters, auth=(client_id,''))
    # Extract JSON data
    json = r.json()

    # Check if the request worked, print out any errors
    if r.status_code == 200:
        data = json['data']
        qt = json['queryTime']
        n = json['totalItemCount']
        print(f"Retrieved {n} items from frost.met.no with query time {qt} s")
    else:
        print('Error! Returned status code %s' % r.status_code)
        print('Message: %s' % json['error']['message'])
        print('Reason: %s' % json['error']['reason'])
        return r.status_code

    print("Converting to df")
    df = pd.json_normalize(json['data'], record_path = 'observations', meta = 'referenceTime')
    df['referenceTime'] = pd.to_datetime(df['referenceTime'])
    df.set_index('referenceTime',drop=True, inplace=True)
    df.sort_index(inplace=True)
    return df

def read_climate_normals(elements,sources,periods):
    """ Reads data from Met.no's frost.met.no.
    Parameters:
    elements: list of element (list)

    Returns:
    Dataframe with requested data
    """
    client_id = os.getenv('FROST_API_CLIENTID')

    sep = ','
    time_sep = '/'
    sources = sep.join(sources)
    elements = sep.join(elements)
    periods = time_sep.join(periods)

    endpoint = 'https://frost.met.no/climatenormals/v0.jsonld'
    parameters = {
                'period': periods,
                'elements': elements,
                'sources': sources
                }

    # Issue an HTTP GET request
    r = requests.get(endpoint, parameters, auth=(client_id,''))
    # Extract JSON data
    json = r.json()

    # Check if the request worked, print out any errors
    if r.status_code == 200:
        data = json['data']
        qt = json['queryTime']
        n = json['totalItemCount']
        print(f"Retrieved {n} items from frost.met.no with query time {qt} s")
    else:
        print('Error! Returned status code %s' % r.status_code)
        print('Message: %s' % json['error']['message'])
        print('Reason: %s' % json['error']['reason'])
        return r.status_code

    print("Converting to df")
    df = pd.json_normalize(json['data'])
    return df


def frost_info(sources, start, end, outdir):
    """
    Script to find out which variable are available for a given station between
    given start and end date

    Parameters:
    sources: list of station number (list)
    elements: list of variables (list)
    start: starttime in formay yyy-mm-dd (str)
    end: endtime in formay yyy-mm-dd (str)
    outdir: path to directory where to store the output data
    """
    client_id = os.getenv('FROST_API_CLIENTID')
    time_sep = '/'
    sep = ','

    if (type(sources)==list) & (len(sources)>1):
        sources = sep.join(sources)
    referencetime = time_sep.join([start,end])

    endpoint = 'https://frost.met.no/observations/availableTimeSeries/v0.jsonld'
    parameters = {
        'sources': sources,
        'referencetime': referencetime,
        }

    # Issue an HTTP GET request
    r = requests.get(endpoint, parameters, auth=(client_id,''))
    json = r.json()
    # data_list = []
    # [data_list.append([d['elementId'],d['validFrom']]) for d in json]
    # data = pd.DataFrame(data_list, names=['elementID','available_from'])
    df = pd.json_normalize(json['data'])
    #outdir = os.path.join('..','data','info')
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for id, data in df.groupby('sourceId'):
        source = id.split(':')[0]
        outfile = f"{source}_{start[:4]}_{end[:4]}.csv"
        outfile = os.path.join(outdir,outfile)
        for col in ['timeSeriesId',
                'exposureCategory','status',
                'uri','codeTable']:
            try: data.drop(col,axis=1,inplace=True)
            except: pass
        data.to_csv(outfile)
    return data

def climate_normals_info(sources):
    endpoint = 'https://frost.met.no/climatenormals/available/v0.jsonld'
    client_id = os.getenv('FROST_API_CLIENTID')
    sep = ','
    sources = sep.join(sources)
    parameters = {
        'sources': sources,
        }
    r = requests.get(endpoint, parameters, auth=(client_id,''))
    json = r.json()
    df = pd.json_normalize(json['data'])
    outdir = os.path.join('..','data','info')
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for id, data in df.groupby('sourceId'):
        source = id.split(':')[0]
        outfile = f"{source}_climate_normals.csv"
        outfile = os.path.join(outdir,outfile)
        data.to_csv(outfile)
    return df