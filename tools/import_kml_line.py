from osgeo import ogr
from shapely.geometry import LineString


def import_kml_line(kml_path):
	# imports a kml line file (ICESat-2 reference ground track) that erroneously resulted in a single point with gpd.read_file(kml_path, driver='KML')	
	# returns a shapely linestring
	
	ds = ogr.Open(kml_path)
	
	track=[] #=LineString([])
	for lyr in ds:
		for feat in lyr:
			geom = feat.GetGeometryRef()
			if geom != None:
				for i in range(0, geom.GetPointCount()):
					t=geom.GetPoint(i)
					track.append(t)
	t=track
	track=LineString(t)		
	return track	
######