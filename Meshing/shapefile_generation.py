
import os
import shapely
from shapely.ops import unary_union
import fiona.crs
import pyproj

UTM_ZONE30 = pyproj.Proj(proj='utm', zone=30, datum='WGS84', units='m', errcheck=True)
schema = {'geometry': 'LineString', 'properties': {'PhysID': 'int'}}
schema_polygon = {'geometry': 'Polygon', 'properties': {'PhysID': 'int'}}
crs = fiona.crs.from_string(UTM_ZONE30.srs)


length = 2000.
width = 1000.

array_length = 500
array_width = 400

circ_radius = 10.


x0, y0, x1, y1 = 0, 0, length, width
x_a0, y_a0, x_a1, y_a1 = 900., 0.5*width - circ_radius, 1500., 0.5*width + circ_radius
circ_x1, circ_x2, circ_y = 900., 1500., 0.5*width

circ_1 = shapely.geometry.Point(circ_x1, circ_y).buffer(circ_radius)
circ_2 = shapely.geometry.Point(circ_x2, circ_y).buffer(circ_radius)
rectangle = shapely.geometry.Polygon([(x_a0, y_a0), (x_a1, y_a0), (x_a1, y_a1), (x_a0, y_a1)])
merged_shape = unary_union([rectangle, circ_1, circ_2])

# Linestring
features = \
    [shapely.geometry.LineString([(x0, y0), (x1, y0)]),
     shapely.geometry.LineString([(x1, y0), (x1, y1)]),
     shapely.geometry.LineString([(x1, y1), (x0, y1)]),
     shapely.geometry.LineString([(x0, y1), (x0, y0)]),
     merged_shape.boundary]

PhysIDs = [0, 1, 2, 3, 4, 4]

folder_name = 'outlines'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)

with fiona.collection("outlines/outlines.shp", "w", "ESRI Shapefile", schema, crs=crs) as output:
    for i in range(len(features)):
        output.write({'geometry': shapely.geometry.mapping(features[i]), 'properties': {'PhysID': PhysIDs[1]}})

print('The domain is rectangular and bound by the following corner coordinates: \n' '(' + str(x0) + ', ' + str(y0) +
      '), (' + str(x1) + ', ' + str(y1) + ')', '\n')

# Array_1
xa0, ya0, xa1, ya1 = 0.5*length, 0.5*width + 100, 0.5*length + array_length, width - 100

array_1 = \
    [shapely.geometry.LineString([(xa0, ya0), (xa1, ya0)]),
     shapely.geometry.LineString([(xa1, ya0), (xa1, ya1)]),
     shapely.geometry.LineString([(xa1, ya1), (xa0, ya1)]),
     shapely.geometry.LineString([(xa0, ya1), (xa0, ya0)]),]



# Array 2
xa0, ya0, xa1, ya1 = 0.5*length, 100, 0.5*length + array_length, 0.5*width - 100

array_2 = \
    [shapely.geometry.LineString([(xa0, ya0), (xa1, ya0)]),
     shapely.geometry.LineString([(xa1, ya0), (xa1, ya1)]),
     shapely.geometry.LineString([(xa1, ya1), (xa0, ya1)]),
     shapely.geometry.LineString([(xa0, ya1), (xa0, ya0)]),]



folder_name = 'subdomains'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)


with fiona.collection("subdomains/array_1.shp", "w", "ESRI Shapefile", schema, crs=crs) as output:
    for i in range(len(array_1)):
        output.write({'geometry': shapely.geometry.mapping(array_1[i]), 'properties': {'PhysID': PhysIDs[2]}})
with fiona.collection("subdomains/array_2.shp", "w", "ESRI Shapefile", schema, crs=crs) as output:
    for i in range(len(array_2)):
        output.write({'geometry': shapely.geometry.mapping(array_2[i]), 'properties': {'PhysID': PhysIDs[3]}})


print('The subdomain is rectangular and bound by the following corner coordinates: \n' '(' + str(x_a0) + ', ' +
      str(y_a0) + '), (' + str(x_a1) + ', ' + str(y_a1) + ')', '\n')
