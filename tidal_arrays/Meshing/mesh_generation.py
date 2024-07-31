

def mesh():

    boundaries = qmesh.vector.Shapes()  # instantiate Shapes class
    boundaries.fromFile('outlines/outlines.shp')  # add the shapes from Shapefile
    loopShapes = qmesh.vector.identifyLoops(boundaries, isGlobal=False, defaultPhysID=1000, fixOpenLoops=True)
    polygonShapes = qmesh.vector.identifyPolygons(loopShapes, smallestNotMeshedArea=100, meshedAreaPhysID=1)

    array_plot = qmesh.vector.Shapes()
    array_plot.fromFile('subdomains/array_1.shp')
    array_shapes = qmesh.vector.identifyLoops(array_plot, isGlobal=False, defaultPhysID=1000, fixOpenLoops=True)
    array_polygon = qmesh.vector.identifyPolygons(array_shapes, smallestNotMeshedArea=100, meshedAreaPhysID=2)

    array_plot_2 = qmesh.vector.Shapes()
    array_plot_2.fromFile('subdomains/array_2.shp')
    array_shapes_2 = qmesh.vector.identifyLoops(array_plot_2, isGlobal=False, defaultPhysID=1000, fixOpenLoops=True)
    array_polygon_2 = qmesh.vector.identifyPolygons(array_shapes_2, smallestNotMeshedArea=100, meshedAreaPhysID=3)

    gradation = []
    gradation.append(qmesh.raster.meshMetricTools.gradationToShapes())
    gradation[-1].setShapes(array_polygon,)
    gradation[-1].setRasterBounds(-1000, 8000., -1000, 3000.0)
    gradation[-1].setRasterResolution(450, 150)
    gradation[-1].setGradationParameters(10.0, 50.0, 500.0, 10)
    gradation[-1].calculateLinearGradation()
    gradation[-1].writeNetCDF('gradation_to_turbine.nc')

    gradation.append(qmesh.raster.meshMetricTools.gradationToShapes())
    gradation[-1].setShapes(array_polygon_2,)
    gradation[-1].setRasterBounds(-1000, 8000., -1000, 3000.0)
    gradation[-1].setRasterResolution(450, 150)
    gradation[-1].setGradationParameters(10.0, 50.0, 500.0, 10)
    gradation[-1].calculateLinearGradation()
    gradation[-1].writeNetCDF('gradation_to_turbine_2.nc')

    if len(gradation) == 1:
        meshMetricRaster = gradation[0]
    else:
        meshMetricRaster = qmesh.raster.meshMetricTools.minimumRaster(gradation)

    meshMetricRaster.writeNetCDF('meshMetric.nc')

    domain = qmesh.mesh.Domain()
    domainLines, domainPolygons = qmesh.vector.insertRegions(loopShapes, polygonShapes, array_shapes, array_polygon)
    domainLines_2, domainPolygons_2 = qmesh.vector.insertRegions(domainLines, domainPolygons, array_shapes_2, array_polygon_2)
    domain.setGeometry(domainLines_2, domainPolygons_2)
    domain.setMeshMetricField(meshMetricRaster)
    domain.setTargetCoordRefSystem('EPSG:32630', fldFillValue=1000.0)

    domain.gmsh(geoFilename='tidal_0.geo',
                fldFilename='tidal_0.fld',
                mshFilename='idealised_channel.msh')

def convertMesh(meshname):

    mesh = qmesh.mesh.Mesh()
    mesh.readGmsh(meshname+'.msh', 'EPSG:32630')
    mesh.writeShapefile(meshname)

if __name__ == '__main__':
    import qmesh
    qmesh.initialise()
    mesh()
    convertMesh('idealised_channel')
