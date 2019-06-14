from mpl_toolkits.basemap import Basemap
import numpy as np
#
#---------function definitions----------
#---------you might want to move along to the main() funtion section. more fun over there---------------
#
#function: assign colour by plate id
def get_colour_by_plateid(plate_id):
    from matplotlib import colors
    converter = colors.ColorConverter()
    plateid_colours = {
        0: list(converter.to_rgba('yellow',alpha=1.0)),
        1: list(converter.to_rgba('aqua',alpha=1.0)),
        2: list(converter.to_rgba('seagreen',alpha=1.0)),
        3: list(converter.to_rgba('fuchsia',alpha=1.0)),
        4: list(converter.to_rgba('slategray',alpha=1.0)),
        5: list(converter.to_rgba('lime',alpha=1.0)),
        6: list(converter.to_rgba('indigo',alpha=1.0)),
        7: list(converter.to_rgba('red',alpha=1.0)),
        8: list(converter.to_rgba('orange',alpha=1.0)),
        9: list(converter.to_rgba('lightsalmon',alpha=1.0)),
        10: list(converter.to_rgba('navy',alpha=1.0)),
    }
    return plateid_colours[plate_id%11]

#function: create mollweide map
def create_mollweide_map():
    m = Basemap(projection='moll', lon_0=0.0, resolution=None)
    m.drawparallels(np.arange(-90.,91.,15.), labels=[True,True,False,False])
    m.drawmeridians(np.arange(-180.,181.,45.), labels=[False,False,False,False])
    return m

#function: create robinson map
def create_robinson_map():
    m = Basemap(projection='robin',lon_0=0.0,resolution=None)
    m.drawparallels(np.arange(-90.,91.,15.), labels=[True,True,False,False])
    m.drawmeridians(np.arange(-180.,181.,45.), labels=[False,False,False,True])
    return m

#function: create rectangular map
def create_rectangular_map():
    lon_0 = 0
    m = Basemap(
            projection='cyl',
            llcrnrlat=-90,
            urcrnrlat=90,
            llcrnrlon=-180.+lon_0,
            urcrnrlon=180.1+lon_0,
            resolution=None)
    m.drawparallels(np.arange(-90.,91.,15.), labels=[True,True,False,False])
    m.drawmeridians(np.arange(-180.,181.,30.), labels=[False,False,False,True])
    return m