import pygplates
import sys

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon

import plate_tectonic_utils, velocity_utils

data_root = 'Data/workshop/'  
coastlines_filename = data_root+'Coastlines/Global_coastlines_2015_v1_low_res.shp'
continental_polygons_filename = data_root+'ContinentalPolygons/Seton_etal_ESR2012_ContinentalPolygons_2012.1.gpmlz'
topology_filenames = []
topology_filenames.append(data_root+'Global_EarthByte_Mesozoic-Cenozoic_plate_boundaries_Matthews_etal.gpml')
topology_filenames.append(data_root+'Global_EarthByte_Paleozoic_plate_boundaries_Matthews_etal.gpml')
fracture_zones_filename = data_root + 'Fracture_zones/Fracture_zones_workshop.shp'
magnetic_picks_filenames = []
magnetic_picks_filenames.append(data_root+'GSFML.IndianOcean/GSFML.IndianOcean.picks_ws.shp')
magnetic_picks_filenames.append(data_root+'GSFML.AtlanticOcean/GSFML.AtlanticOcean.picks_ws.shp')
magnetic_picks_filenames.append(data_root+'GSFML.MarginalBackarcBasins/GSFML.Marginal_BackarcBasins.picks_ws.shp')
magnetic_picks_filenames.append(data_root+'GSFML.PacificOcean/GSFML.PacificOcean.picks_ws.shp')
cob_filename = data_root + 'ContinentOceanBoundaries/Shapefile/Seton_etal_ESR2012_ContinentOceanBoundaries_2012.1_ws.shp'
mineral_filename = data_root + 'mrds/test.shp'
mineral_AU_filename = data_root + 'mrds/AU_deposits.shp'

rotation_filenames = []
rotation_filenames.append(data_root+'Global_EB_250-0Ma_GK07_Matthews_etal.rot')
rotation_filenames.append(data_root+'Global_EB_410-250Ma_GK07_Matthews_etal.rot')

# Base name of ouput file. 
coastlines_output_basename = '/tmp/coastlines'
continental_polygons_output_basename = '/tmp/continental'
topology_output_basename = '/tmp/topology'
fracture_zones_output_basename = '/tmp/fracture'
magnetic_picks_output_basename = '/tmp/magnetic'
cob_output_basename = '/tmp/cob'
mineral_output_basename = '/tmp/mineral'
mineral_AU_output_basename = '/tmp/mineral_au'

class Tutorial(object):
    def __init__(self):
        self.reconstruction_time = 0.
        self.anchor_plate = 0
        self.delta_time = 5.

    def reconstruct_coastlines(self):
        #
        #-----------Use pygplates to carry out the reconstruction 
        #
        print('Reconstructing coastlines...')
        pygplates.reconstruct(
                coastlines_filename, 
                rotation_filenames, 
                coastlines_output_basename+'.shp', 
                self.reconstruction_time, 
                self.anchor_plate)
    
    def reconstruct_continental_polygons(self):
        #
        #-----------Use pygplates to carry out the reconstruction 
        #
        print('Reconstructing continental polygons...')
        pygplates.reconstruct(
                continental_polygons_filename, 
                rotation_filenames, 
                continental_polygons_output_basename+'.shp', 
                self.reconstruction_time, 
                self.anchor_plate)

    def reconstruct_mineral_deposits(self):
        pygplates.reconstruct(
                mineral_filename, 
                rotation_filenames, 
                mineral_output_basename+'.shp', 
                self.reconstruction_time, 
                self.anchor_plate)
        
        pygplates.reconstruct(
                mineral_AU_filename, 
                rotation_filenames, 
                mineral_AU_output_basename+'.shp', 
                self.reconstruction_time, 
                self.anchor_plate)
    
    def reconstruct_topologies(self):
        print('Resolving topologies')
        rotation_model = pygplates.RotationModel(rotation_filenames)
        pygplates.resolve_topologies(
                topology_filenames, 
                rotation_model, 
                topology_output_basename+'.shp', 
                self.reconstruction_time)
    
    
    def reconstruct_fracture_zones(self):
        # Use pygplates to carry out the reconstruction 
        print('Reconstructing fracture zones...')
        pygplates.reconstruct(
                fracture_zones_filename, 
                rotation_filenames, 
                fracture_zones_output_basename+'.shp', 
                self.reconstruction_time, 
                self.anchor_plate)
    
    def reconstruct_continent_ocean_boundaries(self):
        # Use pygplates to carry out the reconstruction 
        pygplates.reconstruct(
                cob_filename, 
                rotation_filenames, 
                cob_output_basename+'.shp', 
                self.reconstruction_time, 
                self.anchor_plate)
    
    def reconstruct_magnetic_picks(self):
        sys.stdout.flush()
        cnt=0
        for filename in magnetic_picks_filenames:
            cnt+=1
            # Use pygplates to carry out the reconstruction
            print(('Reconstructing magnetic picks in {}'.format(filename) ))
            sys.stdout.flush()
            pygplates.reconstruct(
                    filename, 
                    rotation_filenames, 
                    magnetic_picks_output_basename+'_{}.shp'.format(cnt), 
                    self.reconstruction_time, 
                    self.anchor_plate)
    
    def plot_continental_polygons(self, m, facecolor=None, edgecolor='none', alpha=0.1):
        print('Plotting continental polygons...')
        m.readshapefile(continental_polygons_output_basename,'continental',drawbounds=False,color='w')   
        for s in zip(m.continental,m.continental_info):
            poly = Polygon(
                    s[0],
                    facecolor=facecolor, #default
                    edgecolor=edgecolor, #no color
                    alpha=alpha
            )
            plt.gca().add_patch(poly)
    
    def plot_coastlines(self, m, facecolor='default', edgecolor='k', alpha=0.4):
        print('Plotting coastlines...')

        m.readshapefile(coastlines_output_basename,'coastlines',drawbounds=False,color='w')   
        for s in zip(m.coastlines,m.coastlines_info):
            fc = facecolor
            if facecolor == 'default':
                fc = plate_tectonic_utils.get_colour_by_plateid(int(s[1]['PLATEID1']))

            poly = Polygon(
                    s[0],
                    facecolor=fc,
                    edgecolor=edgecolor,
                    alpha=alpha)
            plt.gca().add_patch(poly)
    
    def plot_topologies(self, m, facecolor='default', edgecolor='w', alpha=0.2):
        print('Plotting topologies...')
        m.readshapefile(topology_output_basename,'topologies',drawbounds=False,color='w')  
        for s in zip(m.topologies,m.topologies_info):
            fc = facecolor
            if facecolor == 'default':
                fc = plate_tectonic_utils.get_colour_by_plateid(int(s[1]['PLATEID1']))

            poly = Polygon(
                    s[0],
                    facecolor=fc,
                    edgecolor=edgecolor,
                    alpha=alpha)
            plt.gca().add_patch(poly)
    
    def plot_fracture_zones(self, m, color='default'):
        #plot the fracture zones
        print('Plotting fracture zones...')
        m.readshapefile(fracture_zones_output_basename,'fracture',drawbounds=False,color='b')
        for s in zip(m.fracture,m.fracture_info):
            fc = color
            if color == 'default':
                fc = plate_tectonic_utils.get_colour_by_plateid(int(s[1]['PLATEID1']))
            
            m.plot(
                [x for (x,y) in s[0]],
                [y for (x,y) in s[0]],
                color=fc)
    
    def plot_continent_ocean_boundaries(self, m, color='default'):
        #plot the fracture zones
        m.readshapefile(cob_output_basename,'cob',drawbounds=False,color='b')
        for s in zip(m.cob, m.cob_info):
            fc = color
            if color == 'default':
                fc = plate_tectonic_utils.get_colour_by_plateid(int(s[1]['PLATEID1']))

            m.plot(
                [x for (x,y) in s[0]],
                [y for (x,y) in s[0]],
                color=fc)
    
    def plot_magnetic_picks(self, m, facecolors='none'):
        for i in range(1,5):
            m.readshapefile(magnetic_picks_output_basename+'_{}'.format(i) ,'magnetic_{}'.format(i),drawbounds=False,color='w')
            import hashlib
            colors=[]
            for s in getattr(m,'magnetic_{}_info'.format(i)):
                if all(k in s for k in ['Chron', 'AnomalyEnd']):
                    #not really colouring by plate id. In fact, it is colouring by Chron+AnomalyEnd
                    colors.append(
                        plate_tectonic_utils.get_colour_by_plateid(int(hashlib.sha1(s['Chron']+s['AnomalyEnd']).hexdigest(), 16)))
                else:
                    colors.append(
                        plate_tectonic_utils.get_colour_by_plateid(0))
            m.scatter(
                [x for (x,y) in getattr(m,'magnetic_{}'.format(i))],
                [y for (x,y) in getattr(m,'magnetic_{}'.format(i))],
                facecolors='none',
                edgecolor=colors,
                #s=.1,
                zorder=99)

    def plot_mineral_deposits(self, m, facecolors='none'):
        m.readshapefile(mineral_output_basename,'mineral',drawbounds=False,color='w')
        m.scatter(
                [x for (x,y) in m.mineral],
                [y for (x,y) in m.mineral],
                facecolors=facecolors,
                edgecolors='black',
                s=30,
                zorder=99)
        m.readshapefile(mineral_AU_output_basename,'mineral_au',drawbounds=False,color='w')
        m.scatter(
                [x for (x,y) in m.mineral_au],
                [y for (x,y) in m.mineral_au],
                facecolors=facecolors,
                edgecolors='gold',
                s=30,
                zorder=99)
    
    def plot_velocities(self, m):
        delta_time = 5.
        Xnodes = np.arange(-180,180,5)
        Ynodes = np.arange(-90,90,5)
        Xg,Yg = np.meshgrid(Xnodes,Ynodes)
        Xg = Xg.flatten()
        Yg = Yg.flatten()
        velocity_domain_features = velocity_utils.make_GPML_velocity_feature(Xg,Yg)
        # Load one or more rotation files into a rotation model.
        rotation_model = pygplates.RotationModel(rotation_filenames)
    
        # Load the topological plate polygon features.
        topology_features = []
        for fname in topology_filenames:
            for f in pygplates.FeatureCollection(fname):
                topology_features.append(f)
    
    
        # Call the function we created above to get the velocities
        all_velocities = velocity_utils.Get_Plate_Velocities(velocity_domain_features,
                                              topology_features,
                                              rotation_model,
                                              self.reconstruction_time,
                                              delta_time,
                                              'vector_comp')
    
        uu=[]
        vv=[]
        for vel in all_velocities:
            if not hasattr(vel, 'get_y'): 
                uu.append(vel[1])
                vv.append(vel[0])
            else:
                uu.append(vel.get_y())
                vv.append(vel.get_x())
        u = np.asarray([uu]).reshape((Ynodes.shape[0],Xnodes.shape[0]))
        v = np.asarray([vv]).reshape((Ynodes.shape[0],Xnodes.shape[0]))
    
        # compute native x,y coordinates of grid.
        x, y = m(Xg, Yg)
    
        uproj,vproj,xx,yy = m.transform_vector(u,v,Xnodes,Ynodes,15,15,returnxy=True,masked=True)
        # now plot.
        Q = m.quiver(xx,yy,uproj,vproj,scale=1000,color='grey')
        # make quiver key.
        qk = plt.quiverkey(Q, 0.95, 1.05, 50, '50 mm/yr', labelpos='W')


    
    def create_map(self, name, lat=None, lon=None, width=None, height=None):
        if name == 'mollweide':
            m = Basemap(projection='moll', lon_0=0.0, resolution=None)
            m.drawparallels(np.arange(-90.,91.,15.), labels=[True,True,False,False])
            m.drawmeridians(np.arange(-180.,181.,45.), labels=[False,False,False,False])
        elif name == 'robinson':
            m = Basemap(projection='robin',lon_0=0.0,resolution=None)
            m.drawparallels(np.arange(-90.,91.,15.), labels=[True,True,False,False])
            m.drawmeridians(np.arange(-180.,181.,45.), labels=[False,False,False,True])
        else: #'rectangular' and by default
            if lat:
                m = Basemap(
                    projection='cyl',
                    lat_0=lat,
                    lon_0=lon,
                    width=width,
                    height=height,
                    resolution=None)
            else:
                m = Basemap(
                    projection='cyl',
                    llcrnrlat=-90,
                    urcrnrlat=90,
                    llcrnrlon=-180.,
                    urcrnrlon=180.1,
                    resolution=None)
            m.drawparallels(np.arange(-90.,91.,15.), labels=[True,True,False,False])
            m.drawmeridians(np.arange(-180.,181.,30.), labels=[False,False,False,True])
        return m
    
    def plot_layers(self, layers, m):
        for layer in layers:
            if layer == 'coastlines':
                self.reconstruct_coastlines() 
                self.plot_coastlines(m)
            elif layer == 'continental_polygons':
                self.reconstruct_continental_polygons()
                self.plot_continental_polygons(m)
            elif layer == 'topologies':
                self.reconstruct_topologies()
                self.plot_topologies(m)
            elif layer == 'fracture':
                self.reconstruct_fracture_zones()
                self.plot_fracture_zones(m)
            elif layer == 'magnetic':
                self.reconstruct_magnetic_picks()
                self.plot_magnetic_picks(m)
            elif layer == 'cob':
                self.reconstruct_continent_ocean_boundaries()
                self.plot_continent_ocean_boundaries(m)
            elif layer == 'mineral':
                self.reconstruct_mineral_deposits()
                self.plot_mineral_deposits(m)
            elif layer == 'velocities':
                self.plot_velocities(m)


    def plot_earthquakes(self, m, minmag=0.0, maxmag=100.0):
        earthquakes = pygplates.FeatureCollection('Data/workshop/Earthquakes/earthquakes_new1.shp')

        cm = plt.cm.get_cmap('gnuplot')

        colors = []
        pt_lon = []
        pt_lat = []
        sizes = []
        magnitudes = []
        for q in earthquakes:
            mag = float(q.get_name())

            if mag < minmag:
                continue

            if mag > maxmag:
                continue

            pt = q.get_geometry()
            pt_lon.append(pt.to_lat_lon()[1])
            pt_lat.append(pt.to_lat_lon()[0])
            
            colors.append(mag)
            sizes.append(35)
            magnitudes.append(mag)
    
        x, y = m(pt_lon,pt_lat)
    
        sc = m.scatter(
            x, y,
            c=colors,
            s=sizes,
            vmin=min(magnitudes),
            vmax=max(magnitudes),
            cmap=cm,
            linewidths=0,
            zorder=2
        )

        cbar = plt.colorbar(sc)
        cbar.ax.set_ylabel('Richter magnitude')

