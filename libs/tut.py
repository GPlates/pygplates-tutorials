import pygplates
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon

import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
from shapely.geometry import MultiLineString

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
            print('Reconstructing magnetic picks in {}'.format(filename) )
            sys.stdout.flush()
            pygplates.reconstruct(
                    filename, 
                    rotation_filenames, 
                    magnetic_picks_output_basename+'_{}.shp'.format(cnt), 
                    self.reconstruction_time, 
                    self.anchor_plate)
    
    
    def plot_continental_polygons(self, ax, facecolor=None, edgecolor='none', alpha=0.1):
        print('Plotting continental polygons...')
        shape_feature = ShapelyFeature(Reader(continental_polygons_output_basename).geometries(),
                                ccrs.PlateCarree(), edgecolor=edgecolor)
        ax.add_feature(shape_feature,facecolor=facecolor, alpha=alpha)

        
    
    def plot_coastlines(self, ax, facecolor='default', edgecolor='k', alpha=0.4):
        print('Plotting coastlines...')

        for record in Reader(coastlines_output_basename).records():
            #print(record.attributes)
            fc = facecolor
            if facecolor == 'default':
                fc = plate_tectonic_utils.get_colour_by_plateid(int(record.attributes['PLATEID1']))
            
            shape_feature = ShapelyFeature([record.geometry], ccrs.PlateCarree(), edgecolor=edgecolor)
            ax.add_feature(shape_feature,facecolor=fc, alpha=alpha)
            
    
    def plot_topologies(self, ax, facecolor='default', edgecolor='w', alpha=0.2):
        print('Plotting topologies...')
        
        for record in Reader(topology_output_basename).records():
            #print(record.attributes)
            fc = facecolor
            if facecolor == 'default':
                fc = plate_tectonic_utils.get_colour_by_plateid(int(record.attributes['PLATEID1']))
            
            shape_feature = ShapelyFeature([record.geometry], ccrs.PlateCarree(), edgecolor=edgecolor)
            ax.add_feature(shape_feature,facecolor=fc, alpha=alpha)
        
    
    def plot_fracture_zones(self, ax, color='default'):
        #plot the fracture zones
        print('Plotting fracture zones...')
        
        for record in Reader(fracture_zones_output_basename).records():
            #print(record.geometry)
            c=color
            if color == 'default':
                c = plate_tectonic_utils.get_colour_by_plateid(int(record.attributes['PLATEID1']))
            
            if type(record.geometry) is MultiLineString:
                for line in record.geometry:
                    lon, lat = line.xy
                    ax.plot(lon,lat,transform=ccrs.Geodetic(),color=c)
            else:
                lon, lat = record.geometry.xy
                ax.plot(lon,lat,transform=ccrs.Geodetic(),color=c)
            
                  
    
    def plot_continent_ocean_boundaries(self, ax, color='default'):
        #plot the continent_ocean_boundaries
        
        for record in Reader(cob_output_basename).records():
            #print(record.geometry)
            c = color
            if color == 'default':
                c = plate_tectonic_utils.get_colour_by_plateid(int(record.attributes['PLATEID1']))
            
            if type(record.geometry) is MultiLineString:
                for line in record.geometry:
                    lon, lat = line.xy
                    ax.plot(lon,lat,transform=ccrs.Geodetic(),color=c)
            else:
                lon, lat = record.geometry.xy
                ax.plot(lon,lat,transform=ccrs.Geodetic(),color=c)
            
            
            
    
    def plot_magnetic_picks(self, ax, facecolors='none'):
        for i in range(1,5):
            reader = Reader(f'{magnetic_picks_output_basename}_{i}')
            #m.readshapefile(magnetic_picks_output_basename+'_{}'.format(i) ,'magnetic_{}'.format(i),drawbounds=False,color='w')
            import hashlib
            colors=[]
            lons=[]
            lats=[]
            for s in reader.records():
                if all(k in s.attributes for k in ['Chron', 'AnomalyEnd']):
                    #not really colouring by plate id. In fact, it is colouring by Chron+AnomalyEnd
                    colors.append(
                        plate_tectonic_utils.get_colour_by_plateid(
                            int(hashlib.sha1(
                                (s.attributes['Chron']+s.attributes['AnomalyEnd']).encode('utf-8')).hexdigest(), 16)))
                else:
                    colors.append(
                        plate_tectonic_utils.get_colour_by_plateid(0))
                lons.append(s.geometry.x)
                lats.append(s.geometry.y)
                
            ax.scatter(
                lons,
                lats,
                facecolors='none',
                edgecolor=colors,
                transform=ccrs.PlateCarree(),
                #s=.1,
                zorder=99)
            
            

    def plot_mineral_deposits(self, ax, facecolors='none'):
        reader = Reader(mineral_output_basename)
        
        ax.scatter(
                [point.x for point in reader.geometries()],
                [point.y for point in reader.geometries()],
                facecolors=facecolors,
                transform=ccrs.PlateCarree(),
                edgecolors='black',
                s=30,
                zorder=99)
        
        reader_au = Reader(mineral_AU_output_basename)
       
        ax.scatter(
                [point.x for point in reader_au.geometries()],
                [point.y for point in reader_au.geometries()],
                facecolors=facecolors,
                edgecolors='gold',
                transform=ccrs.PlateCarree(),
                s=30,
                zorder=99)
        
    
    def plot_velocities(self, ax):
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
    
        Q = ax.quiver(Xnodes,Ynodes,u,v,scale=1000,color='grey',transform = ccrs.PlateCarree(), regrid_shape=20)
        # make quiver key.
        qk = plt.quiverkey(Q, 0.95, 1.05, 50, '50 mm/yr', labelpos='W')


    
    def create_map(self, name, lon=0):
        if name == 'mollweide':
            ax = plt.axes(projection=ccrs.Mollweide(central_longitude=lon))
        elif name == 'robinson':
            ax = plt.axes(projection=ccrs.Robinson(central_longitude=lon))
        else: #'rectangular' and by default
            ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=lon))
        ax.gridlines()
        ax.set_global()
        return ax
    
    
    def plot_layers(self, layers, ax):
        for layer in layers:
            if layer == 'coastlines':
                self.reconstruct_coastlines() 
                self.plot_coastlines(ax)
            elif layer == 'continental_polygons':
                self.reconstruct_continental_polygons()
                self.plot_continental_polygons(ax)
            elif layer == 'topologies':
                self.reconstruct_topologies()
                self.plot_topologies(ax)
            elif layer == 'fracture':
                self.reconstruct_fracture_zones()
                self.plot_fracture_zones(ax)
            elif layer == 'magnetic':
                self.reconstruct_magnetic_picks()
                self.plot_magnetic_picks(ax)
            elif layer == 'cob':
                self.reconstruct_continent_ocean_boundaries()
                self.plot_continent_ocean_boundaries(ax)
            elif layer == 'mineral':
                self.reconstruct_mineral_deposits()
                self.plot_mineral_deposits(ax)
            elif layer == 'velocities':
                self.plot_velocities(ax)


    def plot_earthquakes(self, ax, minmag=0.0, maxmag=100.0):
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
    
        #x, y = m(pt_lon,pt_lat)
    
        sc = ax.scatter(
            pt_lon,  
            pt_lat,
            c=colors,
            s=sizes,
            transform=ccrs.PlateCarree(),
            vmin=min(magnitudes),
            vmax=max(magnitudes),
            cmap=cm,
            linewidths=0,
            zorder=2
        )

        cbar = plt.colorbar(sc)
        cbar.ax.set_ylabel('Richter magnitude')

