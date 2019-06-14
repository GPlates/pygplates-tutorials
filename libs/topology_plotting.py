import pygplates
import numpy as np
from mpl_toolkits.basemap import Basemap


def make_GPML_velocity_feature(Long,Lat):
# function to make a velocity mesh nodes at an arbitrary set of points defined in Lat
# Long and Lat are assumed to be 1d arrays. 

    # Add points to a multipoint geometry
    SeedPoints = zip(Lat,Long)
    points = []
    for j in range(0,len(SeedPoints)):
        points.append(SeedPoints[j])
    multi_point = pygplates.MultiPointOnSphere(points)

    # Create a feature containing the multipoint feature, and defined as MeshNode type
    meshnode_feature = pygplates.Feature(pygplates.FeatureType.create_from_qualified_string('gpml:MeshNode'))
    meshnode_feature.set_geometry(multi_point)
    meshnode_feature.set_name('Velocity Mesh Nodes from pygplates')

    output_feature_collection = pygplates.FeatureCollection(meshnode_feature)
    
    # NB: at this point, the feature could be written to a file using
    # output_feature_collection.write('myfilename.gpmlz')
    
    # for use within the notebook, the velocity domain feature is returned from the function
    return output_feature_collection


def get_plate_velocities(velocity_domain_features,
                         topology_features,
                         rotation_model,
                         time,
                         delta_time,
                         rep):

    # Define a function 

    # All domain points and associated (magnitude, azimuth, inclination) velocities for the current time.
    all_domain_points = []
    all_velocities = []

    # Partition our velocity domain features into our topological plate polygons at the current 'time'.
    plate_partitioner = pygplates.PlatePartitioner(topology_features, rotation_model, time)

    for velocity_domain_feature in velocity_domain_features:

        # A velocity domain feature usually has a single geometry but we'll assume it can be any number.
        # Iterate over them all.
        for velocity_domain_geometry in velocity_domain_feature.get_geometries():

            for velocity_domain_point in velocity_domain_geometry.get_points():

                all_domain_points.append(velocity_domain_point)

                partitioning_plate = plate_partitioner.partition_point(velocity_domain_point)
                if partitioning_plate:

                    # We need the newly assigned plate ID to get the equivalent stage rotation of that tectonic plate.
                    partitioning_plate_id = partitioning_plate.get_feature().get_reconstruction_plate_id()

                    # Get the stage rotation of partitioning plate from 'time + delta_time' to 'time'.
                    equivalent_stage_rotation = rotation_model.get_rotation(time, partitioning_plate_id, time + delta_time)

                    # Calculate velocity at the velocity domain point.
                    # This is from 'time + delta_time' to 'time' on the partitioning plate.
                    velocity_vectors = pygplates.calculate_velocities(
                        [velocity_domain_point],
                        equivalent_stage_rotation,
                        delta_time)
                    
                    if rep=='mag_azim':
                        # Convert global 3D velocity vectors to local (magnitude, azimuth, inclination) tuples (one tuple per point).
                        velocities = pygplates.LocalCartesian.convert_from_geocentric_to_magnitude_azimuth_inclination(
                            [velocity_domain_point],
                            velocity_vectors)
                        all_velocities.append(velocities[0])

                    elif rep=='vector_comp':
                        # Convert global 3D velocity vectors to local (magnitude, azimuth, inclination) tuples (one tuple per point).
                        velocities = pygplates.LocalCartesian.convert_from_geocentric_to_north_east_down(
                                [velocity_domain_point],
                                velocity_vectors)
                        all_velocities.append(velocities[0])

                else:
                    all_velocities.append((0,0,0))

    return all_velocities


def plot_velocities_and_topologies(pmap,
                                   topology_features,
                                   rotation_model,
                                   time,
                                   delta_time=5,
                                   res=10,
                                   scale=2000,
                                   lon0=0,
                                   clip_path=None,
                                   alpha=0.4):

    Xnodes = np.arange(-180,180,res)
    Ynodes = np.arange(-90,90,res)
    Xg,Yg = np.meshgrid(Xnodes,Ynodes)

    velocity_domain_features = make_GPML_velocity_feature(Xg.flatten(),Yg.flatten())

    # Call the function we created above to get the velocities
    all_velocities = get_plate_velocities(velocity_domain_features,
                                          topology_features,
                                          rotation_model,
                                          time,
                                          delta_time,
                                          'vector_comp')
    
    # The rest of the cell is for plotting, including rendering resolved topological boundaries to the map
    pt_vel_n=[]
    pt_vel_e=[]
    for vel in all_velocities:
        pt_vel_e.append(vel.get_y())
        pt_vel_n.append(vel.get_x())
    
    u = np.asarray(pt_vel_e).reshape((Ynodes.shape[0],Xnodes.shape[0]))
    v = np.asarray(pt_vel_n).reshape((Ynodes.shape[0],Xnodes.shape[0]))
        
    # Resolve our topological plate polygons (and deforming networks) to the current 'time'.
    # We generate both the resolved topology boundaries and the boundary sections between them.
    resolved_topologies = []
    shared_boundary_sections = []
    pygplates.resolve_topologies(topology_features, rotation_model, resolved_topologies, time, shared_boundary_sections)
    
    # create a dateline wrapper object
    wrapper = pygplates.DateLineWrapper(lon0)
    
    # Iterate over the shared boundary sections.
    for shared_boundary_section in shared_boundary_sections:
    
        # The shared sub-segments contribute either to the ridges or to the subduction zones.
        if shared_boundary_section.get_feature().get_feature_type() == pygplates.FeatureType.create_gpml('MidOceanRidge'):
            # Ignore zero length segments - they don't have a direction.
            for shared_sub_segment in shared_boundary_section.get_shared_sub_segments():
                split_geometry = wrapper.wrap(shared_sub_segment.get_geometry())
                for geometry in split_geometry:
                    X=[]
                    Y=[]
                    for point in geometry.get_points():
                        X.append(point.get_longitude()),Y.append(point.get_latitude())
                    x,y = pmap(X,Y)
                    pmap.plot(x,y,'r',clip_path=clip_path,linewidth=3,alpha=alpha, zorder=1)     
    
        elif shared_boundary_section.get_feature().get_feature_type() == pygplates.FeatureType.create_gpml('SubductionZone'):
            for shared_sub_segment in shared_boundary_section.get_shared_sub_segments():
                split_geometry = wrapper.wrap(shared_sub_segment.get_geometry())
                for geometry in split_geometry:
                    X=[]
                    Y=[]
                    for point in geometry.get_points():
                        X.append(point.get_longitude()),Y.append(point.get_latitude())
                    x,y = pmap(X,Y)
                pmap.plot(x,y,'k',clip_path=clip_path,linewidth=3,alpha=alpha, zorder=1)  
    
        else: #shared_boundary_section.get_feature().get_feature_type() == pygplates.FeatureType.create_gpml('FractureZone'):
            for shared_sub_segment in shared_boundary_section.get_shared_sub_segments():
                split_geometry = wrapper.wrap(shared_sub_segment.get_geometry())
                for geometry in split_geometry:
                    X=[]
                    Y=[]
                    for point in geometry.get_points():
                        X.append(point.get_longitude()),Y.append(point.get_latitude())
                    x,y = pmap(X,Y)
                pmap.plot(x,y,'b',clip_path=clip_path,linewidth=3,alpha=alpha, zorder=1)  
    
    lons, lats = np.meshgrid(Xnodes,Ynodes)
    # compute native x,y coordinates of grid.
    x, y = pmap(lons, lats)
    # define parallels and meridians to draw.
    
    uproj,vproj,xx,yy = \
    pmap.transform_vector(u,v,Xnodes,Ynodes,54,26,returnxy=True,masked=True)
    # now plot.
    Q = pmap.quiver(xx,yy,uproj,vproj,scale=scale,clip_path=clip_path,alpha=alpha)
    #Q2 = pmap.quiver(xx,yy,uproj,vproj,scale=scale,color='')
    # make quiver key.
