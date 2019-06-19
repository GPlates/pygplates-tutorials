
"""
    Copyright (C) 2016 The University of Sydney, Australia
    
    This program is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License, version 2, as published by
    the Free Software Foundation.
    
    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    
    You should have received a copy of the GNU General Public License along
    with this program; if not, write to Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""


##################################################
# Find the convergence rate of subduction zones. #
##################################################


from __future__ import print_function
import argparse
import subprocess
import math
import sys
import pygplates


# Required pygplates version.
PYGPLATES_VERSION_REQUIRED = pygplates.Version(12)

# The default threshold sampling distance along subduction zones.
DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES = 0.5
DEFAULT_THRESHOLD_SAMPLING_DISTANCE_KMS = math.radians(DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES) * pygplates.Earth.equatorial_radius_in_kms

DEFAULT_TIME_RANGE_YOUNG_TIME = 0
DEFAULT_TIME_RANGE_OLD_TIME = 200
DEFAULT_TIME_INCREMENT = 1

DEFAULT_VELOCITY_DELTA_TIME = 1


# Determine the overriding and subducting plates of the subduction shared sub-segment.
def find_overriding_and_subducting_plates(subduction_shared_sub_segment, time):
    
    # Get the subduction polarity of the nearest subducting line.
    subduction_polarity = subduction_shared_sub_segment.get_feature().get_enumeration(pygplates.PropertyName.gpml_subduction_polarity)
    if (not subduction_polarity) or (subduction_polarity == 'Unknown'):
        print('Unable to find the overriding plate of the subducting shared sub-segment "{0}"'.format(
            subduction_shared_sub_segment.get_feature().get_name()), file=sys.stderr)
        print('    subduction zone feature is missing subduction polarity property or it is set to "Unknown".', file=sys.stderr)
        return

    # There should be two sharing topologies - one is the overriding plate and the other the subducting plate.
    sharing_resolved_topologies = subduction_shared_sub_segment.get_sharing_resolved_topologies()
    if len(sharing_resolved_topologies) != 2:
        print('Unable to find the overriding and subducting plates of the subducting shared sub-segment "{0}" at {1}Ma'.format(
            subduction_shared_sub_segment.get_feature().get_name(), time), file=sys.stderr)
        print('    there are not exactly 2 topologies sharing the sub-segment.', file=sys.stderr)
        return

    overriding_plate = None
    subducting_plate = None
    
    geometry_reversal_flags = subduction_shared_sub_segment.get_sharing_resolved_topology_geometry_reversal_flags()
    for index in range(2):

        sharing_resolved_topology = sharing_resolved_topologies[index]
        geometry_reversal_flag = geometry_reversal_flags[index]

        if sharing_resolved_topology.get_resolved_boundary().get_orientation() == pygplates.PolygonOnSphere.Orientation.clockwise:
            # The current topology sharing the subducting line has clockwise orientation (when viewed from above the Earth).
            # If the overriding plate is to the 'left' of the subducting line (when following its vertices in order) and
            # the subducting line is reversed when contributing to the topology then that topology is the overriding plate.
            # A similar test applies to the 'right' but with the subducting line not reversed in the topology.
            if ((subduction_polarity == 'Left' and geometry_reversal_flag) or
                (subduction_polarity == 'Right' and not geometry_reversal_flag)):
                overriding_plate = sharing_resolved_topology
            else:
                subducting_plate = sharing_resolved_topology
        else:
            # The current topology sharing the subducting line has counter-clockwise orientation (when viewed from above the Earth).
            # If the overriding plate is to the 'left' of the subducting line (when following its vertices in order) and
            # the subducting line is not reversed when contributing to the topology then that topology is the overriding plate.
            # A similar test applies to the 'right' but with the subducting line reversed in the topology.
            if ((subduction_polarity == 'Left' and not geometry_reversal_flag) or
                (subduction_polarity == 'Right' and geometry_reversal_flag)):
                overriding_plate = sharing_resolved_topology
            else:
                subducting_plate = sharing_resolved_topology
    
    if overriding_plate is None:
        print('Unable to find the overriding plate of the subducting shared sub-segment "{0}" at {1}Ma'.format(
            subduction_shared_sub_segment.get_feature().get_name(), time), file=sys.stderr)
        print('    both sharing topologies are on subducting side of subducting line.', file=sys.stderr)
        return
    
    if subducting_plate is None:
        print('Unable to find the subducting plate of the subducting shared sub-segment "{0}" at {1}Ma'.format(
            subduction_shared_sub_segment.get_feature().get_name(), time), file=sys.stderr)
        print('    both sharing topologies are on overriding side of subducting line.', file=sys.stderr)
        return
    
    return (overriding_plate, subducting_plate, subduction_polarity)


# Resolves topologies at 'time', tessellates all resolved subduction zones to within
# 'threshold_sampling_distance_radians' radians and the following convergence-related parameters
# at each tessellates point:
#
# - point longitude
# - point latitude
# - subducting convergence (relative to overriding plate) velocity magnitude (in cm/yr)
# - subducting convergence velocity obliquity angle (angle between subduction zone normal vector and convergence velocity vector)
# - subduction zone absolute (relative to anchor plate) velocity magnitude (in cm/yr)
# - subduction zone absolute velocity obliquity angle (angle between subduction zone normal vector and absolute velocity vector)
# - length of arc segment (in degrees) that current point is on
# - subducting arc normal azimuth angle (clockwise starting at North, ie, 0 to 360 degrees) at current point
# - subducting plate ID
# - overriding plate ID
#
# Note that the convergence and absolute velocity magnitudes are negative if the plates are diverging (if convergence obliquity angle
# is greater than 90 or less than -90).

def subduction_convergence(
        # Rotation model or feature collection(s), or list of features, or filename(s)...
        rotation_features_or_model,
        # Topology feature collection(s), or list of features, or filename(s) or any combination of those...
        topology_features,
        # Threshold sampling distance along subduction zones (in radians)...
        threshold_sampling_distance_radians,
        time,
        velocity_delta_time = 1.0,
        anchor_plate_id = 0):
    
    # Turn rotation data into a RotationModel (if not already).
    rotation_model = pygplates.RotationModel(rotation_features_or_model)
    
    # Turn topology data into a list of features (if not already).
    topology_features = pygplates.FeaturesFunctionArgument(topology_features)
    
    # Resolve our topological plate polygons (and deforming networks) to the current 'time'.
    # We generate both the resolved topology boundaries and the boundary sections between them.
    resolved_topologies = []
    shared_boundary_sections = []
    pygplates.resolve_topologies(topology_features.get_features(), rotation_model, resolved_topologies, time, shared_boundary_sections, anchor_plate_id)
    
    # List of tesselated subduction zone shared subsegment points and associated convergence parameters
    # for the current 'time'.
    output_data = []
    
    # Iterate over the shared boundary sections of all resolved topologies.
    for shared_boundary_section in shared_boundary_sections:
    
        # Skip sections that are not subduction zones.
        if shared_boundary_section.get_feature().get_feature_type() != pygplates.FeatureType.gpml_subduction_zone:
            continue
        
        # Iterate over the shared sub-segments of the current subducting line.
        # These are the parts of the subducting line that actually contribute to topological boundaries.
        for shared_sub_segment in shared_boundary_section.get_shared_sub_segments():
        
            # Find the overriding and subducting plates on either side of the shared sub-segment.
            overriding_and_subducting_plates = find_overriding_and_subducting_plates(shared_sub_segment, time)
            if not overriding_and_subducting_plates:
                continue
            overriding_plate, subducting_plate, subduction_polarity = overriding_and_subducting_plates
            overriding_plate_id = overriding_plate.get_feature().get_reconstruction_plate_id()
            subducting_plate_id = subducting_plate.get_feature().get_reconstruction_plate_id()
            
            # The plate ID of the subduction zone line (as opposed to the subducting plate).
            #
            # Update: The plate IDs of the subduction zone line and overriding plate can differ
            # even in a non-deforming model due to smaller plates, not modelled by topologies, moving
            # differently than the larger topological plate being modelled - and the subduction zone line
            # having plate IDs of the smaller plates near them. For that reason we use the plate ID
            # of the subduction zone line whenever we can. Since some subduction zone lines can be
            # topological lines, they might actually be deforming (or intended to be deforming) and
            # hence their plate ID is not meaningful or at least we can't be sure whether it will
            # be zero or the overriding plate (or something else). So if the subduction zone line
            # is a topological line then we'll use the overriding plate ID instead.
            #
            if isinstance(shared_boundary_section.get_topological_section(), pygplates.ResolvedTopologicalLine):
                subduction_zone_plate_id = overriding_plate_id
            else:
                subduction_zone_plate_id = shared_sub_segment.get_feature().get_reconstruction_plate_id()
            
            # Get the rotation of the subducting plate relative to the overriding plate
            # from 'time + velocity_delta_time' to 'time'.
            convergence_relative_stage_rotation = rotation_model.get_rotation(
                    time,
                    subducting_plate_id,
                    time + velocity_delta_time,
                    overriding_plate_id,
                    anchor_plate_id=anchor_plate_id)
            #
            # The subduction zones have been reconstructed using the rotation "R(0->t2,A->M)":
            #
            #   reconstructed_geometry = R(0->t2,A->M) * present_day_geometry
            #
            # We can write "R(0->t2,A->M)" in terms of the stage rotation "R(t1->t2,F->M)" as:
            #
            #   R(0->t2,A->M) = R(0->t2,A->F) * R(0->t2,F->M)
            #                 = R(0->t2,A->F) * R(t1->t2,F->M) * R(0->t1,F->M)
            #                 = R(0->t2,A->F) * stage_rotation * R(0->t1,F->M)
            #
            # ...where 't1' is 't+1' and 't2' is 't' (ie, from 't1' to 't2').
            #
            # So to get the *reconstructed* geometry into the stage rotation reference frame
            # we need to rotate it by "inverse[R(0->t2,A->F)]":
            #
            #   reconstructed_geometry = R(0->t2,A->F) * stage_rotation * R(0->t1,F->M) * present_day_geometry
            #   inverse[R(0->t2,A->F)] * reconstructed_geometry = stage_rotation * R(0->t1,F->M) * present_day_geometry
            #
            # Once we've done that we can calculate the velocities of those geometry points
            # using the stage rotation. Then the velocities need to be rotated back from the
            # stage rotation reference frame using the rotation "R(0->t2,A->F)".
            # 
            from_stage_frame_relative_to_overriding = rotation_model.get_rotation(
                    time,
                    overriding_plate_id,
                    anchor_plate_id=anchor_plate_id)
            to_stage_frame_relative_to_overriding = from_stage_frame_relative_to_overriding.get_inverse()
            
            # Get the rotation of the subduction zone relative to the anchor plate
            # from 'time + velocity_delta_time' to 'time'.
            #
            # Note: We don't need to convert to and from the stage rotation reference frame
            # like the above convergence because this stage rotation is relative to the anchor plate
            # and so the above to/from stage rotation frame conversion "R(0->t2,A->F)" is the
            # identity rotation since the fixed plate (F) is the anchor plate (A).
            subduction_zone_equivalent_stage_rotation = rotation_model.get_rotation(
                    time,
                    subduction_zone_plate_id,
                    time + velocity_delta_time,
                    anchor_plate_id=anchor_plate_id)
            
            # We need to reverse the subducting_normal vector direction if overriding plate is to
            # the right of the subducting line since great circle arc normal is always to the left.
            if subduction_polarity == 'Left':
                subducting_normal_reversal = 1
            else:
                subducting_normal_reversal = -1
            
            # Ensure the shared sub-segment is tessellated to within the threshold sampling distance.
            tessellated_shared_sub_segment_polyline = (
                    shared_sub_segment.get_resolved_geometry().to_tessellated(threshold_sampling_distance_radians))
            
            # Iterate over the great circle arcs of the tessellated polyline to get the
            # arc midpoints, lengths and subducting normals.
            # There is an arc between each adjacent pair of points in the polyline.
            arc_midpoints = []
            arc_lengths = []
            subducting_arc_normals = []
            for arc in tessellated_shared_sub_segment_polyline.get_segments():
                if not arc.is_zero_length():
                    arc_midpoints.append(arc.get_arc_point(0.5))
                    arc_lengths.append(arc.get_arc_length())
                    # The normal to the subduction zone in the direction of subduction (towards overriding plate).
                    subducting_arc_normals.append(subducting_normal_reversal * arc.get_great_circle_normal())
            
            # Shouldn't happen, but just in case the shared sub-segment polyline coincides with a point.
            if not arc_midpoints:
                continue
            
            # The subducting arc normals relative to North (azimuth).
            # Convert global 3D normal vectors to local (magnitude, azimuth, inclination) tuples (one tuple per point).
            subducting_arc_local_normals = pygplates.LocalCartesian.convert_from_geocentric_to_magnitude_azimuth_inclination(
                    arc_midpoints, subducting_arc_normals)
            
            # Calculate the convergence velocities, and subduction zone velocities relative to
            # overriding plate, at the arc midpoints.
            #
            # Note; We need to convert the reconstructed geometry points into the stage rotation
            # reference frame to calculate velocities and then convert the velocities using the
            # reverse transform as mentioned above.
            arc_midpoints_in_stage_frame_relative_to_overriding = [
                    to_stage_frame_relative_to_overriding * arc_midpoint
                            for arc_midpoint in arc_midpoints]
            convergence_velocity_vectors_in_stage_frame_relative_to_overriding = pygplates.calculate_velocities(
                    arc_midpoints_in_stage_frame_relative_to_overriding,
                    convergence_relative_stage_rotation,
                    velocity_delta_time,
                    pygplates.VelocityUnits.cms_per_yr)
            convergence_velocity_vectors = [
                    from_stage_frame_relative_to_overriding * velocity
                            for velocity in convergence_velocity_vectors_in_stage_frame_relative_to_overriding]
            
            # Calculate the absolute velocities at the arc midpoints.
            absolute_velocity_vectors = pygplates.calculate_velocities(
                    arc_midpoints, subduction_zone_equivalent_stage_rotation,
                    velocity_delta_time, pygplates.VelocityUnits.cms_per_yr)
            
            for arc_index in range(len(arc_midpoints)):
                arc_midpoint = arc_midpoints[arc_index]
                arc_length = arc_lengths[arc_index]
                subducting_arc_normal = subducting_arc_normals[arc_index]
                subducting_arc_normal_azimuth = subducting_arc_local_normals[arc_index][1]
                lat, lon = arc_midpoint.to_lat_lon()
                
                # The direction towards which we rotate from the subducting normal in a clockwise fashion.
                clockwise_direction = pygplates.Vector3D.cross(subducting_arc_normal, arc_midpoint.to_xyz())
                
                # Calculate the convergence rate parameters.
                convergence_velocity_vector = convergence_velocity_vectors[arc_index]
                if convergence_velocity_vector.is_zero_magnitude():
                    convergence_velocity_magnitude = 0
                    convergence_obliquity_degrees = 0
                else:
                    convergence_velocity_magnitude = convergence_velocity_vector.get_magnitude()
                    convergence_obliquity_degrees = math.degrees(pygplates.Vector3D.angle_between(
                            convergence_velocity_vector, subducting_arc_normal))
                    # Anti-clockwise direction has range (0, -180) instead of (0, 180).
                    if pygplates.Vector3D.dot(convergence_velocity_vector, clockwise_direction) < 0:
                        convergence_obliquity_degrees = -convergence_obliquity_degrees
                    
                    # See if plates are diverging (moving away from each other).
                    # If plates are diverging (moving away from each other) then make the
                    # velocity magnitude negative to indicate this. This could be inferred from
                    # the obliquity but it seems this is the standard way to output convergence rate.
                    if math.fabs(convergence_obliquity_degrees) > 90:
                        convergence_velocity_magnitude = -convergence_velocity_magnitude
                
                # Calculate the absolute rate parameters.
                absolute_velocity_vector = absolute_velocity_vectors[arc_index]
                if absolute_velocity_vector.is_zero_magnitude():
                    absolute_velocity_magnitude = 0
                    absolute_obliquity_degrees = 0
                else:
                    absolute_velocity_magnitude = absolute_velocity_vector.get_magnitude()
                    absolute_obliquity_degrees = math.degrees(pygplates.Vector3D.angle_between(
                            absolute_velocity_vector, subducting_arc_normal))
                    # Anti-clockwise direction has range (0, -180) instead of (0, 180).
                    if pygplates.Vector3D.dot(absolute_velocity_vector, clockwise_direction) < 0:
                        absolute_obliquity_degrees = -absolute_obliquity_degrees
                    
                    # See if the subduction zone absolute motion is heading in the direction of the
                    # overriding plate. If it is then make the velocity magnitude negative to
                    # indicate this. This could be inferred from the obliquity but it seems this
                    # is the standard way to output convergence rate.
                    #
                    # Note that we are not calculating the motion of the subduction zone
                    # relative to the overriding plate - they are usually attached to each other
                    # and hence wouldn't move relative to each other.
                    if math.fabs(absolute_obliquity_degrees) < 90:
                        absolute_velocity_magnitude = -absolute_velocity_magnitude
                
                # The data will be output in GMT format (ie, lon first, then lat, etc).
                output_data.append((
                        lon,
                        lat,
                        convergence_velocity_magnitude,
                        convergence_obliquity_degrees,
                        absolute_velocity_magnitude,
                        absolute_obliquity_degrees,
                        math.degrees(arc_length),
                        math.degrees(subducting_arc_normal_azimuth),
                        subducting_plate_id,
                        overriding_plate_id))
    
    # Return data sorted since it's easier to compare results (when at least lon/lat is sorted).
    return sorted(output_data)


def write_output_file(output_filename, output_data):
    with open(output_filename, 'w') as output_file:
        for output_line in output_data:
            output_file.write(' '.join(str(item) for item in output_line) + '\n')


def subduction_convergence_over_time(
        output_filename_prefix,
        output_filename_extension,
        rotation_filenames,
        topology_filenames,
        threshold_sampling_distance_radians,
        time_young,
        time_old,
        time_increment,
        velocity_delta_time = 1.0,
        anchor_plate_id = 0):
    
    if time_increment <= 0:
        print('The time increment "{0}" is not positive and non-zero.'.format(time_increment), file=sys.stderr)
        return
    
    if time_young > time_old:
        print('The young time {0} is older (larger) than the old time {1}.'.format(time_young, time_old), file=sys.stderr)
        return
    
    rotation_model = pygplates.RotationModel(rotation_filenames)
    
    # Read/parse the topological features once so we're not doing at each time iteration.
    topology_features = [pygplates.FeatureCollection(topology_filename)
            for topology_filename in topology_filenames]
    
    # Iterate over the time rage.
    time = time_young
    while time <= pygplates.GeoTimeInstant(time_old):
        
        print('Time {0}'.format(time))
        
        # Returns a list of tesselated subduction zone points and associated convergence parameters
        # to write to the output file for the current 'time'.
        output_data = subduction_convergence(
                rotation_model,
                topology_features,
                threshold_sampling_distance_radians,
                time,
                velocity_delta_time,
                anchor_plate_id)
        
        if output_data:
            output_filename = '{0}_{1:0.2f}.{2}'.format(output_filename_prefix, time, output_filename_extension)
            write_output_file(output_filename, output_data)

        # Increment the time further into the past.
        time += time_increment
    
    return 0 # Success


if __name__ == '__main__':
    
    # Check the imported pygplates version.
    if not hasattr(pygplates, 'Version') or pygplates.Version.get_imported_version() < PYGPLATES_VERSION_REQUIRED:
        print('{0}: Error - imported pygplates version {1} but version {2} or greater is required'.format(
                os.path.basename(__file__), pygplates.Version.get_imported_version(), PYGPLATES_VERSION_REQUIRED),
            file=sys.stderr)
        sys.exit(1)
    
    
    __description__ = \
    """Find the convergence rates along subduction zones over time.
    
    For each time (over a range of times) an output xy file is generated containing the resolved subduction zones
    (with point locations as the first two columns x and y) and the following convergence rate parameters in subsequent columns:
    
      - subducting convergence (relative to overriding plate) velocity magnitude (in cm/yr)
      - subducting convergence velocity obliquity angle (angle between subduction zone normal vector and convergence velocity vector)
      - subduction zone absolute (relative to anchor plate) velocity magnitude (in cm/yr)
      - subduction zone absolute velocity obliquity angle (angle between subduction zone normal vector and absolute velocity vector)
      - length of arc segment (in degrees) that current point is on
      - subducting arc normal azimuth angle (clockwise starting at North, ie, 0 to 360 degrees) at current point
      - subducting plate ID
      - overriding plate ID
    
    The obliquity angles are in the range (-180 180). The range (0, 180) goes clockwise (when viewed from above the Earth) from the
    subducting normal direction to the velocity vector. The range (0, -180) goes counter-clockwise.
    You can change the range (-180, 180) to the range (0, 360) by adding 360 to negative angles.
    
    Note that the convergence velocity magnitude is negative if the plates are diverging (if convergence obliquity angle
    is greater than 90 or less than -90). And note that the absolute velocity magnitude is negative if the subduction zone (trench)
    is moving towards the overriding plate (if absolute obliquity angle is less than 90 or greater than -90) - note that this
    ignores the kinematics of the subducting plate.
    
    Each point in the output is the midpoint of a great circle arc between two adjacent points in the subduction zone polyline.
    The subduction zone normal vector used in the obliquity calculations is perpendicular to the great circle arc of each point (arc midpoint)
    and pointing towards the overriding plate (rather than away from it).
    
    Each subduction zone is sampled at approximately uniform intervals along its length (specified via a threshold sampling distance).
    The sampling along the entire length of a subduction zone is not exactly uniform. Each segment along a subduction zone is sampled
    such that the samples have a uniform spacing that is less than or equal to the threshold sampling distance. However each segment
    in a subduction zone might have a slightly different spacing distance (since segment lengths are not integer multiples of
    the threshold sampling distance).
    
    The subducting arc normal (at each arc segment mid-point) always points *towards* the overriding plate.

    NOTE: Separate the positional and optional arguments with '--' (workaround for bug in argparse module).
    For example...

    python %(prog)s -r rotations.rot -m topologies.gpml -t 0 200 -i 1 -v 1 -d 0.5 -e xy -- convergence
     """

    # The command-line parser.
    parser = argparse.ArgumentParser(description = __description__, formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-r', '--rotation_filenames', type=str, nargs='+', required=True,
            metavar='rotation_filename', help='One or more rotation files.')
    parser.add_argument('-m', '--topology_filenames', type=str, nargs='+', required=True,
            metavar='topology_filename', help='One or more topology files to generate resolved subducting lines.')
    parser.add_argument('-a', '--anchor', type=int, default=0,
            dest='anchor_plate_id',
            help='Anchor plate id used for reconstructing. Defaults to zero.')
    
    # Can specify only one of '-i', '-l' or '-t'.
    threshold_sampling_distance_group = parser.add_mutually_exclusive_group()
    threshold_sampling_distance_group.add_argument('-d', '--threshold_sampling_distance_degrees', type=float,
            help='Threshold sampling distance along subduction zones (in degrees). '
                'Defaults to {0} degrees.'.format(DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES))
    threshold_sampling_distance_group.add_argument('-k', '--threshold_sampling_distance_kms', type=float,
            help='Threshold sampling distance along subduction zones (in Kms). '
                'Defaults to {0:.2f} Kms (which is equivalent to {1} degrees).'.format(
                        DEFAULT_THRESHOLD_SAMPLING_DISTANCE_KMS,
                        DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES))

    parser.add_argument('-t', '--time_range', type=float, nargs=2,
            metavar=('young_time', 'old_time'),
            default=[DEFAULT_TIME_RANGE_YOUNG_TIME, DEFAULT_TIME_RANGE_OLD_TIME],
            help='The time range (in Ma) from young time to old time. '
                'Defaults to {0} -> {1} Ma.'.format(
                    DEFAULT_TIME_RANGE_YOUNG_TIME, DEFAULT_TIME_RANGE_OLD_TIME))
    
    def parse_positive_number(value_string):
        try:
            value = float(value_string)
        except ValueError:
            raise argparse.ArgumentTypeError("%s is not a number" % value_string)
        
        if value <= 0:
            raise argparse.ArgumentTypeError("%g is not a positive number" % value)
        
        return value
    
    parser.add_argument('-i', '--time_increment', type=parse_positive_number,
            default=DEFAULT_TIME_INCREMENT,
            help='The time increment in My. Defaults to {0} My.'.format(DEFAULT_TIME_INCREMENT))
    
    parser.add_argument('-v', '--velocity_delta_time', type=parse_positive_number,
            default=DEFAULT_VELOCITY_DELTA_TIME,
            help='The delta time interval used to calculate velocities in My. '
                'Defaults to {0} My.'.format(DEFAULT_VELOCITY_DELTA_TIME))
    
    parser.add_argument('output_filename_prefix', type=str,
            help='The output filename prefix. An output file is created for each geological time in the sequence where '
                'the filename suffix contains the time and the filename extension.')
    parser.add_argument('-e', '--output_filename_extension', type=str, default='xy',
            help='The output xy filename extension. Defaults to "xy".')
    
    # Parse command-line options.
    args = parser.parse_args()
    
    if args.time_range[0] > args.time_range[1]:
        raise argparse.ArgumentTypeError("First (young) value in time range is greater than second (old) value")
    
    # Determine threshold sampling distance.
    if args.threshold_sampling_distance_degrees:
        threshold_sampling_distance_radians = math.radians(args.threshold_sampling_distance_degrees)
    elif args.threshold_sampling_distance_kms:
        threshold_sampling_distance_radians = threshold_sampling_distance_kms / pygplates.Earth.equatorial_radius_in_kms
    else: # default...
        threshold_sampling_distance_radians = math.radians(DEFAULT_THRESHOLD_SAMPLING_DISTANCE_DEGREES)
    
    return_code = subduction_convergence_over_time(
            args.output_filename_prefix,
            args.output_filename_extension,
            args.rotation_filenames,
            args.topology_filenames,
            threshold_sampling_distance_radians,
            args.time_range[0],
            args.time_range[1],
            args.time_increment,
            args.velocity_delta_time,
            args.anchor_plate_id)
    if return_code is None:
        sys.exit(1)
        
    sys.exit(0)
