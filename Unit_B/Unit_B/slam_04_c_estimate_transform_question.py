# For each cylinder in the scan, find its cartesian coordinates,
# in the world coordinate system.
# Find the closest pairs of cylinders from the scanner and cylinders
# from the reference, and the optimal transformation which aligns them.
# Then, output the scanned cylinders, using this transform.
# 04_c_estimate_transform
# Claus Brenner, 14 NOV 2012
from lego_robot import *
from slam_b_library import filter_step
from slam_04_a_project_landmarks import\
     compute_scanner_cylinders, write_cylinders
from math import sqrt

# Given a list of cylinders (points) and reference_cylinders:
# For every cylinder, find the closest reference_cylinder and add
# the index pair (i, j), where i is the index of the cylinder, and
# j is the index of the reference_cylinder, to the result list.
# This is the function developed in slam_04_b_find_cylinder_pairs.
def find_cylinder_pairs(cylinders, reference_cylinders, max_radius):
    cylinder_pairs = []

    # --->>> Enter your code here.
    # Make a loop over all cylinders and reference_cylinders.
    # In the loop, if cylinders[i] is closest to reference_cylinders[j],
    # and their distance is below max_radius, then add the
    # tuple (i,j) to cylinder_pairs, i.e., cylinder_pairs.append( (i,j) ).
    for r in range(len(reference_cylinders)):
        closest_dist = 1000000000000
        closest_cylinder_ind = -1
        for i in range(len(cylinders)):
            dist = ( (cylinders[i][0]-reference_cylinders[r][0])**2 + (cylinders[i][1]-reference_cylinders[r][1])**2 )**0.5
            if dist < closest_dist and dist <= max_radius:
                closest_dist = dist
                closest_cylinder_ind = i
        if closest_cylinder_ind >= 0:
            cylinder_pairs.append((closest_cylinder_ind, r))
    return cylinder_pairs


# Given a point list, return the center of mass.
def compute_center(point_list):
    # Safeguard against empty list.
    if not point_list:
        return (0.0, 0.0)
    # If not empty, sum up and divide.
    sx = sum([p[0] for p in point_list])
    sy = sum([p[1] for p in point_list])
    return (float(sx) / len(point_list), float(sy) / len(point_list))

# Given a left_list of points and a right_list of points, compute
# the parameters of a similarity transform: scale, rotation, translation.
# If fix_scale is True, use the fixed scale of 1.0.
# The returned value is a tuple of:
# (scale, cos(angle), sin(angle), x_translation, y_translation)
# i.e., the rotation angle is not given in radians, but rather in terms
# of the cosine and sine.
def estimate_transform(left_list, right_list, fix_scale = False):
    # If we have 0 or 1 reference points, it's impossible to find
    # a unique transform (only a translation is necessary)
    if len(left_list) <= 1 or len(right_list) <= 1:
        return None   
    
    # Compute left and right center.
    lc = compute_center(left_list)
    rc = compute_center(right_list)

    # --->>> Insert here your code to compute lambda, c, s and tx, ty.
    # Compute reduced coordinates
    reduced_r = []
    reduced_l = []
    for r in right_list:
        reduced_r.append([r[0]-rc[0], r[1]-rc[1]])
    for l in left_list:
        reduced_l.append([l[0]-lc[0], l[1]-lc[1]])
        
    # Compute lambda
    la = 1.0
    if not fix_scale:
        squared_reduced_right_sum = 0.0
        squared_reduced_left_sum = 0.0
        for rr in reduced_r:
            squared_reduced_right_sum += rr[0]**2 + rr[1]**2
        for lr in reduced_l:
            squared_reduced_left_sum += lr[0]**2 + lr[1]**2
        la = (squared_reduced_right_sum/squared_reduced_left_sum)**0.5
    
    # Compute c and s
    cosine_sum = 0.0
    sine_sum = 0.0
    for i in range(len(reduced_r)):
        cosine_sum += reduced_r[i][0]*reduced_l[i][0] + reduced_r[i][1]*reduced_l[i][1]
        sine_sum += -reduced_r[i][0]*reduced_l[i][1] + reduced_r[i][1]*reduced_l[i][0]
    magnitude = ( cosine_sum**2 + sine_sum**2 )**0.5
    c = cosine_sum/magnitude
    s = sine_sum/magnitude
    
    # Compute tx and ty
    t = list(rc[:])
    correction = [c*lc[0] - s*lc[1], s*lc[0] + c*lc[1]]
    correction[0] *= la
    correction[1] *= la
    t[0] -= correction[0]
    t[1] -= correction[1]
    tx = t[0]
    ty = t[1]
    
    return la, c, s, tx, ty

# Given a similarity transformation:
# trafo = (scale, cos(angle), sin(angle), x_translation, y_translation)
# and a point p = (x, y), return the transformed point.
def apply_transform(trafo, p):
    la, c, s, tx, ty = trafo
    lac = la * c
    las = la * s
    x = lac * p[0] - las * p[1] + tx
    y = las * p[0] + lac * p[1] + ty
    return (x, y)


if __name__ == '__main__':
    # The constants we used for the filter_step.
    scanner_displacement = 30.0
    ticks_to_mm = 0.349
    robot_width = 150.0

    # The constants we used for the cylinder detection in our scan.    
    minimum_valid_distance = 20.0
    depth_jump = 100.0
    cylinder_offset = 90.0

    # The maximum distance allowed for cylinder assignment.
    max_cylinder_distance = 300.0

    # The start pose we obtained miraculously.
    pose = (1850.0, 1897.0, 3.717551306747922)

    # Read the logfile which contains all scans.
    logfile = LegoLogfile()
    logfile.read("robot4_motors.txt")
    logfile.read("robot4_scan.txt")

    # Also read the reference cylinders (this is our map).
    logfile.read("robot_arena_landmarks.txt")
    reference_cylinders = [l[1:3] for l in logfile.landmarks]

    out_file = file("estimate_transform.txt", "w")
    for i in xrange(len(logfile.scan_data)):
        # Compute the new pose.
        pose = filter_step(pose, logfile.motor_ticks[i],
                           ticks_to_mm, robot_width,
                           scanner_displacement)

        # Extract cylinders, also convert them to world coordinates.
        cartesian_cylinders = compute_scanner_cylinders(
            logfile.scan_data[i],
            depth_jump, minimum_valid_distance, cylinder_offset)
        world_cylinders = [LegoLogfile.scanner_to_world(pose, c)
                           for c in cartesian_cylinders]

        # For every cylinder, find the closest reference cylinder.
        cylinder_pairs = find_cylinder_pairs(
            world_cylinders, reference_cylinders, max_cylinder_distance)
        
        # Estimate a transformation using the cylinder pairs.
        trafo = estimate_transform(
            [world_cylinders[pair[0]] for pair in cylinder_pairs],
            [reference_cylinders[pair[1]] for pair in cylinder_pairs],
            fix_scale = True)

        # Transform the cylinders using the estimated transform.
        transformed_world_cylinders = []
        if trafo:
            transformed_world_cylinders =\
                [apply_transform(trafo, c) for c in
                 [world_cylinders[pair[0]] for pair in cylinder_pairs]]            

        # Write to file.
        # The pose.
        print >> out_file, "F %f %f %f" % pose
        # The detected cylinders in the scanner's coordinate system.
        write_cylinders(out_file, "D C", cartesian_cylinders)
        # The detected cylinders, transformed using the estimated trafo.
        write_cylinders(out_file, "W C", transformed_world_cylinders)

    out_file.close()
