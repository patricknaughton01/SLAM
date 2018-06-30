import sys
sys.path.append('/home/patricknaughton01/Documents/SLAM/Unit A/CourseFiles/Unit_A')
import matplotlib.pyplot as plt
import time

from math import sin, cos, pi
from lego_robot import *

MIN_VALID_DIST = 20.0
CYLINDER_OFFSET = 90.0
DEPTH_THRESHOLD = 100.0


def main():
    scans = read_scans("scan.txt")
    plt.figure()
    for i in range(278):
        test_scan = scans[i]
        deriv = compute_derivative(test_scan, MIN_VALID_DIST)
        plt.subplot(211)
        plt.plot(test_scan)
        plt.plot(deriv)
        cylinders = scan_derivative(deriv, test_scan, DEPTH_THRESHOLD, MIN_VALID_DIST)
        for cylinder in cylinders:
            plt.plot(cylinder[0], cylinder[1], marker="o")
        plt.subplot(212)
        axes = plt.gca()
        axes.set_xlim([-1000, 5000])
        axes.set_ylim([-1000, 5000])
        plt.plot(0, 0, color='green', marker='x')
        coordinates = get_cylinders_cartesian(cylinders)
        for coor in coordinates:
            plt.plot(coor[0], coor[1], marker="o")
        plt.show(block=False)
        plt.pause(0.1)
        plt.clf()
    

def get_cylinders_cartesian(cylinders):
    """Return a list of the cartesian coordinates of all the 
    cylinders in the robot's frame.
    :param cylinders: A list of cylinders (represented as an
        ordered pair: (scan_index, scan_depth))
    :return: a list of x,y ordered pairs that each represent
        the location of a cylinder in the robot's frame.
    """
    coordinates = []
    log_file = LegoLogfile()        # Helper class from the lego_robot module
    for cylinder in cylinders:
        angle = log_file.beam_index_to_angle(cylinder[0])
        depth = cylinder[1] - CYLINDER_OFFSET       # Cylinder offset was determined experimentally
        coordinates.append((depth*cos(angle), depth*sin(angle)))
    return coordinates


def read_scans(file_name):
    """Read in a list of scans from the file with the
    specified file name.
    :param file_name: The file to read from
    :return: a list of lidar scans where a scan is 
        represented by a list of numbers (distances)
    """
    try:
        in_file = open(file_name, "r")
        lines = in_file.readlines()
        in_file.close()
        scans = []
        for line in lines:
            str_nums = line.split(" ")[3:]
            scan = []
            for str_scan in str_nums:
                scan.append(int(str_scan))
            scans.append(scan)
        return scans
    except Exception:
        print("File " + file_name + " could not be read")
        sys.exit(1)


def scan_derivative(deriv, distances, min_jump, min_depth):
    """Scan deriv (the derivative of a lidar scan) to find
    the scan index of the center of a cylinder. A cylinder
    is marked by a negative derivative of magnitude larger
    than 100 followed by a positive derivative of magnitude
    larger than 100. Multiple negative or positive spikes
    in a row are ignored (they occur when two cylinders
    overlap). Only the most recent spike will be counted.
    :param deriv: A list that represents the derivative of
        the scan function
    :param distances: A list of the original distances
        from which the derivative came
    :param min_jump: The minimum jump in the derivative
        necessary for us to consider this the edge of
        a cylinder
    :param min_depth: Minimum depth that the scanner can
        read. Filters out bad (invalid) scans.
    :return: A list of lists of scan indicies and average
        distances
    """
    on_cylinder = False
    scan_points = []
    scan_depths = []
    cylinders = []
    for i in range(len(deriv)):
        if deriv[i] <= -min_jump:
            on_cylinder = True
            scan_points = []
            scan_depths = []
        if on_cylinder and distances[i] > min_depth:
            scan_points.append(i)
            scan_depths.append(distances[i])
        
        if deriv[i] >= min_jump:
            if on_cylinder:
                cylinders.append([int(average(scan_points)), average(scan_depths)])
            on_cylinder = False
    return cylinders
    
    
def average(arr):
    """Return the average of all the values in arr
    :param arr: an array of values (ints or floats) to 
        take the average of
    :return: float representing the average of all
        the values in arr
    """
    avg = 0.0
    n = len(arr)
    for num in arr:
        avg += num/float(n)
    return avg


def compute_derivative(distances, min_valid):
    """Compute the derivative of the distances found
    by a lidar scan as a function of the scan index.
    :param distances: a list of distances at different
        points in the scan
    :param min_valid: The minimum value a distance can 
        be and still be considered valid (filters
        out some noise)
    :return: a list the same length as distances that 
        represents the derivative of the distances function
    """
    deriv = [0.0]
    for i in range(1, len(distances)-1):
        left = distances[i+1]
        right = distances[i-1]
        # If this is a valid scan
        if left > min_valid and right > min_valid:
            deriv.append((left-right)/2.0)
        else:
            deriv.append(0)
    deriv.append(0.0)
    return deriv
    

if __name__ == "__main__":
    main()

