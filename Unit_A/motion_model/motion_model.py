import matplotlib.pyplot as plt
import numpy as np
from math import sin, cos, pi, fabs

TICKS_TO_DIST = 0.349
WIDTH = 150.0


def main():
    pose = np.matrix(([1850.0], [1897.0], [3.717551]))
    poses = [pose]
    motion_file = open("motion.txt", "r")
    motions = motion_file.readlines()
    for motion in motions:
        motion = motion.split("\t")
        ticks = np.matrix(([int(motion[1])], [int(motion[2])]))
        new_pose = next_pose(pose, ticks)
        pose = new_pose
        poses.append(new_pose)
    motion_file.close()
    output_file = open("poses.out", "w")
    for p in poses:
        output_file.write(str(p[0][0]) + "\t" + str(p[1][0]) + "\t" + str(p[2][0]) + "\n")
    output_file.close()
    plot(poses)
    

def next_pose(pose, ticks):
    """Given the robot's current pose and
    the ticks of the left wheel and right wheel,
    return the next pose.
    :param pose: The current pose of the robot, a
        numpy column vector composed of x, y, heading
    :param ticks: A numpy column vector composed of
        left ticks, right ticks
    :return: The new pose of the robot
    """
    thresh = 0
    new_pose = np.zeros((3, 1), float)
    if fabs(ticks[0][0] - ticks[1][0]) <= thresh:
        new_pose[2][0] = pose[2][0]
        new_pose[0][0] = pose[0][0] + ticks[0][0]*TICKS_TO_DIST*cos(pose[2][0])
        new_pose[1][0] = pose[1][0] + ticks[0][0]*TICKS_TO_DIST*sin(pose[2][0])
    else:
        alpha = TICKS_TO_DIST * (float(ticks[1][0]-ticks[0][0]))/float(WIDTH)
        new_pose[2][0] = (pose[2][0] + alpha + 2*pi) % (2 * pi)
        R = TICKS_TO_DIST * ticks[0][0]/float(alpha)
        dist = float(R) + float(WIDTH)/2
        c = pose[0:2][:] - dist * np.matrix(([sin(pose[2][0])], [-cos(pose[2][0])]), float)
        new_pose[0:2][:] = c + dist*np.matrix(([sin(new_pose[2][0])], [-cos(new_pose[2][0])]), float)
    return new_pose
    

def plot(poses):
    """Plot all the poses as points
    :param poses: List of pose matrices
    :return: None
    """
    x = []
    y = []
    for pose in poses:
        x.append(pose[0][0])
        y.append(pose[1][0])
    plt.plot(x, y, marker="o", ls="-")
    plt.show()
    

if __name__ == "__main__":
    main()

