# Histogram implementation of a bayes filter - combines
# convolution and multiplication of distributions, for the
# movement and measurement steps.
# 06_d_histogram_filter
# Claus Brenner, 28 NOV 2012
from pylab import plot, show, ylim
from distribution import *

def move(distribution, delta):
    """Returns a Distribution that has been moved (x-axis) by the amount of
       delta."""
    return Distribution(distribution.offset + delta, distribution.values)



# --->>> Copy your convolve(a, b) and multiply(a, b) functions here.
def convolve(a, b):
    """Convolve distribution a and b and return the resulting new distribution."""

    # --->>> Put your code here.
    distributions = []
    for i in range(len(b.values)):
        distributions.append(Distribution(a.offset + b.offset+i, [v * b.values[i] for v in a.values]))
    return Distribution.sum(distributions)  # Replace this by your own result.
    
    
def multiply(a, b):
    """Multiply two distributions and return the resulting distribution."""

    # --->>> Put your code here.
    values = []
    if a.start() < b.start():
        offset = a.start()
        if a.stop() > b.stop():
            stop = a.stop()
        else:
            stop = b.stop()
        for i in range(offset, stop):
            if i >= a.start() and i < a.stop() and i >= b.start() and i < b.stop():
                values.append(a.values[i-a.start()] * b.values[i-b.start()])
            else:
                values.append(0.0)
    else:
        offset = b.start()
        if a.stop() > b.stop():
            stop = a.stop()
        else:
            stop = b.stop()
        for i in range(offset, stop):
            if i >= a.start() and i < a.stop() and i >= b.start() and i < b.stop():
                values.append(a.values[i-a.start()] * b.values[i-b.start()])
            else:
                values.append(0.0)
    d = Distribution(offset, values)
    d.normalize()
    return d  # Modify this to return your result.


if __name__ == '__main__':
    arena = (0,2200)

    # Start position. Exactly known - a unit pulse.
    start_position = 10
    position = Distribution.unit_pulse(start_position)
    plot(position.plotlists(*arena)[0], position.plotlists(*arena)[1],
         linestyle='steps')

    # Movement data.
    controls  =    [ 20 ] * 100

    # Measurement data. Assume (for now) that the measurement data
    # is correct. - This code just builds a cumulative list of the controls,
    # plus the start position.
    p = start_position
    measurements = []
    for c in controls:
        p += c
        measurements.append(p)

    # This is the filter loop.
    for i in xrange(len(controls)):
        # Move, by convolution. Also termed "prediction".
        control = Distribution.triangle(controls[i], 10)
        position = convolve(position, control)
        plot(position.plotlists(*arena)[0], position.plotlists(*arena)[1],
             color='b', linestyle='steps')

        # Measure, by multiplication. Also termed "correction".
        measurement = Distribution.triangle(measurements[i], 10)
        position = multiply(position, measurement)
        plot(position.plotlists(*arena)[0], position.plotlists(*arena)[1],
             color='r', linestyle='steps')

    show()
