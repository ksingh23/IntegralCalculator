from math import cos
import Equation as eq
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon


def main():
    formula = "cos(x)"
    func = eq.Expression("cos(x)", ["x"])
    left_bound = 0
    right_bound = np.pi/2
    numRect = 100
    rectWidth = (right_bound - left_bound)/numRect
    change = True
    area = 0
    while change:
        start = 0
        oldArea = area
        area = 0
        for i in range(numRect):
            area += (func(start) * rectWidth)
            start += rectWidth

        # print("Area: " + str(area))
        if abs(area - oldArea) <= 0.001:
            change = False
        else:
            change = True
        numRect *= 2
        rectWidth = (right_bound - left_bound) / numRect
    figure, axis = plt.subplots()
    x = np.linspace(-np.pi, np.pi, 201)
    axis.plot(x, func(x), 'r')
    axis.set_ylim(bottom=0)
    bounds = np.linspace(0, np.pi/2, 201)
    integral_region = func(bounds)
    verts = [(left_bound, 0), *zip(bounds, integral_region), (right_bound, 0)]
    area_under_curve = Polygon(verts, facecolor='0.9', edgecolor='0.2', linewidth=0)
    axis.add_patch(area_under_curve)
    plt.show()

main()