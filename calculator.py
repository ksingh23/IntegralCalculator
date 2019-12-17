import parser
import numpy as np
from math import cos
import Equation as eq
import sys


def main():
    formula = "cos(x)"
    func = eq.Expression("cos(x)", ["x"])
    leftBound = 0
    rightBound = np.pi/2
    numRect = 100
    rectWidth = (rightBound - leftBound)/numRect
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
        if abs(area - oldArea) <= 0.00001:
            change = False
        else:
            change = True
        numRect *= 2
        rectWidth = (rightBound - leftBound) / numRect
    print(area)


main()