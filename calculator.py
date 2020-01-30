from math import cos
import Equation as eq
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr
from sympy.utilities.lambdify import implemented_function, lambdify
import time
from scipy.optimize import fsolve
import math


def integral_with_riemann_sum(func, left_bound, right_bound):
    num_rect = 100
    rect_width = (right_bound - left_bound) / num_rect
    change = True
    area = 0
    while change:
        start = 0
        old_area = area
        area = 0
        for i in range(num_rect):
            area += (func(start) * rect_width)
            start += rect_width

        # print("Area: " + str(area))
        if abs(area - old_area) <= 0.001:
            change = False
        else:
            change = True
        num_rect *= 2
        rect_width = (right_bound - left_bound) / num_rect
    return area


def integral_with_simpsons_rule(func, left_bound, right_bound, deriv):
    error = 0.001
    midpt = (right_bound - left_bound)/2
    return ((right_bound - left_bound)/6) * (func(left_bound) + 4 * func(midpt) + func(right_bound))


def get_step_size(function, left_bound, right_bound, error):
    x = sp.symbols('x')
    second_deriv = sp.diff(function, x, 2)
    third_deriv = sp.diff(second_deriv, x, 1)
    zeros = sp.solve(third_deriv, x)
    zeros += [left_bound, right_bound]
    m = -sys.maxsize - 1
    for z in zeros:
        if z >= left_bound and z <= right_bound:
            value = abs(second_deriv.subs(x, z))
            if value > m:
                m = value
    return math.ceil(math.sqrt((m * (right_bound - left_bound)**3)/(12 * error)))


if __name__ == '__main__':
    formula = "ln(sin(x))"
    x = sp.symbols('x')
    sp.init_printing(use_unicode=True)
    func = eq.Expression(formula, ["x"])
    function = parse_expr(formula)
    lam_f = lambdify(x, function, 'numpy')
    left_bound = sp.pi/6
    right_bound = sp.pi / 2
    start = time.time()
    print(get_step_size(function, left_bound, right_bound, 0.0001))
    print(time.time() - start)
    '''
    
    integral_with_riemann_sum(func, left_bound, right_bound)
    
    integral_with_simpsons_rule(lam_f, left_bound, right_bound, deriv)
    figure, axis = plt.subplots()
    x = np.linspace(-np.pi, np.pi, 201)
    axis.plot(x, func(x), 'r')
    axis.set_ylim(bottom=0)
    bounds = np.linspace(0, np.pi / 2, 201)
    integral_region = func(bounds)
    verts = [(left_bound, 0), *zip(bounds, integral_region), (right_bound, 0)]
    area_under_curve = Polygon(verts, facecolor='0.9', edgecolor='0.2', linewidth=0)
    axis.add_patch(area_under_curve)
    plt.show()

    '''
