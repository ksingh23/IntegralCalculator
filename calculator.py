import Equation as eq
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon
from matplotlib.widgets import AxesWidget
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr
import time
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
            length = sp.N(func.subs(x, start))
            area += (length * rect_width)
            start += rect_width
        if abs(area - old_area) <= 0.001:
            change = False
        else:
            change = True
        num_rect *= 2
        rect_width = (right_bound - left_bound) / num_rect
    return sp.N(area)


def integral_with_simpsons_rule(func, left_bound, right_bound, num_steps):
    step_size = float((right_bound - left_bound)/num_steps)
    def_integral = func.subs(x, (left_bound + sys.float_info.epsilon))
    # Problem is that ln(sin(x)) at 0 is undefined, how do we solve this?
    # Instead of inputting zero, I inputted 0 + epsilon to solve
    i = left_bound + step_size
    itr = 0
    while i < (right_bound - step_size):
        if itr % 2 == 0:
            def_integral += (4 * func.subs(x, i))
        else:
            def_integral += (2 * func.subs(x, i))
        itr += 1
        i += step_size
    # returns zoo = complex infinity
    def_integral += (func.subs(x, i))
    def_integral *= (step_size/3)
    return sp.N(def_integral)


def get_step_size(function, left_bound, right_bound, error):
    x = sp.symbols('x')
    second_deriv = sp.diff(function, x, 2)
    third_deriv = sp.diff(second_deriv, x, 1)
    zeros = sp.solve(third_deriv, x)
    lb = second_deriv.subs(x, left_bound)
    rb = second_deriv.subs(x, right_bound)
    if math.isinf(abs(lb)):
        if math.isinf(abs(rb)):
            m = -sys.maxsize - 1
        else:
            m = abs(rb)
    else:
        if math.isinf(abs(rb)):
            m = abs(lb)
        else:
            m = max(abs(lb), abs(rb))
    for z in zeros:
        if z >= left_bound and z <= right_bound:
            value = abs(second_deriv.subs(x, z))
            if value > m:
                m = value
    num_steps = math.ceil(math.sqrt((m * (right_bound - left_bound)**3)/(12 * error)))
    if num_steps % 2 != 0:
        num_steps += 1
    return num_steps


class DragAxes:
    def __init__(self, line):
        self.line = line
        print(self.line.axes)
        self.x = list(line.get_xdata())
        self.y = list(line.get_ydata())
        # Above instance variables will give us the range of x and y that the lines cover
        self.press = None

    def connect(self):
        self.cidpress = self.line.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cidmotion = self.line.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.cidrelease = self.line.figure.canvas.mpl_connect('button_release_event', self.on_release)

    def on_press(self, event):
        if event.inaxes != self.line.axes:
            return
        contains, attrd = self.line.contains(event)
        if not contains:
            return
        print('event contains', self.x, self.y)
        x0, y0 = self.x, self.y
        self.press = x0, y0, event.xdata, event.ydata

    def on_motion(self, event):
        if not self.press:
            return
        if event.inaxes != self.line.axes:
            return
        x0, y0, xpress, ypress = self.press
        dx = event.xdata - xpress
        new_x = x0 + dx
        self.line.set_xdata([new_x])
        self.x = new_x
        self.line.figure.canvas.draw()

    def on_release(self, event):
        self.press = None
        self.line.figure.canvas.draw()


if __name__ == '__main__':
    formula = "cos(x)"
    x = sp.symbols('x')
    sp.init_printing(use_unicode=True)
    func = eq.Expression(formula, ["x"])
    function = parse_expr(formula)
    left_bound = 0
    right_bound = sp.pi/2
    start = time.time()
    num_steps = get_step_size(function, left_bound, right_bound, 0.0001)
    print(integral_with_simpsons_rule(function, left_bound, right_bound, num_steps))
    print(time.time() - start)
    figure, axis = plt.subplots()
    x = np.linspace(-np.pi, np.pi, 201)
    axis.plot(x, func(x), 'r')
    axis.set_ylim(bottom=0)
    bounds = np.linspace(0, np.pi / 2, 201)
    integral_region = func(bounds)
    verts = [(left_bound, 0), *zip(bounds, integral_region), (right_bound, 0)]
    area_under_curve = Polygon(verts, facecolor='0.9', edgecolor='0.2', linewidth=0)
    axis.add_patch(area_under_curve)
    lb = axis.axvline(x=left_bound)
    left = DragAxes(lb)
    left.connect()
    rb = axis.axvline(x=right_bound)
    right = DragAxes(rb)
    right.connect()
    plt.show()

