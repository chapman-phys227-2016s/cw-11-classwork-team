#! /usr/bin/env python

"""
File: sombrero.py
Copyright (c) 2016 Austin Ayers
License: MIT

Course: PHYS227
Assignment: Plots sombrero and motion of ball in said sombrero.
Date: April 26, 2016
Email: ayers111@mail.chapman.edu
Name: Austin Ayers
Description: Numerically solves Newton's equations for a ball being shaken in a sombrero by sinusoidal driving force.
The pinnacle of groupwork, efficiency,
"""
import math

class Sombrero():
    """
    Class that contains the equations and functions to solve said equations
    Implementation for the 'rk4' function comes from https://rosettacode.org/wiki/Runge-Kutta_method, they're awesome.
    """
    def __init__(self, F, x_0, y_0, w = 1.0, delta = 0.25, m = 1, h = 0.001):
        self.delta = delta
        self.x_0 = x_0
        self.y_0 = y_0
        self.F = F
        self.w = w
        self.m = m
        self.h = h

        rk4_output = rk4(equation_1, equation_2, self.x_0, self.y_0, 10000)
    def equation_1(self, x, y, t):
        """
        x_dot(t) = ...
        """
        return y
    def equation_2(self, x, y, t):
        """
        y_dot(t) = ...
        """
        return ( -1 * self.delta * y + x - x ** 3 + self.F * cos(self.w * t) ) / self.m
    def rk4(self, f1, f2, x0, y0, n, h = self.h):
        """
        Runge Kutta 4 function

        Arguments: f1 (first coupled function being differentiated)
                   f2 (second coupled function being differentiated)
                   x0 (initial x point)
                   y0 (initial y point)
                   n (number of points to be calculated in rk4)
                   h (delta of x)

        Outputs: vx_1 (array of x points, subdivided n times)
                 vy_1 (array of y points, answer array which is output of rk4 - tied to vx_1)
                 vx_2 (array of x points, subdivided n times)
                 vy_2 (array of y points, answer array which is output of rk4 - tied to vx_1)
        """
        vx = [0]*(n + 1)
        vy = [0]*(n + 1)
        vx[0] = x = x0
        vy[0] = y = y0
        t = 0
        for i in range(1, n + 1):
            k1_x = h*f1(x, y, t)
            k1_y = h*f2(x, y, t)

            k2_x = h*f1(x + 0.5*k1_x, y + 0.5*k1_y, t + h / 2.0)
            k2_y = h*f2(x + 0.5*k1_x, y + 0.5*k1_y, t + h / 2.0)

            k3_x = h*f1(x + 0.5*k2_x, y + 0.5*k2_y, t + h / 2.0)
            k3_y = h*f2(x + 0.5*k2_x, y + 0.5*k2_y, t + h / 2.0)

            k4_x = h*f1(x + k3_x, y + k3_y, t + h)
            k4_y = h*f2(x + k3_x, y + k3_y, t + h)


            vx[i] = x = x + (k1_x + k2_x + k2_x + k3_x + k3_x + k4_x)/6
            vy[i] = y = y + (k1_y + k2_y + k2_y + k3_y + k3_y + k4_y)/6

            t += h

        return vx_1, vy_1, vx_2, vy_2

    #vx, vy = rk4(f, 0, 1, 10, 100)
    #for x, y in list(zip(vx, vy))[::10]:
    #    print(x, y, y - (4 + x*x)**2/16)

    def test_linear():
        thisone = Sombrero(0)
        def line1(x):
            return 2*x
        def line2(x):
            return 5*x
        x1, y1, x2, y2 = thisone.rk4(line1, line2, 0, 0, 1000)
        apt = math.fabs(y1[5] == 2 and y2[5] == 5)
        msg = 'Linear test failed.'
        assert apt, msg 
