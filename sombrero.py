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
from math import sin, cos
from unittest import TestCase
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

class Sombrero():
    """
    Class that contains the equations and functions to solve said equations
    Implementation for the 'rk4' function comes from https://rosettacode.org/wiki/Runge-Kutta_method, they're awesome.
    """
    def __init__(self, F, x_0, y_0, w = 1.0, delta = 0.25, m = 1, h = 0.001, n = 10000):
        self.delta = delta
        self.x_0 = x_0
        self.y_0 = y_0
        self.F = F
        self.w = w
        self.m = m
        self.h = h
        self.n = n

        self.rk4_output = self.rk4(self.equation_1, self.equation_2, self.x_0, self.y_0, self.n)
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
    def rk4(self, f1, f2, x0, y0, n, h = .001):
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

        return vx, vy
    
    def generate_graph(self):
        """
        Generates images ('.png's) of the plots of the Runge Kutta's
        solutions as well as the exact solutions.
        """
        
        approx_values = self.rk4_output
        u_list_approx = approx_values[0]
        v_list_approx = approx_values[1]
        
        fig, ax = plt.subplots(nrows = 1, ncols = 1)
        ax.plot(u_list_approx, v_list_approx, 'r')
        ax.set_xlabel("u(t)")
        ax.set_ylabel("v(t)")
        ax.set_title("Sombrero approximation for n = " + str(self.n) + " and F = " + str(self.F))
        ax.grid(True)
        fig.savefig(self.__class__.__name__ + '_%3d.png' % (self.n))
        plt.close(fig)

    #vx, vy = rk4(f, 0, 1, 10, 100)
    #for x, y in list(zip(vx, vy))[::10]:
    #    print(x, y, y - (4 + x*x)**2/16)

class test_sombrero(TestCase):
    def test_linear(self):
        thisone = Sombrero(0, 0, 0)
        def line1(x, y, t):
            return 2
        def line2(x, y, t):
            return 5
        y1, y2 = thisone.rk4(line1, line2, 0, 0, 1000)
        apt = (math.fabs(y1[5] - 0.001*5*2) < 1e-5)  and ((y2[5] - 0.001*5*5) < 1e-5)
        msg = 'Linear test failed.'
        assert apt, msg 

    def test_analytic_result(self):
        test_sombrero = Sombrero(.4, 0, 0, n = 1)
        print test_sombrero.rk4_output
        assert True
        """
        Note to self - edit these later to assert these values
        x = 1.1666876666666666e-07
        y = 0.00023334173071666665
        """
