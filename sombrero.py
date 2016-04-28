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
Description: Numerically solves Newton's equations for a ball being shaken in a sombrero by constant driving force.
"""
from math import sqrt
from math import cos

class Sombrero():
    """
    Class that contains the equations and functions to solve said equations
    Implementation for the 'rk4' function comes from https://rosettacode.org/wiki/Runge-Kutta_method, they're awesome.
    """
    def __init__(self, F, w = 1.0, delta = 0.25, m = 1, h = 0.001):
        self.delta = delta
        self.F = F
        self.w = w
        self.m = m
        self.h = h
    def equation_1(self):
        """
        x_dot(t) = ...
        """
        return y
    def equation_2(self, yt, xt, t):
        """
        y_dot(t) = ...
        """
        return ( -1 * self.delta * yt + xt - xt ** 3 + self.F * cos(self.w * t) ) / self.m
    def rk4(self, f1, f2, x0, y0, n, h = self.h):
        """
        Runge Kutta 4 function

        Arguments: f1 (first coupled function being differentiated)
                   f2 (second coupled function being differentiated)
                   x0 (initial x point)
                   y0 (initial y point)
                   n (number of points to be calculated in rk4)
                   h (delta of x)

        Outputs: vx (array of x points, subdivided n times)
                 vy (array of y points, answer array which is output of rk4)
        """
        vx = [0]*(n + 1)
        vy = [0]*(n + 1)
        vx[0] = x = x0
        vy[0] = y = y0
        for i in range(1, n + 1):
            k1 = h*f(x, y)
            k2 = h*f(x + 0.5*h, y + 0.5*k1)
            k3 = h*f(x + 0.5*h, y + 0.5*k2)
            k4 = h*f(x + h, y + k3)
            vx[i] = x = x0 + i*h
            vy[i] = y = y + (k1 + k2 + k2 + k3 + k3 + k4)/6
        return vx, vy

    #vx, vy = rk4(f, 0, 1, 10, 100)
    #for x, y in list(zip(vx, vy))[::10]:
    #    print(x, y, y - (4 + x*x)**2/16)
