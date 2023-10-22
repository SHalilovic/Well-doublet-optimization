# -*- coding: utf-8 -*-
"""
Definitions of supporting functions

@author: Smajil Halilovic
"""

import numpy as np

def lin_equ(P1, P2):
    """ Line through two points """
    """Line encoded as l=m*x+c."""
    m = (P2[1] - P1[1]) / (P2[0] - P1[0])
    c = (P2[1] - (m * P2[0]))

    # Example Usage:
    # m, c = lin_equ((-40, 30,), (20, 45))
    return m, c

def distance_lines(m, c1, c2):
    """ Compute Euclidean distance between two parallel lines - 2D space """
    # Line 1: y=m*x+c1
    # Line 2: y=m*x+c2
    d = np.abs(c2-c1)/np.sqrt(1 + m**2)
    return d

def rotation(x, y, theta):
    """ Rotate the coordinate system """
    x_r = x*np.cos(theta) + y*np.sin(theta)
    y_r = -x*np.sin(theta) + y*np.cos(theta)
    return x_r, y_r
