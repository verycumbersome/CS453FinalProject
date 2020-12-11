#!/usr/bin/env python
import os
import math

import numpy as np
import numpy.linalg as LA


def linearly_interpolate(a, b, M, m, v):
    """Interpolate between m and M with a range of [a, b] and input v"""
    return((b - a) * ((v - m)/(M - m)) + a)

def classify_singularity(face, x, y):
    dir_vec = []
    for vec in ["vx", "vy"]:
        fx1y1 = face.vertices[face.vertices.index(face.y1)][vec]
        fx2y1 = face.vertices[(face.vertices.index(face.y1) + 1) % 4][vec]
        fx2y2 = face.vertices[(face.vertices.index(face.y1) + 2) % 4][vec]
        fx1y2 = face.vertices[(face.vertices.index(face.y1) + 3) % 4][vec]

        x1, x2, y1, y2 = face.x1.x, face.x2.x, face.y1.y, face.y2.y

        # Get fx
        dir_x = (-1)/(x2 - x1) * (y2 - y)/(y2 - y1) * fx1y1 + \
                (1)/(x2 - x1) * (y2 - y)/(y2 - y1) * fx2y1 + \
                (-1)/(x2 - x1) * (y - y1)/(y2 - y1) * fx1y2 + \
                (1)/(x2 - x1) * (y - y1)/(y2 - y1) * fx2y2

        # Get fy
        dir_y = (x2 - x)/(x2 - x1) * (-1)/(y2 - y1) * fx1y1 + \
                (x - x1)/(x2 - x1) * (-1)/(y2 - y1) * fx2y1 + \
                (x2 - x)/(x2 - x1) * (1)/(y2 - y1) * fx1y2 + \
                (x - x1)/(x2 - x1) * (1)/(y2 - y1) * fx2y2

        # Add fx/gx and fy/gy
        dir_vec.append([dir_x, dir_y])

    dir_vec = np.array(dir_vec)
    eig_vals = LA.eigvals(dir_vec)


    if np.iscomplex(eig_vals[0]) or np.iscomplex(eig_vals[0]):
        return("focus")
    elif np.iscomplex(eig_vals[0]) and np.iscomplex(eig_vals[0]):
        return("center")

    elif eig_vals[0] > 0:
        if eig_vals[1] > 0:
            return("nodal_source")
        else:
            return("saddle_point")
    else:
        if eig_vals[1] > 0:
            return("saddle_point")
        else:
            return("nodal_sink")

    return


def get_singularity(face, v_order):
    # Get vector information at points
    vertices = face.vertices

    fx1y1 = vertices[vertices.index(face.y1)].vx
    fx2y1 = vertices[(vertices.index(face.y1) + 1) % 4].vx
    fx2y2 = vertices[(vertices.index(face.y1) + 2) % 4].vx
    fx1y2 = vertices[(vertices.index(face.y1) + 3) % 4].vx

    gx1y1 = vertices[vertices.index(face.y1)].vy
    gx2y1 = vertices[(vertices.index(face.y1) + 1) % 4].vy
    gx2y2 = vertices[(vertices.index(face.y1) + 2) % 4].vy
    gx1y2 = vertices[(vertices.index(face.y1) + 3) % 4].vy

    """Solve quadratic equation given A, B, and C"""
    # Calculate quadratic to find s1, s2 and t
    a00 = fx1y1
    a10 = fx2y1 - fx1y1
    a01 = fx1y2 - fx1y1
    a11 = fx1y1 - fx2y1 - fx1y2 + fx2y2

    b00 = gx1y1
    b10 = gx2y1 - gx1y1
    b01 = gx1y2 - gx1y1
    b11 = gx1y1 - gx2y1 - gx1y2 + gx2y2

    c00 = a11 * b00 - a00 * b11
    c10 = a11 * b10 - a10 * b11
    c01 = a11 * b01 - a01 * b11

    A = (-a11 * c10)
    B = (-a11 * c00 - a01 * c10 + a10 * c01)
    C = (a00 * c01 - a01 * c00)

    discriminant = ((B ** 2) - 4*A*C)

    if (A == 0):
        quad_solution = []

    elif discriminant > 0:
        s1 = (-B + math.sqrt(discriminant)) / (2*A)
        s2 = (-B - math.sqrt(discriminant)) / (2*A)
        quad_solution = [s1, s2]

    elif (discriminant == 0):
        quad_solution = ((-B + math.sqrt(discriminant)) / (2*A))

    else:
        quad_solution = []


    # Linearly interpolate to find singularity
    singularities = []
    for s in quad_solution:
        t = (-c00 / c01) - (c10 / c01) * s
        if (s > 0 and s < 1) and (t > 0 and t < 1):
            singularities.append((
                linearly_interpolate(0, 1, v_order["x1"].x, v_order["x2"].x, s),
                linearly_interpolate(0, 1, v_order["y1"].y, v_order["y2"].y, t),
                # linearly_interpolate(0, 1, self.y1.z, self.y2.z, 0)
                0,
            ))
    return singularities
