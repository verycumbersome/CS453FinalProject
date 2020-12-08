import os
import pdb
import math
import glfw
import numpy as np
import OpenGL
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *

import shapes

# Load all poly files
ply_filepath = "new_vector_data/"

poly_files = [shapes.read_ply(ply_filepath + x) for x in os.listdir(ply_filepath)]
poly = poly_files[0]


def linearly_interpolate(a, b, M, m, v):
    """Interpolate between m and M with a range of [a, b] and input v"""
    return((b - a) * ((v - m)/(M - m)) + a)


def solve_quadratic(A, B, C):
    """Solve quadratic equation given A, B, and C"""
    if (A == 0):
        return([])

    discriminant = ((B ** 2) - 4*A*C)
    if discriminant > 0:
        s1 = (-B + math.sqrt(discriminant)) / (2*A)
        s2 = (-B - math.sqrt(discriminant)) / (2*A)
        return([s1, s2])

    elif (discriminant == 0):
        return((-B + math.sqrt(discriminant)) / (2*A))

    else:
        return([])


def display_IBFV():
    pass


def render_ply(display_mode, render_streamline):
    """
    Function to render OpenGL shape given function of form T(u, v)
    Input: ply - a PLY file containing a list of vertices anf faces
    Output: renders OpenGl PLY file
    """
    global poly
    poly = poly_files[display_mode]
    singularities = []

    # Find and render all singularities
    poly.calculate_singularities()
    poly.render_singularities()


    # Find and render all streamlines
    if render_streamline:
        poly.render_streamlines()

    # Render all vertices in poly
    glBegin(GL_POINTS)
    for vertex in poly.vertices:
        glColor3f(vertex.rgb["r"], vertex.rgb["g"], vertex.rgb["b"])
        glVertex3fv((vertex.x, vertex.y, vertex.s))
    glEnd()
