import pdb
import math
import glfw
import numpy as np
import OpenGL
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *


def plane(length):
    glBegin(GL_QUADS)
    glNormal3f(0.0, 0.0, (length / 2))

    glTexCoord2d(1, 1)
    glVertex3f((-length / 2), (-length / 2), 0.0)
    glTexCoord2d(1, 0)
    glVertex3f((-length  / 2), (length / 2), 0.0)
    glTexCoord2d(0, 0)
    glVertex3f((length / 2), (length / 2), 0.0)
    glTexCoord2d(0, 1)
    glVertex3f((length / 2), (-length / 2), 0.0)

    glEnd()


def cylinder_dydv(u, v):
    """
    Function to calculate derivative of cylinder function
    """
    radius = 4.0
    height = 10.0
    v = math.radians(v * 360.0)

    return(
        -math.sin(v),
        math.cos(v),
        1
    )

def cylinder_dydu(u, v):
    """
    Function to calculate derivative of cylinder function
    """
    radius = 4.0
    height = 10.0
    v = math.radians(v * 360.0)

    return(
        0,
        0,
        height
    )


def cylinder(u, v):
    """
    Function that returns 3d coordinates for a cylinder from 2d input
    Input: (u, v) -> (height[0-1], theta[0-1])
    Ouput: (x, y, z) coordinate
    """
    radius = 4.0
    height = 10.0
    v = math.radians(v * 360.0)

    return(
        radius * math.cos(v),
        radius * math.sin(v),
        u * height
    )


def vase(u, v):
    """
    Function that returns 3d coordinates for a cylinder from 2d input
    Input: (u, v) -> (height[0-1], theta[0-1])
    Ouput: (x, y, z) coordinate
    """
    radius = 4.0

    v = math.radians(v * 360.0)
    u_rad = math.radians(u * 360.0)

    return(
        radius * math.cos(u_rad) * math.cos(v),
        radius * math.cos(u_rad) * math.sin(v),
        (2 * radius * u) - radius
    )


def sphere(u, v):
    """
    Function that returns 3d coordinates for a cylinder from 2d input
    Input: (u, v) -> (height[0-1], theta[0-1])
    Ouput: (x, y, z) coordinate
    """
    radius = 2.0

    v = math.radians(v * 360.0)
    d = math.sqrt(radius ** 2 - u ** 2)

    return(
        radius * d * math.sin(v),
        radius * d * math.cos(v),
        u * radius
    )


def render_shape(z_res, xy_res, shape_func, multiplier = 1):
    """
    Function to render OpenGL shape given function of form T(u, v)
    Input:  xy_res -> Resolutions of the xy for cylinder
            z_res -> Resolutions of the z for cylinder
    Output: renders OpenGl cylinder
    """
    # Get step size given subdivision
    z_step = (1 / z_res)
    xy_step = (1 / xy_res)

    # Get inputs for cylinder function
    z_res = round(z_res / 2) * multiplier
    z_range = [x * z_step for x in range(-z_res, z_res)]
    xy_range = [x * xy_step for x in range(0, xy_res)]

    for u in z_range:
        for v in xy_range:
            glBegin(GL_QUAD_STRIP)

            # Matrix for storing traingle strip vertex coordinates
            # [1][2] - Top left | Top Right
            # [0][3] - Bottom left | Bottom right
            quad_mat = np.array([
                shape_func(u, v),
                shape_func(u, v + z_step),
                shape_func(u + xy_step, v),
                shape_func(u + xy_step, v + z_step),
            ])

            # Create side vertices
            for i in range(4):
                # p2 - p1
                Vu = quad_mat[(i + 1) % 4] - quad_mat[i]
                # p3 - p1                
                Vv = quad_mat[(i + 2) % 4] - quad_mat[i]

                # normal x = Uy * Vz - Uz * Vy
                normal_x = Vu[1] * Vv[2] - Vu[2] * Vv[1]
                # normal y = Uz * Vx - Ux * Vz
                normal_y = Vu[2] * Vv[0] - Vu[0] * Vv[2]
                # normal z = Ux * Vy - Uy * Vx
                normal_z = Vu[0] * Vv[1] - Vu[1] * Vv[0]

                # glTexCoord2d(0.2 * quad_mat[i][0], 0.2 * quad_mat[i][2])
                glNormal3f(normal_x, normal_y, normal_z)
                # glNormal3f(quad_mat[i][0], quad_mat[i][1], quad_mat[i][2])
                glTexCoord2d(quad_mat[i][0] / xy_res, quad_mat[i][2] / z_res)
                glVertex3fv(quad_mat[i])

            glEnd()

