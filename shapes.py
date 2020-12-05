import pdb
import math
import glfw
import numpy as np
import OpenGL
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *

import utils
poly_files = [
    utils.read_ply("new_vector_data/v1.ply"),
    utils.read_ply("new_vector_data/v3.ply"),
    utils.read_ply("new_vector_data/v4.ply"),
    utils.read_ply("new_vector_data/v5.ply"),
    utils.read_ply("new_vector_data/v6.ply"),
    utils.read_ply("new_vector_data/v8.ply"),
    utils.read_ply("new_vector_data/v9.ply"),
    utils.read_ply("new_vector_data/v10.ply")
]
poly = poly_files[0]


def get_dir(x, y, z):
    dir_vec = []

    for face in poly.faces:
        v = face.get_vert_order()

        if (x <= v["x2"].x and x >= v["x1"].x
            and y <= v["y1"].y and y <= v["y1"].y):
            break

    for vec in ["vx", "vy"]:
        fx1y1 = face.vertices[face.vertices.index(v["y1"])][vec]
        fx2y1 = face.vertices[(face.vertices.index(v["y1"]) + 1) % 4][vec]
        fx2y2 = face.vertices[(face.vertices.index(v["y1"]) + 2) % 4][vec]
        fx1y2 = face.vertices[(face.vertices.index(v["y1"]) + 3) % 4][vec]

        x1, x2, y1, y2 = v["x1"].x, v["x2"].x, v["y1"].y, v["y2"].y

        dir_v = (x2 - x)/(x2 - x1) * (y2 - y)/(y2 - y1) * fx1y1 + \
                (x - x1)/(x2 - x1) * (y2 - y)/(y2 - y1) * fx2y1 + \
                (x2 - x)/(x2 - x1) * (y - y1)/(y2 - y1) * fx1y2 + \
                (x - x1)/(x2 - x1) * (y - y1)/(y2 - y1) * fx2y2

        dir_vec.append(dir_v)

    # Normalize direction
    dir_vec.append(z)
    dir_vec = np.array(dir_vec)
    dir_vec = dir_vec / np.sqrt(np.sum(dir_vec**2))

    return(dir_vec)


def extract_streamline(x, y, z):
    step = 0.1
    count = 200
    max_xyz = [0, 0, 0]
    min_xyz = [0, 0, 0]

    curr = np.array([x, y, z])

    for vertex in poly.vertices:
        if max_xyz[0] < vertex.x:
            max_xyz[0] = vertex.x

        if max_xyz[1] < vertex.y:
            max_xyz[1] = vertex.y

        if max_xyz[2] < vertex.z:
            max_xyz[2] = vertex.z

        if min_xyz[0] > vertex.x:
            min_xyz[0] = vertex.x

        if min_xyz[1] > vertex.y:
            min_xyz[1] = vertex.y

        if min_xyz[2] > vertex.z:
            min_xyz[2] = vertex.z

    if (curr[0] > max_xyz[0] or curr[1] > max_xyz[1] or curr[2] > max_xyz[2]):
        pass

    if (curr[0] < min_xyz[0] or curr[1] < min_xyz[1] or curr[2] < min_xyz[2]):
        pass

    for i in range(count):
        glVertex3fv(tuple(curr))
        curr = curr + get_dir(curr[0], curr[1], curr[2]) * step
        glVertex3fv(tuple(curr))


def render_ply(display_mode):
    """
    Function to render OpenGL shape given function of form T(u, v)
    Input: ply - a PLY file containing a list of vertices anf faces
    Output: renders OpenGl PLY file
    """
    global poly
    poly = poly_files[display_mode]

    # Find and render all streamlines
    glBegin(GL_LINES)
    # extract_streamline(3, 1, 0)
    glEnd()

    # Find and render all singularities
    glPointSize(5.0)
    glBegin(GL_POINTS)

    for face in poly.faces:
        singularity = face.get_singularity()
        for sing in singularity:
            glVertex3fv((sing[0], sing[1], 0))

    glEnd()
    glPointSize(1.0)

    # Render all vertices in poly
    glBegin(GL_POINTS)
    for vertex in poly.vertices:
        glColor3f(vertex.rgb["r"], vertex.rgb["g"], vertex.rgb["b"])
        glVertex3fv((vertex.x, vertex.y, vertex.s))
    glEnd()
