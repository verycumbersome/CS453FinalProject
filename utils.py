import os
import pdb
import math
import random

import glfw
import numpy as np
import OpenGL
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *

from ctypes import *

from PIL import Image

import shapes
import config


# #IBFV stuff
win_width = config.WIN_WIDTH
win_height = config.WIN_HEIGHT
NPN = 64
SCALE = 4.0
Npat = 32
iframe = 0
tmax = win_width / (SCALE*NPN)
dmax = SCALE / win_width

pixels = np.empty(win_width * win_height * 3, "uint8")

DM = float(1.0/(100-1.0))

# Load all poly files
poly = config.POLY_FILES[0]


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


def make_patterns(pixels_reset):
    global pixels

    pixels = pixels_reset

    # pixels = np.empty(win_width * win_height * 3, "uint8")

    # int lut[256]
    # int phase[NPN][NPN]
    # GLubyte pat[NPN][NPN][4]
    # int i, j, k, t

    # for (i = 0 i < 256 i++) lut[i] = i < 127 ? 0 : 255
    # for (i = 0 i < NPN i++)
    # for (j = 0 j < NPN j++) phase[i][j] = rand() % 256

    lut = np.zeros(256)
    phase = np.zeros((NPN, NPN))
    pat = np.zeros((NPN, NPN, 4), "uint8")
    # for i in range(NPN):
        # for j in range(NPN):
            # pat[i][j] = [0,0,0,0]

    i, j, k, t = 0,0,0,0

    for i in range(256):
        lut[i] = (0 if i < 127 else 255)


    for i in range(NPN):
        for j in range(NPN):
            phase[i][j] = random.randrange(0, 10000) % 256


    for k in range(Npat):
        t = k * Npat
        for i in range(NPN):
            for j in range(NPN):
                pat[i][j][0] = lut[int(t + phase[i][j]) % 255]
                pat[i][j][1] = lut[int(t + phase[i][j]) % 255]
                pat[i][j][2] = lut[int(t + phase[i][j]) % 255]
                pat[i][j][3] = int(0.12 * 255)

        glNewList(k + 1, GL_COMPILE)
        glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0, GL_RGBA, GL_UNSIGNED_BYTE, pat)
        glEndList()

        # pixels = pat

        # print(pixels)


    # for (k = 0 k < Npat k++) {
        # t = k * 256 / Npat
        # for (i = 0 i < NPN i++)
        # for (j = 0 j < NPN j++) {
            # pat[i][j][0] =
            # pat[i][j][1] =
            # pat[i][j][2] = lut[(t + phase[i][j]) % 255]
            # pat[i][j][3] = int(0.12 * 255)
            # }
        # glNewList(k + 1, GL_COMPILE)
        # # glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0, GL_RGBA, GL_UNSIGNED_BYTE, pat)
        # glEndList()

def display_IBFV():
    global iframe
    global poly

    glDisable(GL_LIGHTING)
    glDisable(GL_LIGHT0)
    glDisable(GL_LIGHT1)
    glDisable(GL_POLYGON_OFFSET_FILL)
    glDisable(GL_DEPTH_TEST)

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE)

    glEnable(GL_TEXTURE_2D)
    glShadeModel(GL_FLAT)
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

    glClearColor(1.0, 1.0, 1.0, 1.0) #  // background for rendering color coding and lighting
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

    # /*draw the model with using the pixels, using vector field to advert the texture coordinates*/
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, win_width, win_height, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels)
    # glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, img.size[0], img.size[1], 0, GL_RGB, GL_UNSIGNED_BYTE, img_data)

    # double modelview_matrix1[16], projection_matrix1[16]
    modelview_matrix1, projection_matrix1 = np.zeros(16, "double"), np.zeros(16, "double")
    viewport1 = np.zeros(4)

    modelview_matrix1 = glGetDoublev(GL_MODELVIEW_MATRIX)
    projection_matrix1 = glGetDoublev(GL_PROJECTION_MATRIX)
    viewport1 = glGetIntegerv(GL_VIEWPORT)

    # for (int i = 0 i < poly->nquads i++) { //go through all the quads
    for temp_q in poly.faces:
        glBegin(GL_QUADS)

        # for (int j = 0 j < 4 j++) {
        for temp_v in temp_q.vertices:
            x = temp_v.x
            y = temp_v.y

            tx, ty, dummy = 0., 0., 0.

            # gluProject((GLdouble)temp_v->x, (GLdouble)temp_v->y, (GLdouble)temp_v->z,
                # modelview_matrix1, projection_matrix1, viewport1, &tx, &ty, &dummy)
            tx, ty, dummy = gluProject(temp_v.x, temp_v.y, temp_v.z, modelview_matrix1, projection_matrix1, viewport1)

            tx = tx / win_width
            ty = ty / win_height

            # dp = icVector2(temp_v.vx, temp_v.vy)
            # normalize(dp)
            # dx = dp.x
            # dy = dp.y

            dir_vec = [temp_v.vx, temp_v.vy]
            # dir_vec.append(z)
            dir_vec = np.array(dir_vec)
            dir_vec = dir_vec / np.sqrt(np.sum(dir_vec**2))

            dx = dir_vec[0]
            dy = dir_vec[1]

            r = dx * dx + dy * dy
            if (r > dmax*dmax):
                r = math.sqrt(r)
                dx *= dmax / r
                dy *= dmax / r

            px = tx + dx
            py = ty + dy

            glTexCoord2f(px, py)
            glVertex3d(temp_v.x, temp_v.y, temp_v.z)

        glEnd()


    iframe = iframe + 1

    glEnable(GL_BLEND)

    # /*blen2d the drawing with another noise image*/
    glMatrixMode(GL_PROJECTION)
    glPushMatrix()
    glLoadIdentity()


    glMatrixMode(GL_MODELVIEW)
    glPushMatrix()
    glLoadIdentity()

    glTranslatef(-1.0, -1.0, 0.0)
    glScalef(2.0, 2.0, 1.0)

    glCallList(iframe % Npat + 1)

    glBegin(GL_QUAD_STRIP)

    glTexCoord2f(0.0, 0.0)
    glVertex2f(0.0, 0.0)
    glTexCoord2f(0.0, tmax)
    glVertex2f(0.0, 1.0)
    glTexCoord2f(tmax, 0.0)
    glVertex2f(1.0, 0.0)
    glTexCoord2f(tmax, tmax)
    glVertex2f(1.0, 1.0)

    glEnd()
    glDisable(GL_BLEND)

    glMatrixMode(GL_MODELVIEW)
    glPopMatrix()

    glMatrixMode(GL_PROJECTION)
    glPopMatrix()

    glReadPixels(0, 0, win_width, win_height, GL_RGB, GL_UNSIGNED_BYTE, pixels)

    # /*draw the model with using pixels, note the tx and ty do not take the vector on points*/
    glClearColor(0.3, 0.3, 0.3, 1.0)  # background for rendering color coding and lighting
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, win_width, win_height, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels)

    for temp_q in poly.faces:# //go through all the quads
        glBegin(GL_QUADS)
        for temp_v in temp_q.vertices:
            x = temp_v.x
            y = temp_v.y
            tx, ty, dummy = 0.,0.,0.
            # gluproject((gldouble)temp_v.x, (gldouble)temp_v->y, (gldouble)temp_v->z,
                # modelview_matrix1, projection_matrix1, viewport1, &tx, &ty, &dummy)
            tx, ty, dummy = gluProject(temp_v.x, temp_v.y, temp_v.z, modelview_matrix1, projection_matrix1, viewport1)
            tx = tx / win_width
            ty = ty / win_height
            glTexCoord2f(tx, ty)
            glVertex3d(temp_v.x, temp_v.y, temp_v.z)

        glEnd()


    glDisable(GL_TEXTURE_2D)
    glShadeModel(GL_SMOOTH)
    glDisable(GL_BLEND)


def render_ply(display_mode, render_streamline):
    """
    Function to render OpenGL shape given function of form T(u, v)
    Input: ply - a PLY file containing a list of vertices anf faces
    Output: renders OpenGl PLY file
    """
    global poly
    poly = config.POLY_FILES[display_mode]
    singularities = []

    # Find and render all singularities
    poly.calculate_singularities()
    poly.render_singularities()


    # Find and render all streamlines
    if render_streamline:
        poly.render_streamlines()

    else:
        display_IBFV()

    # Render all vertices in poly
    glBegin(GL_POINTS)
    for vertex in poly.vertices:
        glColor3f(vertex.rgb["r"], vertex.rgb["g"], vertex.rgb["b"])
        glVertex3fv((vertex.x, vertex.y, vertex.z))
    glEnd()
