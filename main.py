#!/usr/bin/env python
import os
import pdb
import numpy as np

import glfw
import OpenGL

from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *

from PIL import Image

# Import opengl shapes from file
import utils
import config


display_mode = 1 # Display mode for the vector fields
p_file = 1 # Index of PLY file to load from list of 8 PLY files


def init():
    glClearColor(0.2, 0.2, 0.2, 0.0)
    glShadeModel(GL_FLAT)
    glEnable(GL_DEPTH_TEST)
    glEnable(GL_TEXTURE_2D)
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1)


def keyboard(bkey, x, y):
    global display_mode
    global p_file

    key = bkey.decode("utf-8")
    if key == chr(27):
        os._exit(0)

    # Iterates forward to the next p_file
    if key == "n":
        p_file = (p_file + 1) % 8

    # Iterates backwards to the previous p_file
    if key == "p":
        p_file = (p_file - 1) % 8

    if key == "1":
        display_mode = 1

    if key == "2":
        display_mode = 2

    if key == "3":
        display_mode = 3

    if key == "4":
        display_mode = 4

    if key == "5":
        display_mode = 5

    if key == "q":
        os._exit(0)


def display():
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    glEnable(GL_TEXTURE_2D)
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL)

    # Set Projection window view
    glViewport(0, 0, config.WIN_HEIGHT, config.WIN_WIDTH)
    glMatrixMode(GL_PROJECTION)
    gluPerspective(45, (config.WIN_HEIGHT/config.WIN_WIDTH), 0.1, 50.0)
    glLoadIdentity()
    glOrtho(-20.0, 20.0, -20.0, 20.0, -100000.0, 100000.0)

    # Set Model view
    glMatrixMode(GL_MODELVIEW)

    global display_mode
    global p_file

    if display_mode == 1:
        utils.render_ply(p_file, True)

    if display_mode == 2:
        # glRotatef(0.5, 1, 1, 1)
        utils.render_ply(p_file, False)

    glFlush()


def main():
    glutInit()
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH)
    glutInitWindowSize(config.WIN_HEIGHT, config.WIN_WIDTH)
    glutInitWindowPosition(100, 100)
    wind = glutCreateWindow("Surface Parameterization")

    # lighting()

    init()

    pixels = np.empty(config.WIN_WIDTH * config.WIN_HEIGHT * 3, "uint8")
    utils.make_patterns(pixels)

    glutDisplayFunc(display)
    glutIdleFunc(display)

    glutKeyboardFunc(keyboard)
    glutMainLoop()


if __name__ == "__main__":
    main()
