#!/usr/bin/env python
import os
import math
import numpy as np
from dataclasses import dataclass, field
# from numpy.linalg import eig, eigh, eigvals
import numpy.linalg as LA

import OpenGL
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *

import shapes_helper


@dataclass
class vertex:
    # XYZ coordinates for the vertex
    x: float
    y: float
    z: float

    # Vector field info for these points
    vx: float
    vy: float
    vz: float

    # Scalar value at vertex
    s: float

    # Size 3 dict of rgb: rgb["r"] = red_val, etc..
    rgb: dict = field(repr=False, default_factory=dict)

    def __getitem__(self, index):
        return({
            "x": self.x,
            "y": self.y,
            "z": self.z,
            "vx": self.vx,
            "vy": self.vy,
            "vz": self.vz,
        }[index])


@dataclass
class face:
    # Number of vertices and list of all vertices for the face
    count: int
    vertices: list

    has_singularity: bool = False

    def __post_init__(self):
        self.get_vert_order()

    def get_vert_order(self):
        self.x1 = min(self.vertices, key=lambda k: k.x)
        self.x2 = max(self.vertices, key=lambda k: k.x)
        self.y1 = min(self.vertices, key=lambda k: k.y)
        self.y2 = max(self.vertices, key=lambda k: k.y)
        self.v_order = {
            "x1": self.x1, "x2": self.x2,
            "y1": self.y1, "y2": self.y2,
            }

        return ({
            "x1": self.x1, "x2": self.x2,
            "y1": self.y1, "y2": self.y2,
            })

    def get_singularity(self):
        singularities = shapes_helper.get_singularity(
                self, self.v_order)
        if singularities:
            self.has_singularity = True

        return(singularities)


    def classify_singularity(self, singularity):
        x, y = singularity[0], singularity[1]

        return(shapes_helper.classify_singularity(self, x, y))


@dataclass
class polyline:
    # All coordinates for the streamline
    vertices: list

    rgb: dict = field(repr=False, default_factory=dict)

    def render_arrows(self):
        arrow_size = 1
        for i, v in enumerate(self.vertices):
            if i % round(len(self.vertices) / 15) == 0:
                arrow_size = v.s * 2
                print(v.s)
                dir_vec = [v.vx, v.vy]
                # dir_vec = dir_vec / np.linalg.norm(dir_vec)

                theta = math.atan(dir_vec[1] / dir_vec[0])

                v1 = (v.x, v.y, 0)
                v2 = (v.x - dir_vec[0], v.y - dir_vec[1], 0)

                # v3 = (v.x - dir_vec[0], v.y - dir_vec[1], 0)
                # v4 = (v.x - dir_vec[0], v.y - dir_vec[1], 0)
                v2 = (
                    v.x + math.cos(theta - 0.1) * arrow_size,
                    v.y + math.sin(theta - 0.1) * arrow_size,
                    v.z,
                )
                v3 = (
                    v.x + math.cos(theta + 0.1) * arrow_size,
                    v.y + math.sin(theta + 0.1) * arrow_size,
                    v.z
                )

                # glEnd()
                # glBegin(GL_TRIANGLES)
                # glVertex3fv(v1)
                # glVertex3fv(v2)
                # glVertex3fv(v3)
                # glEnd()
                # glBegin(GL_LINES)

    def render_streamline(self):
        for v in self.vertices:
            glColor3f(v.rgb["r"], v.rgb["g"], v.rgb["b"])
            glVertex3fv((v.x, v.y, v.z))


@dataclass
class poly:
    # All verts and faces for the polygon
    vertices: list
    faces: list

    # List of all singularities 
    singularities: list = field(repr=False, default_factory=list)

    # List of all streamlines
    streamlines: list = field(repr=False, default_factory=list)

    # Dict of all colors by classification
    sing_classes: dict = field(repr=False, default_factory=dict)

    def __post_init__(self):
        # Gets all rgb values for classification
        self.sing_classes = {
                "nodal_source":map(lambda x: x/255.0, (225, 0, 0)),
                "nodal_sink":map(lambda x: x/255.0, (0, 225, 0)),
                "saddle_point":map(lambda x: x/255.0, (0, 0, 225)),
                "center":map(lambda x: x/255.0, (104, 134, 197)),
                "focus":map(lambda x: x/255.0, (90, 164, 105)),
            }
        # Maps tuple for each color calue
        self.sing_classes = {k:tuple(v) for k, v in self.sing_classes.items()}

        # Get max and min vertices for poly
        max_xyz = [0, 0, 0]
        min_xyz = [0, 0, 0]

        for v in self.vertices:
            if max_xyz[0] < v.x:
                max_xyz[0] = v.x
            if max_xyz[1] < v.y:
                max_xyz[1] = v.y
            if max_xyz[2] < v.z:
                max_xyz[2] = v.z
            if min_xyz[0] > v.x:
                min_xyz[0] = v.x
            if min_xyz[1] > v.y:
                min_xyz[1] = v.y
            if min_xyz[2] > v.z:
                min_xyz[2] = v.z
        self.max_xyz, self.min_xyz = max_xyz, min_xyz

        self.calculate_singularities()
        if not self.streamlines:
            self.streamlines = []
            for s in self.singularities:
                # Calculates streamline for singularity
                self.calculate_streamline(s)

        # for i in range(-10, 10, 5):
            # for j in range(-10, 10, 5):
                # print(i, j)
                # self.calculate_streamline((i , j, 0))


    def get_dir(self, coords):
        """
        Gets the direction of point on the poly
        Output: (dir_x, dir_y, dir_z) - normalized vec in direction at coords
        """
        x, y, z = coords[0], coords[1], coords[2]
        dir_vec = []

        # Iterate through faces until x is found within face vertices
        for face in self.faces:
            v = face.get_vert_order()

            if (x <= v["x2"].x and x >= v["x1"].x
                and y <= v["y2"].y and y >= v["y1"].y):
                break

        # Linearly interpolate between face vertices to find vx vy at x, y, z
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

    def calculate_streamline(self, vert):
        """Calculate the streamline for a poly given coordinates"""
        x, y, z = vert["coordinates"]
        step = 0.01
        count = 2000

        streamline = []

        for i in [1, -1]:
            curr = np.array([x, y, z])

            for j in range(count):
                # If calculated value is out of bounds return nothing
                if (curr[0] > self.max_xyz[0] or curr[1] > self.max_xyz[1]
                                              or curr[2] > self.max_xyz[2]):
                    break

                if (curr[0] < self.min_xyz[0] or curr[1] < self.min_xyz[1]
                                              or curr[2] < self.min_xyz[2]):
                    break

                # Get vx vy for coordinates
                dir_vec = self.get_dir(curr)
                # Add vertices to polyline
                streamline.append(vertex(
                    curr[0],
                    curr[1],
                    curr[2],
                    dir_vec[0],
                    dir_vec[1],
                    dir_vec[2],
                    0,
                    rgb={"r": dir_vec[0], "g": dir_vec[1], "b": 1},
                ))

                curr = curr + dir_vec * i * step # i indicates -1 or 1
                streamline.append(vertex(
                    curr[0],
                    curr[1],
                    curr[2],
                    dir_vec[0],
                    dir_vec[1],
                    dir_vec[2],
                    0,
                    rgb={"r": dir_vec[0], "g": dir_vec[1], "b": 1},
                ))

        self.streamlines.append(polyline(streamline))

    def calculate_singularities(self):
        """Calculate and find all the singularities for a given poly"""
        if not self.singularities:
            self.singularities = []

            for face in self.faces:
                for s in face.get_singularity():
                    # Classification for singularity
                    c = face.classify_singularity(s)

                    self.singularities.append({
                        "type":c,
                        "coordinates":s,
                        "rgb":self.sing_classes[c],
                        "s":face.vertices[0].s,
                    })

    def render_singularities(self):
        glPointSize(10.0)
        glBegin(GL_POINTS)
        for s in self.singularities:
            # Colors singularity based on class
            glColor3f(s["rgb"][0], s["rgb"][1], s["rgb"][2])
            # Renders singularity vertex
            glVertex3fv(s["coordinates"])
        glEnd()
        glPointSize(2.0)

    def render_streamlines(self):
        glPointSize(2.0)
        glBegin(GL_LINES)

        for stream in self.streamlines:
            stream.render_streamline()
            # stream.render_arrows()
        glEnd()


def read_ply(filename):
    """
    Read and parse a PLY file
    Output: PLY(vertices, faces) - vertices and faces parsed from file
    """
    with open(filename, "r") as ply_file:
        header = ply_file.read().split("end_header")
        header, data = header[0], header[1][1:].split("\n")

        n_verts = [x for x in header.split("\n") if "element vertex" in x][0]
        n_faces = [x for x in header.split("\n") if "element face" in x][0]

        n_verts = int(n_verts.split()[-1])
        n_faces = int(n_faces.split()[-1])

        vertices = []
        for i in range(n_verts):
            v_info = data.pop(0).split()

            v = vertex(
                    int(v_info[0]),    # x
                    int(v_info[1]),    # y
                    int(v_info[2]),    # z
                    float(v_info[3]),  # vx
                    float(v_info[4]),  # vy
                    float(v_info[5]),  # vz
                    float(v_info[6]),  # s
                    {"r": 1.0, "g": 1.0, "b": 1.0},  #rgb
                )
            vertices.append(v)

        faces = []
        for i in range(n_faces):
            f_info = data.pop(0).split()
            f = face(
                    int(f_info[0]),
                    [
                        vertices[int(f_info[1])], # v1
                        vertices[int(f_info[2])], # v2
                        vertices[int(f_info[3])], # v3
                        vertices[int(f_info[4])], # v4
                    ]
                )
            faces.append(f)

        return(poly(vertices, faces))
