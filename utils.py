import math
import numpy as np
from dataclasses import dataclass, field

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

    def get_vert_order(self):
        self.x1 = min(self.vertices, key = lambda k: k.x)
        self.x2 = max(self.vertices, key = lambda k: k.x)
        self.y1 = min(self.vertices, key = lambda k: k.y)
        self.y2 = max(self.vertices, key = lambda k: k.y)

        return ({
            "x1":self.x1, "x2":self.x2,
            "y1":self.y1, "y2":self.y2,
            })

    def get_singularity(self):
        fx1y1 = self.vertices[self.vertices.index(self.y1)].vx
        fx2y1 = self.vertices[(self.vertices.index(self.y1) + 1) % 4].vx
        fx2y2 = self.vertices[(self.vertices.index(self.y1) + 2) % 4].vx
        fx1y2 = self.vertices[(self.vertices.index(self.y1) + 3) % 4].vx

        gx1y1 = self.vertices[self.vertices.index(self.y1)].vy
        gx2y1 = self.vertices[(self.vertices.index(self.y1) + 1) % 4].vy
        gx2y2 = self.vertices[(self.vertices.index(self.y1) + 2) % 4].vy
        gx1y2 = self.vertices[(self.vertices.index(self.y1) + 3) % 4].vy

        print(fx1y1, fx2y1, fx1y2, fx2y2)
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

        print(a11)

        A = (-a11 * c10)
        B = (-a11 * c00 -a01 * c10 + a10 * c01)
        C = (a00 * c01 - a01 * c00)

        solve_quadratic(A, B, C)

@dataclass
class poly:
    # All verts and faces for the polygon
    vertices: list
    faces: list


def solve_quadratic(A, B, C):
    discriminant = ((B ** 2) - 4*A*C)
    if discriminant > 0:
        x = (-B + math.sqrt(discriminant)) / (2*A)

    return(x)



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
                    int(v_info[0]),   # x
                    int(v_info[1]),   # y
                    int(v_info[2]),   # z
                    float(v_info[3]), # vx
                    float(v_info[4]), # vy
                    float(v_info[5]), # vz
                    float(v_info[6]), # s
                    {"r": 1.0, "g": 1.0, "b": 1.0}, #rgb
                )
            vertices.append(v)

        faces = []
        for i in range(n_faces):
            f_info = data.pop(0).split()
            f = face(
                    int(f_info[0]),
                    [
                        vertices[int(f_info[1])],
                        vertices[int(f_info[2])],
                        vertices[int(f_info[3])],
                        vertices[int(f_info[4])],
                    ]
                )
            faces.append(f)

        return(poly(vertices, faces))
