from dataclasses import dataclass, field
import numpy as np

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
    # Number of vertices and list of vertices for the face
    count: int
    vertices: []

    def get_dir(self):
        x1 = min(self.vertices, key = lambda k: k.x)
        x2 = max(self.vertices, key = lambda k: k.x)
        y1 = min(self.vertices, key = lambda k: k.y)
        y2 = max(self.vertices, key = lambda k: k.y)

        return ({
            "x1":x1, "x2":x2,
            "y1":y1, "y2":y2,
            })


@dataclass
class poly:
    # All verts and faces for the polygon
    vertices: list
    faces: list


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
