import os
import numpy as np
import shapes

WIN_WIDTH = 1500
WIN_HEIGHT = 1500

# Load all poly files
ply_filepath = "new_vector_data/"
POLY_FILES = [shapes.read_ply(ply_filepath + x) for x in os.listdir(ply_filepath)]
