# -*- coding: utf-8 -*-
"""
Script for finding if a vector from an instrument intersects with any part 
of the spacecraft

INPUT:
  - model path: a character array with the path of the spacecraft model, 
       including the extension                                            [.obj or .stl]
       e.g. 'C:/Users/Darth_Vader/Documents/Cassini_3D_model_example/Cassini_NASA_model.obj'
  - [x,y,z]: the instrument's location on the model                       [model units / scale]
  - [vector_matrix]: n x 3 matrix with the vector of interest, e.g.
      instrument-->Sun                                                    [model units / scale]

OUTPUT:
  -  index: index of whether the given vector intersects the spacecraft
      or not; 1 if it does, 0 if it doesn't                               [0 or 1]
  -  phi: the phi angle of each of the vectors                            [rad] -> [-pi pi]
  -  theta: the theta angle of each of the vectors                        [rad] -> [-pi/2 pi/2]

------------------------------------
NOTES:
- As it takes quite some time to find the contour for a model it is to
calculate the contour once and save it. Then every time you want to use
the contour load the saved one.
- The units for the instruments x,y,z should be in the same units as the
model. e.g. if the 1 unit of the model corresponds to 1 metre, the
instrument's x,y,z "offset" from the s/c origin should also follow this.
- The units of the vector matrix should follow the units of the model.

------------------------------------

This is part of the open-access published work "A simple spacecraft - vector 
intersection methodology and applications", DOI:10.1093/rasti/rzae012
https://academic.oup.com/rasti/article/3/1/166/7634356

Author: George Xystouris (george.xystouris@gmail.com)
Originally Created for MATLAB in August 2020
Translated for Python in November 2025

v1.0
"""

import numpy as np
from shapely.geometry import Point
from shapely.prepared import prep

from svi_instr_field_of_view import svi_instr_field_of_view


# %%

def svi_instr_sc_vector_intersection(model_path, instr_x, instr_y, instr_z, process_faces, vector_matrix):

    # --- Create the instrument's FOV ---
    fov = svi_instr_field_of_view(model_path, instr_x, instr_y, instr_z, process_faces)
    unified_fov = fov["unified"]

    # --- Convert the vectors from cartesian to spherical coordinates ---
    x, y, z = vector_matrix[:, 0], vector_matrix[:, 1], vector_matrix[:, 2]
    vec_phi = np.arctan2(y, x)  # azimuth
    vec_theta = np.arctan2(z, np.sqrt(x**2 + y**2))  # elevation

    # Convert negative azimuths to positive (0–2π)
    vec_phi[vec_phi < 0] += 2 * np.pi

    # --- Prepare unified geometry for fast lookups ---
    prepared_fov = prep(unified_fov)

    # --- Vectorized check for interior points ---
    # Create shapely Points for all vectors at once
    points = [Point(phi, theta) for phi, theta in zip(vec_phi, vec_theta)]
    index = np.fromiter((prepared_fov.contains(p) for p in points), dtype=bool)

    # Convert boolean array to integers (1=inside, 0=outside)
    index = index.astype(int)

    return index, vec_phi, vec_theta













