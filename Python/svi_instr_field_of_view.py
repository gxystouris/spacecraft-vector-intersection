# -*- coding: utf-8 -*-
"""

Script for finding an instrument's field-of-view (FOV)

INPUT:
  - spacecraft model path: a character array with the path of the 
       spacecraft model including the extension                             [.obj or .stl]
       e.g. 'C:/Users/Darth_Vader/Documents/Cassini_3D_model_example/Cassini_NASA_model.obj'
  - [x,y,z]: the instrument's location on the model                         [model units / scale]
  - process_faces: flag on how to process the "problematic" faces           [case insensitive string]
       "none":    no process is taking place
       'remove':  all the faces of interest are removed
       'split':   all the faces of interest are splitted into faces in 
                  negative y and in postivie y. Those faces are added 
                  to the model.

OUTPUT:
  - fov structure:
      .faces:     the individual faces                                      [Polygon]
      .unified:   the contour of all the faces                              [MultiPolygon]
          ***SHAPE LIMITS: azimuth: [0 2pi]; elevation: [-pi/2 pi/2]
          
          

------------------------------------
NOTES:
- As it takes quite some time to find the contour for a model it is to
calculate the contour once and save it. Then every time you want to use
the contour load the saved one.

- The units for the instruments x,y,z should be in the same units as the
model. e.g. if the 1 unit of the model corresponds to 1 metre, the
instrument's x,y,z "offset" from the s/c origin should also follow this.

- Problematic faces are faces that cross the 2pi->0 line, i.e. going from
negative to positive y -or vice versa- on the model. The issue is that 
in those cases, e.g. a face going crossing the 2pi->0 line, the programme 
cannot understand that the face must be cut into two parts: one for the
negative y and one for the positive y. Instead it stretches the face
from negative y to positive y the "other" way around.
  E.g. if the face has an azimuth of 20 degrees and crosses the 2pi->0 
  line, MATLAB will create a face of 340 degrees, i.e. what's missing 
  from a circle.
Therefore, here we are splitting the face in 3, as follow:
    Let the face has vertices A, B, and C. Vertex A is the "lone" on 
    the "other side" of B and C. P1 is the point on line AB where y=0, 
    and P2 is the point on line AC where y=0. Therefore, the face ABC is 
    splitted into faces AP1P2, BP1C, and CP1P2.
    
Regarding the splitting: in we included the rare cases that a face will span 
over the x=0 and y=0 planes. In this case we first split the face on the 
x plane, and then we take the face(s) and split them over the y plane.

------------------------------------

This is part of the open-access published work "A simple spacecraft - vector 
intersection methodology and applications", DOI:10.1093/rasti/rzae012
https://academic.oup.com/rasti/article/3/1/166/7634356

Author: George Xystouris (george.xystouris@gmail.com)
Originally Created for MATLAB in August 2020
Translated for Python in November 2025

v1.0
"""

# %% Imports and environment setup

import os

import numpy as np    
import trimesh # might need to be installed. Install by typing    pip install trimesh shapely    at the console. If needed to upgrade type       !{sys.executable} -m pip install --upgrade pip
from shapely.geometry import Polygon
from shapely.ops import unary_union


# %% Definition of functions

def cart2sph(x, y, z):
    """ Convert Cartesian coordinates to Spherical coordinates.
    Returns: azimuth: [0  2pi]; elevation: [-pi/2  pi/2]"""
    r = np.sqrt(x**2 + y**2 + z**2)
    az = np.arctan2(y, x)
    el = np.arcsin(z / r)

    # shift azimuth to go from 0 to 2pi
    az = np.where(az < 0, az + 2*np.pi, az)
    return az, el, r


# Parameterized line intersection with x=0 plane
def line_zero_x(p1, p2):
    t = -p1[0] / (p2[0] - p1[0])
    return p1 + t * (p2 - p1)


# Parameterized line intersection with y=0 plane
def line_zero_y(p1, p2):
    t = -p1[1] / (p2[1] - p1[1])
    return p1 + t * (p2 - p1)


# %% MAIN SCRIPT

def svi_instr_field_of_view(model_path, instr_x, instr_y, instr_z, process_faces='none'):
    
    # --- File check ---
    # Check that file is in a valid extension (.obj or .stl). If it's not, the script stops.
    ext = os.path.splitext(model_path)[1].lower()
    if ext not in [".obj", ".stl"]:
        raise ValueError("Error in opening the model file. Fast checks:\n - File name: must end in .obj or .stl)\n - File location")


    # --- Load the model ---
    mesh = trimesh.load(model_path, force='mesh')
    vertices = mesh.vertices.copy()
    faces = mesh.faces.copy()
    
    
    # --- Move instrument location to orgin ---
    vertices[:, 0] -= instr_x
    vertices[:, 1] -= instr_y
    vertices[:, 2] -= instr_z
    
    # --- Check if anything is EXACTLY on y=0 or x=0. If so move the vertex slightly to the right
    ind_zero = np.where(vertices[:,0] == 0)[0]
    if len(ind_zero) > 0:
        vertices[ind_zero,1] = vertices[ind_zero,0] + 10**-10
    ind_zero = np.where(vertices[:,1] == 0)[0]
    if len(ind_zero) > 0:
        vertices[ind_zero,1] = vertices[ind_zero,1] + 10**-10

    # --- Convert to spherical coordinates ---
    Vaz, Vel, _ = cart2sph(vertices[:, 0], vertices[:, 1], vertices[:, 2])

    # --- Process faces ---
    face_polygons = []

    for f in faces:
        i1, i2, i3 = f
        
        # Vertex coordinates in spherical form
        az = [Vaz[i1], Vaz[i2], Vaz[i3]]
        el = [Vel[i1], Vel[i2], Vel[i3]]

        # Cartesian coordinates for splitting process
        xyz = vertices[f]

        # Determine if crossing y=0 plane
        # If a face crosses the y=0 plane, the y's will have different signs. 
        # Regarding the x's: all or some might be x>0
       
        x_signs = np.sign(xyz[:, 0])
        crosses_x = np.any(xyz[:,0] > 0)  # at least one vertex is at x>0
        
        y_signs = np.sign(xyz[:, 1])
        crosses_y = (not np.all(y_signs == y_signs[0]))  # at least two different signs
        
        
        # --- Manipulate faces based on the input processing way
        if process_faces.lower() == 'remove':
            if crosses_x and crosses_y:
                continue  # skip this face entirely

        elif process_faces.lower() == 'split' and crosses_x and crosses_y:
            # Check if need to split on x-axis
            if np.all(x_signs > 0):     # if all 3 vertices are on x>0, no splitting required
                split_faces = [xyz]

            else:
                # --- Split on x-axis
                # Count how many times + and - occur, to define if the "lone" face is in the negatives or positives
                if np.sum(x_signs > 0) == 1:
                    idx_a = np.where(x_signs > 0)[0][0]   # the lone positive
                elif np.sum(x_signs < 0) == 1:
                    idx_a = np.where(x_signs < 0)[0][0]   # the lone negative
    
                # Get the other two vertices
                idx_bc = [i for i in range(3) if i != idx_a]
    
                # Gather all vertices 
                A = xyz[idx_a]          # the vertex on the "other" side
                B = xyz[idx_bc[0]]
                C = xyz[idx_bc[1]]
    
    
                # The +/-10**-15 value on the y-component of each face is there for
                # Python to give correct results: if a point lays EXACTLY on 0, 
                # Python sometimes tends to plot wrongly (e.g. if it "comes" from
                # the positive side, sometimes Python thinks it's on negative). 
                # We are working around this by adding a TINY positive or negative 
                # value.
                #
                # P1 and P2 are for the faces in the positive side of the x axis, 
                # and P1_neg and P2_neg are for the faces in the negative side
                P1 = line_zero_x(A, B)
                P2 = line_zero_x(A, C)
                P1 = np.float64([10**-15, P1[1], P1[2]])
                P2 = np.float64([10**-15, P2[1], P2[2]])            
                P1_neg = np.float64([-10**-15, P1[1], P1[2]])
                P2_neg = np.float64([-10**-15, P2[1], P2[2]])  
                
                # Gather all vertices - it's easier to convert them to sph. coord.
                new_vertices = np.array([A, B, C,
                                P1, P2,
                                P1_neg, P2_neg])
                
                az_new, el_new, _ = cart2sph(new_vertices[:, 0], new_vertices[:, 1], new_vertices[:, 2])
                
                # --- Construct the new faces ---
                # If A is on y<0, A gets the P1/2_neg, otherwise B and C get the negatives
                #
                # Add the faces that are on the x < 0 to the rest of the model, and prepare
                # those at x>0 for further splitting
                if np.sum(x_signs > 0) == 1:    # lone positive -> 1 triangle for splitting
                    polys = [
                        Polygon([(az_new[5], el_new[5]), (az_new[1], el_new[1]), (az_new[2], el_new[2])]),
                             Polygon([(az_new[5], el_new[5]), (az_new[2], el_new[2]), (az_new[6], el_new[6])])
                             ]
                    split_faces = [np.array([A, P1, P2])]
                else:                            # lone negative -> 2 triangles for splitting
                    polys = [ 
                        Polygon([(az_new[0], el_new[0]), (az_new[5], el_new[5]), (az_new[6], el_new[6])])
                        ]
                    split_faces = [np.array([B, P1, C]), np.array([C, P1, P2])]
                face_polygons.extend(polys)
            
            
            # Loop through all faces need splitting in split_faces
            for sf in split_faces:
                # --- Split on y-axis
                # Same process as splitting on x-axis
                xyz = sf
                y_signs = np.sign(xyz[:, 1])
                crosses_y = (not np.all(y_signs == y_signs[0]))
            
                if crosses_y:
                    # --- same logic as your current y-split ---
                    if np.sum(y_signs > 0) == 1:
                        idx_a = np.where(y_signs > 0)[0][0]
                    elif np.sum(y_signs < 0) == 1:
                        idx_a = np.where(y_signs < 0)[0][0]
            
                    idx_bc = [i for i in range(3) if i != idx_a]
            
                    A = xyz[idx_a]
                    B = xyz[idx_bc[0]]
                    C = xyz[idx_bc[1]]
            
                    P1 = line_zero_y(A, B)
                    P2 = line_zero_y(A, C)
                    P1 = np.float64([P1[0], 10**-15, P1[2]])
                    P2 = np.float64([P2[0], 10**-15, P2[2]])
                    P1_neg = np.float64([P1[0], -10**-15, P1[2]])
                    P2_neg = np.float64([P2[0], -10**-15, P2[2]])
            
                    new_vertices = np.array([A, B, C, P1, P2, P1_neg, P2_neg])
                    az_new, el_new, _ = cart2sph(new_vertices[:, 0], new_vertices[:, 1], new_vertices[:, 2])
            
                    if np.sum(y_signs > 0) == 1:
                        polys = [
                            Polygon([(az_new[0], el_new[0]), (az_new[3], el_new[3]), (az_new[4], el_new[4])]),
                            Polygon([(az_new[5], el_new[5]), (az_new[1], el_new[1]), (az_new[2], el_new[2])]),
                            Polygon([(az_new[5], el_new[5]), (az_new[2], el_new[2]), (az_new[6], el_new[6])])
                        ]
                    else:
                        polys = [
                            Polygon([(az_new[0], el_new[0]), (az_new[5], el_new[5]), (az_new[6], el_new[6])]),
                            Polygon([(az_new[3], el_new[3]), (az_new[1], el_new[1]), (az_new[2], el_new[2])]),
                            Polygon([(az_new[3], el_new[3]), (az_new[2], el_new[2]), (az_new[4], el_new[4])])
                        ]
                    face_polygons.extend(polys)
            
        else:
        # no manipulation of the face in any way (removal or splitting)
            az, el, _ = cart2sph(xyz[:,0], xyz[:,1], xyz[:,2])
            face_polygons.append(Polygon([(az[0], el[0]), (az[1], el[1]), (az[2], el[2])]))


    # -------------------------------------
    # MERGE POLYGONS INTO A SINGLE SHAPE
    # -------------------------------------
    fov_faces = face_polygons
    fov_unified = unary_union(face_polygons)


    # -------------------------------------
    # OUTPUTS
    # -------------------------------------
    return {
        "faces": fov_faces,
        "unified": fov_unified
    }











