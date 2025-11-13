# -*- coding: utf-8 -*-
"""

This script retrieves the faces and vertices of a 3D model, calculates 
the normal of each face, and plots the model (in 3D) with the normals.
The model can either be an official one (e.g. the given Cassini model)
or a custom-built one.

INPUT: 
    model_path: a character array with the path of the spacecraft model 
        including the extension                            [.obj or .stl]
        e.g. 'C:/Users/Darth_Vader/Documents/Cassini_3D_model_example/Cassini_NASA_model.obj'
    normals: flag whether to plot the normals of the model [True / False]
        default: False
    viewer: viewer type to visualise the model             ["matplotlib", 'pyvista', or 'mayavi']

OUTPUT: 
    sc_model structure:
       .v: n x 3 array containing three coordinates for each vertex.
       .f: m x 3 array containing three indices of the triangle vertices.
       .n: n x 3 array containing the normal vector for each vertex.
       .c: n x 3 array containing the centroid of each polygon.

How to use:
    Firs run this script (F5 in windows), and then you can call it as
    sc_model = model_visualisation(model_path, normals=True, viewer="pyvista")

----------------------------------------------------
NOTES:
- The model mesh should consist of triangles
----------------------------------------------------

Author: George Xystouris (george.xystouris@gmail.com)
Originally Created for MATLAB in August 2020
Translated for Python in August 2025

v1.0

"""

# %% Import and environment setup

import os, sys
import trimesh # might need to be installed. Install by typing    pip install trimesh shapely    at the console. If needed to upgrade type       !{sys.executable} -m pip install --upgrade pip
import numpy as np


# %% User inputs

model_path = "give-the-path-to-the-model-here--See-input-section-above"


# %% Main script

def model_visualisation(model_path, normals=False, viewer="matplotlib"):

    # %% File check    
    # Check that file is in a valid extension (.obj or .stl). If it's not, the script stops.
    ext = os.path.splitext(model_path)[1].lower()
    if ext not in [".obj", ".stl"]:
        raise ValueError("Error in opening the model file. Fast checks:\n - File name: must end in .obj or .stl)\n - File location")
    
    
    # %% Load model 
    
    # If no error is found, load the mesh
    mesh = trimesh.load(model_path, force='mesh')
    sc_model = {
        "v": mesh.vertices,
        "f": mesh.faces,
        "n": mesh.face_normals,
        "c": mesh.triangles_center}
    
    
    # %% Visualise
    
    # select viewer
    
    if viewer.lower() == "matplotlib":
       _viewer_matplotlib(sc_model, normals)
    elif viewer.lower() == "pyvista":
       _viewer_pyvista(sc_model, normals)
    elif viewer.lower() == "mayavi":
       _viewer_mayavi(sc_model, normals)
    else:
        raise ValueError("Viewer must be 'matplotlib', 'pyvista', or 'mayavi'")

    return sc_model
    
    
    
# %% --- Viewers ---

def _viewer_matplotlib(sc, normals):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    
    v, f, n, c = sc["v"], sc["f"], sc["n"], sc["c"]
    
    fig = plt.figure(figsize=(15,5))
        
    # List of projections: (title, elev, azim)
    projections = [('YZ', 0, 0), ('XZ', 0, -90), ('XY', 90, -90)]
        
    for i, (title, elev, azim) in enumerate(projections, start=1):
        ax = fig.add_subplot(1, 3, i, projection='3d')
        
        # Add mesh
        mesh_collection = Poly3DCollection(v[f],
                          facecolors=(0.95, 0.69, 0.06, 1.0),
                          edgecolors='k',
                          linewidths=0.3,
                          alpha=1.0)
        ax.add_collection3d(mesh_collection)
        
        # Set view
        ax.view_init(elev=elev, azim=azim)
        
        # Hide perpendicular axis
        if title == 'XY':
            ax.set_zlabel('')
            ax.set_zticks([])
            ax.zaxis.line.set_color((0,0,0,0))
        elif title == 'XZ':
            ax.set_ylabel('')
            ax.set_yticks([])
            ax.yaxis.line.set_color((0,0,0,0))
        elif title == 'YZ':
            ax.set_xlabel('')
            ax.set_xticks([])
            ax.xaxis.line.set_color((0,0,0,0))
        
        # Set remaining axes labels
        if title != 'YZ': ax.set_xlabel('X')
        if title != 'XZ': ax.set_ylabel('Y')
        if title != 'XY': ax.set_zlabel('Z')
        
        # Title closer to plot
        ax.set_title(f'{title} Projection', y=0.95)
        
        # Optionally draw normals
        if normals:
            length_scale = np.mean(np.linalg.norm(v, axis=1)) * 0.05
            ax.quiver(c[:,0], c[:,1], c[:,2],
                      n[:,0], n[:,1], n[:,2],
                      color='r', length=length_scale, normalize=True, linewidth=0.25)
    
    plt.tight_layout()
    plt.show()



def _viewer_pyvista(sc, normals):
    
    # check if pyvista is installed - if not install it
    try:
        import pyvista as pv
    except ImportError:
        print("PyVista not found. Attempting to install...")
        import subprocess
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install", "pyvista"])
            print("PyVista installed successfully.")
            import pyvista as pv  # try importing again
        except Exception as e:
            print("Could not install PyVista. Error:", e)
            print("Falling back to matplotlib viewer.")
            return _viewer_matplotlib(sc, normals)

    import numpy as np  # ensure np exists here
    
    v, f, n, c = sc["v"], sc["f"], sc["n"], sc["c"]
    
    # PyVista expects faces as: [3, v1, v2, v3, 3, v1, v2, v3, ...]
    faces = np.hstack([np.full((f.shape[0], 1), 3), f]).astype(np.int32)
    
    mesh = pv.PolyData(v, faces)
    mesh.compute_normals(inplace=True)
    
    p = pv.Plotter()
    p.add_mesh(mesh, color="lightgrey", show_edges=True)
    
    if normals:
        # Scale arrows to 5% of model size:
        model_size = np.linalg.norm(v.ptp(axis=0))   # bounding box diagonal length
        arrow_scale = model_size * 0.5
    
        p.add_arrows(c, n, mag=arrow_scale)
    
    p.show()




def _viewer_mayavi(sc, normals):

    # Try importing mayavi; install if missing
    try:
        from mayavi import mlab
    except ImportError:
        print("Mayavi not found. Attempting to install...")
        import subprocess
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install", "mayavi"])
            print("Mayavi installed successfully.")
            from mayavi import mlab
        except Exception as e:
            print("Could not install Mayavi. Falling back to Matplotlib. Error:", e)
            return _viewer_matplotlib(sc, normals)
    
    v, f, n, c = sc["v"], sc["f"], sc["n"], sc["c"]
    
    mlab.triangular_mesh(v[:, 0], v[:, 1], v[:, 2], f, color=(0.8, 0.8, 0.8))
    
    if normals:
        mlab.quiver3d(c[:, 0], c[:, 1], c[:, 2],
                      n[:, 0], n[:, 1], n[:, 2],
                      scale_factor=5,
                      color=(1, 0, 0))
    
    mlab.show()





























