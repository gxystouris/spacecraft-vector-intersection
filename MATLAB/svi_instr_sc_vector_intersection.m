%
% DESCRIPTION
% Script for finding if a vector (with one end to an instrument) intesects
% any part of the spacecraft
%
% INPUT:
%   - model path: a character array with the path of the spacecraft model, 
%        including the extension                                            [.obj or .stl]
%        e.g. 'C:\Users\Darth_Vader\Documents\Cassini_3D_model_example\Cassini_NASA_model.obj'
%   - [x,y,z]: the instrument's location on the model                       [model units / scale]
%   - [vector_matrix]: n x 3 matrix with the vector of interest, e.g.
%       instrument-->Sun                                                    [model units / scale]
%
% OUTPUT:
%   -  index: index of whether the given vector intersects the spacecraft
%       or not; 1 if it does, 0 if it doesn't                               [0 or 1]
%   -  phi: the phi angle of each of the vectors                            [rad] -> [-pi pi]
%   -  theta: the theta angle of each of the vectors                        [rad] -> [-pi/2 pi/2]
%
% ------------------------------------
% NOTES:
% - As it takes quite some time to find the contour for a model it is to
% calculate the contour once and save it. Then every time you want to use
% the contour load the saved one.
% - The units for the instruments x,y,z should be in the same units as the
% model. e.g. if the 1 unit of the model corresponds to 1 metre, the
% instrument's x,y,z "offset" from the s/c origin should also follow this.
% - The units of the vector matrix should follow the units of the model.
%
% ------------------------------------
% Author: George Xystouris (23 March 2023)
% ------------------------------------
% Credits: George Xystouris
%          Oleg Shebanits
%          Chris Arridge
%
% This is part of the open-access published work "A simple spacecraft - vector 
% intersection methodology and applications", DOI:10.1093/rasti/rzae012
% https://academic.oup.com/rasti/article/3/1/166/7634356
% ------------------------------------
% v1
% v1.1 - replaced big part of the script with function "svi_instr_fov"

function [index, vec_phi, vec_theta ] = svi_instr_sc_vector_intersection(model_path, instr_x, instr_y, instr_z, process_faces, vector_matrix)


% CREATE THE INSTRUMENT'S FOV
% ----------------------------------
% Call the function to calculate the instrument's FOV
[sc_contour] = svi_instr_field_of_view(model_path, instr_x, instr_y, instr_z, process_faces);


% CHECK FOR VECTOR INTERSECTIONS
% ----------------------------------
% Create the output variables
index = zeros(size(vector_matrix,1),1);

% Convert the vector from cartesian to spherical
[vec_phi, vec_theta, ~] = cart2sph(vector_matrix(:,1),vector_matrix(:,2),vector_matrix(:,3));
% We want the angle increasing from 0 all the way to the point. MATLAB 
% gives the negative angle from -pi to zero. It is fixed by adding 2pi 
% to the negative angle (e.g. the point -pi/2 becomes 3pi/2 )
vec_phi(vec_phi < 0) = 2.*pi + vec_phi(vec_phi < 0);

% Check whether there is an intersection
wake_ind = isinterior(sc_contour.unified,vec_phi,vec_theta); % this gives a logical

% If there an intersection replace the 0 with 1 on the index variable
index(wake_ind) = 1;


end
