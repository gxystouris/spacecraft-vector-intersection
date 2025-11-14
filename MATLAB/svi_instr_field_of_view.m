%
% DESCRIPTION
% Script for finding an instrument's field-of-view (FOV).
%
% INPUT:
%   - spacecraft model path: a character array with the path of the 
%        spacecraft model including the extension [.obj or .stl]
%        e.g. 'C:\Users\Darth_Vader\Documents\Cassini_3D_model_example\Cassini_NASA_model.obj'
%   - [x,y,z]: the instrument's location on the model     [model units / scale]
%   - process_faces: flag on how to process the "problematic" faces [case insensitive string]
%        'none':    no process is taking place
%        'remove':  all the faces of interest are removed
%        'split':   all the faces of interest are splitted into faces in 
%                   negative y and in postivie y. Those faces are added 
%                   to the model.
%
% OUTPUT:
%   - fov structure:
%       .faces:     the individual faces as polyshapes
%       .unified:   the contour of all the faces
%           ***SHAPE LIMITS: azimuth: [0 2pi]; elevation: [-pi/2 pi/2]
%
%------------------------------------
% NOTES:
% - As it takes quite some time to find the contour for a model it is to
% calculate the contour once and save it. Then every time you want to use
% the contour load the saved one.
% - The units for the instruments x,y,z should be in the same units as the
% model. e.g. if the 1 unit of the model corresponds to 1 metre, the
% instrument's x,y,z "offset" from the s/c origin should also follow this.
% - Problematic faces are faces that cross the 2pi->0 line, i.e. going from
% negative to positive y -or vice versa- on the model. The issue is that 
% in those cases, e.g. a face going crossing the 2pi->0 line MATLAB cannot 
% understand that the face must be cut into two parts: one for the
% negative y and one for the positive y. Instead it stretches the face
% from negative y to positive y the "other" way around.
%   E.g. if the face has an azimuth of 20 degrees and crosses the 2pi->0 
%   line, MATLAB will create a face of 340 degrees, i.e. what's missing 
%   from a circle.
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
% v1.1 - added feature for manipulating "problematic" vertices

function [fov] = svi_instr_field_of_view(model_path, instr_x, instr_y, instr_z, process_faces)

% -----------------
% CHECKS:
% -----------------

% Check that file is in a valid file name (.obj or .stl). If it's not, the script stops.
if ~strcmp(model_path(end-3:end),'.obj') && ~strcmp(model_path(end-3:end),'.stl')
    fprintf('Error in opening the model file. Fast checks:\n -File name (the extension .obj or .stl must be part of the file)\n -File location\n');
    return
end


% -----------------
% MAIN SCRIPT:
% -----------------

% Create the model (face, vertices, and normals)
if strcmp(model_path(end-3:end),'.obj')
    [sc_model.v, sc_model.f, ~, ~, sc_model.n,~] = readOBJ(model_path);
else strcmp(model_path(end-3:end),'.stl');
    [sc_model.f, sc_model.v, sc_model.n]  = stlread(model_path);
end

% Shift the model to have the instrument at the origin (0,0,0)
sc_shifted = sc_model;
sc_shifted.v = [sc_shifted.v(:,1)-instr_x  sc_shifted.v(:,2)-instr_y  sc_shifted.v(:,3)-instr_z];

% Convert from cartesian to spherical coordinates
[Vaz, Vel, ~] = cart2sph(sc_shifted.v(:,1), sc_shifted.v(:,2), sc_shifted.v(:,3));
% We want the angle from 0 all the way to the point. MATLAB gives the minus
% from -pi to zero. It is fixed by adding 2pi to the negative angle (e.g.
% -pi/2 -> 3pi/2)
Vaz(sc_shifted.v(:,2) < 0) = 2.*pi + Vaz(sc_shifted.v(:,2) < 0); 


% Assign the vertices for each face in the new coordinates
face = [];
for i_face = 1:length(sc_shifted.f)
    %
    % REMOVE CROSSING FACES
    if strcmpi(process_faces,'remove')
        % You proceed only if the condition does not apply. If it does, you
        % do nothing.
        if ~((sc_shifted.v(sc_shifted.f(i_face,1),1) > 0 && sc_shifted.v(sc_shifted.f(i_face,2),1) > 0 && sc_shifted.v(sc_shifted.f(i_face,3),1) > 0 )  & ... % x>0 condition
                ~isequal( sign(sc_shifted.v(sc_shifted.f(i_face,1),2)), sign(sc_shifted.v(sc_shifted.f(i_face,2),2)), sign(sc_shifted.v(sc_shifted.f(i_face,3),2)))  ) % different sign for y condition
            face = [face; polyshape( [Vaz(sc_shifted.f(i_face,1)), Vaz(sc_shifted.f(i_face,2)), Vaz(sc_shifted.f(i_face,3))] , ...
                [Vel(sc_shifted.f(i_face,1)), Vel(sc_shifted.f(i_face,2)), Vel(sc_shifted.f(i_face,3))] ) ];
        end
    %
    % SPLIT CROSSING FACES
    elseif strcmpi(process_faces,'split')
        if (sc_shifted.v(sc_shifted.f(i_face,1),1) > 0 && sc_shifted.v(sc_shifted.f(i_face,2),1) > 0 && sc_shifted.v(sc_shifted.f(i_face,3),1) > 0 )  & ... % x>0 condition
                ~isequal( sign(sc_shifted.v(sc_shifted.f(i_face,1),2)), sign(sc_shifted.v(sc_shifted.f(i_face,2),2)), sign(sc_shifted.v(sc_shifted.f(i_face,3),2))) % different sign for y condition
            % Each face is a triangle. We're rearranging the vertices: the subscript a is the
            % vertice that's on the other side. b and c are on the same side
            %  Find which face is on the "other" side
            if isequal( sign(sc_shifted.v(sc_shifted.f(i_face,1),2)), sign(sc_shifted.v(sc_shifted.f(i_face,2),2)) ) % if 1st and 2nd vertices have the same sign, the 3rd is on the other side
                x_a = sc_shifted.v(sc_shifted.f(i_face,3),1);
                y_a = sc_shifted.v(sc_shifted.f(i_face,3),2);
                z_a = sc_shifted.v(sc_shifted.f(i_face,3),3);
                %
                x_b = sc_shifted.v(sc_shifted.f(i_face,1),1);
                y_b = sc_shifted.v(sc_shifted.f(i_face,1),2);
                z_b = sc_shifted.v(sc_shifted.f(i_face,1),3);
                %
                x_c = sc_shifted.v(sc_shifted.f(i_face,2),1);
                y_c = sc_shifted.v(sc_shifted.f(i_face,2),2);
                z_c = sc_shifted.v(sc_shifted.f(i_face,2),3);
            %
            elseif isequal( sign(sc_shifted.v(sc_shifted.f(i_face,1),2)), sign(sc_shifted.v(sc_shifted.f(i_face,3),2)) ) % if 1st and 3rd vertices have the same sign, the 2rd is on the other side
                x_a = sc_shifted.v(sc_shifted.f(i_face,2),1);
                y_a = sc_shifted.v(sc_shifted.f(i_face,2),2);
                z_a = sc_shifted.v(sc_shifted.f(i_face,2),3);
                %
                x_b = sc_shifted.v(sc_shifted.f(i_face,1),1);
                y_b = sc_shifted.v(sc_shifted.f(i_face,1),2);
                z_b = sc_shifted.v(sc_shifted.f(i_face,1),3);
                %
                x_c = sc_shifted.v(sc_shifted.f(i_face,3),1);
                y_c = sc_shifted.v(sc_shifted.f(i_face,3),2);
                z_c = sc_shifted.v(sc_shifted.f(i_face,3),3);
            %
            else % if none of the above then the 1st is on the other side
                x_a = sc_shifted.v(sc_shifted.f(i_face,1),1);
                y_a = sc_shifted.v(sc_shifted.f(i_face,1),2);
                z_a = sc_shifted.v(sc_shifted.f(i_face,1),3);
                %
                x_b = sc_shifted.v(sc_shifted.f(i_face,2),1);
                y_b = sc_shifted.v(sc_shifted.f(i_face,2),2);
                z_b = sc_shifted.v(sc_shifted.f(i_face,2),3);
                %
                x_c = sc_shifted.v(sc_shifted.f(i_face,3),1);
                y_c = sc_shifted.v(sc_shifted.f(i_face,3),2);
                z_c = sc_shifted.v(sc_shifted.f(i_face,3),3);
            end
            % We're parametrising the line equations from A->B and A->C 
            % introducing a parameter t, i.e. x(t) = (x2-x1)t + x1 etc. We
            % find t(y=0)=0, i.e. we're finding the t where y=0
            % Find the crossing points: c1 for A->B and c2 for A->C
            t_c1 = roots([(y_b-y_a),y_a]);
            t_c2 = roots([(y_c-y_a),y_a]);
            % Find the coordinates for x and z at c1 and c2. y = 0 for both
            x_c1 = (x_b - x_a)*t_c1 + x_a;
            z_c1 = (z_b - z_a)*t_c1 + z_a;
            x_c2 = (x_c - x_a)*t_c2 + x_a;
            z_c2 = (z_c - z_a)*t_c2 + z_a;
            y_c1 = 0; y_c2 = 0;
            % Define the three new triangles: one from the lone point to
            % the x-line and two on the quadrilateral on the other side.
            % We are defining it with the same way as before: picking
            % vertices for each face from a vertice matrix. This way we
            % don't get duplicate vertices.
            % For faces that approach from the 2pi: instead of having the y
            % to zero we it very very small negative number - otherwise we
            % would face the problem mentioned in the notes
            new_vertices = [x_a y_a z_a; x_b y_b z_b; x_c y_c z_c; ...
                            x_c1 y_c1 z_c1; x_c2 y_c2 z_c2;...      % this is for vertices that approach the zero from positive y
                            x_c1 -10e-15 z_c1; x_c2 -10e-15 z_c2;]; % this is for vertices that approach the zero from negative y
            [new_vertices_Vaz, new_vertices_Vel, ~] = cart2sph(new_vertices(:,1), new_vertices(:,2), new_vertices(:,3));
            new_vertices_Vaz(new_vertices(:,2) < 0) = 2.*pi + new_vertices_Vaz(new_vertices(:,2) < 0);
            % add the new faces to the faces matrix. For the faces that 
            % are in negative y, the faces are ending at 2pi instead of 0
            if y_a > 0
                face = [face; ...
                    polyshape( [new_vertices_Vaz(1), new_vertices_Vaz(4), new_vertices_Vaz(5),], [new_vertices_Vel(1), new_vertices_Vel(4), new_vertices_Vel(5)] ); ... % A-C1-C2
                    polyshape( [new_vertices_Vaz(6), new_vertices_Vaz(2), new_vertices_Vaz(3),], [new_vertices_Vel(6), new_vertices_Vel(2), new_vertices_Vel(3)] ); ... % C1-B-C
                    polyshape( [new_vertices_Vaz(6), new_vertices_Vaz(3), new_vertices_Vaz(7),], [new_vertices_Vel(6), new_vertices_Vel(3), new_vertices_Vel(7)] )];    % C1-C-C2
            else
                face = [face; ...
                    polyshape( [new_vertices_Vaz(1), new_vertices_Vaz(6), new_vertices_Vaz(7),], [new_vertices_Vel(1), new_vertices_Vel(6), new_vertices_Vel(7)] ); ... % A-C1-C2
                    polyshape( [new_vertices_Vaz(4), new_vertices_Vaz(2), new_vertices_Vaz(3),], [new_vertices_Vel(4), new_vertices_Vel(2), new_vertices_Vel(3)] ); ... % C1-B-C
                    polyshape( [new_vertices_Vaz(4), new_vertices_Vaz(3), new_vertices_Vaz(5),], [new_vertices_Vel(4), new_vertices_Vel(3), new_vertices_Vel(5)] )];    % C1-C-C2
            end
        else   % If the split condition is not satisfied (no faces are crossing the y=0 line)
            face = [face; polyshape( [Vaz(sc_shifted.f(i_face,1)), Vaz(sc_shifted.f(i_face,2)), Vaz(sc_shifted.f(i_face,3))] , ...
                [Vel(sc_shifted.f(i_face,1)), Vel(sc_shifted.f(i_face,2)), Vel(sc_shifted.f(i_face,3))] ) ];
        end
    %
    % NO SPLITTING OR REMOVING
    else
        face = [face; polyshape( [Vaz(sc_shifted.f(i_face,1)), Vaz(sc_shifted.f(i_face,2)), Vaz(sc_shifted.f(i_face,3))] , ...
            [Vel(sc_shifted.f(i_face,1)), Vel(sc_shifted.f(i_face,2)), Vel(sc_shifted.f(i_face,3))] ) ];
    end
end


% -----------------
% SAVE VARIABLES:
% -----------------
fov.faces = face; % save all the multiple-polygons plot
fov.unified = union(face); % save the contour of the multiple-polygons plot


end





