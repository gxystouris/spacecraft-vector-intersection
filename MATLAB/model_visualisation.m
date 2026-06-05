%
%   This script retrieves the faces and vertices of a 3D model, calculates 
%   the normal of each face, and plots the model (in 3D) with the normals.
%   The model can either be an official one (e.g. the given Cassini model)
%   or a custom-built one.
%
%   INPUT: 
%     model path: a character array with the path of the spacecraft model 
%        including the extension                            [.obj or .stl]
%        e.g. 'C:\Users\Darth_Vader\Documents\Cassini_3D_model_example\Cassini_NASA_model.obj'
%     normals: flag whether to plot the normals of the model [char - case insensitive]
%        'normals' or 'on' will plot the normals
%
%   OUTPUT: 
%     sc_model structure:
%       .v: n x 3 array containing three coordinates for each vertex.
%       .f: m x 3 array containing three indices of the triangle vertices.
%       .n: n x 3 array containing the normal vector for each vertex.
%       .c: n x 3 array containing the centroid of each polygon.
%   
% ----------------------------------------------------
%   NOTES:
%   - The model should consist of triangles
% ----------------------------------------------------
%
%   Author: George Xystouris (20 August 2020) 
%   based on cassiniModel.m by Chris Arridge
%
% ----------------------------------------------------
% v1
% v1.1 - added option for showing normals or not

function sc_model = model_visualisation(model_path, normals)

% Check that file is in a valid file name (.obj or .stl). If it's not, the script stops.
if ~strcmp(model_path(end-3:end),'.obj') && ~strcmp(model_path(end-3:end),'.stl')
    fprintf('Error in opening the model file. Fast checks:\n -File name (the extension .obj or .stl must be part of the file)\n -File location\n')
    return
end


% CREATE THE MODEL (face, vertices, and normals)
% ----------------------------------
if strcmp(model_path(end-3:end),'.obj')
    [sc_model.v, sc_model.f] = readOBJ(model_path);
% Uncomment and use this "else" with a custom stlread function if using earlier versions than R2018b
%  [sc_model.f, sc_model.v, sc_model.n]  = stlread(model_path);
% -------------
% Comment the next four lines if using the custom stlread function
    sc_model_tr  = stlread(model_path);
    sc_model.v = sc_model_tr.Points;
    sc_model.f = sc_model_tr.ConnectivityList;
    sc_model.n = faceNormal(sc_model_tr);
end    


% Compute normal and centroid of each face
% Each centroid is computed in order to draw the corresponding normal on it
faces_num = size(sc_model.f,1); % number of faces

normal = zeros(faces_num,3);
for i=1:faces_num
	edge1 = [sc_model.v(sc_model.f(i,2),1)-sc_model.v(sc_model.f(i,1),1), ...
        sc_model.v(sc_model.f(i,2),2)-sc_model.v(sc_model.f(i,1),2), ...
        sc_model.v(sc_model.f(i,2),3)-sc_model.v(sc_model.f(i,1),3)];
    edge2 = [sc_model.v(sc_model.f(i,3),1)-sc_model.v(sc_model.f(i,1),1), ...
        sc_model.v(sc_model.f(i,3),2)-sc_model.v(sc_model.f(i,1),2), ...
        sc_model.v(sc_model.f(i,3),3)-sc_model.v(sc_model.f(i,1),3)];
    normal(i,:) = cross(edge1,edge2);
    normal(i,:) = normal(i,:)./sqrt(sum(normal(i,:).^2));
end
sc_model.n = normal;

centroid = zeros(faces_num,3);
for i=1:faces_num
	centroid(i,1) = mean(sc_model.v(sc_model.f(i,:),1));
    centroid(i,2) = mean(sc_model.v(sc_model.f(i,:),2));
    centroid(i,3) = mean(sc_model.v(sc_model.f(i,:),3));
end
sc_model.c = centroid;


% PLOT THE MODEL
% ----------------------------------
figure;
patch('Faces',sc_model.f,'Vertices',sc_model.v,'FaceColor',[0.95 0.69 0.06],'EdgeColor','black')
hold on;
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal; view(90,180);
if strcmpi(normals,'normals') | strcmpi(normals,'on')
    % add the normalised vectors on each face
    quiver3(sc_model.c(:,1), sc_model.c(:,2), sc_model.c(:,3), sc_model.n(:,1), sc_model.n(:,2), sc_model.n(:,3), 'r')
end

end




