classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) SimpleTriMesh3D < TriMesh3D
% Simple class for 3D triangular mesh with limited edition possiblities.
%
%   MESH = SimpleTriMesh3D(V, F)
%
%   Example
%   SimpleTriMesh3D
%
%   See also
%     Mesh3D
 
% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2019-02-07,    using Matlab 9.4.0.813654 (R2018a)
% Copyright 2018 INRAE - BIS - BIBS.

properties
    % Coordinates of vertices, as a NV-by-3 array.
    Vertices;
    
    % Vertex indices for each face, as a NF-by-3 array.
    Faces;
    
end

%% Constructor
methods
    function obj = SimpleTriMesh3D(varargin)
        % Constructor for the SimpleTriMesh3D class.
        
        var1 = varargin{1};
        if isnumeric(var1)
            obj.Vertices = varargin{1};
            obj.Faces = varargin{2};
            
        elseif nargin == 1 && isa(var1, 'SimpleTriMesh3D')
            % Copy constructor from another SimpleTriMesh3D instance.
            obj.Vertices = var1.Vertices;
            obj.Faces = var1.Faces;
            
        elseif nargin == 1 && isa(var1, 'SimpleTriMesh3D')
            % Copy constructor from another SimpleTriMesh3D instance.
            obj.Vertices = vertexPositions(var1);
            obj.Faces = faceVertexIndices(var1);
            
        elseif isstruct(var1)
            % Copy constructor from a structure.
            obj.Vertices = var1.vertices;
            obj.Faces = var1.faces;
        end
        
    end
end


%% Global procesing of mesh
methods
end


%% Vertex management methods
methods
    function nv = vertexNumber(obj)
        % Get the number of vertices in the mesh.
        nv = size(obj.Vertices, 1);
    end
end

%% Edge management methods
methods
end


%% Face management methods
methods
    function nf = faceNumber(obj)
        % Get the number of faces in the mesh.
        nf = size(obj.Faces, 1);
    end

    % TODO: return 3D polygon
    function poly = facePolygon(obj, ind)
        poly = obj.Vertices(obj.Faces(ind, :), :);
    end
end


%% Methods implementing Geometry3D
methods

end % end methods


%% Access methods
methods
    function v = vertexPositions(obj)
        % Return the Nv-by-3 array of vertex coordinates.
        v = obj.Vertices;
    end
    
    function f = faceVertexIndices(obj)
        % Return the Nf-by-3 array of face vertex indices.
        f = obj.Faces;
    end
end

%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'SimpleTriMesh3D', ...
            'Vertices', obj.Vertices, ...
            'Faces', obj.Faces);
    end
end
methods (Static)
    function mesh = fromStruct(str)
        % Create a new instance from a structure.
        if ~(isfield(str, 'Vertices') && isfield(str, 'Faces'))
            error('Requires fields vertices and faces');
        end
        if size(str.Faces, 2) ~= 3
            error('Requires a triangular face array');
        end
        mesh = SimpleTriMesh3D(str);
    end
end

end