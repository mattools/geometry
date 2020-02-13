classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) MutableTriMesh3D < TriMesh3D
% Representation of a 3D triangular mesh keeping topological information.
%
%   MESH = MutableTriMesh3D(V, F)
%
%   Example
%   MutableTriMesh3D
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
    
    % Vertex indices for each edge, as a NE-by-2 array (optional).
    % Can be empty.
    Edges = [];
    
    % Vertex indices for each face, as a NF-by-3 array.
    Faces;
    
    % Mapping of faces associated to each edge (optional).
    % updated with method 'computeEdgeFaces'
    EdgeFaces = [];
end

%% Constructor
methods
    function obj = MutableTriMesh3D(varargin)
        % Constructor for the MutableTriMesh3D class.
        
        var1 = varargin{1};
        if isnumeric(var1)
            obj.Vertices = varargin{1};
            obj.Faces = varargin{2};
            
        elseif nargin == 1 && isa(var1, 'MutableTriMesh3D')
            % Copy constructor from another MutableTriMesh3D instance.
            obj.Vertices = var1.Vertices;
            obj.Edges = var1.Edges;
            obj.Faces = var1.Faces;
            obj.EdgeFaces = var1.EdgeFaces;
            
        elseif nargin == 1 && isa(var1, 'TriMesh3D')
            % Copy constructor from another TriMesh3D instance.
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

    
    function res = subdivide(obj, n)
        % Create a finer version of the mesh by subdividing each face.
        
        % compute the edge array
        computeEdges(obj);
        nEdges = size(obj.Edges, 1);
        
        % index of edges around each face
%         faceEdgeIndices = meshFaceEdges(obj.Vertices, obj.Edges, obj.Faces);
        faceEdgeIndices = TriMesh3D.computeFaceEdgeList(edges, faces);
        
        % Create new vertices on edges
        
        % several interpolated positions
        t = linspace(0, 1, n + 1)';
        coef2 = t(2:end-1);
        coef1 = 1 - t(2:end-1);
        
        % initialise the array of new vertices
        vertices2 = obj.Vertices;
        
        % keep an array containing index of new vertices for each original edge
        edgeNewVertexIndices = zeros(nEdges, n-1);
        
        % create new vertices on each edge
        for iEdge = 1:nEdges
            % extract each extremity as a point
            v1 = obj.Vertices(obj.Edges(iEdge, 1), :);
            v2 = obj.Vertices(obj.Edges(iEdge, 2), :);
            
            % compute new points
            newPoints = coef1 * v1 + coef2 * v2;
            
            % add new vertices, and keep their indices
            edgeNewVertexIndices(iEdge,:) = size(vertices2, 1) + (1:n-1);
            vertices2 = [vertices2 ; newPoints]; %#ok<AGROW>
        end
        
        
        % create array
        faces2 = zeros(0, 3);
        
        % iterate on faces of initial mesh
        nFaces = size(obj.Faces, 1);
        for iFace = 1:nFaces
            % compute index of each corner vertex
            face = obj.Faces(iFace, :);
            iv1 = face(1);
            iv2 = face(2);
            iv3 = face(3);
            
            % compute index of each edge
            faceEdges = faceEdgeIndices{iFace};
            ie1 = faceEdges(1);
            ie2 = faceEdges(2);
            ie3 = faceEdges(3);
            
            % indices of new vertices on edges
            edge1NewVertexIndices = edgeNewVertexIndices(ie1, :);
            edge2NewVertexIndices = edgeNewVertexIndices(ie2, :);
            edge3NewVertexIndices = edgeNewVertexIndices(ie3, :);
            
            % keep vertex 1 as reference for edges 1 and 3
            if obj.Edges(ie1, 1) ~= iv1
                edge1NewVertexIndices = edge1NewVertexIndices(end:-1:1);
            end
            if obj.Edges(ie3, 1) ~= iv1
                edge3NewVertexIndices = edge3NewVertexIndices(end:-1:1);
            end
            
            % create the first new face, on 'top' of the original face
            topVertexInds = [edge1NewVertexIndices(1) edge3NewVertexIndices(1)];
            newFace = [iv1 topVertexInds];
            faces2 = [faces2; newFace]; %#ok<AGROW>
            
            % iterate over middle strips
            for iStrip = 2:n-1
                % index of extreme vertices of current row
                ivr1 = edge1NewVertexIndices(iStrip);
                ivr2 = edge3NewVertexIndices(iStrip);
                
                % extreme vertices as points
                v1 = vertices2(ivr1, :);
                v2 = vertices2(ivr2, :);
                
                % create additional vertices within the bottom row of the strip
                t = linspace(0, 1, iStrip+1)';
                coef2 = t(2:end-1);
                coef1 = 1 - t(2:end-1);
                newPoints = coef1 * v1 + coef2 * v2;
                
                % compute indices of new vertices in result array
                newInds = size(vertices2, 1) + (1:iStrip-1);
                botVertexInds = [ivr1 newInds ivr2];
                
                % add new vertices
                vertices2 = [vertices2 ; newPoints]; %#ok<AGROW>
                
                % create top faces of current strip
                for k = 1:iStrip-1
                    newFace = [topVertexInds(k) botVertexInds(k+1) topVertexInds(k+1)];
                    faces2 = [faces2; newFace]; %#ok<AGROW>
                end
                
                % create bottom faces of current strip
                for k = 1:iStrip
                    newFace = [topVertexInds(k) botVertexInds(k) botVertexInds(k+1)];
                    faces2 = [faces2; newFace]; %#ok<AGROW>
                end
                
                % bottom vertices of current strip are top vertices of next strip
                topVertexInds = botVertexInds;
            end
            
            % for edge 2, keep vertex 2 of the current face as reference
            if obj.Edges(ie2, 1) ~= iv2
                edge2NewVertexIndices = edge2NewVertexIndices(end:-1:1);
            end
            
            % consider new vertices together with extremities
            botVertexInds = [iv2 edge2NewVertexIndices iv3];
            
            % create top faces for last strip
            for k = 1:n-1
                newFace = [topVertexInds(k) botVertexInds(k+1) topVertexInds(k+1)];
                faces2 = [faces2; newFace]; %#ok<AGROW>
            end
            
            % create bottom faces for last strip
            for k = 1:n
                newFace = [topVertexInds(k) botVertexInds(k) botVertexInds(k+1)];
                faces2 = [faces2; newFace]; %#ok<AGROW>
            end
        end

        % create the resulting data structure
        res = MutableTriMesh3D(vertices2, faces2);
    end
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
    function ne = edgeNumber(obj)
        % Get the number of edges in the mesh.
        
        % ne = edgeNumber(mesh)
        computeEdges(obj);
        ne = size(obj.Edges, 1);
    end
        
    function edgeList = edges(obj)
        % edgeList = edges(mesh);
        if isempty(obj.Edges)
            computeEdges(obj);
        end
        edgeList = obj.Edges;
    end
end

methods (Access = private)
    function computeEdges(obj)
        % Update the property Edges.
        
        % compute total number of edges
        % (3 edges per face)
        nFaces  = size(obj.Faces, 1);
        nEdges  = nFaces * 3;
        
        % create vertex indices for all edges (including duplicates)
        edges = zeros(nEdges, 2);
        for i = 1:nFaces
            f = obj.Faces(i, :);
            edges(((i-1)*3+1):i*3, :) = [f' f([2:end 1])'];
        end
        
        % remove duplicate edges, and sort the result
        obj.Edges = sortrows(unique(sort(edges, 2), 'rows'));
    end
    
    function edgeFaces = computeEdgeFaces(obj)
        % Update the property EdgeFaces.
        
        % ensure edge array is computed
        if isempty(obj.Edges)
            computeEdges(obj);
        end
        edges = obj.Edges;
        
        % allocate memory for result
        nEdges = size(obj.Edges, 1);
        obj.EdgeFaces = zeros(nEdges, 2);

        % iterate on faces
        nFaces = size(obj.Faces, 1);
        for iFace = 1:nFaces
            face = obj.Faces(iFace, :);
            
            % iterate on edges of current face
            for j = 1:length(face)
                % build edge: array of vertices
                j2 = mod(j, length(face)) + 1;
                
                % do not process edges with same vertices
                if face(j) == face(j2)
                    continue;
                end
                
                % vertex indices of current edge
                currentEdge = [face(j) face(j2)];
                
                % find index of current edge, assuming face is left-located
                b1 = ismember(obj.Edges, currentEdge, 'rows');
                indEdge = find(b1);
                if ~isempty(indEdge)
                    if obj.EdgeFaces(indEdge, 1) ~= 0
                        error('MutableTriMesh3D:subdivide:IllegalTopology', ...
                            'Two faces were found on left side of edge %d ', indEdge);
                    end
                    
                    obj.EdgeFaces(indEdge, 1) = iFace;
                    continue;
                end
                
                % otherwise, assume the face is right-located
                b2 = ismember(edges, currentEdge([2 1]), 'rows');
                indEdge = find(b2);
                if ~isempty(indEdge)
                    if obj.EdgeFaces(indEdge, 2) ~= 0
                        error('MutableTriMesh3D:subdivide:IllegalTopology', ...
                            'Two faces were found on left side of edge %d ', indEdge);
                    end
                    
                    obj.EdgeFaces(indEdge, 2) = iFace;
                    continue;
                end
                
                % If face was neither left nor right, error
                warning('MutableTriMesh3D:subdivide:IllegalTopology', ...
                    'Edge %d of face %d was not found in edge array', ...
                    j, iFace);
                continue;
            end
        end
        
        edgeFaces = obj.EdgeFaces;
    end
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
        str = struct('Type', 'MutableTriMesh3D', ...
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
        mesh = MutableTriMesh3D(str);
    end
end

end