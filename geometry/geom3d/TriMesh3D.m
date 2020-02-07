classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) TriMesh3D < Mesh3D
% Abstract class for representing a 3D triangular mesh.
%
%   Parent class for all implementations representing 3D triangular meshes.
%
%   Example
%   TriMesh3D
%
%   See also
%     Mesh3D
 
% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2019-02-07,    using Matlab 9.4.0.813654 (R2018a)
% Copyright 2018 INRAE - BIS - BIBS.

properties
end

%% Static factories
methods (Static)
    function mesh = create(vertices, faces)
        % Create a new trimesh from vertex and face arrays.
        % 
        % Usage:
        %   mesh = TriMesh3D.create(vertices, faces);
        %
        mesh = SimpleTriMesh3D(vertices, faces);
    end
end

%% Constructor
methods (Access = protected)
    function obj = TriMesh3D(varargin)
        % Constructor for the TriMesh3D class.
        % (called by derived class constructors)
    end
end

%% Drawing functions
methods
    function h = drawFaceNormals(obj, varargin)
        % Draw the normal of each face.
        %
        % h = drawFaceNormals(mesh);
        pts = faceCentroids(obj);
        pos = pts.Coords;
        vn = faceNormals(obj);
        h = quiver3(pos(:, 1), pos(:, 2), pos(:, 3), ...
            vn(:, 1), vn(:, 2), vn(:, 3), 0, varargin{:});
    end
end

%% Global procesing of mesh
methods
    function res = smooth(obj, varargin)
        %SMOOTH Smooth a mesh.
        %
        % mesh2 = smooth(mesh, nIter);
        % Smoothes the mesh by replacing each vertex by the average of its
        % neighbors.
        
        % determine number of iterations
        nIter = 1;
        if ~isempty(varargin)
            nIter = varargin{1};
        end
        
        % compute adjacency matrix,
        % result is a Nv-by-Nv matrix with zeros on the diagonal
        adj = vertexAdjacencyMatrix(obj);
        
        % Add "self adjacencies"
        nv = vertexNumber(obj);
        adj = adj + speye(nv);
        
        % weight each vertex by the number of its neighbors
        w = spdiags(full(sum(adj, 2).^(-1)), 0, nv, nv);
        adj = w * adj;
        
        % do averaging to smooth the field
        v2 = vertexPositions(obj);
        for k = 1:nIter
            v2 = adj * v2;
        end
        
        % return new TriMesh
        res = TriMesh3D.create(v2, obj.Faces);
    end
    
    function res = subdivide(obj, varargin)
        % Create a finer version of the mesh by subdividing each face.
        
        % determine number of subdivisions
        n = 2;
        if ~isempty(varargin)
            n = varargin{1};
        end

        % Extract vertex and faces mesh information
        vertices  = vertexPositions(obj);
        faces = faceVertexIndices(obj);

        % compute the edge array
        edges = edgeVertexIndices(obj);
        nEdges = size(edges, 1);
        
        % index of edges around each face
        faceEdgeIndices = TriMesh3D.computeFaceEdgeList(edges, faces);
        
        
        % Create new vertices on edges
        
        % several interpolated positions
        t = linspace(0, 1, n + 1)';
        coef2 = t(2:end-1);
        coef1 = 1 - t(2:end-1);
        
        % initialise the array of new vertices
        vertices2 = vertices;
        
        % keep an array containing index of new vertices for each original edge
        edgeNewVertexIndices = zeros(nEdges, n-1);
        
        % create new vertices on each edge
        for iEdge = 1:nEdges
            % extract each extremity as a point
            v1 = vertices(edges(iEdge, 1), :);
            v2 = vertices(edges(iEdge, 2), :);
            
            % compute new points
            newPoints = coef1 * v1 + coef2 * v2;
            
            % add new vertices, and keep their indices
            edgeNewVertexIndices(iEdge,:) = size(vertices2, 1) + (1:n-1);
            vertices2 = [vertices2 ; newPoints]; %#ok<AGROW>
        end
        
        
        % create new face array
        faces2 = zeros(0, 3);
        
        % iterate on faces of initial mesh
        nFaces = size(faces, 1);
        for iFace = 1:nFaces
            % compute index of each corner vertex
            face = faces(iFace, :);
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
            if edges(ie1, 1) ~= iv1
                edge1NewVertexIndices = edge1NewVertexIndices(end:-1:1);
            end
            if edges(ie3, 1) ~= iv1
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
            if edges(ie2, 1) ~= iv2
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
        res = TriMesh3D.create(vertices2, faces2);
    end
end

methods (Static, Access = protected)
    function FE = computeFaceEdgeList(edges, faces)
        % Computes list of edge indices for each face.
        %
        %   FE = computeFaceEdgeList(E, F)
        %   Returns a 1-by-NF cell array containing for each face, the set of edge
        %   indices corresponding to adjacent edges.
        %
        %   Example
        %   meshFaceEdges
        %
        %   See also
        %     meshes3d, meshEdgeFaces
        
        % ------
        % Author: David Legland
        % e-mail: david.legland@inra.fr
        % Created: 2013-08-22,    using Matlab 7.9.0.529 (R2009b)
        % Copyright 2013 INRA - Cepia Software Platform.
        
        nFaces = size(faces, 1);
        
        FE = cell(nFaces, 1);
        
        % impose ordering of edge indices
        edges = sort(edges, 2);
        
        for iFace = 1:nFaces
            % extract vertex indices of current face
            face = meshFace(faces, iFace);
            nv = length(face);
            
            % for each couple of adjacent vertices, find the index of the matching
            % row in the edges array
            fei = zeros(1, nv);
            for iEdge = 1:nv
                % compute index of each edge vertex
                edge = sort([face(iEdge) face(mod(iEdge, nv) + 1)]);
                v1 = edge(1);
                v2 = edge(2);
                
                % find the matching row
                ind = find(edges(:,1) == v1 & edges(:,2) == v2);
                fei(iEdge) = ind;
                
            end
            FE{iFace} = fei;
        end
    end
end

%% Geometric information about mesh
methods
    function vol = volume(obj)
        % (signed) volume enclosed by this mesh.
        %
        % See Also
        %   surfaceArea

        % get mesh data
        vertices = vertexPositions(obj);
        faces = faceVertexIndices(obj);

        % initialize an array of volume
        nFaces = size(faces, 1);
        vols = zeros(nFaces, 1);

        % Shift all vertices to the mesh centroid
        centroid = mean(vertices, 1);
        
        % compute volume of each tetraedron
        for iFace = 1:nFaces
            % consider the tetrahedron formed by face and mesh centroid
            tetra = vertices(faces(iFace, :), :);
            tetra = bsxfun(@minus, tetra, centroid);
            
            % volume of current tetrahedron
            vols(iFace) = det(tetra) / 6;
        end
        
        vol = sum(vols);
    end
    
    function area = surfaceArea(obj)
        % Surface area of this mesh, obtained by summing face areas.
        %
        % See Also
        %   volume
        
        % get mesh data
        vertices = vertexPositions(obj);
        faces = faceVertexIndices(obj);

        % compute two direction vectors of each trinagular face, using the
        % first vertex of each face as origin
        v1 = vertices(faces(:, 2), :) - vertices(faces(:, 1), :);
        v2 = vertices(faces(:, 3), :) - vertices(faces(:, 1), :);
        
        % area of each triangle is half the cross product norm
        % see also crossProduct3d in MatGeom
        vn = zeros(size(v1));
        vn(:) = bsxfun(@times, v1(:,[2 3 1],:), v2(:,[3 1 2],:)) - ...
                bsxfun(@times, v2(:,[2 3 1],:), v1(:,[3 1 2],:));
        vn = sqrt(sum(vn .* vn, 2));
        
        % sum up and normalize
        area = sum(vn) / 2;
    end
    
%     function mb = meanBreadth(obj)
%         % Mean breadth of this mesh
%         % Mean breadth is proportionnal to the integral of mean curvature
%         %
%         % See Also
%         %   trimeshMeanBreadth
%         
%         mb = trimeshMeanBreadth(obj.Vertices, obj.Faces);
%     end
end


%% Vertex management methods
methods
    
    function verts = vertices(obj)
        % Return vertices in the mesh as a MultiPoint3D.
        verts = MultiPoint3D(vertexPositions(obj));
    end
    
    function adj = vertexAdjacencyMatrix(obj)
        % Get the adjacency matrix of mesh vertices.
        
        % forces faces to be floating point array, for sparse function
        faces = faceVertexIndices(obj);
        if ~isfloat(faces)
            faces = double(faces);
        end
        
        % populate a sparse matrix
        adj = sparse(...
            [faces(:,1); faces(:,1); faces(:,2); faces(:,2); faces(:,3); faces(:,3)], ...
            [faces(:,3); faces(:,2); faces(:,1); faces(:,3); faces(:,2); faces(:,1)], ...
            1.0);
        
        % remove double adjacencies
        adj = min(adj, 1);
        
        % ensure the size of the matrix is Nv-by-Nv
        % (this may happen if some vertices are not referenced)
        nv = vertexNumber(obj);
        if size(adj, 1) < nv
            adj(nv, nv) = 0;
        end
    end
end

%% Edge management methods
methods
    function ne = edgeNumber(obj)
        % Get the number of edges in the mesh.
        %
        % See Also
        %   edgeVertexIndices
        ne = size(edgeVertexIndices(obj), 1);
    end
    
    function edges = edgeVertexIndices(obj)
        % Compute edge vertex indices.
        
        % get mesh face vertex indices
        faces = faceVertexIndices(obj);
        
        % create vertex indices for all edges (including duplicates)
        edges = [...
            faces(:,1) faces(:,2) ; ...
            faces(:,2) faces(:,3) ; ...
            faces(:,3) faces(:,1)];
        
        % remove duplicate edges, and sort the result
        edges = sortrows(unique(sort(edges, 2), 'rows'));
    end
end


%% Face management methods
methods
    function normals = faceNormals(obj, inds)
        % Compute the normal vector to each face.
        %
        % vn = faceNormals(mesh);
        
        vertices = vertexPositions(obj);
        faces = faceVertexIndices(obj);
        
        nf = size(faces, 1);
        if nargin == 1
            inds = 1:nf;
        end

        % compute vector of each edge
        v1 = vertices(faces(inds,2),:) - vertices(faces(inds,1),:);
        v2 = vertices(faces(inds,3),:) - vertices(faces(inds,1),:);

        % compute normals using cross product (vectors have same size)
        normals = cross(v1, v2, 2);
    end
    
    function pts = faceCentroids(obj, inds)
        % Compute the centroid of each face.
        %
        % pts = faceCentroids(mesh);
        
        % determine indices of face to process
        nf = faceNumber(obj);
        if nargin == 1
            inds = 1:nf;
        end
        
        % get mesh data
        vertices = vertexPositions(obj);
        faces = faceVertexIndices(obj);
        
        % process centroids in parallel
        pts = zeros(length(inds), 3);
        for i = 1:3
            pts = pts + vertices(faces(inds,i),:) / 3;
        end
        
        % create Geometry object for storing result
        pts = MultiPoint3D(pts);
    end

    function poly = facePolygon(obj, ind)
        poly = obj.Vertices(obj.Faces(ind, :), :);
    end
end


%% Methods implementing Geometry3D
methods
    function box = boundingBox(obj)
        % Return the bounding box of this mesh.
        mini = min(vertexPositions(obj));
        maxi = max(vertexPositions(obj));
        box = Box3D([mini(1) maxi(1) mini(2) maxi(2) mini(3) maxi(3)]);
    end
    
    function h = draw(varargin)
        %DRAW Draw the current mesh, eventually specifying the style.
        
        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});

        % add default drawing options
        options = {'FaceColor', [.7 .7 .7]};

        % extract optional drawing options
        if length(varargin) > 1 && ischar(varargin{1})
            options = [options varargin];
        end
        
        h = patch('Parent', ax, ...
            'vertices', vertexPositions(obj), ...
            'faces', faceVertexIndices(obj), ...
            options{:} );

        % optionnally add style processing
        if ~isempty(style)
            apply(style, hh);
        end
                
        if nargout > 0
            h = hh;
        end
    end
    
    function res = transform(obj, transfo)
        % Apply a transform to this mesh.
        vt = transformPoint(transfo, vertexPositions(obj));
        res = TriMesh3D.create(vt, obj.Faces);
    end
    
    function res = scale(obj, varargin)
        % Return a scaled version of this mesh.
        factor = varargin{1};
        res = TriMesh3D.create(vertexPositions(obj) * factor, obj.Faces);
    end
    
    function res = translate(obj, varargin)
        % Return a translated version of this mesh.
        shift = varargin{1};
        res = TriMesh3D.create(bsxfun(@plus, vertexPositions(obj), shift), obj.Faces);
    end
    
end % end methods

%% Access methods
% Declares new methods.
methods (Abstract)
    % Return the Nv-by-3 array of vertex coordinates.
    v = vertexPositions(obj);
    
    % Return the Nf-by-3 array of face vertex indices.
    f = faceVertexIndices(obj);
end


%% Serialization methods
% Uses default implementation that creates / save only vertices and faces.
% Specialized implementation may manage more data.

methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'TriMesh3D', ...
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
        mesh = TriMesh3D.create(str.Vertices, str.Faces);
    end
end

end