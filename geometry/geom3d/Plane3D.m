classdef Plane3D < Geometry3D
% A plane in 3-dimensional space.
%
%   Class Plane3D
%
%   Example
%   Plane3D
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2021-11-03,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE - BIA-BIBS.



%% Static factories
methods (Static)
    function plane = medianPlane(p1, p2)
        % Create the median plane between two points.
        %
        % PLANE = Plane3D.medianPlane(PT1, PT2);
        
        % middle point
        p0 = Point3D.centroid(p1, p2);
        % normal to plane
        n = normalize(Vector3D(p1, p2));
        % create plane from point and normal
        plane = Plane3D(p0, n);
    end
    
    function plane = XY()
        % The 3D plane through the origin, containing the Ox and Oy axes.
        plane = Plane3D([0 0 0], [1 0 0], [0 1 0]);
    end
    
    function plane = XZ()
        % The 3D plane through the origin, containing the Ox and Oz axes.
        plane = Plane3D([0 0 0], [1 0 0], [0 0 1]);
    end
    
    function plane = YZ()
        % The 3D plane through the origin, containing the Oy and Oz axes.
        plane = Plane3D([0 0 0], [0 1 0], [0 0 1]);
    end
end


%% Properties
properties
    % The origin of the plane, as a 1-by-3 row vector.
    Origin;
    % The first direction vector of the plane, as a 1-by-3 row vector.
    Direction1;
    % The second direction vector of the plane, as a 1-by-3 row vector.
    Direction2;
    
end % end properties


%% Constructor
methods
    function obj = Plane3D(varargin)
        % Constructor for Plane3D class.
        
        if nargin == 0
            % Empty constructor
            
        elseif nargin == 2
            % Build from position and normal vector
            p0 = varargin{1};
            n = varargin{2};
            if isa(p0, 'Point3D') && isa(n, 'Vector3D')
                % arguments given as geometry instances
                v0 = Vector3D([1 0 0]);
                if norm(crossProduct(n, v0)) < 1e-12
                    v0 = Vector3D([0 1 0]);
                end
                
                % create direction vectors
                v1 = normalize(crossProduct(n, v0));
                v2 = -normalize(crossProduct(v1, n));
                
                % store result in inner properties
                obj.Origin = [p0.X p0.Y p0.Z];
                obj.Direction1 = [v1.X v1.Y v1.Z];
                obj.Direction2 = [v2.X v2.Y v2.Z];
                
            elseif isnumeric(p0) && size(p0, 2) == 3 && isnumeric(n) && size(n, 2) == 3 
                % arguments given as numeric row vectors
                n = normalize(Vector3D(n));
                v0 = Vector3D([1 0 0]);
                if norm(crossProduct(n, v0)) < 1e-12
                    v0 = Vector3D([0 1 0]);
                end
                
                % create direction vectors
                v1 = normalize(crossProduct(n, v0));
                v2 = -normalize(crossProduct(v1, n));
                
                % store result in inner properties
                obj.Origin = p0;
                obj.Direction1 = [v1.X v1.Y v1.Z];
                obj.Direction2 = [v2.X v2.Y v2.Z];
                
            else
                error('Wrong arguments');
            end
            
        elseif nargin == 3
            % Build from origin and two direction vectors
            var1 = varargin{1};
            var2 = varargin{2};
            var3 = varargin{3};
            
            if isa(var1, 'Point3D') && isa(var2, 'Vector3D') && isa(var3, 'Vector3D')
                % Create a plane from origin and two directions
                p0 = var1;
                v1 = var2;
                v2 = var3;
                
                % store result in inner properties
                obj.Origin = [p0.X p0.Y p0.Z];
                obj.Direction1 = [v1.X v1.Y v1.Z];
                obj.Direction2 = [v2.X v2.Y v2.Z];
                
            elseif isnumeric(var1) && size(var1, 2) == 3 && isnumeric(var2) && ...
                    size(var2, 2) == 3 && isnumeric(var3) && size(var3, 2) == 3
                % Create a plane from origin and two directions, given as
                % numeric row vectors
                
                % store result in inner properties
                obj.Origin = var1;
                obj.Direction1 = var2;
                obj.Direction2 = var3;
                
            elseif isa(var1, 'Point3D') && isa(var2, 'Point3D') && isa(var3, 'Point3D')
                % Create the plane passing through three points
                p0 = var1;
                v1 = Vector3D(var1, var2);
                v2 = Vector3D(var1, var3);
                
                % store result in inner properties
                obj.Origin = [p0.X p0.Y p0.Z];
                obj.Direction1 = [v1.X v1.Y v1.Z];
                obj.Direction2 = [v2.X v2.Y v2.Z];
                obj = normalize(obj);

            else
                error('Can not parse input arguments in plane creation');
        
            end
            
        else
            error('Wrong number of arguments in plane creation.');
        end

    end

end % end constructors


%% Methods specific to plane
methods
    function o = origin(obj)
        % Return the origin point of this plane.
        o = Point3D(obj.Origin);
    end
    
    function v1 = direction1(obj)
        % Return the first direction vector of this plane.
        v1 = Vector3D(obj.Direction1);
    end
    
    function v2 = direction2(obj)
        % Return the second direction vector of this plane.
        v2 = Vector3D(obj.Direction2);
    end
    
    function n = normal(obj)
        % Return the normal vector to this plane.
        %
        % The normal of the plane is equal to the cross product of the two
        % direction vectors.
        n = crossProduct(Vector3D(obj.Direction1), Vector3D(obj.Direction2));
    end
    
    function res = normalize(obj)
        % Return a normalized parametric representation of this plane.
        %
        %   PLANE2 = normalize(PLANE1);
        %   Transforms the plane PLANE1 into the pane PLANE2 such that:
        %   - the origin of the plane is the point of the plane the closest
        %       to the global origin
        %   - the first direction vector has norm equal to 1
        %   - the second direction vector has norm equal to 1 and is
        %       orthogonal to first direction vector
        %
        
        % compute first direction vector
        d1  = normalize(Vector3D(obj.Direction1));
        
        % compute second direction vector
        n   = normalize(normal(obj));
        d2  = normalize(crossProduct(n, d1));
        
%         % compute origin point of the plane
%         origin0 = origin(obj);
%         p0 = projPointOnPlane(origins, [plane1(:,1:3) d1 d2]);
        
        % create the resulting plane
        res = Plane3D(origin(obj), d1, d2);
    end
    
    function poly = clip(obj, bounds)
        % Clip this plane with a Bounds3D and return a 3D polygon.
        
        % retrieve min/max coords
        xmin = bounds.XMin;
        xmax = bounds.XMax;
        ymin = bounds.YMin;
        ymax = bounds.YMax;
        zmin = bounds.ZMin;
        zmax = bounds.ZMax;
        
        % create lines corresponding to the edges of the bounding box
        vx = Vector3D(1, 0, 0);
        lineX00 = StraightLine3D(Point3D(xmin, ymin, zmin), vx);
        lineX01 = StraightLine3D(Point3D(xmin, ymin, zmax), vx);
        lineX10 = StraightLine3D(Point3D(xmin, ymax, zmin), vx);
        lineX11 = StraightLine3D(Point3D(xmin, ymax, zmax), vx);
        vy = Vector3D(0, 1, 0);
        lineY00 = StraightLine3D(Point3D(xmin, ymin, zmin), vy);
        lineY01 = StraightLine3D(Point3D(xmin, ymin, zmax), vy);
        lineY10 = StraightLine3D(Point3D(xmax, ymin, zmin), vy);
        lineY11 = StraightLine3D(Point3D(xmax, ymin, zmax), vy);
        vz = Vector3D(0, 0, 1);
        lineZ00 = StraightLine3D(Point3D(xmin, ymin, zmin), vz);
        lineZ01 = StraightLine3D(Point3D(xmin, ymax, zmin), vz);
        lineZ10 = StraightLine3D(Point3D(xmax, ymin, zmin), vz);
        lineZ11 = StraightLine3D(Point3D(xmax, ymax, zmin), vz);

        % compute intersection points with each line
        lines = [...
            lineX00; lineX10; lineX01; lineX11; ...
            lineY00; lineY10; lineY01; lineY11; ...
            lineZ00; lineZ10; lineZ01; lineZ11];
        [points, valids] = intersectLine(obj, lines);
        
        % Keep only points within bounding box, using some tolerance
        acc = 1e-8;
        ivx = [points.X] >= xmin-acc & [points.X] <= xmax+acc;
        ivy = [points.Y] >= ymin-acc & [points.Y] <= ymax+acc;
        ivz = [points.Z] >= zmin-acc & [points.Z] <= zmax+acc;
        valids = valids & ivx & ivy & ivz;
        points = points(valids);
        
        % if there is no intersection point, escape.
        if length(points) < 3
            disp('plane is outside the bounds');
            poly = [];
            return;
        end
        
        % the two spanning lines of the plane
        line1 = StraightLine3D(origin(obj), direction1(obj));
        line2 = StraightLine3D(origin(obj), direction2(obj));
                
        % position of intersection points in plane coordinates
        u1 = position(line1, points);
        u2 = position(line2, points);
        
        % reorder vertices in the correct order
        inds = convhull(u1, u2);
        inds = inds(1:end-1);
        points = points(inds);
        
        % convert to polygon
        poly = LinearRing3D(points);
    end
    
    function [point, valid] = intersectLine(obj, line, varargin)
        % Intersection point of the plane with the given line.
        %
        % Usage:
        %   P = intersectLine(PLANE, LINE);
        %   [P, VALID] = intersectLine(PLANE, LINE);
        %  If line and plane are parallel, return Point with [NaN NaN NaN]
        %  coordinates.
        
        % extract tolerance if needed
        tol = 1e-14;
        if nargin > 2
            tol = varargin{1};
        end
        
        np = numel(obj);
        nl = numel(line);
        point(np, nl) = Point3D();
        valid = false(np, nl);
        
        % iterate over (plane, line) pairs
        for iPlane = 1:np
            plane_i = obj(iPlane);
            
            for iLine = 1:nl
                % current line
                line_i = line(iLine);
                
                % plane normal
                n = normal(obj(iPlane));
                
                % difference between origins of plane and line
                dp = Vector3D(origin(line_i), origin(plane_i));
                
                lineVect = direction(line_i);
                denom = dotProduct(n, lineVect);
                
                % relative position of intersection point on line (can be infinite
                % in case of a line parallel to the plane)
                t = dotProduct(n, dp) / denom;
                
                % compute coord of intersection point
                point(iPlane, iLine) = origin(line_i) + lineVect * t;
                
                % check validity of intersection
                valid(iPlane, iLine) = abs(denom) > tol;
                if ~valid(iPlane, iLine)
                    point(iPlane, iLine) = Point3D([NaN NaN NaN]);
                end
            end
        end
    end
    
    function below = isBelow(obj, point)
        % Test whether a point is below or above a plane.
        %
        % See also
        %   normal
        
        
        % compute position of point projected on 3D line corresponding to
        % plane normal, and returns true for points locatd below the plane
        % (pos<=0).
        normalLine = StraightLine3D(origin(obj), normal(obj));
        below = position(normalLine, point) <= 0;
    end
end % end methods


%% Methods implementing the Geometry3D interface
methods
    function bnd = bounds(obj) %#ok<MANU>
        % Return infinite bounds in each direction.
        bnd = Bounds3D([-inf inf -inf inf -inf inf]);
    end
    
    function h = draw(varargin)
        %DRAW Draw this plane, eventually specifying the style.
        %
        % draw(PLANE)
        % draw(PLANE, COLOR)
        % H = draw(...)
        
        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        
        % default style for drawing planes
        if length(varargin) ~= 1
            varargin = [{'m'}, varargin];
        end

        % extract bounding box of the current axis
        bounds = Bounds3D.fromAxis(ax);
        
        % clip lines with current axis bounds
        poly = clip(obj, bounds);
        
        if ~isempty(poly)
            % draw current clipped line
            hh = fill(ax, poly, varargin{:});
        
            if ~isempty(style)
                apply(style, hh);
            end
        else
            hh = -1;
        end
        
        if nargout > 0
            h = hh;
        end
    end
    
    function res = transform(obj, transform)
        % Apply a geometric transform to this plane.
        %
        % PLANE2 = transform(PLANE, TRANSFO);
        origin = transformPoint(transform, obj.Origin);
        dir1 = transformVector(transform, obj.Direction1);
        dir2 = transformVector(transform, obj.Direction2);
        res = Plane3D(origin, dir1, dir2);
    end
end


%% Methods implementing the Geometry3D interface (more)
methods
    function res = scale(obj, varargin)
        % Return a scaled version of this geometry.
        factor = varargin{1};
        if ~isscalar(factor)
            error('Requires scaling factor to be a scalar');
        end
        origin = obj.Origin * factor;
        dir1 = obj.Direction1 * factor;
        dir2 = obj.Direction2 * factor;
        res = Plane3D(origin, dir1, dir2);
    end
    
    function res = translate(obj, varargin)
        % Return a translated version of this geometry.
        shift = varargin{1};
        newOrigin = obj.Origin + [shift.X shift.Y shift.Z];
        res = Plane3D(newOrigin, obj.Direction1, obj.Direction2);
    end    
end % end methods


%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'Plane3D', ...
            'Origin', obj.Origin, ...
            'Direction1', obj.Direction1, ...
            'Direction2', obj.Direction2);
    end
end
methods (Static)
    function plane = fromStruct(str)
        % Create a new instance from a structure.
        p0 = Point3D(str.Origin);
        v1 = Vector3D(str.Direction1);
        v2 = Vector3D(str.Direction2);
        plane = Plane3D(p0, v1, v2);
    end
end

end % end classdef

