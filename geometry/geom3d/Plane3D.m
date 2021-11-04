classdef Plane3D < handle
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
end


%% Properties
properties
    Origin;
    Vector1;
    Vector2;
end % end properties


%% Constructor
methods
    function obj = Plane3D(varargin)
        % Constructor for Plane3D class.
        
        if nargin == 0
            % Empty constructor
            
        elseif nargin == 2
            if isa(varargin{1}, 'Point3D') && isa(varargin{2}, 'Vector3D')
                % Build from position and normal vector
                p0 = varargin{1};
                n = normalize(varargin{2});
                v0 = Vector3D([1 0 0]);
                if norm(crossProduct(n, v0)) < 1e-12
                    v0 = Vector3D([0 1 0]);
                end
                
                % create direction vectors
                v1 = normalize(crossProduct(n, v0));
                v2 = -normalize(crossProduct(v1, n));
                
                % store result in inner properties
                obj.Origin = [p0.X p0.Y p0.Z];
                obj.Vector1 = [v1.X v1.Y v1.Z];
                obj.Vector2 = [v2.X v2.Y v2.Z];

            else
                error('Wrong arguments');
            end
            
        elseif nargin == 3
            % Build from origin and two direction vectors
            var1 = varargin{1};
            var2 = varargin{2};
            var3 = varargin{3};
            
            if isa(var1, 'Point3D') && isa(var2, 'Point3D') && isa(var3, 'Point3D')
                p0 = var1;
                v1 = Vector3D(var1, var2);
                v2 = Vector3D(var1, var3);
                
                % store result in inner properties
                obj.Origin = [p0.X p0.Y p0.Z];
                obj.Vector1 = [v1.X v1.Y v1.Z];
                obj.Vector2 = [v2.X v2.Y v2.Z];
                obj = normalize(obj);

            elseif isa(var1, 'Point3D') && isa(var2, 'Vector3D') && isa(var3, 'Vector3D')
                p0 = var1;
                v1 = var2;
                v2 = var3;
                
                % store result in inner properties
                obj.Origin = [p0.X p0.Y p0.Z];
                obj.Vector1 = [v1.X v1.Y v1.Z];
                obj.Vector2 = [v2.X v2.Y v2.Z];

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
        v1 = Vector3D(obj.Vector1);
    end
    
    function v2 = direction2(obj)
        % Return the second direction vector of this plane.
        v2 = Vector3D(obj.Vector2);
    end
    
    function n = normal(obj)
        % Return the normal vector to this plane.
        %
        % The normal of the plane is equal to the cross product of the two
        % direction vectors.
        n = crossProduct(Vector3D(obj.Vector1), Vector3D(obj.Vector2));
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
        d1  = normalize(Vector3D(obj.Vector1));
        
        % compute second direction vector
        n   = normalize(normal(obj));
        d2  = normalize(crossProduct(n, d1));
        
%         % compute origin point of the plane
%         origin0 = origin(obj);
%         p0 = projPointOnPlane(origins, [plane1(:,1:3) d1 d2]);
        
        % create the resulting plane
        res = Plane3D(origin(obj), d1, d2);
    end
    
    function [point, valid] = intersectLine(obj, line, varargin)
        % Intersection point of the plane with the given line.
        %  If line and plane are parallel, return Point with [NaN NaN NaN]
        %  coordinates.
        
        % extract tolerance if needed
        tol = 1e-14;
        if nargin > 2
            tol = varargin{1};
        end
        
        % plane normal
        n = normal(obj);
        
        % difference between origins of plane and line
        dp = Vector3D(origin(line), origin(obj));
        
        lineVect = direction(line);
        denom = dotProduct(n, lineVect);
        
        % relative position of intersection point on line (can be infinite
        % in case of a line parallel to the plane)
        t = dotProduct(n, dp) / denom;
        
        % compute coord of intersection point
        point = origin(line) + lineVect * t;
        
        valid = abs(denom) > tol;
        if ~valid
            point.X = NaN;
            point.Y = NaN;
            point.Z = NaN;
        end
    end
    
end % end methods

end % end classdef

