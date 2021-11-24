classdef Sphere3D < Geometry3D
% A sphere, defined by a center and a radius.
%
%   Class Sphere3D
%
%   Example
%   Sphere3D
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2020-03-06,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2020 INRAE - BIA-BIBS.

%% Static methods
methods (Static)
    function alpha = sphericalAngle(p1, p2, p3)
        % Compute angle between points on the sphere.
        %
        %   ALPHA = sphericalAngle(P1, P2, P3)
        %   Computes angle (P1, P2, P3), i.e. the angle, measured at point P2,
        %   between the direction (P2, P1) and the direction (P2, P3).
        %   The result is given in radians, between 0 and 2*PI.
        %
        
        % convert points to (normalized) vectors
        p1 = Point3D(normalize(Vector3D(p1)));
        p2 = Point3D(normalize(Vector3D(p2)));
        p3 = Point3D(normalize(Vector3D(p3)));
        
        % create the plane tangent to the unit sphere and containing central point
        plane = Plane3D(p2, Vector3D(p2));
        
        % project the two other points on the plane
        pp1 = position(plane, projection(plane, p1));
        pp3 = position(plane, projection(plane, p3));
        
        % compute angle on the tangent plane
        pp2 = zeros(max(size(pp1, 1), size(pp3,1)), 2);
        alpha = angle3Points(pp1, pp2, pp3);
    end
end


%% Properties
properties
    % The coordinates of the center.
    Center = [0 0 0];
    
    % The radius of the sphere. Default is 1.
    Radius  = 1;
    
end % end properties


%% Constructor
methods
    function obj = Sphere3D(varargin)
        % Constructor for Sphere3D class
        if nargin == 1
            var1 = varargin{1};
            if isa(var1, 'Sphere3D')
                % copy constructor
                obj.Center = var1.Center;
                obj.Radius  = var1.Radius;
                
            elseif isnumeric(var1)
                % initialize with a 1-by-4 row vector
                obj.Center = var1(1, 1:3);
                obj.Radius  = var1(4);
            end
        elseif nargin == 2
            % center + radius
            
            % initialize center
            var1 = varargin{1};
            if isa(var1, 'Point3D')
                obj.Center = [var1.X var1.Y var1.Z];
            elseif isnumeric(var1)
                obj.Center = var1(1, 1:3);
            else
                error('Can not interpret first argument');
            end
            
            % initialize radius
            obj.Radius = varargin{2};
            
        elseif nargin > 2
            error('Can not initialize sphere');
        end
    end

end % end constructors


%% Methods specific to Sphere3D
methods
    function s = surfaceArea(obj)
        s = reshape([obj.Radius], size(obj)) .^4 * 4 * pi;
    end
    
    function [point, valid] = intersectLine(obj, lines, varargin)
        % Intersection points between a line and a sphere.
        
        % check if user-defined tolerance is given
        tol = 1e-14;
        if ~isempty(varargin)
            tol = varargin{1};
        end
        
        dim1 = size(obj);
        dim2 = size(lines);
        dim = max(dim1, dim2);
        nd = length(dim);
        
        % allocate memory for output
        dims2 = [dim 2];
        point(prod(dims2)) = Point3D;
        point = reshape(point, dims2);
        valid = true(dim);
        
        % iterate over the pairs of inputs
        for i = 1:prod(dim)
            % convert linear index of output to sub indices
            inds = cell(nd,1);
            [inds{:}] = ind2sub(dim, i);
%             inds = [inds{:}];
            
            % convert sub indices to linear index of inputs
            tmp = num2cell(min([inds{:}], dim1)); ind1 = sub2ind(dim1, tmp{:});
            tmp = num2cell(min([inds{:}], dim2)); ind2 = sub2ind(dim2, tmp{:});
            
            sphere = obj(ind1);
            line = lines(ind2);
            
            % difference between centers
            dc = Vector3D(line.origin(), Point3D(sphere.Center));
            
            % equation coefficients
            dl = line.direction();
            a = dotProduct(dl, dl);
            b = 2 * dotProduct(dc, dl);
            c = dotProduct(dc, dc) - sphere.Radius^2;
            
            % solve equation
            delta = b.*b - 4*a.*c;
            
            if delta > tol
                % process couples with two intersection points
                
                % delta positive: find two roots of second order equation
                u1 = (-b -sqrt(delta)) / 2 / a;
                u2 = (-b +sqrt(delta)) / 2 / a;
                
                % convert into 3D coordinate
                origin = coordinates(line.origin());
                point1 = Point3D(origin + coordinates(dl * u1));
                point2 = Point3D(origin + coordinates(dl * u2));
                
                % convert to linear indices in result array
                inds2 = [inds ; {1}]; indR1 = sub2ind(dims2, inds2{:});
                point(indR1) = point1;
                inds2 = [inds ; {2}]; indR2 = sub2ind(dims2, inds2{:});
                point(indR2) = point2;
                
            elseif abs(delta) <= tol
                % process couples with one intersection point
                
                % delta around zero: find unique root, and convert to 3D coord.
                u = -b / 2 / a;
                point = Point3D(coordinates(origin(line)) + coordinates(dl * u));
                
                inds2 = [inds ; {1}]; indR1 = sub2ind(dims2, inds2{:});
                inds2 = [inds ; {2}]; indR2 = sub2ind(dims2, inds2{:});
                point([indR1 indR2]) = [point point];

            else
                valid(i) = false;
            end
        end
        
        if nd == 2 && dim(2) == 1
            point = squeeze(point)';
        end
    end
end


%% Methods implementing the Geometry3D interface
methods
    function res = transform(obj, transform) %#ok<STOUT,INUSD>
        % Apply a geometric transform to this geometry.
        error('Method not implemented');
    end
    
    function box = bounds(obj)
        % Return the bounding box of this shape.
        
        box = Bounds3D([...
            obj.Center(1) - obj.Radius obj.Center(1) + obj.Radius ...
            obj.Center(2) - obj.Radius obj.Center(2) + obj.Radius ...
            obj.Center(3) - obj.Radius obj.Center(3) + obj.Radius]);
    end
    
    function h = draw(varargin)
        %DRAW Draw this sphere, eventually specifying the style.
        
        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        
        % default values for drawing
        nPhi    = 32;
        nTheta  = 16;
        
        % process input options: when a string is found, assumes this is the
        % beginning of options
        options = {'FaceColor', 'g', 'LineStyle', 'none'};
        if length(varargin) == 1
            options = {'FaceColor', varargin{1}, 'LineStyle', 'none'};
        else
            options = [options varargin];
        end
        
        % compute spherical coordinates
        theta   = linspace(0, pi, nTheta+1);
        phi     = linspace(-pi, +pi, nPhi+1);
        
        % convert to cartesian coordinates
        sintheta = sin(theta);
        x = obj.Center(1) + cos(phi') * sintheta * obj.Radius;
        y = obj.Center(2) + sin(phi') * sintheta * obj.Radius;
        z = obj.Center(3) + ones(length(phi),1) * cos(theta) * obj.Radius;
        
        % draw the geometric primitive
        hh = surf(ax, x, y, z, options{:});
        
        % optionnally add style processing
        if ~isempty(style)
            apply(style, hh);
        end
        
        % format output argument
        if nargout > 0
            h = hh;
        end
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
        center = obj.Center * factor;
        radius = obj.Radius * factor;
        res = Sphere3D(center, radius);
    end
    
    function res = translate(obj, varargin)
        % Return a translated version of this geometry.
        shift = varargin{1};
        center = obj.Center + shift;
        res = Sphere3D(center, obj.Radius);
    end    
end % end methods


%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'Sphere3D', ...
            'Center', obj.Center, ...
            'Radius', obj.Radius);
    end
end
methods (Static)
    function sphere = fromStruct(str)
        % Create a new instance from a structure.
        sphere = Sphere3D(str.Center, str.Radius);
    end
end

end % end classdef

