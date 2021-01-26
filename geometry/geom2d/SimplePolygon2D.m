classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) SimplePolygon2D < Polygon2D
% A simply-connected polygon in the plane.
%
%   Represents a polygon defined by a series of Coords. Its boundary is a
%   linear ring.
%
%   Data are stored as a NV-by-2 array.
%
%   Example
%     data = [10 10; 30 10; 30 20; 20 20; 20 30; 10 30];
%     poly = SimplePolygon2D(data);
%     figure; axis equal; axis([0 50 0 50]); hold on;
%     draw(poly, 'DrawVertices', true);
% 
%   See also
%     Polygon2D, LinearRing2D

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2018-08-14,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2013 INRAE - Cepia Software Platform.


%% Properties
properties
    % The set of vertex coordinates, as a N-by-2 array of double.
    Coords;
    
end % end properties


%% Constructor
methods
    function obj = SimplePolygon2D(varargin)
    % Constructor for SimplePolygon2D class.
    
        if nargin == 0
            % empty constructor
            obj.Coords = [];
            
        elseif nargin == 1
            var1 = varargin{1};
            if isa(var1, 'SimplePolygon2D')
                % copy constructor
                obj.Coords = var1.Coords;
                
            elseif isa(var1, 'LinearRing2D') || isa(var1, 'LineString2D')
                % initialization constructor from another geometry
                obj.Coords = var1.Coords;
                
            elseif isnumeric(var1) 
                % initialization from numerical input
                if size(var1, 2) == 2
                    obj.Coords = var1;
                else
                    error('Creating a polygon requires an array with two columns, not %d', size(var1, 2));
                end
            else
                error('Unable to interpret input argument');
            end

        elseif nargin == 2
            % initialize from two arrays.
            var1 = varargin{1};
            var2 = varargin{2};
            obj.Coords = [var1(:) var2(:)];
        end
    end

end % end constructors

%% Methods specific to SimplePolygon2D
methods
    function centro = centroid(obj)
        % Compute the centroid of this polygon.
        %
        % Return the result as a point.
        
        % isolate coordinates
        px = obj.Coords(:,1);
        py = obj.Coords(:,2);

        % indices of next vertices
        N = length(obj.Coords);
        iNext = [2:N 1];
        
        % compute cross products
        common = px .* py(iNext) - px(iNext) .* py;
        sx = sum((px + px(iNext)) .* common);
        sy = sum((py + py(iNext)) .* common);
        
        % area and centroid
        area = sum(common) / 2;
        centro = Point2D([sx sy] / 6 / area);
    end
    
    function a = area(obj)
        % Compute the area of this polygon.
        
        % isolate coordinates
        px = obj.Coords(:,1);
        py = obj.Coords(:,2);

        % indices of next vertices
        N = length(obj.Coords);
        iNext = [2:N 1];
        
        % compute area (vectorized version)
        a = sum(px .* py(iNext) - px(iNext) .* py) / 2;
    end
    
    function p = perimeter(obj)
        % Compute the perimeter (boundary length) of this polygon.
        dp = diff(obj.Coords([1:end 1], :), 1, 1);
        p = sum(hypot(dp(:, 1), dp(:, 2)));
    end
    

end


%% Implementation of the Polygon2D interface
methods
    function verts = vertices(obj)
        verts = MultiPoint2D(obj.Coords);
    end
    
    function nv = vertexCount(obj)
        nv = size(obj.Coords, 1);
    end
    
    function coords = vertexCoordinates(obj)
        coords = obj.Coords;
    end
end


%% Methods considering the Polygon as a domain
methods
    function bnd = boundary(obj)
        % Return the boundary of this polygon as a LinearRing2D.
        bnd = LinearRing2D(obj.Coords);
    end
end


%% Methods implementing the Geometry2D interface
methods
    function res = transform(obj, transfo)
        % Apply a geometric transform to this polygon.
        res = SimplePolygon2D(transformPoint(transfo, obj.Coords));
    end
    
    function box = boundingBox(obj)
        % Return the bounding box of this polygon.
        mini = min(obj.Coords);
        maxi = max(obj.Coords);
        box = Box2D([mini(1) maxi(1) mini(2) maxi(2)]);
    end
    
    function h = draw(varargin)
        %DRAW Draw the polygon, eventually specifying the style.
        
        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        holdState = ishold(ax);
        hold(ax, 'on');

        % default options
        fillInside = false;
        drawContour = true;
        drawVertices = false;
        if ~isempty(style)
            fillInside = style.FillVisible;
            drawContour  = style.LineVisible;
            drawVertices = style.MarkerVisible;
        end
        
        % parse some options
        inds = strcmpi(varargin, 'fillInside');
        if any(inds)
            inds = find(inds(1));
            fillInside = varargin{inds+1};
            varargin([inds inds+1]) = [];
        end
        inds = strcmpi(varargin, 'drawVertices');
        if any(inds)
            inds = find(inds(1));
            drawVertices = varargin{inds+1};
            varargin([inds inds+1]) = [];
        end
        
        % extract data
        xdata = obj.Coords(:,1);
        ydata = obj.Coords(:,2);

        % draw outline
        h0 = [];
        if fillInside
            options = {'LineStyle', 'none'};
            h0 = fill(ax, xdata, ydata, 'y', options{:});
            if ~isempty(style)
                apply(style, h0);
            end
        end
        
        % draw outline
        h1 = [];
        if drawContour
            if isempty(varargin)
                varargin = {'Color', 'b', 'LineStyle', '-'};
            end
            h1 = plot(ax, xdata([1:end 1]), ydata([1:end 1]), varargin{:});
            if ~isempty(style)
                apply(style, h1);
            end
        end
        
        % optionnally draw markers
        h2 = [];
        if drawVertices
            options = {'Marker', 's', 'Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', 'w'};
            h2 = plot(ax, xdata, ydata, options{:});
            if ~isempty(style)
                apply(style, h2);
            end
        end
        
        if holdState
            hold(ax, 'on');
        else
            hold(ax, 'off');
        end
        
        if nargout > 0
            h = [h0 h1 h2];
        end
    end
end

%% Methods implementing the Geometry2D interface (more)
methods
    function res = scale(obj, varargin)
        % Return a scaled version of this geometry.
        factor = varargin{1};
        res = SimplePolygon2D(obj.Coords * factor);
    end
    
    function res = translate(obj, varargin)
        % Return a translated version of this geometry.
        shift = varargin{1};
        res = SimplePolygon2D(bsxfun(@plus, obj.Coords, shift));
    end
    
    function res = rotate(obj, angle, varargin)
        % Return a rotated version of this polygon.
        %
        % POLY2 = rotate(POLY, THETA)
        % POLY2 = rotate(POLY, THETA, CENTER)
        % THETA is given in degrees, in counter-clockwise order.
        
        origin = [0 0];
        if ~isempty(varargin)
            origin = varargin{1};
        end
        
        rot = createRotation(origin, deg2rad(angle));
        verts = transformPoint(obj.Coords, rot);
        
        res = SimplePolygon2D(verts);
    end
end % end methods


%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'SimplePolygon2D', 'Coordinates', obj.Coords);
    end
end
methods (Static)
    function poly = fromStruct(str)
        % Create a new instance from a structure.
        if isfield(str, 'Coordinates')
            poly = SimplePolygon2D(str.Coordinates);
        elseif isfield(str, 'coordinates')
            poly = SimplePolygon2D(str.coordinates);
        else
            error('Field <Coordinates> of SimplePolygon2D is not defined');
        end
    end
end

end % end classdef

