classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) LineString2D < Curve2D
% An open polyline composed of several line segments. 
%
%   Represents a polyline defined be a series of vertex coordinates. 
%
%   Data are represented by a NV-by-2 array.
%
%   Example
%     poly1 = LineString2D([10 10;20 10;20 20;10 20;10 30;20 30]);
%     figure; axis equal; axis([0 50 0 50]);hold on;
%     draw(poly1, 'b');
%     poly2 = poly1.translate([-10 -20]).scale(2).translate([25 25]);
%     draw(poly2, 'm');
%
%   See also
%   Geometry2d, Polygon2D

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-08-14,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2013 INRA - Cepia Software Platform.


%% Properties
properties
    % the set of Coords, given as a N-by-2 array of double.
    Coords;
    
end % end properties


%% Constructor
methods
    function obj = LineString2D(varargin)
    % Constructor for LineString2D class.
    
        if ~isempty(varargin)
            var1 = varargin{1};
            if size(var1, 2) ~= 2
                error('Creating a linestring requires an array with two columns, not %d', size(var1, 2));
            end
            obj.Coords = var1;

        else
            obj.Coords = [];
        end
    end

end % end constructors


%% Methods specific to LineString2D
methods
    function res = smooth(obj, M)
        % Smooth a polyline using local averaging.

        % create convolution vector
        v2 = ones(M, 1) / M;
        
        % allocate memory for result
        res = zeros(size(obj.Coords));
        
        % iterate over dimensions
        for d = 1:2
            v0 = obj.Coords(1, d);
            v1 = obj.Coords(end, d);
            vals = [v0(ones(M, 1)) ; obj.Coords(:,d) ; v1(ones(M, 1))];
            resd = conv(vals, v2, 'same');
            res(:,d) = resd(M+1:end-M);
        end
        
        % convert result to LineString object
        res = LineString2D(res);
    end

    function poly2 = resample(obj, n)
        % RESAMPLE Resample this polyline with a given number of vertices.
        %
        %   Syntax:  POLY2 = resample(POLY, N);

        % compute arc length along each ertex
        s = verticesArcLength(obj);
        
        % distribute N points equally spaced
        Lmax = s(end);
        pos = linspace(0, Lmax, n);

        poly2 = zeros(n, size(obj.Coords, 2));
        for i = 1:n
            % index of surrounding vertices before and after
            ind0 = find(s <= pos(i), 1, 'last');
            ind1 = find(s >= pos(i), 1, 'first');
            
            if ind0 == ind1
                % get position of a vertex in input polyline
                poly2(i, :) = obj.Coords(ind0, :);
                continue;
            end
            
            % position of surrounding vertices
            pt0 = obj.Coords(ind0, :);
            pt1 = obj.Coords(ind1, :);
            
            % weights associated to each neighbor
            l0 = pos(i) - s(ind0);
            l1 = s(ind1) - pos(i);
            
            % linear interpolation of neighbor positions
            if (l0 + l1) > Lmax * 1e-12
                poly2(i, :) = (pt0 * l1 + pt1 * l0) / (l0 + l1);
            else
                % if neighbors are too close, do not use interpolation
                poly2(i, :) = pt0;
            end
        end
        
        % Convert result to a LineString2D instance
        poly2 = LineString2D(poly2);
    end
    
    function al = verticesArcLength(obj)
        % Return the arc length at each vertex of the polyline.
        
        % compute the cumulative  sum of the length of each line segment,
        % and add 0 for the first vertex.
        al = [0 ; cumsum(sqrt(sum(diff(obj.Coords).^2, 2)))];
    end

    function l = length(obj)
        % Return the length of this polyline.

        % compute the sum of the length of each line segment
        l = sum(sqrt(sum(diff(obj.Coords).^2, 2)));
    end
    
    function centro = centroid(obj)
        % Return the centroid of this polyline.
        
        % compute center and length of each line segment
        centers = (obj.Coords(1:end-1,:) + obj.Coords(2:end,:))/2;
        lengths = sqrt(sum(diff(obj.Coords).^2, 2));
        
        % centroid of edge centers weighted by edge lengths
        centro = Point2D(sum(bsxfun(@times, centers, lengths), 1) / sum(lengths));
    end
    
    function nv = vertexNumber(obj)
        % Get the number of vertices.
        nv = size(obj.Coords, 1);
    end
end

%% Methods
methods
    function res = transform(obj, transform)
        % Apply a geometric transform to this geometry.
        res = LineString2D(transformCoords(transform, obj.Coords));
    end
    
    function box = boundingBox(obj)
        % Return the bounding box of this polyline.
        mini = min(obj.Coords);
        maxi = max(obj.Coords);
        box = Box2D([mini(1) maxi(1) mini(2) maxi(2)]);
    end
    
    function verts = vertices(obj)
        % Get the vertices as a new instance of MultiPoint2D.
        verts = MultiPoint2D(obj.Coords);
    end
    
    function h = drawVertices(varargin)
        % Draw vertices of this polyline, with optional drawing options.
        
        % extract drawing options
        [ax, obj, style, varargin] = parseDrawOptions(varargin{:});
        holdState = ishold(ax);
        hold(ax, 'on');
        
        % default options
        if isempty(varargin)
            varargin = {'Marker', 's', 'Color', 'k', 'LineStyle', 'none'};
        end
        
        % extract data
        xdata = obj.Coords(:,1);
        ydata = obj.Coords(:,2);
        
        hh = plot(ax, xdata, ydata, varargin{:});
        if ~isempty(style)
            apply(style, hh);
        end
        
        if holdState
            hold(ax, 'on');
        else
            hold(ax, 'off');
        end
        
        if nargout > 0
            h = hh;
        end
    end
    
    function h = draw(varargin)
        %DRAW Draw the current geometry, eventually specifying the style.
        
        % extract drawing options
        [ax, obj, style, varargin] = parseDrawOptions(varargin{:});
        holdState = ishold(ax);
        hold(ax, 'on');

        % default options
        drawLines = true;
        drawVertices = false;
        if ~isempty(style)
            drawLines = style.LineVisible;
            drawVertices = style.MarkerVisible;
        end
        
        % parse some options
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
        h1 = [];
        if drawLines
            if isempty(varargin)
                varargin = {'Color', 'b', 'LineStyle', '-'};
            end
            h1 = plot(ax, xdata, ydata, varargin{:});
            if ~isempty(style)
                apply(style, h1);
            end
        end
        
        % optionnally draw markers
        h2 = [];
        if drawVertices
            options = {'Marker', 's', 'Color', 'k', 'LineStyle', 'none'};
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
            h = [h1 h2];
        end
    end
    
    function res = scale(obj, varargin)
        % Return a scaled version of this geometry.
        factor = varargin{1};
        res = LineString2D(obj.Coords * factor);
    end
    
    function res = translate(obj, varargin)
        % Return a translated version of this geometry.
        shift = varargin{1};
        res = LineString2D(bsxfun(@plus, obj.Coords, shift));
    end
    
    function res = rotate(obj, angle, varargin)
        % Return a rotated version of this polyline.
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
        
        res = LineString2D(verts);
    end
end % end methods

%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'LineString2D', 'Coordinates', obj.Coords);
    end
end
methods (Static)
    function poly = fromStruct(str)
        % Create a new instance from a structure.
        if isfield(str, 'Coordinates')
            poly = LineString2D(str.Coordinates);
        elseif isfield(str, 'coordinates')
            poly = LineString2D(str.coordinates);
        else
            error('Field <Coordinates> of LineString2D is not defined');
        end
    end
end

end % end classdef

