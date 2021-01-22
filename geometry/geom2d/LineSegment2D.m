classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) LineSegment2D < Curve2D
% A line segment defined by its two extremities.
%
%   Usage
%   L = LineSegment2D(P1, P2);
%   Creates the line segment by specifying its two extremities. P1 and P2
%   must be either instances of Point2D, or 1-by-2 numeric arrays
%   containing coordinates of extremities. 
%
%   Example
%     P1 = Point2D(20, 10);
%     P2 = Point2D(40, 20);
%     L = LineSegment2D(P1, P2);
%     draw(L);
%
%   See also
%     Point2D
%

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2019-10-15,    using Matlab 9.7.0.1190202 (R2019b)
% Copyright 2019 INRA - BIA-BIBS.


%% Properties
properties
    % The coordinate of the source point, as a 1-by-2 numeric array.
    P1 = [0 0];

    % The coordinate of the target point, as a 1-by-2 numeric array.
    P2 = [1 0];
    
end % end properties


%% Constructor
methods
    function obj = LineSegment2D(varargin)
        % Constructor for the LineSegment2D class.
        
        if nargin == 0
            % Default constructor: unit line segment
            
        elseif nargin == 1
            % Copy constructor
            if ~isa(varargin{1}, 'LineSegment2D')
                error('Requires a LineSegment2D as input');
            end
            var1 = varargin{1};
            obj.P1 = var1.P1;
            obj.P2 = var1.P2;
            
        elseif nargin == 2
            p1 = varargin{1};
            if isa(p1, 'Point2D')
                obj.P1 = [p1.X p1.Y];
            elseif isnumeric(p1) && all(size(p1) == [1 2])
                obj.P1 = p1;
            else
                error('Unable to interpret first input arument');
            end
              
            p2 = varargin{2};
            if isa(p2, 'Point2D')
                obj.P2 = [p2.X p2.Y];
            elseif isnumeric(p2) && all(size(p2) == [1 2])
                obj.P2 = p2;
            else
                error('Unable to interpret second input arument');
            end
        end

    end

end % end constructors


%% Methods specific to LineSegment2D
methods
    function [dist, pos] = distancePoint(obj, point)
        % Minimum distance between point(s) and this line segment.
        
        % direction vector of the line segment
        vx = obj.P2(1) - obj.P1(1);
        vy = obj.P2(2) - obj.P1(2);
        
        % squared length of edges, with a check of validity
        delta = vx .* vx + vy .* vy;
        
        % difference of coordinates between points and first point
        dx  = point(:, 1) - obj.P1(1);
        dy  = point(:, 2) - obj.P1(2);
        
        % compute position of points projected on the supporting line,
        % by using normalised dot product (NP-by-NE array)
        if delta > eps
            % ensure projected point is located on the edge
            pos = min(max((dx * vx + dy * vy) / delta, 0), 1);
        else
            % consider point1 is the closest egde point
            pos = 0;
        end
        
        % compute distance between point and its projection on the edge
        dist = hypot(pos * vx - dx, pos * vy - dy);
    end
    
    function pm = middlePoint(obj)
        % Return the middle point of this line segment.
        xm = mean(obj.P1);
        ym = mean(obj.P2);
        pm = Point2D(xm, ym);
    end
end

%% Methods generic to curve objects
methods
    function l = length(obj)
        % Return the length of this line segment.
        dx = obj.P2(1) - obj.P1(1);
        dy = obj.P2(2) - obj.P1(2);
        l = hypot(dx, dy);
    end
    
    function res = reverse(obj)
        % Reverse this line segment.
        res = LineSegment2D(obj.P2, obj.P1);
    end
    
    function p1 = firstPoint(obj)
        % Return the first point of this line segment.
        p1 = Point2D(obj.P1);
    end
    
    function p2 = lastPoint(obj)
        % Return the last point of this line segment.
        p2 = Point2D(obj.P2);
    end
end


%% Methods implementing the Geometry2D interface
methods
    function [dist, pos] = distance(obj, point)
        % Distance between point(s) and this line segment.
        %
        % Example:
        %   seg = LineSegment2D(Point2D(10, 10), Point2D(70, 55));
        %   p = Point2D(56, 32);
        %   distance(seg, p)
        %   ans =
        %       10
        %
        % see also
        %   distancePoint
        if isa(point, 'Point2D')
            [dist, pos] = distancePoint(obj, [point.X point.Y]);
        end
    end
    
    function res = transform(obj, transfo)
        % Apply a geometric transform to this line segment.
        p1t = transformPoint(transfo, [obj.X1 obj.Y1]);
        p2t = transformPoint(transfo, [obj.X2 obj.Y2]);
        res = LineSegment2D(p1t, p2t);
    end
    
    function box = boundingBox(obj)
        % Returns the bounding box of this geometry.
        x = sort([obj.P1(1) obj.P2(1)]);
        y = sort([obj.P1(2) obj.P2(2)]);
        box = Box2D([x y]);
    end
    
    function h = draw(varargin)
        % Draws the current geometry, eventually specifying the style.

        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        
        % plot line segment
        xdata = [obj.P1(1) obj.P2(1)];
        ydata = [obj.P1(2) obj.P2(2)];
        hh = plot(ax, xdata, ydata, varargin{:});
        
        if ~isempty(style)
            apply(style, hh);
        end
        
        if nargout > 0
            h = hh;
        end
    end
    
    function res = scale(obj, factor)
        % Returns a scaled version of this line segment.
        res = LineSegment2D(obj.P1 * factor, obj.P2 * factor);
    end
    
    function res = translate(obj, shift)
        % Returns a translated version of this line segment.
        res = LineSegment3D(obj.P1 + shift, obj.P2 + shift);
    end
    
    function res = rotate(obj, varargin)
        % Returns a rotated version of this line segment.
        origin = [0 0];
        if ~isempty(varargin)
            origin = varargin{1};
        end
        
        rot = createRotation(origin, deg2rad(angle));
        p1t = transformPoint(rot, obj.P1);
        p2t = transformPoint(rot, obj.P2);
        
        res =  LineSegment2D(p1t, p2t);
    end
    
end % end methods


%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'LineSegment2D', ...
            'P1', obj.P1, ...
            'P2', obj.P2);
    end
end
methods (Static)
    function line = fromStruct(str)
        % Create a new instance from a structure.
        line = LineSegment2D(str.P1, str.P2);
    end
end

end % end classdef

