classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) OrientedBox2D < Geometry2D
% One-line description here, please.
%
%   Class OrientedBox2D
%
%   Example
%   OrientedBox2D
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2021-01-22,    using Matlab 9.8.0.1323502 (R2020a)
% Copyright 2021 INRAE - BIA-BIBS.


%% Properties
properties
    CenterX = 0;
    CenterY = 0;
    Size1 = 1;
    Size2 = 1;
    % The orientation in degrees, CCW
    Orientation = 0;
end % end properties


%% Constructor
methods
    function obj = OrientedBox2D(varargin)
        % Constructor for OrientedBox2D class.

        switch nargin
            case 0
                % nothing to do
            case 1
                var1 = varargin{1};
                if size(var1, 2) ~= 5
                    error('Creating an oriented box requires an array with five columns, not %d', size(var1, 2));
                end
                obj.CenterX = var1(1);
                obj.CenterY = var1(2);
                obj.Size1 = var1(3);
                obj.Size2 = var1(4);
                obj.Orientation = var1(5);
        end
    end

end % end constructors

%% Methods specific to Ellipse2D
methods
    function center = center(obj)
        % Returns the center of this circle as a Point2D.
        center = Point2D(obj.CenterX, obj.CenterY);
    end
    
    function poly = asPolyline(obj, varargin)
        % Converts this box into a (closed) polyline.
        %
        % POLY = asPolyline(OBJ);
        % Returns the result as an instance of LinearRing2D.
        
        % easier to compute with w and h divided by 2
        w1 = obj.Size1 / 2;
        w2 = obj.Size2 / 2;
        
        % generate array of vertex coordinates, centered on (0,0)
        M = bsxfun (@times, [-1 1; 1 1; 1 -1; -1 -1], [w1 w2]);
        
        % apply rotation and translation
        v = [cosd(obj.Orientation); sind(obj.Orientation)];
        tx  = obj.CenterX + M * v;
        ty  = obj.CenterY + M(4:-1:1,[2 1]) * v;
        
        % convert to polyon
        poly = LinearRing2D.create([tx ty]);
    end
end

%% Methods implementing the Geometry2D interface
methods
    function res = transform(obj, transform) 
        % Apply a geometric transform to this geometry.
        res = transform(asPolyline(obj), transform);
    end
    
    function box = bounds(obj)
        % Return the bounding box of this geometry.
        extX = [obj.CenterX - obj.Radius obj.CenterX + obj.Radius];
        extY = [obj.CenterY - obj.Radius obj.CenterY + obj.Radius];
        box = Bounds2D([extX extY]);
    end
    
    function h = draw(varargin)
        %DRAW Draw the current geometry, eventually specifying the style.
        
        % parse arguments using protected method implemented in Geometry
        [ax, obj, style] = parseDrawInputArguments(varargin{:});
        
        hh = draw(ax, asPolyline(obj));
        if ~isempty(style)
            apply(style, hh);
        end
        
        if nargout > 0
            h = hh;
        end
    end
    
    function res = scale(obj, factor)
        % Return a scaled version of this geometry.
        res = OrientedBox2D([[obj.CenterX obj.CenterY obj.Size1 obj.Size2] * factor obj.Orientation]);
    end
    
    function res = translate(obj, shift)
        % Return a translated version of this geometry.
        res = OrientedBox2D([obj.CenterX+shift(1) obj.CenterY+shift(2) obj.Size2 obj.Radius2 obj.Orientation]);
    end
    
    function res = rotate(obj, angle, varargin)
        % Return a rotated version of this ellipse.
        %
        % POLY2 = rotate(POLY, THETA)
        % POLY2 = rotate(POLY, THETA, CENTER)
        % THETA is given in degrees, in counter-clockwise order.
        
        
        transfo = AffineTransform2D.createRotation(angle, varargin{:});
        center2 = transfo.transformPoint([obj.CenterX obj.CenterY]);
        orientation2 = obj.Orientation + rad2deg(angle);
        res = OrientedBox2D([center2  obj.Radius1 obj.Radius2 orientation2]);
    end
end % end methods


%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'OrientedBox2D', ...
            'CenterX', obj.CenterX, ...
            'CenterY', obj.CenterY, ...
            'Size1', obj.Size1, ...
            'Size2', obj.Size2, ...
            'Orientation', obj.Orientation);
    end
end
methods (Static)
    function circ = fromStruct(str)
        % Create a new instance from a structure.
        circ = OrientedBox2D([str.CenterX str.CenterY str.Size1 str.Size2 str.Orientation]);
    end
end

end % end classdef

