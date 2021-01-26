classdef  (InferiorClasses = {?matlab.graphics.axis.Axes}) Polyline2D < Curve2D
% Abstract class for polyline geometries.
%
%
%   The Polyline2D is an abstract class for all geometries that can be
%   represented by a continuous set of 2D line segments. The two main
%   implementations are LineString2D (for open polylines) and LinearRing2D
%   (for closed polylines).
%
%   The Polyline2D class serves both as tagging class for subclasses, and
%   as a stub for providing implementations for shared functionnalities.
%
%   Example
%     poly1 = Polyline2D.create([10 10;20 10;20 20;10 20;10 30;20 30]);
%     figure; axis equal; axis([0 50 0 50]);hold on;
%     draw(poly1, 'b');
%     poly2 = poly1.translate([-10 -20]).scale(2).translate([25 25]);
%     draw(poly2, 'm');
%
%   See also
%     LineString2D, LinearRing2D 

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2021-01-24,    using Matlab 9.8.0.1323502 (R2020a)
% Copyright 2021 INRAE - BIA-BIBS.


methods (Static)
    function poly = create(coords, varargin)
        % Create a polyline from an array of vertex coordinates.
        
        % parse input argument
        closed = false;
        if length(varargin) > 1
            pname = varargin(1);
            if strcmpi(pname, 'closed')
                closed = varargin{2};
            else
                error('Unnown input argument');
            end
        end
        
        % create polyline
        if closed
            poly = LinearRing2D(coords);
        else
            poly = LineString2D(coords);
        end
    end
end

%% Abstract methods
methods (Abstract)
    coords = vertexCoordinates(obj);
    c = isClosed(obj);
end


%% Constructor
methods
    function obj = Polyline2D(varargin)
        % Constructor for Polyline2D class.

    end

end % end constructors


%% Methods
methods
        
    function box = boundingBox(obj)
        % Return the bounding box of this polyline.
        coords = vertexCoordinates(obj);
        mini = min(coords);
        maxi = max(coords);
        box = Box2D([mini(1) maxi(1) mini(2) maxi(2)]);
    end
    
    function verts = vertices(obj)
        % Get the vertices as a new instance of MultiPoint2D.
        verts = MultiPoint2D(vertexCoordinates(obj));
    end
    
    function nv = vertexCount(obj)
        % Get the number of vertices.
        % (Default implementation that can be overloaded)
        nv = size(vertexCoordinates(obj), 1);
    end
    
    function h = drawVertices(varargin)
        % Draw vertices of this polyline, with optional drawing options.
        
        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        holdState = ishold(ax);
        hold(ax, 'on');
        
        % default options
        if isempty(varargin)
            varargin = {'Marker', 's', 'Color', 'k', 'LineStyle', 'none'};
        end
        
        % extract data
        coords = vertexCoordinates(obj);
        xdata = coords(:,1);
        ydata = coords(:,2);
        
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
end % end methods

end % end classdef

