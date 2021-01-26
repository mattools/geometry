classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) Polygon2D < Geometry2D
% Parent class for representing polygons in the plane.
%
%   The Polygon2D class is an abstract class for all geometries whose
%   boundary can be represented by a set of 2D line segments forming one or
%   more closed polylines.
%
%   The class SimplePolygon2D is the main implementation of Polygon2D.
%
%
%   Example
%     % Create and draw a simple polygon
%     poly = Polygon2D.create([0 10; 10 0; 20 10;10 20]);
%     figure; draw(poly); axis equal;
%
%   See also
%     Geometry2D, Polyline2D, LinearRing2D

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2018-08-14,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2013 INRA - Cepia Software Platform.


%% Static factory methods
methods (Static)
    function poly = create(coords)
        % Create a polygon from an array of vertex coordinates.
        % 
        % Example:
        %   poly = Polygon2D.create([0 10; 10 0; 20 10;10 20]);
        %
        poly = SimplePolygon2D(coords);
    end
end


%% Abstract methods
methods (Abstract)
    % Return vertices as an instance of MultiPoint2D.
    verts = vertices(obj);
    % The number of vertices in this polygon.
    nv = vertexCount(obj);
    % Return the coordinates of vertices as a numeric array.
    coords = vertexCoordinates(obj);
    
    % Return the boundary of this polygon.
    bnd = boundary(obj);
end


%% Constructor
methods
    function obj = Polygon2D(varargin)
        % Constructor for Polygon2D class.
    end
end % end constructors

end % end classdef

