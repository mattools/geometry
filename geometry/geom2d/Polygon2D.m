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
    
    function [hull, inds] = convexHull(varargin)
        % Compute the convex hull of the set of input points.
        %
        % This function is mainly a wrapper to the convhull function, that
        % format the result to a Polygon2D.
        %
        % Example
        %   % Generate a point set, and compute its convex hull
        %   coords = [10 30;10 60;30 20; 30 40; 30 70; 40 80; 50 40; 50 60; 70 30; 70 50; 80 70];
        %   pts = MultiPoint2D(coords);
        %   figure; draw(pts, 'bo'); axis([0 100 0 100]); hold on;
        %   hull = Polygon2D.convexHull(pts);
        %   draw(hull, 'b');
        
        % get coordinates of points
        if nargin == 1
            var1 = varargin{1};
            if isa(var1, 'MultiPoint2D')
                points = var1.Coords;
            elseif isnumeric(var1) && size(var1, 2) == 2
                points = var1;
            else
                error('Input argument must be either numeric or a MultiPoint2D');
            end
            
        elseif nargin > 1
            np = length(varargin);
            points = zeros(np, 2);
            
            for i = 1:np
                arg = varargin{i};
                if isnumeric(arg)
                    points(i, :) = arg;
                elseif isa(arg, 'Point2D')
                    points(i, :) = [arg.X arg.Y];
                else
                    error('Can not interpret input argument');
                end
            end
        else
            error('Requires at least one input argument');
        end
        
        % checkup on array size
        if size(points, 1) < 3
            hull = points;
            inds = 1:size(points, 1);
            return;
        end
        
        % parse simplify option
        simplify = true;
%         if nargin > 2 && strcmpi(varargin{1}, 'simplify')
%             simplify = varargin{2};
%         end
        
        % compute convex hull by calling the 'convhull' function
        inds = convhull(points(:,1), points(:,2), 'simplify', simplify);
        
        % avoid considering the same vertex at the beginning and at the end
        if inds(1) == inds(end)
            inds(end) = [];
        end
        
        % create polygon
        hull = SimplePolygon2D(points(inds, :));
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

