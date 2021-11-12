classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) LineSegment3D < Geometry3D
% A 3D line segment defined by its two extremities.
%
%   Class LineSegment3D
%
%   Example
%     P1 = Point3D(30, 20, 10);
%     P2 = Point3D(50, 40, 20);
%     L = LineSegment3D(P1, P2);
%     draw(L);
%
%   See also
%     Point3D, LineSegment2D

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2020-02-06,    using Matlab 9.7.0.1190202 (R2019b)
% Copyright 2019 INRA - BIA-BIBS.


%% Properties
properties
    % The x-coordinate of the source point.
    X1 = 0;
    % The y-coordinate of the source point.
    Y1 = 0;
    % The z-coordinate of the source point.
    Z1 = 0;

    % The x-coordinate of the target point.
    X2 = 1;
    % The y-coordinate of the target point.
    Y2 = 0;
    % The z-coordinate of the target point.
    Z2 = 0;
    
end % end properties


%% Constructor
methods
    function obj = LineSegment3D(varargin)
        % Constructor for LineSegment3D class.
        
        if nargin == 0
            % Default constructor: unit line segment
            return;
            
        elseif nargin == 1
            % Copy constructor
            if ~isa(varargin{1}, 'LineSegment3D')
                error('Requires a LineSegment3D as input');
            end
            var1 = varargin{1};
            obj.X1 = var1.X1;
            obj.Y1 = var1.Y1;
            obj.Z1 = var1.Z1;
            obj.X2 = var1.X2;
            obj.Y2 = var1.Y2;
            obj.Z2 = var1.Z2;
            
        elseif nargin == 2
            % Creation from two points, either as Point3D or as numeric
            % arrays
            
            % check input validity
            p1 = varargin{1};
            p2 = varargin{2};
            if any(size(p1) ~= size(p2))
                error('Both input arguments must have same size');
            end
                
            if isa(p1, 'Point3D') && isa(p2, 'Point3D')
                n1 = size(p1, 1);
                n2 = size(p1, 2);
                obj(n1, n2) = LineSegment3D();
                for i = 1:n1
                    for j = 1:n2
                        obj(i,j).X1 = p1(i,j).X;
                        obj(i,j).Y1 = p1(i,j).Y;
                        obj(i,j).Z1 = p1(i,j).Z;
                        obj(i,j).X2 = p2(i,j).X;
                        obj(i,j).Y2 = p2(i,j).Y;
                        obj(i,j).Z2 = p2(i,j).Z;
                    end
                end
                
            elseif isnumeric(p1) && size(p1, 2) == 3 && isnumeric(p2) && size(p2, 2) == 3
                n1 = size(p1, 1);
                obj(n1, 1) = LineSegment3D();
                for i = 1:n1
                        obj(i).X1 = p1(i, 1);
                        obj(i).Y1 = p1(i, 2);
                        obj(i).Z1 = p1(i, 3);
                        obj(i).X2 = p2(i, 1);
                        obj(i).Y2 = p2(i, 2);
                        obj(i).Z2 = p2(i, 3);
                end
                
            end
        end

    end

end % end constructors

%% Geometry methods
methods
    function point = planeIntersection(obj, plane, varargin)
        % Return intersection point between a plane and this line segment.
        %
        % INTERS = planeIntersection(LINESEG, PLANE);
        %
        % Support array of line segments.
        
        % extract tolerance for determination of parallel edges and planes
        tol = 1e-12;
        if ~isempty(varargin)
            tol = varargin{1};
        end
        
        % initialize empty arrays
        nEdges = numel(obj);
        nPlanes = numel(plane);
        if nPlanes > 1
            error('multiple planes not supported');
        end
        
        % initialize empty arrays
        point(nEdges, nPlanes) = Point3D();
        t = zeros(nEdges, 1);
        
        % plane normal
        n = normal(plane);
        
        % origin and direction of edge supporting line
        p1 = firstPoint(obj);
        lineDir = Vector3D(p1, lastPoint(obj));
        
% %         p2 = reshape(Point3D([[obj.X2]' [obj.Y2]' [obj.Z2]']), size(obj));
%         line = StraightLine3D(firstPoint(obj), lastPoint(obj));
%         lineDir = direction(line);
        
        % get indices of edge and plane which are parallel
        par = abs(dotProduct(n, lineDir)) < tol;
        if any(par)
            point(par).X = NaN;
            point(par).Y = NaN;
            point(par).Z = NaN;
            t(par) = NaN;
        end
        
        % difference between origins of plane and edge
        dp = Vector3D(p1, repmat(origin(plane), size(p1)));
        
        % relative position of intersection on line
        %t = dot(n(~par,:), dp(~par,:), 2)./dot(n(~par,:), line(~par,4:6), 2);
        t(~par) = dotProduct(n, dp(~par)) ./ dotProduct(n, lineDir(~par));
        
        % compute coord of intersection point
        %point(~par, :) = line(~par,1:3) + repmat(t,1,3).*line(~par,4:6);
        point(~par) = p1(~par) + lineDir(~par) * t(~par);
        
        % set points outside of edge to [NaN NaN NaN]
        if any(t < 0)
            point(t<0).X = NaN;
            point(t<0).Y = NaN;
            point(t<0).Z = NaN;
        end
        if any(t > 1)
            point(t>1).X = NaN;
            point(t>1).Y = NaN;
            point(t>1).Z = NaN;
        end
    end
end

%% Methods generic to curve objects
methods
    function l = length(obj)
        dx = [obj.X2] - [obj.X1];
        dy = [obj.Y2] - [obj.Y1];
        dz = [obj.Z2] - [obj.Z1];
        l = sqrt(dx.*dx + dy.*dy + dz.*dz);
        l = reshape(l, size(obj));
    end
    
    function res = reverse(obj)
        res = LineSegment3D(lastPoint(obj), firstPoint(obj));
    end
    
    function p1 = firstPoint(obj)
        p1 = reshape(Point3D([[obj.X1]' [obj.Y1]' [obj.Z1]']), size(obj));
    end
    
    function p2 = lastPoint(obj)
        p2 = reshape(Point3D([[obj.X2]' [obj.Y2]' [obj.Z2]']), size(obj));
    end
end


%% Methods implementing the Geometry2D interface
methods
    function res = transform(obj, transfo)
        % Apply a geometric transform to this line segment.
        p1t = transformPoint(transfo, [obj.X1 obj.Y1 obj.Z1]);
        p2t = transformPoint(transfo, [obj.X2 obj.Y2 obj.Z2]);
        res = LineSegment3D(p1t, p2t);
    end
    
    function box = bounds(obj)
        % Returns the bounding box of this geometry.
        x = sort([obj.X1 obj.X2]);
        y = sort([obj.Y1 obj.Y2]);
        z = sort([obj.Z1 obj.Z2]);
        box = Bounds3D([x y z]);
    end
    
    function h = draw(varargin)
        % Draws the current geometry, eventually specifying the style.

        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        
        % plot line segment
        xdata = [obj.X1 obj.X2];
        ydata = [obj.Y1 obj.Y2];
        zdata = [obj.Z1 obj.Z2];
        hh = plot3(ax, xdata, ydata, zdata, varargin{:});
        
        if ~isempty(style)
            apply(style, hh);
        end
        
        if nargout > 0
            h = hh;
        end
    end
    
    function res = scale(obj, factor)
        % Returns a scaled version of this geometry.
        res = LineSegment3D([obj.X1 obj.Y1 obj.Z2] * factor, [obj.X2 obj.Y2 obj.ZZ] * factor);
    end
   
    function res = translate(obj, shift)
        % Returns a translated version of this geometry.       
        res = LineSegment3D([obj.X1 obj.Y1 obj.Z2] + shift, [obj.X2 obj.Y2 obj.ZZ] + shift);
    end
    
end % end methods


%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'LineSegment3D', ...
            'X1', obj.X1, 'Y1', obj.Y1, 'Z1', obj.Z1, ...
            'X2', obj.X2, 'Y2', obj.Y2, 'Z2', obj.Z2);
    end
end
methods (Static)
    function line = fromStruct(str)
        % Create a new instance from a structure.
        p1 = [str.X1 str.Y1 str.Z1];
        p2 = [str.X2 str.Y2 str.Z2];
        line = LineSegment3D(p1, p2);
    end
end

end % end classdef

