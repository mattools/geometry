classdef StraightLine3D < Geometry3D
% A 3D straight line, unbounded in each direction.
%
%   Class StraightLine3D
%
%   Example
%     P1 = Point3D([40 10 20]);
%     P2 = Point3D([20 40 30]);
%     L = StraightLine3D(P1, P2);
%     figure; axis equal;axis([0 50 0 50 0 50]; hold on;
%     draw(L); view(3);
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2021-11-04,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE - BIA-BIBS.


%% Properties
properties
    % The origin of this straight line, as a 1-by-3 numeric array.
    Origin = [0 0 0];
    
    % The direction vector of this straight line, as a 1-by-3 numeric array.
    Direction = [1 0 0];
end % end properties


%% Constructor
methods
    function obj = StraightLine3D(varargin)
        % Constructor for StraightLine3D class.
        
        if nargin == 0
            % empty constructor -> already initialized to origin
            
        elseif nargin == 1
            var1 = varargin{1};
            if isa(var1, 'StraightLine3D')
                % copy constructor
                var1 = varargin{1};
                obj.Origin    = var1.Origin;
                obj.Direction = var1.Direction;
                
            elseif isnumeric(var1) && all(size(var1) == [1 6])
                % numeric input consider [X0 Y0 Z0 DX DY DZ].
                obj.Origin    = [var1(1) var1(2) var1(3)];
                obj.Direction = [var1(4) var1(5) var1(6)];
            else
                error('Can not parse input for StraightLine3D');
            end
            
        elseif nargin == 2
            % initialisation from two arguments
            
            % Check both inputs have same size
            p1 = varargin{1};
            p2 = varargin{2};
            if any(size(p1) ~= size(p2))
                error('Both input arguments must have same size');
            end
            
            if isa(p1, 'Point3D') && isa(p2, 'Point3D')
                n1 = size(p1, 1);
                n2 = size(p1, 2);
                obj(n1, n2) = StraightLine3D();
                for i = 1:numel(p1)
                    obj(i).Origin = coordinates(p1(i));
                    obj(i).Direction = coordinates(p2(i)) - obj(i).Origin;
                end
                
            elseif isa(p1, 'Point3D') && isa(p2, 'Vector3D')
                n1 = size(p1, 1);
                n2 = size(p1, 2);
                obj(n1, n2) = StraightLine3D();
                for i = 1:numel(p1)
                    obj(i).Origin = coordinates(p1(i));
                    obj(i).Direction = coordinates(p2(i));
                end
                
            elseif isnumeric(p1) && size(p1, 3) == 3 && isnumeric(p2) && size(p2, 3) == 3
                % when  inputs are numeric, consider second one corresponds
                % to another point  
                n1 = size(p1, 1);
                obj(n1, 1) = StraightLine3D();
                for i = 1:numel(p1)
                    obj(i).Origin = p1(i,:);
                    obj(i).Direction = p1(2,:) - p1(i,:);
                end
                
%             end
%             % parse origin
%             var1 = varargin{1};
%             if isa(var1, 'Point3D')
%                 obj.Origin = [var1.X var1.Y var1.Z];
%             elseif isnumeric(var1)
%                 obj.Origin = var1(1, 1:3);
%             else
%                 error('Can not interpret first argument');
%             end
            
%             % second argument can be either another point, or the direction
%             var2 = varargin{2};
%             if isa(var2, 'Point3D')
%                 obj.Direction = [var2.X var2.Y var2.Z] - obj.Origin;
%             elseif isa(var2, 'Vector3D')
%                 obj.Direction = [var2.X var2.Y var2.Z];
%             elseif isnumeric(var2)
%                 % numeric input consider another point as default.
%                 obj.Direction = var2 - obj.Origin;
            else
                error('Can not interpret input arguments');
            end
            
        else
            error('Wrong number of input arguments.');
        end
        
    end

end % end constructors


%% Information Methods
methods
    function p = origin(obj)
        nl = numel(obj);
        p = Point3D(reshape([obj.Origin], [3 nl])');
        p = reshape(p, size(obj));
    end
    
    function v = direction(obj)
        nl = numel(obj);
        v = Vector3D(reshape([obj.Direction], [3 nl])');
        v = reshape(v, size(obj));
    end
end

%% General processing methods
methods
    function pos = position(obj, point)
        % Compute position of a point on this line.
        %
        %  POS = line.position(PT);
        
        % vector from line origin to point
        ptCoords = [[point.X]' [point.Y]' [point.Z]'];
        dp = bsxfun(@minus, ptCoords, obj.Origin);
        
        % direction vector of the line
        vl = obj.Direction;
        
        % precompute and check validity of denominator
        denom = sum(vl.^2, 2);
        if sqrt(denom) < eps
            error('degenerated line');
        end
        
        % compute position using dot product normalized with norm of line vector.
        pos = bsxfun(@rdivide, sum(bsxfun(@times, dp, vl), 2), denom);
    end
    
    function lineSeg = clip(obj, bounds)
        % Clip this line with a Bounds3D object and return a line segment.
        
        % extract limits of the box
        xmin = bounds.XMin; xmax = bounds.XMax;
        ymin = bounds.YMin; ymax = bounds.YMax;
        zmin = bounds.ZMin; zmax = bounds.ZMax;
        
        % extreme corners of the box
        p000 = Point3D([xmin ymin zmin]);
        p111 = Point3D([xmax ymax zmax]);
        
        % main direction vectors
        ex   = Vector3D([1 0 0]);
        ey   = Vector3D([0 1 0]);
        ez   = Vector3D([0 0 1]);
        
        % box faces parallel to Oxy
        planeZ0 = Plane3D(p000, ex, ey);
        planeZ1 = Plane3D(p111, ex, ey);
        
        % box faces parallel to Oxz
        planeY0 = Plane3D(p000, ex, ez);
        planeY1 = Plane3D(p111, ex, ez);
        
        % box faces parallel to Oyz
        planeX0 = Plane3D(p000, ey, ez);
        planeX1 = Plane3D(p111, ey, ez);


        % compute intersection point with each plane
        [points(1), valid(1)] = intersectLine(planeZ0, obj);
        [points(2), valid(2)] = intersectLine(planeZ1, obj);
        [points(3), valid(3)] = intersectLine(planeY0, obj);
        [points(4), valid(4)] = intersectLine(planeY1, obj);
        [points(5), valid(5)] = intersectLine(planeX0, obj);
        [points(6), valid(6)] = intersectLine(planeX1, obj);
        
        % keep only valid points
        points = points(valid);
        
        % sort points with respect to their position
        pos = position(obj, points);
        [pos, inds] = sort(pos); %#ok<ASGLU>
        points  = points(inds);

        % keep median points wrt to position. These points define the limit
        % of the clipped edge.
        nv = length(inds)/2;
        
        % create resulting edge.
        lineSeg = LineSegment3D(points(nv), points(nv+1));
        
        % check that middle point of the edge is contained in the box
        midPoint = Point3D.centroid(points(nv), points(nv+1));


        % check that middle point of the edge is contained in the box
        xOk  = xmin <= midPoint.X & midPoint.X <= xmax;
        yOk  = ymin <= midPoint.Y & midPoint.Y <= ymax;
        zOk  = zmin <= midPoint.Z & midPoint.Z <= zmax;

        % if one of the bounding condition is not met, set edge to NaN
        inds = ~(xOk & yOk & zOk);
        if any(inds)
            lineSeg(inds).X1 = NaN;
        end
        
    end

end


%% Methods implementing the Geometry3D interface
methods
    function bnd = bounds(obj) %#ok<MANU>
        % Returns infinite bounds in each direction.
        bnd = Bounds3D([-inf inf -inf inf -inf inf]);
    end
    
    function h = draw(varargin)
        %DRAW Draw this line, eventually specifying the style.
        
        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        
        % default style for drawing lines
        if length(varargin) ~= 1
            varargin = [{'color', 'b'}, varargin];
        end

        % extract bounding box of the current axis
        xlim = get(ax, 'xlim');
        ylim = get(ax, 'ylim');
        zlim = get(ax, 'zlim');
        
        % clip lines with current axis box
        edge = clip(obj, Bounds3D([xlim ylim zlim]));
        
        if ~isnan(edge.X1)
            % draw current clipped line
            hh = draw(ax, edge, varargin{:});
        
            if ~isempty(style)
                apply(style, hh);
            end
        else
            hh = -1;
        end
        
        if nargout > 0
            h = hh;
        end
    end
    
    function res = transform(obj, transform)
        % Apply a geometric transform to this geometry.
        origin = transformPoint(transform, obj.Origin);
        direction = transformVector(transform, obj.Direction);
        res = StraightLine3D(origin, direction);
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
        origin2 = origin(obj) * factor;
        direction2 = direction(obj) * factor;
        res = StraightLine3D(origin2, direction2);
    end
    
    function res = translate(obj, varargin)
        % Return a translated version of this geometry.
        shift = varargin{1};
        origin2 = origin(obj) + shift;
        res = StraightLine3D(origin2, direction(obj));
    end    
end % end methods



%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'StraightLine3D', ...
            'Origin', obj.Origin, ...
            'Direction', obj.Direction);
    end
end
methods (Static)
    function line = fromStruct(str)
        % Create a new instance from a structure.
        line = StraightLine3D(str.Origin, str.Origin+str.Direction);
    end
end

end % end classdef

