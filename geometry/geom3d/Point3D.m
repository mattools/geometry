classdef (InferiorClasses = {?matlab.graphics.axis.Axes, ?Vector3D}) Point3D < Geometry3D
% A point in the 3-dimensional space.
%
%   Usage:
%   P = Point3D(COORDS)
%   where COORDS is a 1-by-3 array of numeric values
%   P = Point3D(X, Y, Z)
%   where each of X, Y and Z are numeric scalars
%   P = Point3D(PT)
%   where PT is another instance of Point3D
%
%   Example
%   Point3D
%
%   See also
%     Geometry3D, MultiPoint3D, Point2D

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-02-07,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2019 INRA - BIA-BIBS.

%% Static factories
methods (Static)
    function p = centroid(varargin)
        % Compute the centroid of several points.
        xc = 0;
        yc = 0;
        zc = 0;
        for i = 1:length(varargin)
            var_i = varargin{i};
            if ~isa(var_i, 'Point3D')
                error('Requires all input to be Point3D');
            end
            xc = xc + var_i.X;
            yc = yc + var_i.Y;
            zc = zc + var_i.Z;
        end
        p = Point3D([xc yc zc] / length(varargin));
    end
end


%% Properties
properties
    X = 0;
    Y = 0;
    Z = 0;
end % end properties


%% Constructor
methods
    function obj = Point3D(varargin)
    % Constructor for Point3D class.

        % empty constructor -> initialize to origin
        if isempty(varargin)
            return;
        end
        
        % copy constructor
        if isa(varargin{1}, 'Point3D')
            that = varargin{1};
            obj.X = that.X;
            obj.Y = that.Y;
            obj.Z = that.Z;
            return;
        end
        
        % initialisation constructor with one argument
        if nargin == 1
            var1 = varargin{1};
            if isnumeric(var1) && size(var1, 2) == 3
                np = size(var1, 1);
                obj = Point3D.empty(np, 0);
                for ip = 1:np
                    obj(ip).X = var1(ip,1);
                    obj(ip).Y = var1(ip,2);
                    obj(ip).Z = var1(ip,3);
                end
            else
                error('Can not parse input for Point3D');
            end
        end
        
        % initialisation from three scalar numeric arguments
        if nargin == 3
            var1 = varargin{1};
            var2 = varargin{2};
            var3 = varargin{3};
            if isnumeric(var1) && isnumeric(var2) && isnumeric(var3) && isscalar(var1) && isscalar(var2) && isscalar(var3)
                obj.X = var1;
                obj.Y = var2;
                obj.Z = var3;
            else
                error('Can not parse inputs for Point3D');
            end
        end
    end

end % end constructors

%% Methods specific to Point3D
methods
    function res = plus(obj, obj2)
        % Implement plus operator for points and vectors.
        if ~isa(obj2, 'Vector3D')
            error('Requires second argument to be a Vector3D');
        end
        res = Point3D(obj.X+obj2.X, obj.Y+obj2.Y, obj.Z+obj2.Z);
    end

end

%% Methods implementing the Geometry3D interface
methods
    function res = transform(obj, transform)
        % Apply a geometric transform to this geometry.
        res = Point3D(transformPoint(transform, [obj.X obj.Y obj.Z]));
    end
    
    function box = bounds(obj)
        % Return the bounding box of this shape.
        box = Bounds3D([obj.X obj.X obj.Y obj.Y obj.Z obj.Z]);
    end
    
    function h = draw(varargin)
        %DRAW Draw this point, eventually specifying the style.

        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        
        % draw the geometric primitive
        hh = plot3(ax, obj.X, obj.Y, obj.Z, varargin{:});

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
        res = Point3D([obj.X obj.Y obj.Z] * factor);
    end
    
    function res = translate(obj, varargin)
        % Return a translated version of this geometry.
        shift = varargin{1};
        res = Point3D(bsxfun(@plus, [obj.X obj.Y obj.Z], shift));
    end    
end % end methods


%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'Point3D', 'X', obj.X, 'Y', obj.Y, 'Z', obj.Z);
    end
end
methods (Static)
    function point = fromStruct(str)
        % Create a new instance from a structure.
        point = Point3D([str.X str.Y str.Z]);
    end
end

end % end classdef
