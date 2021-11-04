classdef Vector3D < handle
% A vector in 3-dimensional space.
%
%   Class Vector3D
%
%   Example
%   Vector3D
%
%   See also
%     Point3D, StraightLine3D

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2021-11-04,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE - BIA-BIBS.


%% Properties
properties
    % The X-coordinate of this vector.
    X = 0;
    
    % The Y-coordinate of this vector.
    Y = 0;
    
    % The Z-coordinate of this vector.
    Z = 0;
    
end % end properties


%% Constructor
methods
    function obj = Vector3D(varargin)
        % Constructor for Vector3D class.
        
        if nargin == 0
            % empty constructor -> already initialized to origin
            
        elseif nargin == 1
            var1 = varargin{1};
            if isa(var1, 'Point3D') || isa(var1, 'Vector3D')
                % copy constructor, or initialize from a point
                var1 = varargin{1};
                obj.X = var1.X;
                obj.Y = var1.Y;
                obj.Z = var1.Z;
                
            elseif isnumeric(var1) && ~any(size(var1) ~= [1 3])
                % Initialize from a 1-by-3 numeric array
                obj.X = var1(1);
                obj.Y = var1(2);
                obj.Z = var1(3);
            else
                error('Can not parse input for Vector3D');
            end
            
        elseif nargin == 2
            var1 = varargin{1};
            var2 = varargin{2};
            if isa(var1, 'Point3D') && isa(var2, 'Point3D')
                % create vector as the difference between two points
                obj.X = var2.X - var1.X;
                obj.Y = var2.Y - var1.Y;
                obj.Z = var2.Z - var1.Z;
            else
                error('Can not parse input for Vector3D');
            end
            
        elseif nargin == 3
            % initialisation from three scalar numeric argument
            for i = 1:3
                var_i = varargin{i};
                if ~isnumeric(var_i) || ~isscalar(var_i)
                    error('Requires the three input argument to be numeric scalar values');
                end
            end
            
            obj.X = varargin{1};
            obj.Y = varargin{2};
            obj.Z = varargin{3};
        else
            error('Wrong number of input arguments.');
        end
        
    end

end % end constructors


%% Methods specific to Vector2D
methods
    function p = dotProduct(v1, v2)
        % Dot product of two 3D vectors.
        p = v1.X * v2.X + v1.Y * v2.Y + v1.Z * v2.Z;
    end
    
    function res = crossProduct(vect1, vect2)
        % Dot product of two 3D vectors.
        % 
        % Example
        %   v1 = Vector3D([2 0 0]);
        %   v2 = Vector3D([0 3 0]);
        %   v  = crossProduct(v1, v2)
        %   v = 
        %     Vector3D with properties:
        %       X: 0
        %       Y: 0
        %       Z: 6

        v1 = [vect1.X vect1.Y vect1.Z];
        v2 = [vect2.X vect2.Y vect2.Z];
        v = v1([2 3 1]) .* v2([3 1 2]) - v2([2 3 1]) .* v1([3 1 2]);
        res = Vector3D(v);
    end
    
    function n = norm(obj)
        % Return the norm of this vector.
        %
        % Example
        %   v = Vector3D([5 4 3]);
        %   norm(v)
        %   ans = 
        %       5
        n = hypot(hypot(obj.X, obj.Y), obj.Z);
    end
    
    function vn = normalize(obj)
        % Return the unit norm vector with same direction.
        %
        % Example
        %   v = Vector3D([5 4 3]);
        %   vn = normalize(v);
        %   norm(vn)
        %   ans = 
        %       1
        
        n = hypot(hypot(obj.X, obj.Y), obj.Z);
        vn = Vector3D(obj.X / n, obj.Y / n, obj.Z / n);
    end
    
end % end methods

%% Linear algebra methods for vectors
methods
    function res = plus(obj, obj2)
        % Implement plus operator for Vector3D objects.
        res = Vector3D(obj.X+obj2.X, obj.Y+obj2.Y, obj.Z+obj2.Z);
    end
    
    function res = minus(obj, obj2)
        % Implement minus operator for Vector3D objects.
        res = Vector3D(obj.X-obj2.X, obj.Y-obj2.Y, obj.Z-obj2.Z);
    end
    
    function res = mtimes(obj, k)
        % Implement times operator for Vector3D objects.
        res = Vector3D(obj.X * k, obj.Y * k, obj.Z * k);
    end
    
    function res = mrdivide(obj, k)
        % Implement divides operator for Vector3D objects.
        res = Vector3D(obj.X / k, obj.Y / k, obj.Z / k);
    end
    
    function res = uminus(obj)
        res = Vector3D(-obj.X, -obj.Y, -obj.Z);
    end
end


%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('type', 'Vector3D', 'X', obj.X, 'Y', obj.Y, 'Z', obj.Z);
    end
end
methods (Static)
    function vect = fromStruct(str)
        % Create a new instance from a structure.
        vect = Vector3D(str.X, str.Y, str.Z);
    end
end

end % end classdef

