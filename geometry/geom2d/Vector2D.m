classdef Vector2D < handle
% A Vector in 2D plane.
%
%   Class Vector2D
%
%   Example
%   Vector2D
%
%   See also
%     Point2D
%

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2020-03-22,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2020 INRAE - BIA-BIBS.


%% Properties
properties
    % The X-coordinate of this vector.
    X = 0;
    
    % The Y-coordinate of this vector.
    Y = 0;
    
end % end properties


%% Constructor
methods
    function obj = Vector2D(varargin)
        % Constructor for Vector2D class
        
        if nargin == 0
            % empty constructor -> already initialized to origin
            
        elseif nargin == 1
            var1 = varargin{1};
            if isa(var1, 'Point2D') || isa(var1, 'Vector2D')
                % copy constructor, or initialize from a point
                var1 = varargin{1};
                obj.X = var1.X;
                obj.Y = var1.Y;
                
            elseif isnumeric(var1) && ~any(size(var1) ~= [1 2])
                % Initialize from a 1-by-2 numeric array
                obj.X = var1(1);
                obj.Y = var1(2);
            else
                error('Can not parse input for Vector');
            end
            
        elseif nargin == 2
            var1 = varargin{1};
            var2 = varargin{2};
            if isnumeric(var1) && isnumeric(var2) && isscalar(var1) && isscalar(var2)
                % initialisation from two scalar numeric argument
                obj.X = var1;
                obj.Y = var2;
                
            elseif isa(var1, 'Point2D') && isa(var2, 'Point2D')
                % create vector as the difference between two points
                obj.X = var2.X - var1.X;
                obj.Y = var2.Y - var1.Y;
                
            else
                error('Can not parse inputs for Vector2D');
            end
            
        else
            error('Wrong number of input arguments.');
        end
        
    end
    
end % end constructors


%% Methods specific to Vector2D
methods
    function n = norm(obj)
        % Return the norm of this vector.
        %
        % Example
        %   v = Vector2D([4 3]);
        %   norm(v)
        %   ans = 
        %       5
        n = hypot(obj.X, obj.Y);
    end
    
    function vn = normalize(obj)
        % Return the unit norm vector with same direction.
        %
        % Example
        %   v = Vector2D([4 3]);
        %   vn = normalize(v);
        %   norm(vn)
        %   ans = 
        %       1
        n = hypot(obj.X, obj.Y);
        vn = Vector2D(obj.X / n, obj.Y / n);
    end
    
end % end methods


%% Methods implementing the Geometry2D interface (more)
methods
    function res = scale(obj, varargin)
        % Return a scaled version of this geometry.
        factor = varargin{1};
        res = Vector2D([obj.X obj.Y] * factor);
    end
    
    function res = translate(obj, varargin)
        % Return a translated version of this geometry.
        shift = varargin{1};
        res = Point2D(bsxfun(@plus, [obj.X obj.Y], shift));
    end
    
    function res = rotate(obj, angle)
        % Return a rotated version of this vector.
        %
        % V2 = rotate(V, THETA)
        % THETA is given in radians, in counter clockwise order.
        
        % precomputes angles
        cot = cos(angle);
        sit = sin(angle);
        
        % compute rotated coordinates
        vx = cot * obj.X - sit * obj.Y;
        vy = sit * obj.X + cot * obj.Y;
        
        res = Vector2D([vx vy]);
    end
    
    function res = rotate90(obj, n)
        % Apply rotation(s) by 90 degrees to this vector.
        %
        %   V2 = rotate90(V);
        %   V2 = rotate90(V, NROT);
        
        if nargin == 1
            n = 1;
        end
        vx = obj.X;
        vy = obj.Y;
        for i = 1:n
            tmp = vx;
            vx = -vy;
            vy = tmp;
        end
        res = Vector2D([vx vy]);
    end
    
    function res = plus(obj, obj2)
        % Implement plus operator for Vector2D objects.
        res = Vector2D(obj.X+obj2.X, obj.Y+obj2.Y);
    end
    
    function res = minus(obj, obj2)
        % Implement minus operator for Vector2D objects.
        res = Vector2D(obj.X-obj2.X, obj.Y-obj2.Y);
    end
    
    function res = mtimes(obj, k)
        % Implement times operator for Vector2D objects.
        res = Vector2D(obj.X * k, obj.Y * k);
    end
    
    function res = mrdivide(obj, k)
        % Implement divides operator for Vector2D objects.
        res = Vector2D(obj.X / k, obj.Y / k);
    end
    
    function res = uminus(obj)
        res = Vector2D(-obj.X, -obj.Y);
    end
end % end methods

%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('type', 'Vector2D', 'X', obj.X, 'Y', obj.Y);
    end
end
methods (Static)
    function vect = fromStruct(str)
        % Create a new instance from a structure.
        vect = Vector2D(str.X, str.Y);
    end
end

end % end classdef

