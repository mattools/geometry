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
                n1 = size(var1, 1);
                n2 = size(var1, 2);
                obj(n1, n2) = Vector3D();
                for i = 1:numel(var1)
                    obj(i).X = var1(i).X;
                    obj(i).Y = var1(i).Y;
                    obj(i).Z = var1(i).Z;
                end
                
            elseif isnumeric(var1) && size(var1, 2) == 3
                % Initialize from a 1-by-3 numeric array
                nv = size(var1, 1);
                obj(nv, 1) = Vector3D();
                for i = 1:nv
                    obj(i).X = var1(i, 1);
                    obj(i).Y = var1(i, 2);
                    obj(i).Z = var1(i, 3);
                end
                
            else
                error('Can not parse input for Vector3D');
            end
            
        elseif nargin == 2
            var1 = varargin{1};
            var2 = varargin{2};
            if isa(var1, 'Point3D') && isa(var2, 'Point3D')
                % create vector as the difference between two points
                
                % compute size of output array
                dim1 = size(var1);
                dim2 = size(var2);
                [dim, nd] = computeOutputArraySize(dim1, dim2);
                
                % create output array
                obj(prod(dim)) = Vector3D();
                obj = reshape(obj, dim);
                
                % process each item in output array
                for i = 1:prod(dim)
                    % convert linear index of output to sub indices
                    inds = cell(nd,1);
                    [inds{:}] = ind2sub(dim, i);
                    inds = [inds{:}];
                    
                    % convert sub indices to linear index of inputs
                    ind1 = linearIndex(dim1, inds);
                    ind2 = linearIndex(dim2, inds);
                    
                    % store result in inner properties
                    obj(i).X = var2(ind2).X - var1(ind1).X;
                    obj(i).Y = var2(ind2).Y - var1(ind1).Y;
                    obj(i).Z = var2(ind2).Z - var1(ind1).Z;
                end
            else
                error('Can not parse input for Vector3D');
            end
            
        elseif nargin == 3
            % initialisation from three numeric arguments
            for i = 1:3
                var_i = varargin{i};
                if ~isnumeric(var_i) || size(var_i, 2) ~= 1
                    error('Requires the three input argument to be numeric with only one column');
                end
            end
            var1 = varargin{1};
            var2 = varargin{2};
            var3 = varargin{3};
            if size(var1, 1) ~= size(var2, 1) || size(var1, 1) ~= size(var3, 1)
                error('Requires the three inputs to have same number of rows');
            end
            
            nv = size(var1, 1);
            obj(nv, 1) = Vector3D();
            for i = 1:nv
                obj(i).X = var1(i);
                obj(i).Y = var2(i);
                obj(i).Z = var3(i);
            end
        else
            error('Wrong number of input arguments.');
        end
        
    end

end % end constructors


%% Methods specific to Vector2D
methods
    function p = dotProduct(v1, v2)
        % Dot product of two 3D vectors.
        p = reshape([v1.X], size(v1)) .* reshape([v2.X], size(v2)) + ...
            reshape([v1.Y], size(v1)) .* reshape([v2.Y], size(v2)) + ...
            reshape([v1.Z], size(v1)) .* reshape([v2.Z], size(v2));
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

        v1 = [vertcat(vect1.X) vertcat(vect1.Y) vertcat(vect1.Z)];
        v2 = [vertcat(vect2.X) vertcat(vect2.Y) vertcat(vect2.Z)];
        v = v1(:, [2 3 1]) .* v2(:, [3 1 2]) - v2(:, [2 3 1]) .* v1(:, [3 1 2]);
        res = reshape(Vector3D(v), size(vect1));
    end
    
    function n = norm(obj)
        % Return the norm of this vector.
        %
        % Example
        %   v = Vector3D([5 4 3]);
        %   norm(v)
        %   ans = 
        %       5
        n = reshape(hypot(hypot([obj.X]', [obj.Y]'), [obj.Z]'), size(obj));
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
        
        vx = reshape(vertcat(obj.X), size(obj));
        vy = reshape(vertcat(obj.Y), size(obj));
        vz = reshape(vertcat(obj.Z), size(obj));
        n = hypot(hypot(vx, vy), vz);
        vn = Vector3D(vx ./ n, vy ./ n, vz ./ n);
    end
    
end % end methods


%% Access methods
methods
    function coords = coordinates(obj)
        % Get the coordinates of this vector as a numeric array.
        coords = [[obj.X]' [obj.Y]' [obj.Z]'];
    end
    
end

%% Linear algebra methods for vectors
methods
    function res = plus(obj, obj2)
        % Implement plus operator for Vector3D objects.
        res = Vector3D(...
            reshape([obj.X], size(obj)) + reshape([obj2.X], size(obj2)), ...
            reshape([obj.Y], size(obj)) + reshape([obj2.Y], size(obj2)), ...
            reshape([obj.Z], size(obj)) + reshape([obj2.Z], size(obj2)));
    end
    
    function res = minus(obj, obj2)
        % Implement minus operator for Vector3D objects.
        res = Vector3D(...
            reshape([obj.X], size(obj)) - reshape([obj2.X], size(obj2)), ...
            reshape([obj.Y], size(obj)) - reshape([obj2.Y], size(obj2)), ...
            reshape([obj.Z], size(obj)) - reshape([obj2.Z], size(obj2)));
    end
    
    function res = mtimes(obj, k)
        % Implement times operator for Vector3D objects.
        
        if isa(k, 'Vector3D') && isnumeric(obj)
            % swap inputs
            tmp = k; k = obj; obj = tmp;
        end
        
        res = reshape(Vector3D([[obj.X]' [obj.Y]' [obj.Z]'] .* k), size(obj));
    end
    
    function res = mrdivide(obj, k)
        % Implement divides operator for Vector3D objects.
        res = reshape(Vector3D([[obj.X]' [obj.Y]' [obj.Z]'] ./ k), size(obj));
    end
    
    function res = uminus(obj)
        res = reshape(Vector3D(-[[obj.X]' [obj.Y]' [obj.Z]']), size(obj));
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


%% Utility functions

function [dim, nd] = computeOutputArraySize(dim1, dim2)
% Compute size of output array from size of inputs.
nd1 = length(dim1);
nd2 = length(dim2);
nd = max(nd1, nd2);
dim = ones(1, nd);
dim(1:nd1) = dim1;
dim(1:nd2) = max(dim(1:nd2), dim2);
end

function ind = linearIndex(dim, inds)
nd = length(dim);
tmp = num2cell(min(inds(1:nd), dim));
ind = sub2ind(dim, tmp{:});
end
