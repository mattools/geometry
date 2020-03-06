classdef Sphere3D < Geometry3D
% A sphere, defined by a center and a radius.
%
%   Class Sphere3D
%
%   Example
%   Sphere3D
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2020-03-06,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2020 INRAE - BIA-BIBS.


%% Properties
properties
    XC = 0;
    YC = 0;
    ZC = 0;
    R  = 1;
end % end properties


%% Constructor
methods
    function obj = Sphere3D(varargin)
        % Constructor for Sphere3D class
        if nargin == 1
            var1 = varargin{1};
            if isa(var1, 'Sphere3D')
                % copy constructor
                obj.XC = var1.XC;
                obj.YC = var1.YC;
                obj.ZC = var1.ZC;
                obj.R  = var1.R;
                
            elseif isnumeric(var1)
                % initialize with a 1-by-4 row vector
                obj.XC = var1(1);
                obj.YC = var1(2);
                obj.ZC = var1(3);
                obj.R  = var1(4);
            end
        elseif nargin == 2
            % center + radius
            
            % initialize center
            var1 = varargin{1};
            if isa(var1, 'Point3D')
                obj.XC = var1.X;
                obj.YC = var1.Y;
                obj.ZC = var1.Z;
            elseif isnumeric(var1)
                obj.XC = var1(1);
                obj.YC = var1(2);
                obj.ZC = var1(3);
            else
                error('Can not interpret first argument');
            end
            
            % initialize radius
            obj.R = varargin{2};
            
        elseif nargin > 2
            error('Can not initialize sphere');
        end
    end

end % end constructors



%% Methods implementing the Geometry3D interface
methods
    function res = transform(obj, transform) %#ok<STOUT,INUSD>
        % Apply a geometric transform to this geometry.
        error('Method not implemented');
    end
    
    function box = boundingBox(obj)
        % Return the bounding box of this shape.
        
        box = Box3D([...
            obj.XC - obj.R obj.XC + obj.R ...
            obj.YC - obj.R obj.YC + obj.R ...
            obj.ZC - obj.R obj.ZC + obj.R]);
    end
    
    function h = draw(varargin)
        %DRAW Draw this point, eventually specifying the style.

        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        
        % defualt values for drawing
        nPhi    = 32;
        nTheta  = 16;

        % process input options: when a string is found, assumes this is the
        % beginning of options
        options = {'FaceColor', 'g', 'LineStyle', 'none'};
        if length(varargin) == 1
            options = {'FaceColor', varargin{1}, 'LineStyle', 'none'};
            varargin= {};
        else
            options = [options varargin];
        end

        % compute spherical coordinates
        theta   = linspace(0, pi, nTheta+1);
        phi     = linspace(0, 2*pi, nPhi+1);
        
        % convert to cartesian coordinates
        sintheta = sin(theta);
        x = obj.XC + cos(phi') * sintheta * obj.R;
        y = obj.YC + sin(phi') * sintheta * obj.R;
        z = obj.ZC + ones(length(phi),1) * cos(theta) * obj.R;

        % draw the geometric primitive
        hh = surf(ax, x, y, z, options{:});

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
        center = [obj.XC obj.YC obj.ZC] * factor;
        res = Sphere3D(center, obj.R * factor);
    end
    
    function res = translate(obj, varargin)
        % Return a translated version of this geometry.
        shift = varargin{1};
        center = [obj.XC obj.YC obj.ZC] + shift;
        res = Sphere3D(center, obj.R);
    end    
end % end methods


%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'Sphere3D', 'XC', obj.XC, 'YC', obj.YC, 'ZC', obj.ZC, 'R', obj.R);
    end
end
methods (Static)
    function point = fromStruct(str)
        % Create a new instance from a structure.
        point = Sphere3D([str.XC str.YC str.ZC str.R]);
    end
end

end % end classdef

