classdef Box3D < handle
% Bounding box of a 3D shape.
%
%   Class Box3D
%   Defined by max extent in each dimension:
%   * XMin, XMax, YMin, YMax, ZMin, ZMax
%
%   Example
%   box = Box3D([0, 10, 0, 10, 0, 10])
%
%   See also
%     Geometry3D, Box2D
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-02-06,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2019 INRA - BIA-BIBS.


%% Properties
properties
    XMin;
    XMax;
    YMin;
    YMax;
    ZMin;
    ZMax;

end % end properties


%% Constructor
methods
    function obj = Box3D(varargin)
        % Constructor for Box3D class.
        
        if ~isempty(varargin)
            var1 = varargin{1};
            if size(var1, 1) ~= 1
                error('Creating a box requires an array with one row, not %d', size(var1, 1));
            end
            if size(var1, 2) ~= 6
                error('Creating a box requires an array with six columns, not %d', size(var1, 2));
            end
            data = var1;
        else
            % default box is unit square, with origin as lower-left corner.
            data = [0 1 0 1 0 1];
        end
        
        obj.XMin = data(1);
        obj.XMax = data(2);
        obj.YMin = data(3);
        obj.YMax = data(4);
        obj.ZMin = data(5);
        obj.ZMax = data(6);
    end

end % end constructors


%% Methods
methods
    function box = boundingBox(obj)
        % Return the bounding box of obj shape.
        box = Box2D([obj.XMin obj.XMax obj.YMin obj.YMax obj.ZMin obj.ZMax]);
    end
    
    function varargout = draw(obj, varargin)
        %DRAW Draw the box on the current axis.
        
        % extract style agument if present
        style = [];
        if nargin > 1 && isa(varargin{1}, 'Style')
            style = varargin{1};
            varargin(1) = [];
        end
        
        % draw the box
        h = drawBox3d([obj.XMin obj.XMax obj.YMin obj.YMax obj.ZMin obj.ZMax], varargin{:});
        
        % eventually apply style
        if ~isempty(style)
            apply(style, h);
        end
        
        % return handle if requested
        if nargout > 0
            varargout = {h};
        end
    end
    
    function res = scale(obj, varargin)
        % Return a scaled version of this box.
        factor = varargin{1};
        res = Box3D([obj.XMin obj.XMax obj.YMin obj.YMax obj.ZMin obj.ZMax] * factor);
    end
    
    function res = translate(obj, varargin)
        % Return a translated version of this box.
        shift = varargin{1};
        data2 = [obj.XMin obj.XMax obj.YMin obj.YMax obj.ZMin obj.ZMax] + shift(1, [1 1 2 2 3 3]);
        res = Box3D(data2);
    end
end % end methods


%% Serialization methods
methods
    function write(obj, fileName, varargin)
        %WRITE Write box representation into a JSON file.
        %
        % Requires implementation of the "toStruct" method.
        if exist('savejson', 'file') == 0
            error('Requires the ''jsonlab'' library');
        end
        savejson('', toStruct(obj), 'FileName', fileName, varargin{:});
    end
    
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('type', 'Box3D', ...
            'XMin', obj.XMin, 'XMax', obj.XMax, ...
            'YMin', obj.YMin, 'YMax', obj.YMax, ...
            'ZMin', obj.ZMin, 'ZMax', obj.ZMax);
    end
end

methods (Static)
    function box = read(fileName)
        % Reads box information from a file in JSON format
        if exist('loadjson', 'file') == 0
            error('Requires the ''jsonlab'' library');
        end
        box = Box3D.fromStruct(loadjson(fileName));
    end
    
    function box = fromStruct(str)
        % Creates a new instance from a structure
        box = Box3D([str.XMin str.XMax str.YMin str.YMax str.ZMin str.ZMax]);
    end
end

end % end classdef

