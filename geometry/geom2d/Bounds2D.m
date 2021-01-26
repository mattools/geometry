classdef Bounds2D < handle
% The bounds of a planar shape in each direction.
%
%   Class Bounds2D
%   Defined by max extent in each dimension:
%   * XMin, XMax, YMin, YMax.
%
%   Example
%     box = Bounds2D([5 15 6 14]);
%     figure; axis([0 20 0 20]); hold on
%     draw(box, 'b')
%
%   See also
%     Geometry2D, Box3D

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-08-14,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2013 INRA - Cepia Software Platform.


%% Properties
properties
    XMin;
    XMax;
    YMin;
    YMax;
    
end % end properties


%% Constructor
methods
    function obj = Bounds2D(varargin)
        % Constructor for Bounds2D class.
    
        if ~isempty(varargin)
            var1 = varargin{1};
            if size(var1, 1) ~= 1
                error('Creating Bounds2D requires an array with one row, not %d', size(var1, 1));
            end
            if size(var1, 2) ~= 4
                error('Creating Bounds2D requires an array with four columns, not %d', size(var1, 2));
            end
            data = var1;
        else
            % default box is unit square, with origin as lower-left corner.
            data = [0 1 0 1];
        end
        
        obj.XMin = data(1);
        obj.XMax = data(2);
        obj.YMin = data(3);
        obj.YMax = data(4);
    end

end % end constructors


%% Methods
methods
    function b = isFinite(obj)
        b = all(isfinite([obj.XMin obj.XMax  obj.YMin obj.YMax]));
    end
    
    function varargout = draw(obj, varargin)
        %DRAW Draw the current geometry, eventually specifying the style.
        
        % extract style agument if present
        style = [];
        if nargin > 1 && isa(varargin{1}, 'Style')
            style = varargin{1};
            varargin(1) = [];
        end
        
        % draw the box
        tx = [obj.XMin obj.XMax obj.XMax obj.XMin obj.XMin];
        ty = [obj.YMin obj.YMin obj.YMax obj.YMax obj.YMin];
        h = plot(tx, ty, varargin{:});
        
        % eventually apply style
        if ~isempty(style)
            apply(style, h);
        end
        
        % return handle if requested
        if nargout > 0
            varargout = {h};
        end
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
        str = struct('type', 'Bounds2D', ...
            'XMin', obj.XMin, 'XMax', obj.XMax, ...
            'YMin', obj.YMin, 'YMax', obj.YMax);
    end
end

methods (Static)
    function box = read(fileName)
        %READ Read box information from a file in JSON format.
        if exist('loadjson', 'file') == 0
            error('Requires the ''jsonlab'' library');
        end
        box = Bounds2D.fromStruct(loadjson(fileName));
    end
    
    function box = fromStruct(str)
        % Create a new instance from a structure.
        box = Bounds2D([str.XMin str.XMax str.YMin str.YMax]);
    end
end

end % end classdef

