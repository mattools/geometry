classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) Bounds3D < handle
% Bounding box of a 3D shape.
%
%   Class Bounds3D
%   Defined by min and max extent in each dimension:
%   * XMin, XMax, YMin, YMax, ZMin, ZMax
%
%   Example
%     bnd = Bounds3D([0, 10, 0, 10, 0, 10])
%
%   See also
%     Geometry3D, Bounds2D
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-02-06,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2019 INRA - BIA-BIBS.

%% Static factories

methods (Static)
    function res = fromAxis(varargin)
        % Create a Bounds3D from the current axis bounds.
        %
        % BNDS = Bounds3D.fromAxis;
        % BNDS = Bounds3D.fromAxis(AX);
        %
        % Example:
        %   % create 3D bounds for axis object of an empty figure
        %   close all;
        %   bounds = Bounds3D.fromAxis(gca)
        %   bounds = 
        %     Bounds3D with properties:
        %       XMin: 0
        %       XMax: 1
        %       YMin: 0
        %       YMax: 1
        %       ZMin: 0
        %       ZMax: 1
        
        if ~isempty(varargin)
            hAx = varargin{1};
        else
            hAx = gca;
        end
        
        % extract axis bounds to crop plane
        xlim = get(hAx, 'xlim');
        ylim = get(hAx, 'ylim');
        zlim = get(hAx, 'zlim');
        res = Bounds3D([xlim ylim zlim]);
    end
    
end

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
    function obj = Bounds3D(varargin)
        % Constructor for Bounds3D class.
        
        if ~isempty(varargin)
            var1 = varargin{1};
            if size(var1, 1) ~= 1
                error('Creating a Bounds3D requires an array with one row, not %d', size(var1, 1));
            end
            if size(var1, 2) ~= 6
                error('Creating a Bounds3D requires an array with six columns, not %d', size(var1, 2));
            end
            data = var1;
        else
            % default box is unit square, with origin as lower-left corner.
            data = [0 1  0 1  0 1];
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
    function b = isFinite(obj)
        b = all(isfinite([obj.XMin obj.XMax  obj.YMin obj.YMax  obj.ZMin obj.ZMax]));
    end
end

%% Methods
methods
    function box = boundingBox(obj)
        % Return the bounding box of this shape.
        box = Box2D([obj.XMin obj.XMax obj.YMin obj.YMax obj.ZMin obj.ZMax]);
    end
    
    function varargout = draw(varargin)
        %DRAW Draw the box on the current axis.
        
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        if isempty(varargin)
            varargin = {'color', 'k'};
        end
        
        % create corner points
        p000 = Point3D(obj.XMin, obj.YMin, obj.ZMin);
        p100 = Point3D(obj.XMax, obj.YMin, obj.ZMin);
        p010 = Point3D(obj.XMin, obj.YMax, obj.ZMin);
        p110 = Point3D(obj.XMax, obj.YMax, obj.ZMin);
        p001 = Point3D(obj.XMin, obj.YMin, obj.ZMax);
        p101 = Point3D(obj.XMax, obj.YMin, obj.ZMax);
        p011 = Point3D(obj.XMin, obj.YMax, obj.ZMax);
        p111 = Point3D(obj.XMax, obj.YMax, obj.ZMax);
        
        % draw the box, by drawing each of the 12 edges as line segments
        % first the lower face (z=zmin)
        he(1) = draw(ax, LineSegment3D(p000, p100), varargin{:});
        he(2) = draw(ax, LineSegment3D(p000, p010), varargin{:});
        he(3) = draw(ax, LineSegment3D(p100, p110), varargin{:});
        he(4) = draw(ax, LineSegment3D(p010, p110), varargin{:});

        % left face (x=ymin)
        he(5) = draw(ax, LineSegment3D(p000, p001), varargin{:});
        he(6) = draw(ax, LineSegment3D(p010, p011), varargin{:});
        he(7) = draw(ax, LineSegment3D(p001, p011), varargin{:});

        % front face (y=xmin)
        he(8) = draw(ax, LineSegment3D(p100, p101), varargin{:});
        he(9) = draw(ax, LineSegment3D(p001, p101), varargin{:});
        
        % the last 3 remaining edges
        he(10) = draw(ax, LineSegment3D(p110, p111), varargin{:});
        he(11) = draw(ax, LineSegment3D(p101, p111), varargin{:});
        he(12) = draw(ax, LineSegment3D(p011, p111), varargin{:});
        
        gh = hggroup;
        set(he, 'Parent', gh);
        
        % eventually apply style
        if ~isempty(style)
            apply(style, he);
        end
        
        % return handle if requested
        if nargout > 0
            varargout = {ge};
        end
    end
    
    function res = scale(obj, varargin)
        % Return a scaled version of this box.
        factor = varargin{1};
        res = Bounds3D([obj.XMin obj.XMax obj.YMin obj.YMax obj.ZMin obj.ZMax] * factor);
    end
    
    function res = translate(obj, varargin)
        % Return a translated version of this box.
        shift = varargin{1};
        data2 = [obj.XMin obj.XMax obj.YMin obj.YMax obj.ZMin obj.ZMax] + shift(1, [1 1 2 2 3 3]);
        res = Bounds3D(data2);
    end
end % end methods

methods (Access = protected)
    function [ax, obj, style, varargin] = parseDrawInputArguments(varargin)
        % Return the different elements necessary to draw the object.
        %
        % Usage:
        %   [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        %
        % Returns the following:
        % - 'ax' is the handle of the axis to draw in, that can be used as
        %     first input for plot functions. If not specified, the current
        %     axis is returned.
        % - 'obj' is the instance of the object, that should be a subclass
        %     of Geometry
        % - 'style' is an optional 'Style', that can be used to update the
        %     drawing style of the graphical object.
        % - 'varargin' are the remaining input arguments.
        %
        % See Also
        %   parseDrawInputArguments in Geometry class
        
        % identify the variable corresponding to class instance
        ind = cellfun(@(x)isa(x, 'Bounds3D'), varargin);
        obj = varargin{ind};
        varargin(ind) = [];
        
        % extract handle of axis to draw on
        if ~isempty(varargin) && isscalar(varargin{1}) && ishghandle(varargin{1}) && strcmpi(get(varargin{1}, 'type'), 'axes')
            ax = varargin{1};
            varargin(1) = [];
        else
            ax = gca;
        end
        
        % parse optional style info
        style = [];
        ind = cellfun(@(x) isa(x, 'Style'), varargin);
        if any(ind)
            style = varargin{ind};
            varargin(ind) = [];
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
        str = struct('type', 'Bounds3D', ...
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
        box = Bounds3D.fromStruct(loadjson(fileName));
    end
    
    function box = fromStruct(str)
        % Creates a new instance from a structure
        box = Bounds3D([str.XMin str.XMax str.YMin str.YMax str.ZMin str.ZMax]);
    end
end

end % end classdef

