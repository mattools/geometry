classdef LineSegment2D < Geometry2D
% A line segment defined by its two extremities.
%
%   Class LineSegment2D
%
%   Example
%     P1 = Point2D(20, 10);
%     P2 = Point2D(40, 20);
%     L = LineSegment2D(P1, P2);
%     draw(L);
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-10-15,    using Matlab 9.7.0.1190202 (R2019b)
% Copyright 2019 INRA - BIA-BIBS.


%% Properties
properties
    % the source point
    P1;
    % the target point
    P2;
    
end % end properties


%% Constructor
methods
    function obj = LineSegment2D(varargin)
        % Constructor for LineSegment2D class
        
        if nargin == 0
            % Default constructor: unit line segment
            obj.P1 = Point2D(0,0);
            obj.P2 = Point2D(1,0);
            
        elseif nargin == 1
            % Copy constructor
            if ~isa(varargin{1}, 'LineSegment2D')
                error('Requires a LineSegment2D as input');
            end
            var1 = varargin{1};
            obj.P1 = var1.P1;
            obj.P2 = var1.P2;
            
        elseif nargin == 2
            obj.P1 = varargin{1};
            obj.P2 = varargin{2};
        end

    end

end % end constructors


%% Methods implementing the Curve2D interface
methods
end


%% Methods implementing the Geometry2D interface
methods
    function box = boundingBox(obj)
        % Returns the bounding box of this geometry
        x = sort([obj.P1.X obj.P2.X]);
        y = sort([obj.P1.Y obj.P2.Y]);
        box = Box2D([x y]);
    end
    
    function varargout = draw(obj, varargin)
        % Draws the current geometry, eventually specifying the style

        % parse style info
        style = [];
        ind = cellfun(@(x)isa(x, 'Style'), varargin);
        if any(ind)
            style = varargin{ind};
            varargin(ind) = [];
        end
        
        % plot line segment
        xdata = [obj.P1.X obj.P2.X];
        ydata = [obj.P1.Y obj.P2.Y];
        h = plot(xdata, ydata, varargin{:});
        
        if ~isempty(style)
            apply(style, h);
        end
        
        if nargout > 0
            varargout = {h};
        end
    end
    
    function res = scale(obj, factor)
        % Returns a scaled version of this geometry
        res = LineSegment2D(scale(obj.P1, factor), scale(obj.P2, factor));
    end
    
    function res = translate(obj, shift)
        % Returns a translated version of this geometry       
        res = LineSegment2D(translate(obj.P1, shift), translate(obj.P2, shift));
    end
    
    function res = rotate(obj, varargin)
        % Returns a rotated version of this geometry
        origin = [0 0];
        if ~isempty(varargin)
            origin = varargin{1};
        end
        
        rot = createRotation(origin, deg2rad(angle));
        p1t = transformPoint(obj.P1, rot);
        p2t = transformPoint(obj.P2, rot);
        
        res =  LineSegment2D(p1t, p2t);
    end
    
end % end methods


%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization
        str = struct('Type', 'LineSegment2D', 'P1', toStruct(obj.P1), 'P2', toStruct(obj.P2));
    end
end
methods (Static)
    function line = fromStruct(str)
        % Create a new instance from a structure
        p1 = Point2D.fromStruct(str.P1);
        p2 = Point2D.fromStruct(str.P2);
        line = LineSegment2D(p1, p2);
    end
end

end % end classdef

