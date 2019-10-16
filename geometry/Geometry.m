classdef Geometry < handle
% Abstract class for a geometry of any dimensionality.
%
%   Class Geometry
%
%   Example
%   Geometry
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-08-13,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2018 INRA - BIA-BIBS.


%% Properties
properties
end % end properties


%% Constructor
methods
end % end constructors


%% Protected methods for general use
methods
    function [ax, obj, style, varargin] = parseDrawOptions(varargin)
        % Return the different elements necessary to draw the object.
        %
        % Usage:
        % [ax, obj, style, otherOptions] = parseDrawOptions(varargin{:});
        %
        
        % extract handle of axis to draw on
        arg = varargin{1};
        if isscalar(arg) && ishandle(arg) && strcmpi(get(arg, 'type'), 'axes')
            ax = varargin{1};
            varargin(1) = [];
        else
            ax = gca;
        end

        % parse optional style info
        style = [];
        ind = cellfun(@(x)isa(x, 'Style'), varargin);
        if any(ind)
            style = varargin{ind};
            varargin(ind) = [];
        end
        
        % assumes obj is first argument of the remaining ones
        obj = varargin{1};
        
        % update varargin array that can contains additional options
        varargin(1) = [];
    end
    
end % end methods

%% Serialization methods
methods
    function write(obj, fileName, varargin)
        % Write geometry into a JSON file.
        % 
        % Requires implementation of the "toStruct" method.
        
        if exist('savejson', 'file') == 0
            error('Requires the ''jsonlab'' library');
        end
        if ~ismethod(obj, 'toStruct')
            error('Requires implementation of the ''toStruct'' method');
        end
        
        savejson('', toStruct(obj), 'FileName', fileName, varargin{:});
    end
end

methods (Static)
    function geom = fromStruct(str)
        % Create a new transform instance from a structure.
        
        % check existence of 'type' field
        if isfield(str, 'Type')
            type = str.Type;
        elseif isfield(str, 'type')
            type = str.type;
        else
            error('Requires a field with name "Type"');
        end

        % parse transform
        try
            geom = eval([type '.fromStruct(str)']);
        catch ME
            error(['Unable to parse Geometry with type: ' type]);
        end
    end
    
    function geom = read(fileName)
        %READ Read a geometry from a file in JSON format.
        if exist('loadjson', 'file') == 0
            error('Requires the ''jsonlab'' library');
        end
        geom = Geometry.fromStruct(loadjson(fileName));
    end
end

end % end classdef

