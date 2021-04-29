classdef Ellipsoid3D < Geometry3D
% An ellipsoid defined by center, size and orientation.
%
%   Class Ellipsoid3D
%
%   Example
%   Ellipsoid3D
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2020-03-12,    using Matlab 9.7.0.1247435 (R2019b) Update 2
% Copyright 2020 INRAE - BIA-BIBS.


%% Properties
properties
    % The coordinates of the center.
    Center = [0 0 0];

    % The three radiusses (semi-axis lengths), from largest to smallest.
    Radius = [1 1 1];
    
    % The three Euler angle, in degrees.
    % The correspond to the azimut of the main axis, to the elevation of
    % the main axis, and to the rotation around the main axis, in that
    % order.
    EulerAngles = [0 0 0];
    
end % end properties


%% Static factories
methods (Static)
    function elli = equivalentEllipsoid(points)
        % Compute the Equivalent ellipsoid of a set of 3D points.
        %
        % ELL = equivalentEllipsoid(PTS)
        % Compute the equivalent ellipsoid of the set of points PTS.
        %
        % Example
        %     pts = MultiPoint3D(randn(300, 3));
        %     pts = transform(pts, AffineTransform3D.createScaling([6 4 2]));
        %     pts = transform(pts, AffineTransform3D.createRotationOx(pi/6));
        %     pts = transform(pts, AffineTransform3D.createRotationOy(pi/4));
        %     pts = transform(pts, AffineTransform3D.createRotationOz(pi/3));
        %     pts = transform(pts, AffineTransform3D.createTranslation([5 4 3]));
        %     elli = Ellipsoid3D.equivalentEllipsoid(pts);
        %     figure; hold on; axis equal;
        %     draw(pts); draw(elli);
        %     draw(elli, ...
        %         'drawEllipses', true, 'EllipseColor', 'b', 'EllipseWidth', 3);
        %
        %   See also
        %     Sphere3D, Ellipse2D.equivalentEllipse, draw
        %     rotationMatrixToEulerAngles 
        
        % ensure input is a numeric array
        if isa(points, 'MultiPoint3D')
            points = points.Coords;
        end
        
        % number of points
        n = size(points, 1);
        
        % compute centroid
        center = mean(points);
        
        % compute the covariance matrix
        covPts = cov(points)/n;
        
        % perform a principal component analysis with 3 variables,
        % to extract equivalent axes
        [U, S] = svd(covPts);
        
        % extract length of each semi axis
        radii = sqrt(5) * sqrt(diag(S)*n)';
        
        % sort axes from greater to lower
        [radii, ind] = sort(radii, 'descend');
        
        % format U to ensure first axis points to positive x direction
        U = U(ind, :);
        if U(1,1) < 0
            U = -U;
            % keep matrix determinant positive
            U(:,3) = -U(:,3);
        end
        
        % convert axes rotation matrix to Euler angles
        angles = rotationMatrixToEulerAngles(U);
        
        % concatenate result to form an ellipsoid object
        elli = Ellipsoid3D([center, radii, angles]);
    end
end


%% Constructor
methods
    function obj = Ellipsoid3D(varargin)
        % Constructor for Ellipsoid3D class.
        %
        %   ELLI = Ellipsoid3D();
        %   Empty constructor, initialized to an ellipsoid centered at
        %   (0,0,0) and with all radiusses equal to 1.
        %
        %   ELLI = Ellipsoid3D(DATA);
        %   Initialize from a 1-by-9 row vector in the format: 
        %   ELLI = [CX CY CZ R1 R2 R3 PHI THETA PSI]
        %   
        %   ELLI = Ellipsoid3D(ELL0);
        %   Copy constructor from another Ellipsoid3D object. ELL0 can also
        %   be  sphere object.
        %
        
        if nargin == 0
            % Empty contructor
            return;
        end
        
        if nargin == 1
            var1 = varargin{1};
            if isa(var1, 'Ellipsoid3D')
                % copy constructor
                obj.Center = var1.Center;
                obj.Radius = var1.Radius;
                obj.EulerAngles = var1.EulerAngles;
                
            elseif isa(var1, 'Sphere3D')
                % copy constructor from a sphere.
                obj.Center = var1.Center;
                obj.Radius = [var1.Radius var1.Radius var1.Radius];
                
            elseif isnumeric(var1)
                % initialize with a 1-by-9 row vector
                obj.Center = var1(1:3);
                obj.Radius = var1(4:6);
                obj.EulerAngles = var1(7:9);
            end
            
        elseif nargin == 2
            % center + radius
            
            % initialize center
            var1 = varargin{1};
            if isa(var1, 'Point3D')
                obj.Center = [var1.X var1.Z var1.Y];
            elseif isnumeric(var1)
                obj.Center = var1(1, 1:3);
            else
                error('Can not interpret first argument');
            end
            
            % initialize radius
            var2 = varargin{2};
            if ~isnumeric(var2) || any(size(var2) == [1 3])
                error('Second argument must be a 1-by-3 row vector');
            end
            obj.Radius = var2(1, 1:3);
            
        elseif nargin > 2
            error('Can not initialize 3D ellipsoid');
        end
    end

end % end constructors


%% Methods implementing the Geometry3D interface
methods
    function res = transform(obj, transform) %#ok<STOUT,INUSD>
        % Apply a geometric transform to this geometry.
        error('Method not implemented');
    end
    
    function bnd = bounds(obj)
        % Return the (approximated) bounds of this ellipsoid.

        nPhi    = 64;
        nTheta  = 32;
        [x, y, z] = surfaceVertices(obj, nPhi, nTheta);
        bnd = Bounds3D([min(x(:)) max(x(:)) min(y(:)) max(y(:)) min(z(:)) max(z(:))]);
    end
    
    function h = draw(varargin)
        %DRAW Draw the ellipsoid, eventually specifying the style.
        
        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        
        % default values for drawing
        nPhi    = 32;
        nTheta  = 16;
        
        % process input options: when a string is found, assumes this is the
        % beginning of options
        options = {'FaceColor', 'g', 'LineStyle', 'none'};
        if length(varargin) == 1
            options = {'FaceColor', varargin{1}, 'LineStyle', 'none'};
        else
            options = [options varargin];
        end
        
        [x, y, z] = surfaceVertices(obj, nPhi, nTheta);
     
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
    
    function [x, y, z] = surfaceVertices(obj, nPhi, nTheta)
        % Compute coordinates of vertices used to draw ellipsoid surface.
        %
        % [X, Y, Z] = surfaceVertices(ELL, NPHI, NTHETA);
        %
        
        % convert unit basis to ellipsoid basis
        sca  = AffineTransform3D.createScaling(obj.Radius);
        rotZ = AffineTransform3D.createRotationOz(obj.EulerAngles(1) * pi / 180);
        rotY = AffineTransform3D.createRotationOy(obj.EulerAngles(2) * pi / 180);
        rotX = AffineTransform3D.createRotationOx(obj.EulerAngles(3) * pi / 180);
        tra  = AffineTransform3D.createTranslation(obj.Center);
        
        % concatenate transforms
        trans = tra * rotZ * rotY * rotX * sca;
        
        % parametrisation of ellipsoid in spherical coordinates
        theta   = linspace(0, pi, nTheta+1);
        phi     = linspace(0, 2*pi, nPhi+1);
        
        % convert to cartesian coordinates
        sintheta = sin(theta);
        x = cos(phi') * sintheta;
        y = sin(phi') * sintheta;
        z = ones(length(phi),1) * cos(theta);
        
        % transform mesh vertices
        pts2 = transformPoint(trans, [x(:) y(:) z(:)]);
        x(:) = pts2(:, 1);
        y(:) = pts2(:, 2);
        z(:) = pts2(:, 3);
    end
end


%% Methods implementing the Geometry3D interface (more)
methods
    function res = scale(obj, varargin)
        % Return a scaled version of this ellipsoid.
        factor = varargin{1};
        if ~isscalar(factor)
            error('Requires scaling factor to be a scalar');
        end
        center = obj.Center * factor;
        radius = obj.Radius * factor;
        data = [center radius obj.EulerAngles];
        res = Ellipsoid3D(data);
    end
    
    function res = translate(obj, varargin)
        % Return a translated version of this ellipsoid.
        shift = varargin{1};
        center = obj.Center + shift;
        data = [center obj.Radius obj.EulerAngles];
        res = Ellipsoid3D(data);
    end    
end % end methods


%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'Ellipsoid3D', ...
            'Center', obj.Center, ...
            'Radius', obj.Radius, ...
            'EulerAngles', obj.EulerAngles);
    end
end

methods (Static)
    function elli = fromStruct(str)
        % Create a new instance from a structure.
        elli = Ellipsoid3D([str.Center str.Radius str.EulerAngles]);
    end
end

end % end classdef

