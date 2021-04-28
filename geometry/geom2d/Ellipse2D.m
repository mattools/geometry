classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) Ellipse2D < Curve2D
% A planar ellipse.
%
%   An ellipse is defined by five parameters:
%   * CenterX     the x-coordinate of the center
%   * CenterY     the y-coordinate of the center
%   * Radius1     the length of the major semi-axis
%   * Radius2     the length of the minor semi-axis
%   * Orientation the orientation of the major axis
%
%   
%   Example
%     % Create and display a basic ellipse
%     ELLI = Ellipse2D([50 50  40 20  30]);
%     figure; axis([0 100 0 100]); hold on; axis square;
%     draw(ELLI, 'b');
%     draw(center(ELLI), 'bo');
%
%   See also
%     Circle2D, fit, asPolyline, draw
%

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2019-05-17,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2019 INRA - BIA-BIBS.

%% Static factories
methods (Static)
    function elli = fit(points)
        % FIT Fit the best ellipse from a set of points.
        %
        %  ELLI = Ellipse2D.fit(PTS)
        %  PTS can be either a N-by-2 array of point coordinates, or an
        %  instance of MultiPoint2D.
        %
        % Example
        %     % Fit an ellipse from jittered points around a ref ellipse
        %     elli = Ellipse2D([50 40 30 20 10]);
        %     figure; draw(elli, 'k');
        %     hold on; axis equal; axis([0 100 0 80]);
        %     poly = asPolyline(elli, 100);
        %     pts = poly.Coords + randn(100, 2);
        %     draw(MultiPoint2D(pts), 'b+');
        %     elli2 = Ellipse2D.fit(pts);
        %     draw(elli2, 'm');
        %     legend({'Original ellipse', 'Points', 'Fitted ellipse'});
        % 
        % Reference
        % https://fr.mathworks.com/matlabcentral/fileexchange/3215-fit_ellipse
        %
        % See also
        %   fromCartesianCoeffs
        
        % parse input
        if isa(points, 'MultiPoint2D')
            points = points.Coords;
        end
        
        % tolerance for assuming orientation zero
        tol = 1e-3;

        % recenter coordinates to improve numerical accuracy
        center = mean(points);
        xi = points(:,1) - center(1);
        yi = points(:,2) - center(2);
        
        % build coefficients
        X = [xi.*xi  xi.*yi  yi.*yi  xi  yi];
        coeffs = sum(X) / (X'*X);
        
        % extract parameters from the conic equation
        A = coeffs(1);
        B = coeffs(2);
        C = coeffs(3);
        D = coeffs(4);
        E = coeffs(5);
        
        % compute coefficients of ellipse oriented with main axis
        if min(abs(B / A), abs(B / C)) > tol
            % compute orientation of ellipse
            theta = 0.5 * atan(B / (A - C));
            cost = cos(theta);
            sint = sin(theta);
            
            % compute coefficients of axis-oriented ellipse
            [A, C, D, E] = deal(...
                A * cost^2 + B * cost * sint + C * sint^2, ...
                A * sint^2 - B * cost * sint + C * cost^2, ...
                D * cost + E * sint, ...
                D * sint - E * cost);
            
            % rotation matrix from rotated space to original space
            R = [cost -sint; sint cost];
            % new center after rotation
            center = center * R; % same as: center = (R' * center')'
        else
            theta = 0;
            R = [1 0;0 1];
        end
        
        % checkup
        if A * C <= 0
            error('Ellipse2D.fit: points seem to correspond to another conic type');
        end
        
        % ensure A coefficient is positive
        if A < 0
            [A, C, D, E] = deal(-A, -C, -D, -E);
        end
        
        % adjust position of center (in axis-aligned basis)
        center = center - [ D/(2*A) E/(2*C) ];
        
        % rotate the axes backwards to find the center point of the original TILTED ellipse
        center = center * R'; % same as: center = (R * center')'
        
        % identify ellipse radiusses
        F   = 1 + (D^2) / (4 * A) + (E^2) / (4 * C);
        R1  = sqrt(F / A);
        R2  = sqrt(F / C);
        if R1 < R2
            [R1, R2] = deal(R2, R1);
        end
        
        % create ellipse instance to encapsulate parameters
        elli = Ellipse2D([center R1 R2 rad2deg(theta)]);
    end
    
    function elli = fromCartesianCoeffs(coeffs)
        % Identify an ellipse from its cartesian coefficients.
        %
        % References
        % https://en.wikipedia.org/wiki/Ellipse#Standard_parametric_representation
        %
        % See Also
        %   cartesianCoefficients, fit

        % get coefficients by their usual names
        A = coeffs(1);
        B = coeffs(2);
        C = coeffs(3);
        D = coeffs(4);
        E = coeffs(5);
        F = coeffs(6);
        
        % retrieve center
        delta = B * B - 4 * A * C;
        xc = (2 * C * D - B * E) / delta;
        yc = (2 * A * E - B * D) / delta;
        
        % find orientation
        theta = 0.5 * atan(B / (A - C));
        
        % retrieve length of semi-axes
        common = 2 * (A * E^2 + C * D^2 - B * D * E + delta * F);
        root = sqrt((A - C)^2 + B^2);
        a1 = -sqrt(common * ((A + C) + root)) / delta;
        a2 = -sqrt(common * ((A + C) - root)) / delta;
        
        elli = Ellipse2D([xc yc a1 a2 rad2deg(theta)]);
    end
end

%% Properties
properties
    CenterX = 0;
    CenterY = 0;
    Radius1 = 1;
    Radius2 = 1;
    % The orientation in degrees, CCW
    Orientation = 0;
end % end properties


%% Constructor
methods
    function obj = Ellipse2D(varargin)
    % Constructor for Ellipse2D class.
    % 
    % Syntax
    %   ELLI = Ellipse2D([CX CY  A B  THETA]);
    %   Creates a new Ellipse2D object, with cente rgiven by CX and CY,
    %   length of semi-axes given by A and B, and orientation given by
    %   THETA (in degrees, counter-clockwise).
    %
    % Example
    %   ELLI = Ellipse2D([50 50  40 20  30]);
    %   figure; axis([0 100 0 100]); hold on; axis square;
    %   draw(ELLI, 'b');
    %   draw(center(ELLI), 'bo');
    %

        switch nargin
            case 0
                % nothing to do
            case 1
                var1 = varargin{1};
                if size(var1, 2) ~= 5
                    error('Creating an ellipse requires an array with five columns, not %d', size(var1, 2));
                end
                obj.CenterX = var1(1);
                obj.CenterY = var1(2);
                obj.Radius1 = var1(3);
                obj.Radius2 = var1(4);
                obj.Orientation = var1(5);
        end
    end

end % end constructors


%% Methods specific to Ellipse2D
methods
    function a = area(obj)
        % The area of the ellipse.
        a = pi * obj.Radius1 * obj.Radius2;
    end
    
    function p = perimeter(obj, varargin)
        % An approximation of the perimeter length of this ellipse.
        
        % relative tolerance
        tol = 1e-10;
        if ~isempty(varargin)
            tol = varargin{1};
        end
        
        ra = obj.Radius1;
        rb = obj.Radius2;
        
        % function to integrate
        f = @(t) sqrt(ra .^ 2 .* cos(t) .^ 2 + rb .^ 2 .* sin(t) .^ 2) ;
        
        % absolute tolerance from relative tolerance
        eps = tol * max(ra, rb);
        
        % integrate on first quadrant
        % (Requires Matlab 2012a)
        p = 4 * integral(f, 0, pi/2, 'AbsTol', eps);
    end
    
    function center = center(obj)
        % Return the ellipse center as a Point2D.
        %
        % See Also
        %   Point2D, area, bounds
        center = Point2D(obj.CenterX, obj.CenterY);
    end
    
    function poly = asPolyline(obj, varargin)
        % Convert this ellipse into a (closed) polyline.
        %
        % POLY = asPolyline(OBJ);
        % POLY = asPolyline(OBJ, NPTS);
        % Returns the result as an instance of LinearRing2D. Can specify
        % the number of vertices of the polyline as second argument.
        %
        % See Also
        %   Polyline2D, perimeter
        
        % determines number of points
        N = 72;
        if ~isempty(varargin)
            N = varargin{1};
        end
        
        % create time basis
        t = linspace(0, 2*pi, N+1)';
        t(end) = [];
        
        % get ellipse parameters
        xc = obj.CenterX;
        yc = obj.CenterY;
        r1 = obj.Radius1;
        r2 = obj.Radius2;
        theta = obj.Orientation;
        
        % pre-compute trig functions (angles is in degrees)
        cot = cosd(theta);
        sit = sind(theta);
        
        % position of points
        x = xc + r1 * cos(t) * cot - r2 * sin(t) * sit;
        y = yc + r1 * cos(t) * sit + r2 * sin(t) * cot;
        
        poly = LinearRing2D([x y]);
    end
    
    function params = cartesianCoefficients(obj)
        % Return the six parameters of the cartesian representation.
        %
        % References
        % https://en.wikipedia.org/wiki/Ellipse#Standard_parametric_representation
        %
        % See Also
        %   fromCartesianCoeffs
        
        % get ellipse parameters
        xc = obj.CenterX;
        yc = obj.CenterY;
        a2 = obj.Radius1 * obj.Radius1;
        b2 = obj.Radius2 * obj.Radius2;

        % pre-compute trig functions (angles is in degrees)
        cot = cosd(obj.Orientation);
        sit = sind(obj.Orientation);
        
        % identification of each parameter
        A = a2 * sit * sit + b2 * cot * cot;
        B = 2 * (b2 - a2) * sit * cot;
        C = a2 * cot * cot + b2 * sit * sit;
        D = - 2 * A * xc - B * yc;
        E = - B * xc - 2 * C * yc;
        F = A * xc * xc + B * xc * yc + C * yc * yc - a2 * b2;
        
        % concatenate into a single row vector
        params = [A B C D E F];
    end
end

%% Methods implementing the Geometry2D interface
methods
    function res = transform(obj, transform)
        % Apply a geometric transform to this ellipse.
        %
        % Example
        %     % apply an arbitrary transform to a simple ellipse
        %     elli = Ellipse2D([5 4 3 2 0]);
        %     rot = AffineTransform2D.createRotation(pi/6);
        %     sca = AffineTransform2D.createScaling([2.5 1.5]);
        %     tra = AffineTransform2D.createTranslation([4 3]);
        %     transfo = sca * rot * tra;
        %     elli2 = transform(elli, transfo);
        %     % display original and transformed ellipses
        %     figure; hold on; axis equal; axis([0 20 0 20]);
        %     draw(elli, 'k');
        %     draw(elli2, 'b');
        %     % Compare with transform on polygonal approximation
        %     draw(transform(asPolyline(elli, 100), transfo), 'm');
        %
        % Reference
        % https://math.stackexchange.com/questions/3076317/what-is-the-equation-for-an-ellipse-in-standard-form-after-an-arbitrary-matrix-t
        
        % first extract coefficients of cartesian representation
        coeffs = cartesianCoefficients(obj);
        
        % writing the matrix of the general conic equation x^t * Q * x = 0
        Q = [...
             coeffs(1)   coeffs(2)/2  coeffs(4)/2; ...
            coeffs(2)/2   coeffs(3)   coeffs(5)/2; ...
            coeffs(4)/2  coeffs(5)/2   coeffs(6)];
        
        % compute the matrix form of the transformed ellipse
        Minv = inv(affineMatrix(transform));
        Q2 = Minv' * Q * Minv;
        
        coeffs2 = [Q2(1,1) 2*Q2(1,2) Q2(2,2) 2*Q2(1,3) 2*Q2(2,3) Q2(3,3)];
        res = Ellipse2D.fromCartesianCoeffs(coeffs2);
    end
    
    function box = bounds(obj)
        % Return the bounds of this ellipse.
        %
        % See Also
        %   Bounds2D, center
        
        extX = [obj.CenterX - obj.Radius obj.CenterX + obj.Radius];
        extY = [obj.CenterY - obj.Radius obj.CenterY + obj.Radius];
        box = Bounds2D([extX extY]);
    end
    
    function h = draw(varargin)
        %DRAW Draw the ellipse, eventually specifying the style.
        %
        % See also
        %   drawAxes, center
        
        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        
        % default drawing argument
        if isempty(varargin)
            varargin = {'b-'};
        end
        
        % get ellipse parameters
        xc = obj.CenterX;
        yc = obj.CenterY;
        r1 = obj.Radius1;
        r2 = obj.Radius2;
        theta = obj.Orientation;
        
        % pre-compute trig functions (angles is in degrees)
        cot = cosd(theta);
        sit = sind(theta);
        
        % compute position of several points along ellipse
        t = linspace(0, 2*pi, 145);
        xt = xc + r1 * cos(t) * cot - r2 * sin(t) * sit;
        yt = yc + r1 * cos(t) * sit + r2 * sin(t) * cot;
        
        % stores handle to graphic object
        hh = plot(ax, xt, yt, varargin{:});
        
        if ~isempty(style)
            apply(style, hh);
        end
        
        if nargout > 0
            h = hh;
        end
    end
    
    function h = drawAxes(varargin)
        % Draw the axes of ellipse, eventually specifying the style.
        %
        % See also
        %   draw
       
        % parse arguments using protected method implemented in Geometry
        [ax, obj, style, varargin] = parseDrawInputArguments(varargin{:});
        
        % default drawing argument
        if isempty(varargin)
            varargin = {'b-'};
        end
        
        % get ellipse parameters
        xc = obj.CenterX;
        yc = obj.CenterY;
        r1 = obj.Radius1;
        r2 = obj.Radius2;
        theta = obj.Orientation;
        
        % pre-compute trig functions (angles is in degrees)
        cot = cosd(theta);
        sit = sind(theta);
        
        % compute position of several points along ellipse
        t = [0 0.5 1.0 1.5] * pi;
        xt = xc + r1 * cos(t) * cot - r2 * sin(t) * sit;
        yt = yc + r1 * cos(t) * sit + r2 * sin(t) * cot;
        
        % stores handle to graphic object
        h1 = plot(ax, xt([3 1]), yt([3 1]), varargin{:});
        h2 = plot(ax, xt([4 2]), yt([4 2]), varargin{:});
        hh = [h1 h2];
        
        if ~isempty(style)
            apply(style, hh);
        end
        
        if nargout > 0
            h = hh;
        end
    end
    
    function res = scale(obj, factor)
        % Return a scaled version of this ellipse.
        res = Ellipse2D([[obj.CenterX obj.CenterY obj.Radius1 obj.Radius2] * factor obj.Orientation]);
    end
    
    function res = translate(obj, shift)
        % Return a translated version of this ellipse.
        res = Ellipse2D([obj.CenterX+shift(1) obj.CenterY+shift(2) obj.Radius1 obj.Radius2 obj.Orientation]);
    end
    
    function res = rotate(obj, angle, varargin)
        % Return a rotated version of this ellipse.
        %
        % POLY2 = rotate(POLY, THETA)
        % POLY2 = rotate(POLY, THETA, CENTER)
        % THETA is given in radians, in counter-clockwise order.
        %
        % Example
        %   elli = Ellipse2D([40 20 25 15 0]);
        %   elli2 = elli.rotate(pi/4);
        %   figure; hold on; axis equal; axis([0 80 0 60]);
        %   draw(elli, 'b'); draw(elli.center(), 'bo');
        %   draw(elli2, 'm'); draw(elli2.center(), 'mo');
        
        transfo = AffineTransform2D.createRotation(angle, varargin{:});
        center2 = transfo.transformPoint([obj.CenterX obj.CenterY]);
        orientation2 = obj.Orientation + rad2deg(angle);
        res = Ellipse2D([center2  obj.Radius1 obj.Radius2  orientation2]);
    end
end % end methods


%% Serialization methods
methods
    function str = toStruct(obj)
        % Convert to a structure to facilitate serialization.
        str = struct('Type', 'Ellipse2D', ...
            'CenterX', obj.CenterX, ...
            'CenterY', obj.CenterY, ...
            'Radius1', obj.Radius1, ...
            'Radius2', obj.Radius2, ...
            'Orientation', obj.Orientation);
    end
end
methods (Static)
    function elli = fromStruct(str)
        % Create a new instance from a structure.
        elli = Ellipse2D([str.CenterX str.CenterY str.Radius1 str.Radius2 str.Orientation]);
    end
end


end % end classdef

