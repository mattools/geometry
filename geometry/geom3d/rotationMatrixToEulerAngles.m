function varargout = rotationMatrixToEulerAngles(mat)
% Convert a 4-by-4 rotation matrix into three ZYX Euler angles.
%
%   [PHI, THETA, PSI] = rotationMatrixToEulerAngles(MAT)
%   Computes Euler angles PHI, THETA and PSI (in degrees) from a 3D 4-by-4
%   or 3-by-3 rotation matrix.
%
%   ANGLES = rotationMatrixToEulerAngles(MAT)
%   Concatenates the angles into  a single 1-by-3 row vector. This format
%   is used for representing 3D geometries such as ellipsoids.
%
%
%   Example
%     rotation3dToEulerAngles
%
%   References
%   Code from '1994 - Shoemake - Graphics Gems IV: Euler Angle Conversion:
%   http://webdocs.cs.ualberta.ca/~graphics/books/GraphicsGems/gemsiv/euler_angle/EulerAngles.c
%   (see also rotm2eul, that is part of MATLAB's Robotics System Toolbox)
%   Modified using explanations in:
%   http://www.gregslabaugh.net/publications/euler.pdf
%   https://www.geometrictools.com/Documentation/EulerAngles.pdf
%
%   See also
%     AffineTransform3D, Ellipsoid3D, rotation3dToEulerAngles

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% INRAE - BIA Research Unit - BIBS Platform (Nantes)
% Created: 2021-04-29,    using Matlab 9.8.0.1323502 (R2020a)
% Copyright 2021 INRAE.

% extract |cos(theta)|
cy = hypot(mat(1,1), mat(2,1));

% avoid dividing by 0
if cy > 16*eps
    % normal case: theta <> 0
    phi   = atan2( mat(2,1), mat(1,1));
    theta = atan2(-mat(3,1), cy);
    psi   = atan2( mat(3,2), mat(3,3));
else
    phi   = 0;
    theta = atan2(-mat(3,1), cy);
    psi   = atan2(-mat(2,3), mat(2,2));
end

% format output arguments
if nargout <= 1
    % one array
    varargout{1} = rad2deg([phi theta psi]);
else
    % three separate arrays
    varargout = cellfun(@rad2deg, {phi theta psi}, 'uni', 0);
end
