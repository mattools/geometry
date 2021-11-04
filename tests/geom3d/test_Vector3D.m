function tests = test_Vector3D
% Test suite for the file Vector3D.
%
%   Test suite for the file Vector3D
%
%   Example
%   test_Vector3D
%
%   See also
%     Vector3D

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2021-11-04,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE - BIA-BIBS.

tests = functiontests(localfunctions);

function test_crossProduct(testCase) %#ok<*DEFNU>
% Test call of function without argument.

v1 = Vector3D([2 0 0]);
v2 = Vector3D([0 3 0]);

v = crossProduct(v1, v2);

assertEqual(testCase, v.X, 0);
assertEqual(testCase, v.Y, 0);
assertEqual(testCase, v.Z, 6);


