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


function test_constructor_twoPoints(testCase) %#ok<*DEFNU>

p1 = Point3D([10 20 30]);
p2 = Point3D([12 23 34]);

v = Vector3D(p1, p2);

exp = Vector3D([2 3 4]);
assertEqual(testCase, norm(v - exp), 0.0, 'AbsTol', 0.01);


function test_constructor_twoPointArrays(testCase) %#ok<*DEFNU>

p1 = repmat(Point3D([10 20 30]), [3 1]);
p2 = repmat(Point3D([12 23 34]), [1 4]);

v = Vector3D(p1, p2);

assertEqual(testCase, size(v), [3 4]);


function test_constructor_twoPointArraysDifferentDimensions(testCase) %#ok<*DEFNU>

p1 = repmat(Point3D([10 20 30]), [3 1]);
p2 = repmat(Point3D([12 23 34]), [1 4 2]);

v = Vector3D(p1, p2);

assertEqual(testCase, size(v), [3 4 2]);


function test_crossProduct(testCase) %#ok<*DEFNU>
% Test call of function without argument.

v1 = Vector3D([2 0 0]);
v2 = Vector3D([0 3 0]);

v = crossProduct(v1, v2);

assertEqual(testCase, v.X, 0);
assertEqual(testCase, v.Y, 0);
assertEqual(testCase, v.Z, 6);


