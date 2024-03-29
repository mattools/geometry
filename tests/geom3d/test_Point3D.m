function tests = test_Point3D(varargin)
%TEST_POINT3D  Test case for the file Point3D
%
%   Test case for the file Point3D

%   Example
%   test_Point3D
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-02-07,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2019 INRA - Cepia Software Platform.

tests = functiontests(localfunctions);

function test_Constructor_Single(testCase) %#ok<*DEFNU>
% Test call of function without argument

p = Point3D([5 4 3]);
assertEqual(testCase, 5, p.X);
assertEqual(testCase, 4, p.Y);
assertEqual(testCase, 3, p.Z);


function test_plusVector(testCase) %#ok<*DEFNU>
% Test call of function without argument

p1 = Point3D([10 10 10]);
v = Vector3D([1 2 3]);

res = p1 + v;

assertEqual(testCase, res.X, 11);
assertEqual(testCase, res.Y, 12);
assertEqual(testCase, res.Z, 13);


function test_distance(testCase) %#ok<*DEFNU>

p1 = Point3D([10 10 10]);
p2 = Point3D([12 13 16]);

d = distance(p1, p2);

assertEqual(testCase, d, 7.0, 'AbsTol', 0.01);


function test_distance_p2Array(testCase) %#ok<*DEFNU>

p1 = Point3D([10 10 10]);
p2 = Point3D([12 13 16;10 10 10;20 20 20]);

d = distance(p1, p2);

assertEqual(testCase, size(d), [3 1]);
assertEqual(testCase, d, [7.0 ; 0.0; sqrt(300)], 'AbsTol', 0.01);


function test_distance_p1p2Array(testCase) %#ok<*DEFNU>

p1 = Point3D([10 10 10;0 0 0]);
p2 = Point3D([12 13 16;10 10 10;20 20 20])';

d = distance(p1, p2);

assertEqual(testCase, size(d), [2 3]);


function test_Serialize(testCase)
% Test call of function without argument

p = Point3D([5 4 3]);
str = toStruct(p);
p2 = Point3D.fromStruct(str);
assertEqual(testCase, 5, p2.X);
assertEqual(testCase, 4, p2.Y);
assertEqual(testCase, 3, p2.Z);




