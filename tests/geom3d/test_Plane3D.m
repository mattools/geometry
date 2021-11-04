function tests = test_Plane3D
% Test suite for the file Plane3D.
%
%   Test suite for the file Plane3D
%
%   Example
%   test_Plane3D
%
%   See also
%     Plane3D

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2021-11-04,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE - BIA-BIBS.

tests = functiontests(localfunctions);

function test_CreateFrom3Points(testCase) %#ok<*DEFNU>
% Test call of function without argument.

p1 = Point3D([1 0 0]);
p2 = Point3D([0 1 0]);
p3 = Point3D([0 0 1]);

plane = Plane3D(p1, p2, p3);

assertTrue(testCase, isa(plane, 'Plane3D'));



function test_medianPlane_twoPoints(testCase) %#ok<*DEFNU>
% Test call of function without argument.

p1 = Point3D([0 0 0]);
p2 = Point3D([2 4 6]);

plane = Plane3D.medianPlane(p1, p2);

assertTrue(testCase, isa(plane, 'Plane3D'));


function test_intersectLine(testCase) %#ok<*DEFNU>

p1 = Point3D([0 0 0]);
p2 = Point3D([10 0 0]);
med12 = Plane3D.medianPlane(p1, p2);
line12 = StraightLine3D(p1, p2);

pinter = intersectLine(med12, line12);

assertTrue(testCase, isa(pinter, 'Point3D'));

