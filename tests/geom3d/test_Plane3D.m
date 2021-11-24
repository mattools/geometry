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


function test_clipPlaneXY(testCase) %#ok<*DEFNU>

plane = Plane3D(Point3D([5 5 5]), Vector3D([1 0 0]), Vector3D([0 1 0]));
bounds = Bounds3D([0 10 0 10 0 10]);

poly = clip(plane, bounds);

assertTrue(testCase, isa(poly, 'LinearRing3D'));
assertEqual(testCase, 4, vertexCount(poly));
assertTrue(testCase, ismember([ 0  0 5], poly.Coords, 'rows'));
assertTrue(testCase, ismember([10  0 5], poly.Coords, 'rows'));
assertTrue(testCase, ismember([ 0 10 5], poly.Coords, 'rows'));
assertTrue(testCase, ismember([10 10 5], poly.Coords, 'rows'));


function test_clipPlaneXZ(testCase) %#ok<*DEFNU>

plane = Plane3D(Point3D([5 5 5]), Vector3D([1 0 0]), Vector3D([0 0 1]));
bounds = Bounds3D([0 10 0 10 0 10]);

poly = clip(plane, bounds);

assertTrue(testCase, isa(poly, 'LinearRing3D'));
assertEqual(testCase, 4, vertexCount(poly));
assertTrue(testCase, ismember([ 0 5  0], poly.Coords, 'rows'));
assertTrue(testCase, ismember([10 5  0], poly.Coords, 'rows'));
assertTrue(testCase, ismember([ 0 5 10], poly.Coords, 'rows'));
assertTrue(testCase, ismember([10 5 10], poly.Coords, 'rows'));


function test_clipPlaneYZ(testCase) %#ok<*DEFNU>

plane = Plane3D(Point3D([5 5 5]), Vector3D([0 1 0]), Vector3D([0 0 1]));
bounds = Bounds3D([0 10 0 10 0 10]);

poly = clip(plane, bounds);

assertTrue(testCase, isa(poly, 'LinearRing3D'));
assertEqual(testCase, 4, vertexCount(poly));
assertTrue(testCase, ismember([5  0  0], poly.Coords, 'rows'));
assertTrue(testCase, ismember([5 10  0], poly.Coords, 'rows'));
assertTrue(testCase, ismember([5  0 10], poly.Coords, 'rows'));
assertTrue(testCase, ismember([5 10 10], poly.Coords, 'rows'));


function test_intersectPlane_XYvsYZ(testCase) %#ok<*DEFNU>

Oxy = Plane3D.XY;
Oyz = Plane3D.YZ;
lineOy = intersectPlane(Oxy, Oyz);

assertEqual(testCase, distance(origin(lineOy), Point3D([0 0 0])), 0, 'AbsTol', 0.01);


function test_projection_XY(testCase)

planeXY = Plane3D.XY;
point = Point3D([5 4 3]);

proj = planeXY.projection(point);

exp = Point3D([5 4 0]);
assertEqual(testCase, distance(proj, exp), 0.0, 'AbsTol', 0.01);

