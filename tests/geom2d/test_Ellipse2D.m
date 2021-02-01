function tests = test_Ellipse2D
% Test suite for the file Ellipse2D.
%
%   Test suite for the file Ellipse2D
%
%   Example
%   test_Ellipse2D
%
%   See also
%     Ellipse2D

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2021-02-01,    using Matlab 9.8.0.1323502 (R2020a)
% Copyright 2021 INRAE - BIA-BIBS.

tests = functiontests(localfunctions);


function test_Simple(testCase) %#ok<*DEFNU>
% Test call of function without argument.

elli = Ellipse2D([50 40 30 20 10]);

assertTrue(testCase, isa(elli, 'Ellipse2D'));


function test_asPolyline(testCase) %#ok<*DEFNU>
% Test call of function without argument.

elli = Ellipse2D([50 40 30 20 10]);

poly = asPolyline(elli, 20);

assertTrue(testCase, isa(poly, 'Polyline2D'));
assertEqual(testCase, vertexCount(poly), 20);


function test_area(testCase) %#ok<*DEFNU>

elli = Ellipse2D([5 4 3 2 10]);

a = area(elli);

assertEqual(testCase, a, pi * 2 * 3, 'AbsTol', 0.01);


function test_transform_rotation(testCase)

elli = Ellipse2D([5 4 3 2 0]);
transfo = AffineTransform2D.createRotation(pi/6);

elli2 = transform(elli, transfo);
centerT = transform(center(elli), transfo);

assertTrue(testCase, isa(elli2, 'Ellipse2D'));
assertEqual(testCase, elli2.CenterX, centerT.X, 'AbsTol', 0.01);
assertEqual(testCase, elli2.CenterY, centerT.Y, 'AbsTol', 0.01);
assertEqual(testCase, elli2.Radius1, elli.Radius1, 'AbsTol', 0.01);
assertEqual(testCase, elli2.Radius2, elli.Radius2, 'AbsTol', 0.01);
assertEqual(testCase, elli2.Orientation, 30, 'AbsTol', 0.01);


function test_transform_scaling(testCase)

elli = Ellipse2D([5 4 3 2 0]);
transfo = AffineTransform2D.createScaling([3 2]);

elli2 = transform(elli, transfo);
centerT = transform(center(elli), transfo);

assertTrue(testCase, isa(elli2, 'Ellipse2D'));
assertEqual(testCase, elli2.CenterX, centerT.X, 'AbsTol', 0.01);
assertEqual(testCase, elli2.CenterY, centerT.Y, 'AbsTol', 0.01);
assertEqual(testCase, elli2.Radius1, 9, 'AbsTol', 0.01);
assertEqual(testCase, elli2.Radius2, 4, 'AbsTol', 0.01);
assertEqual(testCase, elli2.Orientation, 0, 'AbsTol', 0.01);


function test_fromCartesianCoefficients(testCase) %#ok<*DEFNU>
% Test call of function without argument.

elli = Ellipse2D([50 40 30 20 10]);

coeffs = cartesianCoefficients(elli);
elli2 = Ellipse2D.fromCartesianCoeffs(coeffs);

assertTrue(testCase, isa(elli2, 'Ellipse2D'));
assertEqual(testCase, elli2.CenterX, elli.CenterX, 'AbsTol', 0.01);
assertEqual(testCase, elli2.CenterY, elli.CenterY, 'AbsTol', 0.01);
assertEqual(testCase, elli2.Radius1, elli.Radius1, 'AbsTol', 0.01);
assertEqual(testCase, elli2.Radius2, elli.Radius2, 'AbsTol', 0.01);
assertEqual(testCase, elli2.Orientation, elli.Orientation, 'AbsTol', 0.01);

