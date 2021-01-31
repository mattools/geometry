function tests = test_AffineTransform2D(varargin)
%TEST_AFFINETRANSFORM2D  Test case for the file AffineTransform2D
%
%   Test case for the file AffineTransform2D

%   Example
%   test_AffineTransform2D
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-04-20,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2019 INRA - Cepia Software Platform.

tests = functiontests(localfunctions);

function test_EmptyConstructor(testCase) %#ok<*DEFNU>
% Test call of function without argument
trans = AffineTransform2D();
assertTrue(testCase, isa(trans, 'AffineTransform2D'));

function test_createTranslation(testCase) %#ok<*DEFNU>

trans = AffineTransform2D.createTranslation([20 10]);
assertTrue(testCase, isa(trans, 'AffineTransform2D'));


function test_createRotation(testCase)
trans = AffineTransform2D.createRotation(pi/3);
assertTrue(testCase, isa(trans, 'AffineTransform2D'));


function test_createRotation90(testCase)

p = Point2D([20 0]);
trans = AffineTransform2D.createRotation90(1);

p2 = transform(p, trans);

assertTrue(testCase, isa(p2, 'Point2D'));
assertEqual(testCase, p2.X,  0, 'AbsTol', 0.01);
assertEqual(testCase, p2.Y, 20, 'AbsTol', 0.01);


function test_createRotation90_centerCoords(testCase)

p = Point2D([15 20]);
trans = AffineTransform2D.createRotation90(1, [10 20]);

p2 = transform(p, trans);

assertTrue(testCase, isa(p2, 'Point2D'));
assertEqual(testCase, p2.X, 10, 'AbsTol', 0.01);
assertEqual(testCase, p2.Y, 25, 'AbsTol', 0.01);


function test_createRotation90_centerPoint(testCase)

p = Point2D([15 20]);
trans = AffineTransform2D.createRotation90(1, Point2D(10, 20));

p2 = transform(p, trans);

assertTrue(testCase, isa(p2, 'Point2D'));
assertEqual(testCase, p2.X, 10, 'AbsTol', 0.01);
assertEqual(testCase, p2.Y, 25, 'AbsTol', 0.01);



function test_isIdentity_true(testCase) 
trans = AffineTransform2D();
assertTrue(testCase, isIdentity(trans));

function test_isIdentity_false(testCase) 
trans = AffineTransform2D.createRotation(pi/3);
assertFalse(testCase, isIdentity(trans));


function test_toAndFromStruct(testCase)
trans = AffineTransform2D.createTranslation([30 20]);
str = toStruct(trans);
trans2 = AffineTransform2D.fromStruct(str);
assertEqual(testCase, trans.Coeffs, trans2.Coeffs);
