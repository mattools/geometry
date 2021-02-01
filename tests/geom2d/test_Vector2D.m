function tests = test_Vector2D
% Test suite for the file Vector2D.
%
%   Test suite for the file Vector2D
%
%   Example
%   test_Vector2D
%
%   See also
%     Vector2D

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2021-02-01,    using Matlab 9.8.0.1323502 (R2020a)
% Copyright 2021 INRAE - BIA-BIBS.

tests = functiontests(localfunctions);

function test_Simple(testCase) %#ok<*DEFNU>
% Test call of function without argument.

vect = Vector2D(3, 2);

assertTrue(testCase, isa(vect, 'Vector2D'));


function test_plus(testCase) %#ok<*DEFNU>
% Test call of function without argument.

v1 = Vector2D(5, 3);
v2 = Vector2D(2, 1);

res = v1 + v2;

assertTrue(testCase, isa(res, 'Vector2D'));
assertEqual(testCase, res.X, 7, 'AbsTol', 0.01);
assertEqual(testCase, res.Y, 4, 'AbsTol', 0.01);


function test_minus(testCase) %#ok<*DEFNU>
% Test call of function without argument.

v1 = Vector2D(5, 3);
v2 = Vector2D(2, 1);

res = v1 - v2;

assertTrue(testCase, isa(res, 'Vector2D'));
assertEqual(testCase, res.X, 3, 'AbsTol', 0.01);
assertEqual(testCase, res.Y, 2, 'AbsTol', 0.01);


function test_mtimes(testCase) %#ok<*DEFNU>
% Test call of function without argument.

v1 = Vector2D(4, 3);

res = v1 * 2;

assertTrue(testCase, isa(res, 'Vector2D'));
assertEqual(testCase, res.X, 8, 'AbsTol', 0.01);
assertEqual(testCase, res.Y, 6, 'AbsTol', 0.01);


function test_mrdivide(testCase) %#ok<*DEFNU>
% Test call of function without argument.

v1 = Vector2D(4, 2);

res = v1 / 2;

assertTrue(testCase, isa(res, 'Vector2D'));
assertEqual(testCase, res.X, 2, 'AbsTol', 0.01);
assertEqual(testCase, res.Y, 1, 'AbsTol', 0.01);


function test_norm(testCase) %#ok<*DEFNU>
% Test call of function without argument.

v1 = Vector2D(4, 3);

n = norm(v1);

assertEqual(testCase, n, 5, 'AbsTol', 0.01);


function test_normalize(testCase) %#ok<*DEFNU>
% Test call of function without argument.

v1 = Vector2D(4, 3);

res = normalize(v1);

assertTrue(testCase, isa(res, 'Vector2D'));
assertEqual(testCase, res.X, 4/5, 'AbsTol', 0.01);
assertEqual(testCase, res.Y, 3/5, 'AbsTol', 0.01);

