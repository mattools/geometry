function tests = test_Polygon2D(varargin)
%TEST_POLYGON2D  Test case for the file Polygon2D
%
%   Test case for the file Polygon2D

%   Example
%   test_Polygon2D
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@nantes.inra.fr
% Created: 2018-09-01,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2018 INRA - Cepia Software Platform.

tests = functiontests(localfunctions);


function test_create(testCase) %#ok<*DEFNU>
% Test 'create' static factory method

vertices = [10 10; 20 10; 20 20; 10 20];

poly = Polygon2D.create(vertices);

assertTrue(testCase, isa(poly, 'Polygon2D'));
assertEqual(testCase, size(vertexCoordinates(poly), 1), 4);
