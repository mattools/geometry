function tests = test_Polyline2D
% Test suite for the file Polyline2D.
%
%   Test suite for the file Polyline2D
%
%   Example
%   test_Polyline2D
%
%   See also
%     Polyline2D

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2021-01-26,    using Matlab 9.8.0.1323502 (R2020a)
% Copyright 2021 INRAE - BIA-BIBS.

tests = functiontests(localfunctions);

function test_CreateOpen(testCase) %#ok<*DEFNU>
% Test static method create.

coords = [0 0; 10 0; 10 10;0 10];

poly = Polyline2D.create(coords);

assertTrue(testCase, isa(poly, 'Polyline2D'));
assertEqual(testCase, vertexCount(poly), 4);
assertFalse(testCase, isClosed(poly));


function test_CreateClosed(testCase) %#ok<*DEFNU>
% Test static method create.

coords = [0 0; 10 0; 10 10;0 10];

poly = Polyline2D.create(coords, 'Closed', true);

assertTrue(testCase, isa(poly, 'Polyline2D'));
assertEqual(testCase, vertexCount(poly), 4);
assertTrue(testCase, isClosed(poly));
