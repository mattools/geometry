function tests = test_MultiPoint2D(varargin)
%TEST_POINT3D  Test case for the file MultiPoint2D
%
%   Test case for the file MultiPoint2D

%   Example
%   test_MultiPoint2D
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

v = [0 0;1 0;0 1;1 1];
pts = MultiPoint2D(v); 
assertTrue(testCase, isa(pts, 'MultiPoint2D'));


function test_Serialize(testCase)
% Test call of function without argument

v = [0 0;1 0;0 1;1 1];
pts = MultiPoint2D(v); 

str = toStruct(pts);

pts2 = MultiPoint2D.fromStruct(str);
assertEqual(testCase, pts.Coords, pts2.Coords);


function test_draw_simple(testCase) %#ok<*DEFNU>
% Test call of function without argument

v = [0 0;1 0;0 1;1 1];
pts = MultiPoint2D(v); 

f = figure;
h = draw(pts);
assertTrue(testCase, ishandle(h));
close(f);


function test_draw_drawOptions(testCase) %#ok<*DEFNU>
% Test call of function without argument

v = [0 0;1 0;0 1;1 1];
pts = MultiPoint2D(v); 

f = figure;
h = draw(pts, 'Color', 'b', 'Marker', '+');
assertTrue(testCase, ishandle(h));
assertEqual(testCase, h.Color, [0 0 1]);
assertEqual(testCase, h.Marker, '+');
close(f);

function test_draw_style(testCase) %#ok<*DEFNU>
% Test call of function without argument

v = [0 0;1 0;0 1;1 1];
pts = MultiPoint2D(v); 
style = Style('MarkerColor', 'b', 'MarkerStyle', '+');

f = figure;
h = draw(pts, style);

assertTrue(testCase, ishandle(h));
assertEqual(testCase, h.Color, [0 0 1]);
assertEqual(testCase, h.Marker, '+');

close(f);
