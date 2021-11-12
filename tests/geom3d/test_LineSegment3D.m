function tests = test_LineSegment3D
%TEST_LineSegment3D  Test case for the file LineSegment3D
%
%   Test case for the file LineSegment3D
%
%   Example
%   test_LineSegment3D
%
%   See also
%   LineSegment3D

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2019-10-15,    using Matlab 9.7.0.1190202 (R2019b)
% Copyright 2019 INRA - Cepia Software Platform.

tests = functiontests(localfunctions);

function test_Constructor_TwoPoints(testCase) %#ok<*DEFNU>
% Test call of function without argument

p1 = Point3D(30, 20, 10);
p2 = Point3D(60, 40, 20);

seg = LineSegment3D(p1, p2);

assertTrue(testCase, isa(seg, 'LineSegment3D'));

function test_Constructor_TwoPointsNumeric(testCase) %#ok<*DEFNU>
% Test call of function without argument

p1 = [30, 20, 10];
p2 = [60, 40, 20];

seg = LineSegment3D(p1, p2);

assertTrue(testCase, isa(seg, 'LineSegment3D'));

function test_Constructor_Empty(testCase) %#ok<*DEFNU>
% Test call of function without argument

seg = LineSegment3D();

assertTrue(testCase, isa(seg, 'LineSegment3D'));


function test_Constructor_Copy(testCase) %#ok<*DEFNU>
% Test call of function without argument

p1 = Point3D(30, 20, 10);
p2 = Point3D(60, 40, 20);

seg1 = LineSegment3D(p1, p2);
seg2 = LineSegment3D(seg1);

assertTrue(testCase, isa(seg2, 'LineSegment3D'));


function test_Constructor_Point3DArray(testCase) %#ok<*DEFNU>
% Create an array of line segment from tow arrays of points

p1 = Point3D([0 0 0;30 20 10]);
p2 = Point3D([10 10 10;60 40 20]);

segArray = LineSegment3D(p1, p2);

assertTrue(testCase, isa(segArray, 'LineSegment3D'));
assertEqual(testCase, size(segArray), [2 1]);


function test_planeIntersection_simple(testCase)

p1 = Point3D(30, 20, -10);
p2 = Point3D(30, 20, +10);
seg = LineSegment3D(p1, p2);
plane = Plane3D.XY;

point = seg.planeIntersection(plane);

assertEqual(testCase, 30, point.X);
assertEqual(testCase, 20, point.Y);
assertEqual(testCase,  0, point.Z);


function test_draw_simple(testCase) %#ok<*DEFNU>
% Test call of function without argument

p1 = Point3D(30, 20, 10);
p2 = Point3D(60, 40, 20);
seg = LineSegment3D(p1, p2);

f = figure;
h = draw(seg);
assertTrue(testCase, ishandle(h));
close(f);


function test_draw_drawOptions(testCase) %#ok<*DEFNU>
% Test call of function without argument

p1 = Point3D(30, 20, 10);
p2 = Point3D(60, 40, 20);
seg = LineSegment3D(p1, p2);

f = figure;
h = draw(seg, 'Color', 'b', 'LineWidth', 2);
assertTrue(testCase, ishandle(h));
assertEqual(testCase, h.Color, [0 0 1]);
assertEqual(testCase, h.LineWidth, 2);
close(f);

function test_draw_axis(testCase) %#ok<*DEFNU>
% Test call of function without argument

p1 = Point3D(30, 20, 10);
p2 = Point3D(60, 40, 20);
seg = LineSegment3D(p1, p2);

f = figure;
ax = gca;

h = draw(ax, seg, 'Color', 'b', 'LineWidth', 2);
assertTrue(testCase, ishandle(h));
assertEqual(testCase, h.Color, [0 0 1]);
assertEqual(testCase, h.LineWidth, 2);
close(f);

function test_draw_style(testCase) %#ok<*DEFNU>
% Test call of function without argument

p1 = Point3D(30, 20, 10);
p2 = Point3D(60, 40, 20);
seg = LineSegment3D(p1, p2);
style = Style('LineColor', 'b', 'LineWidth', 2);

f = figure;
h = draw(seg, style);

assertTrue(testCase, ishandle(h));
assertEqual(testCase, h.Color, [0 0 1]);
assertEqual(testCase, h.LineWidth, 2);

close(f);
