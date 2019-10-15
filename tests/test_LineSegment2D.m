function tests = test_LineSegment2D
%TEST_LINESEGMENT2D  Test case for the file LineSegment2D
%
%   Test case for the file LineSegment2D
%
%   Example
%   test_LineSegment2D
%
%   See also
%   LineSegment2D

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-10-15,    using Matlab 9.7.0.1190202 (R2019b)
% Copyright 2019 INRA - Cepia Software Platform.

tests = functiontests(localfunctions);

function test_Constructor_TwoPoints(testCase) %#ok<*DEFNU>
% Test call of function without argument

p1 = Point2D(20, 10);
p2 = Point2D(40, 20);

seg = LineSegment2D(p1, p2);

assertTrue(testCase, isa(seg, 'LineSegment2D'));

function test_Constructor_Empty(testCase) %#ok<*DEFNU>
% Test call of function without argument

seg = LineSegment2D();

assertTrue(testCase, isa(seg, 'LineSegment2D'));

function test_Constructor_Copy(testCase) %#ok<*DEFNU>
% Test call of function without argument

p1 = Point2D(20, 10);
p2 = Point2D(40, 20);

seg1 = LineSegment2D(p1, p2);
seg2 = LineSegment2D(seg1);

assertTrue(testCase, isa(seg2, 'LineSegment2D'));
