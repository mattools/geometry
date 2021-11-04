function tests = test_StraightLine3D
% Test suite for the file StraightLine3D.
%
%   Test suite for the file StraightLine3D
%
%   Example
%   test_StraightLine3D
%
%   See also
%     StraightLine3D

% ------
% Author: David Legland
% e-mail: david.legland@inrae.fr
% Created: 2021-11-04,    using Matlab 9.10.0.1684407 (R2021a) Update 3
% Copyright 2021 INRAE - BIA-BIBS.

tests = functiontests(localfunctions);

function test_position(testCase) %#ok<*DEFNU>
% Test call of function without argument.

p1 = Point3D([10 10 10]);
p2 = Point3D([50 10 10]);
line = StraightLine3D(p1, p2);
% choose a point located at 3/4 of the line direction
pt = Point3D([40 20 20]);

pos = position(line, pt);

assertEqual(testCase, pos, 0.75);


function test_clip_Bounds_horizLine(testCase)

p1 = Point3D([10 10 10]);
p2 = Point3D([50 10 10]);
line = StraightLine3D(p1, p2);
bounds = Bounds3D([0 100 0 100 0 100]);

seg = clip(line, bounds);

assertTrue(testCase, isa(seg, 'LineSegment3D'));
