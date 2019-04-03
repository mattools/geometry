function tests = test_Scene(varargin)
%TEST_SCENE  Test case for the file Scene
%
%   Test case for the file Scene

%   Example
%   test_Scene
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-09-20,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2018 INRA - Cepia Software Platform.

tests = functiontests(localfunctions);

function test_Creation(testCase) %#ok<*DEFNU>
% Test call of function without argument
scene = Scene();
assertTrue(testCase, isa(scene, 'Scene'));

function test_read(testCase)

scene = Scene.read('kandinsky.scene');
assertTrue(isa(scene, 'Scene'));
assertEqual(3, length(scene.Shapes));
assertEqual(testCase, [0 90], scene.XAxis.Limits, 'AbsTol', .01);


function test_viewBox(testCase)
% test conversion viewbox to scene axis limits
box = [1 2 3 4 5 6];
scene = Scene();
setViewBox(scene, box);
box2 = viewBox(scene);
assertEqual(testCase, box, box2, 'AbsTol', .01);

