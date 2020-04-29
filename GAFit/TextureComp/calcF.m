function[F] = calcF(euler)
%fname is the filename of a set of euler angles, probably tex00XX.out
% if nargin>

F = calcT(euler(:,1), euler(:,2), euler(:,3));