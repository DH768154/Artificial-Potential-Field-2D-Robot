function T = Rot2d(theta,varargin)
T = [cos(theta) -sin(theta);sin(theta) cos(theta)];
if nargin==2
    T = [T,varargin{1};0,0,1];
end
end