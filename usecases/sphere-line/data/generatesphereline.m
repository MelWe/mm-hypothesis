function [sample, spheresample, linesample, samplesize] = generatesphereline(sampleparam,rgparam,varargin)
%==========================================================================
% Syntax 
%   [sample, spheresample, linesample, samplesize] = generatesphereline(sampleparam,rgparam)
%   [sample, spheresample, linesample, samplesize] = generatesphereline(sampleparam,rgparam,option)
%
% generatesphereline - GENERATE SPHERE and LINE  
%   This function generates a random sample from the figure consisting of a
%   sphere and two line segments emanating from the sphere. 
%   The sphere has radius 1/2 and goes through the origin: 
%   x = rsin(theta)cos(phi)
%   y = rsin(theta)sin(phi)
%   z = rcos(theta)
%   theta in [0,pi]
%   phi in [0,2*pi]
% 
% Example of inputs:
% sampleparam = 1000; rgparam = 5; 
%==========================================================================
% Inputs: 
%       sampleparam - 
%       rgparam     - 
%       option      - 1 to display data set, 0 otherwise.
% 
% Outputs:   
%       sample  -   Coordinates of figure
%       spheresample - Coordinates of sphere
%       linesample   - Coordinates of line
%       samplesize   - 
%==========================================================================
% Reference : N/A
% Author   	: KYD
% Created	: Jul 18, 2017 at L. Ness
% Revised	: Aug 11, 2017 by KYD (mainly clean up)
%==========================================================================

% Initialize parameters
rng(rgparam);
c = length(varargin);
if c == 1
	option = varargin{1};
else
    option = 0;
end 

% Compute sphere coordinates
theta  = pi*rand(floor(pi*sampleparam),1);
phi = 2*pi*rand(floor(pi*sampleparam),1);
spheresample = [.5*sin(theta).*cos(phi),.5*sin(theta).*sin(phi),.5*cos(theta)];

% Compute line coordinates
sz = floor(.5*sampleparam);
r = rand(sz,1);
r = 2*(r - .5*ones(sz,1));
ind = find(r >= .5 | r <= -.5);
linesample = zeros(length(ind),3);
linesample(:,1) = r(ind);
sample = [spheresample;linesample];

% Add two points on the intersection of the line and the spere
% These points should have intrinisic dimension 3

sample = [sample;.5,0,0;-.5,0,0];
samplesize = length(sample);
if option == 1
    x = sample(:,1); y = sample(:,2); z = sample(:,3);
    figure; scatter3(x,y,z)
end

end
