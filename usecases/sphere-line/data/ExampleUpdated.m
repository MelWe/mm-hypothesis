function [ sample,spheresample,linesample,samplesize ] = ExampleUpdated( sampleparam,rgparam)
%This function generates samples from a figure consisting of a sphere of
%radius 1/2 and a through the origin and then computes the SSV Energy
%function for the figure.
%  
%x = rsin(theta)cos(phi)
%y = rsin(theta)sin(phi)
%z = rcos(theta)
%theta in [0,pi]
%phi in [0,2*pi]
rng(rgparam);
theta  = pi*rand(floor(pi*sampleparam),1);
phi = 2*pi*rand(floor(pi*sampleparam),1);
spheresample = [.5*sin(theta).*cos(phi),.5*sin(theta).*sin(phi),.5*cos(theta)];
sz = floor(.5*sampleparam);
r = rand(sz,1);
r = 2*(r - .5*ones(sz,1));
ind = find(r >= .5 | r <= -.5);
linesample = zeros(length(ind),3);
linesample(:,1) = r(ind);
sample = [spheresample;linesample];
%add two points on the intersection of the line and the spere
%these points should have intrinisic dimension 3
sample = [sample;.5,0,0;-.5,0,0];
samplesize = length(sample);
x = sample(:,1); y = sample(:,2); z = sample(:,3);
scatter3(x,y,z)

end

