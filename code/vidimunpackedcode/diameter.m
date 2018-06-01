function [ diam ] = diameter( D)
%This computes the diameter of a data set D 
%rows of D are the data points
[m,n] = size(D);
diam = 0;
for i = 1:m
    for j = 1:(i-1)
        n = norm(D(i,:)- D(j,:));
        diam = max(n,diam);
    end
end
end

