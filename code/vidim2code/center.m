function [ centereddata ] = center( data )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
[m,n] = size(data);
EV = (1/m)*sum(data);
centereddata = data - ones(m,1)*EV;
end

