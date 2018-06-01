function [ idimsigsummary] = idimsig( idimstats )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[m,n] = size(idimstats.ballcount);
idim = idimstats(1).idim;
count = idimstats(1).ballcount;
scale = idimstats(1).firstscaleindex;
sig = zeros(m,1);
M = zeros(m,3);
M(:,2) = idim; 
for i = 1:m
    s = scale(i,1);
    d = idim(i,1);
    scalecount = count(i,s);
    %sig = scalecount >= d*log(d)
    sig = scalecount >= (1.2)*d;
    M(i,1) = sig;
    M(i,3) = scalecount;
end
u = unique(M,'rows');
l = length(u);
S = zeros(l,4);
for i = 1:l
    r = u(i,:);
    ind = find(M(:,1) == r(1) & M(:,2) == r(2) & M(:,3) == r(3));
    S(i,1:3) = u(i,:);
    S(i,4) = length(ind);
end
%display('sig idim scalecount number')
idimsigsummary = S;

    

