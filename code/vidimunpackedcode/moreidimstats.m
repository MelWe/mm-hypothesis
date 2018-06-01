function [ idimstats ] = moreidimstats( idimstats )
%three new fields are computed for the idimstats structure
%idimstats(1).sigratio is the ratio c/idim where c is the count of the
%number of points in the largest neighborhood to witness the intrinsic
%dimension; it is >= 1 if idim is defined (i.e. > 0) so it indicates
%the significance of the linear relationship revealed by idim
%idimstats(1).lexsorted is a 3 column matrix idim,first scale index,
%sigratio, lexicographically sorted in increasing order for idim,
%decreasing order for first scale index and decreasing order for sigratio
%idimstats(1).lexsortedindex = J where J transforms the original order
%to the lexicographically sorted order
[m,n] = size(idimstats.ballcount);
M = zeros(m,3);
M(:,1) = idimstats.idim;
M(:,2) = idimstats.firstscaleindex;
sigratio = zeros(m,1);
for i = 1:m
    f = idimstats.firstscaleindex(i,1);
    c = idimstats.ballcount(i,f);
    d = idimstats.idim(i,1);
    if d > 0
        sigratio(i,1) = c/d;
    else
        %d = -1
        sigratio(i,1) = -1;
    end 
end
M(:,3) = sigratio;
idimstats(1).sigratio = sigratio;
N = zeros(m,3);
N(:,1) = M(:,1);
N(:,2) = -1*M(:,2);
u = unique(N,'rows');
l = length(u);
Msorted = zeros(0,3);
J = zeros(0,1);
for i = 1:l
    ind = find(N(:,1) == u(i,1) & N(:,2) == u(i,2));
    A = M(ind,3);
    [B,I] = sort(A,'descend');
    C = M(ind,:);
    D = C(I,:);
    Msorted = [Msorted;D];
    J = [J; ind(I,1)];
end
idimstats(1).lexsorted = Msorted;
idimstats(1).lexsortedindex = J;
end

