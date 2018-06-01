function [ SSVs,dense,ballcount ] = ComputeSSVs2( data,radii,cutoff )
% data is an m x n matrix, whose rows are the data
%r is a 1 x q list of decreasing radii (e.g 2^0, 2^(-1), 2^(-2), ....);
%This function computes the squared singular values for each point p and each
%subset d(p,r) of data points at distance <= r from data point p 
%Results are returned in the m x q*maxdim matrix SSVs, maxdim = min(m,n);
%For the ith data point p = data(i,:) 
%SSVs for p and d(p,r(j)) are in row i,  columns [1:maxdim] + maxdim*(j-1)
%SSVs are not computed if |d(p,r)| < cutoff, instead -1's are stored in SSVs \
%dense is an m x q matrix of 0's and 1's
%dense(i,j)  = |d(p,r(j))| < cutoff

[m,n] = size(data);
q = length(radii);
%maxdim is the maximum number of non-zero singular values
maxdim = min(m,n);
%dense(i,j) = 1 if |d(p,r)| > cutoff else 0, p = data(i,:)
dense = ones(m,q);
SSVs = -1*ones(m,q*maxdim);
ballcount = zeros(m,q);

%for each point 
for i = 1:m
    if mod(i,1000) == 0
        i
    end
    p = data(i,:);
    ballmembership = zeros(m,q);
    for j = 1:q
        r = radii(j);
        [b,B] = ball(data,p,r);
        ballmembership(:,j) = B;
        %compute EV of b and translate by -EV so EV(c) = 0;
        c = center(b);
        cardc = size(c,1);
        if cardc > cutoff
            s = svd(c);
            sq = s.*s;
            inds2 = [1:maxdim] + maxdim*(j-1);
            sqextended = zeros(size(inds2));
            sqextended(1:length(sq)) = sq;
            SSVs(i,inds2) = sqextended; 
        else
            dense(i,j) = 0;
        end      
    end
    ballcount(i,:) = sum(ballmembership);
end
end

