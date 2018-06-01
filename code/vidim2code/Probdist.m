function [ Pr,S,Cum,ind ] = Probdist(M,epszero)
%M is a m x n matrix 
%S is an m x 1 matrix; Cum and Pr are m x n matrices. 
%S,Cum,Pr are initialized to -1;
%for the rows of M which are non-negative with positive sum (> epszero)
%S is the sum of the row
%Pr is the discrete probability distribution determined by the row
%Cum is the cumulative distribution determined by the row
%Example row = [1,2,3]; S = 6, Pr = [1/6,2/6,3/6], Cum = [1/6,1/2,1]
%For the other rows of M, the entries of S, Pr, and Cu are -1.
%ind is the index of the non-negative rows with positive sum


[m,n] = size(M);
S = -1*ones(m,1);
Pr = -1*ones(m,n);
Cum = -1*ones(m,n);
nonneg = M >= 0;
possum = sum(M,2);
N = sum(nonneg,2);
ind = find(N == n & possum > epszero);
S(ind,:)  = sum(M(ind,:),2);
k = length(ind);
if k > 0
    Tmp = zeros(k,1);
    for i = 1:n
        Pr(ind,i) = M(ind,i)./S(ind,1);
        Tmp = Tmp + Pr(ind,i);
        Cum(ind,i) = Tmp;
    end
end


