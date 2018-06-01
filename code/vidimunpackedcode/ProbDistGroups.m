function [ PrG,SG,CumG,indG ] = ProbDistGroups(M,q,maxdim,dense,epszero)
%M is an m x q*maxdim matrix
%dense is an m x q matrix
%epszero very small and positive and used to test for zero e.g 10^(-12)
%SG is m x q; PrG and SG are m x q*maxdim; indG is m x q
%For the jth group of maxdim columns, for each row i, if dense(i,j) == 1 and
%if the entries are >=0  with positive sum (> epszero) , indG(i,j) = 1 and
%the rowsum,probability distribution, and cumlative probability distribution
%of the row are computed and stored in the matrices SG, PrG and CumG.
%for the other rows in the group indG = 0 and the entries of 
%SG, PrG and CumG are -1's.
%if M does not have q*maxdim columns PrG,SG,CumG have -1 entries and inG = 0

[m,n] = size(M);
PrG = -1*ones(m,n);
CumG = -1*ones(m,n);
SG = -1*ones(m,q);
indG = zeros(m,q);
if n == q*maxdim
    for j = 1:q
        inddense = find(dense(:,j) == 1);
        indj = [1:maxdim] + maxdim*(j-1);
        [ Pr,S,Cum,ind ] = Probdist(M(inddense,indj),epszero);
        PrG(inddense,indj) = Pr;
        SG(inddense,j) = S;
        CumG(inddense,indj) = Cum;
        indG(inddense(ind),j) = 1;
    end
end
end




