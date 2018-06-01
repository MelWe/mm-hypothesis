function [ idimGsig,scalessig,totalsigscales ] = significantidim(idimG,ballcount )
%idimG and ballcount are m x q matrices
%each row of the matrices contain q local statistics for a data point
%the statistics are computed for the data in  q neigborhoods of the point 
%idimG rows contain local intrinsic dimensions or -1 if idim not computed
%because there were too few points in the neighborhood.
%ballcount rows contain the number of data points in the q local
%neighborhoods of p
%the function computes idimGsig, an updated version of local intrinsic dimension
%for each neighborhood which resets non-significant values to -1 and two
%idimstat functions: the scales  indicator function which indicates
%which neighborhood indice have significant idim values and
%totalsigscales -- the number of significant scales for each data point 
idimGsig = idimG;
scalessig = zeros(size(idimG));
m = size(idimG,1);
for i = 1:m
    d = idimG(i,:);
    ct = ballcount(i,:);
    sig = d > 0 & ct >= d.*log(d);
    %if sig(1,j) = 0, the intrinsic dimension isn't small enough relative to the
    %number of data points in the neighborhood, so that local value of idim will
    %be reset to -1, so it isn't used to compute the value of idim for the
    %set of neighborhoods of the point
    ind2 = find(sig == 0);
    idimGsig(i,ind2) = -1;
    ind3 = find(sig == 1);
    scalessig(i,ind3) = 1;
end
totalsigscales = sum(scalessig,2);
    
    
    



 


end

