function [ idim,EV,W,totalvar] = idimlocal( data,varthresh )
%This function computes the variance based intrinsic dimension for 
%a data set by finding the smallest i, such sum(first i svd
%variances)/totalvariance < varthresh
%It also computes the expected value EV of the data, and the right singular
%vectors (columns of  W ) which span the best fitting subspace if dim idim
%andd the total variance. It is not necessary to computethe sum of the squared distances
%of the data points to the best fittign subspace of the computed intrinsic
%dimension because d2 < (1-varthresh)*totalvar
%thus an upper bound for the expected value of the squared distance d2 
%to the idimensional subspace is 
%given by  1/m * d2 < 1/m*(1-varthresh)*totalvar. 
[m,n] = size(data);
if m > 1
    EV = (1/m)*sum(data);
elseif m == 1
    EV = data;
end
centereddata = data - repmat(EV,m,1);
[U,s,V] = svd(centereddata);
if size(s,1) > 1
    s2 = diag(s).*diag(s);
    totalvar = sum(s2);
    l = length(diag(s));
elseif size(s,1) == 1
    diags = s(1,1);
    s2 = diags*diags;
    totalvar = s2;
    l = 1;
end
initialvar = 0;
i = 0;
p = 0;
while p < varthresh & i < l
    i = i+1;
    initialvar = initialvar + s2(i);
    p = initialvar/totalvar;
end
idim = i;
W = V(:,1:idim);
end

