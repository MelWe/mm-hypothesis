function [ d2 ] = distancetosubspace2( data,p,B)
%Given an m x n matrix of data (rows are data, columns are measurements)
%compute the sum of the squared distances of the data set to the linear
%space containing the point p and spanned by the columns of B;
%It is assumed that the columns of B are orthnormal;

 

[m,n] = size(data);
data2 = data - repmat(p,m,1);
S = svd(data2);
%sum of square lengths of rows of data2 = sum of squares of singular values
x = sum(S.*S);
%project rows of data2 onto columns of B
projdata2 = data2*B;
S2 = svd(projdata2);
%sum of square lengths of rows of projdata2 = sum of squares of sing values
y = sum(S2.*S2);
d2 = x - y;







end

