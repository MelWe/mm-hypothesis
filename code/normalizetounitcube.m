function [ normalizeddata ] = normalizetounitcube( data,eps )
%data is an m x n matrix; the rows are the data points
%normalized data has coordinates in the the n-dimensional unit cube
%eps is used to test for zero -- e.g. intialize eps = 10^(-14);
[m,n] = size(data);
normalizeddata = zeros(m,n);
for i = 1:n
   c = data(:,i);
   if max(c) - min(c) < eps
       cnorm = .5*ones(m,1);
   else
       cnorm = (1/(max(c) - min(c)))*(c - min(c));
   end
   normalizeddata(:,i) = cnorm;
end
end

