function [ stats,pathindex, ptindex,cubes,M ] = IntDimPCtopdown5( data,SVDbound,idimbound, dim, maxscale,K, stats,pathindex, ptindex,cubes, currentcube,varthresh,M )
%The difference between this version 5 and version 3 is an additional
%input parameter SVDbound and an additional input parameter idimbound.
%The use of the parameter SVDbound is explained in NOTE 2 below.
%The parameter idimbound is used to stop the recursion if the
%computed intrinsic dimension idim for the cube is <= idimbound.

%data is an m x n matrix
%the entries must in the interval [0,1], so each row is a set of coordinates
%for a point in a unit cube of dimension n
%dim is the current dimension for dyadic splitting
%This function computes PCs for the counting measure top down
%The binary tree strucure divides each of the dimension intervals in half
%and then repeats that process on the subintervals
%outputs: stats, pathindex, ptindex, and M all have the same number of rows
%there is one row for each node in the tree with non-zero count.


%stats is a 10 column matrix:dim,scale,count,product coefficient,
%idim,rowparent,rowleftchild,rowrightchild,minidim,TV

%pathindex has maxscale*n columns; the positive entries in 
%each row have value 0 or 1; the remaining entries are -1; the positive
%entries are the labels of the path from the root to the node in the tree 
%the first row is all -1's and represents the root node
%ptindex has m columns; the entries are 0's and 1's; each row
%of ptindex is the index for the subset of the data points whose 
%coordinates all lie in the "quadrant" determined by the row's tree node

%currentcube is an nx2 matrix whose ith row is endpoints of the interval
%for the ith dimension of the current cube.
%cubes is an n x 2*(rows in stats) matrix; it is a concatenation 
%of the current cube matrices.


%M is a matrix with n*(n+1) columns; the first n columns is the Expected
%Value of the points associated with the node (given by pt index); the next
%idim groups of n columns are the first idim right singular vectors; the
%other groups of n entries (if any) have value -2; M specifies the
%"multi-scale multi-manifold" approximating the data set, by specifying the linearspace
%associated with each node. Think of this as a multi-manifold or 
%each particular scale s(i.e depth n*s in the tree). Each linearspace
%should be constrained to be in the cube associated with the node.

%NOTE: for large n, this matrix will be too large, so perhaps only a
%restricted number R of right singular vectors should be reported, or perhaps 
%effectively giving complete information only for idim < a specified bound

%NOTE2: In this version, idim and M and TV are not computed for a cube unless the number
%of points in the cube <= SVDbound This is to avoid SVD computations that
%are too large. Let rowcube denote the row for such a cube. 
%stats(rowcube,5) = -2 stats(rowcube,10)  -2 and M(row,cube) = zeros(1,n*(n+1);

%initialization:
%dim = 1; stats = [0,0,m,-2,-2,-2,-2,-2,2*n,-2];pathindex = -1*ones(1,maxscale*n); 
%ptindex = ones(1,m); currentcube = [zeros(1,n);ones(1,n)]'; 
%cubes = currentcube; varthresh = .95; idimancestormin = n;
%M = -2**ones(0,n*(n+1));
%NOTE: the depth of the tree is maxscale*n not just maxscale; 

if min(data) < 0 | max(data) > 1
    display('line30: PCtopdown2.m: data matrix does not have entries in [0,1]')
end

[m,n] = size(data);
row = size(stats,1);
row2 = size(ptindex,1);
if ~(row == row2)
    display('error')
    row
    row2
end
currentptindex = find(ptindex(row,:) == 1);
currentstatsrow = stats(row,:);
count = currentstatsrow(1,3);
currentpathindexrow = pathindex(row,:);

indfornextdepth = find(currentpathindexrow == -1);
if length(indfornextdepth) == 0
    atmaxdepth = 1;
else
    atmaxdepth = 0;
    nextdepth = min(find(currentpathindexrow == -1));
    nextscale = floor((nextdepth - 1)/n) + 1;
end
inputrow = row;

%compute rows for left and right


if length(currentptindex) > K*log(K)
    l = currentcube(dim,1);
    r = currentcube(dim,2);
    mid = (l+r)/2;
    leftinterval = [l,mid];
    rightinterval = [mid,r];
    d = data(currentptindex,dim);
    indleft = find(d < mid & d >= l);
    leftcount = length(indleft);
    if r < 1
        indright = find(d >= mid & d < r);
    elseif r == 1
        indright = find(d >= mid & d <= r);
    end
    rightcount = length(indright);
    pc = (leftcount - rightcount)/length(currentptindex);
    stats(row,4) = pc;
    if length(currentptindex) <= SVDbound
        [idim,EV,V,totalvar]  = idimlocal( data(currentptindex,:),varthresh );
        stats(row,5) = idim;
        stats(row,10) = totalvar;
        idimmin = min(idim,stats(row,9));
        stats(row,9) = idimmin;
        U = -2*ones(0,n*(n+1));
        if idim >= 1
            W = reshape(V,[1,idim*n]);
            U(1,1:(idim*n + n)) = [EV,W];
            M = [M;U];
        else
            M = [M;U];
        end
    else
        U = -2*ones(1,n*(n+1));
        M = [M;U];
        idim = 2*n;
        idimmin = min(idim,stats(row,9));
        stats(row,9) = idimmin;
    end
        
    if atmaxdepth == 0
        leftstatsrow = [dim,nextscale,leftcount,-2,-2,row,-2,-2,idimmin,-2];
        rightstatsrow = [dim,nextscale,rightcount,-2,-2,row,-2,-2,idimmin,-2];
        leftptindex = currentptindex(indleft);
        leftptindexrow = zeros(1,m);
        leftptindexrow(1,leftptindex) = 1;
        rightptindex = currentptindex(indright);
        rightptindexrow = zeros(1,m);
        rightptindexrow(1,rightptindex) = 1;
        leftpathindexrow = currentpathindexrow;
        leftpathindexrow(1,nextdepth) = 0;
        rightpathindexrow = currentpathindexrow;
        rightpathindexrow(1,nextdepth) = 1;

        if dim == n
            nextdim = 1;
        else 
            nextdim = dim + 1;
        end
        
        
        leftstats = [stats;leftstatsrow];
        leftpathindex = [pathindex;leftpathindexrow];
        leftptindex = [ptindex; leftptindexrow];
        leftcube = currentcube;
        leftcube(dim,:) = leftinterval;
        %cubes = [cubes,leftcube];
        if (leftcount > K*log(K)) & idimmin > idimbound
             leftstats(row,7) = row + 1;
            cubes = [cubes,leftcube];
            [ stats,pathindex, ptindex,cubes,M ] = IntDimPCtopdown5( data,SVDbound,idimbound,nextdim, maxscale,K, leftstats,leftpathindex, leftptindex,cubes, leftcube,varthresh,M );
        end
        
        
        rightstats = [stats;rightstatsrow];
        rightpathindex = [pathindex;rightpathindexrow];
        rightptindex = [ptindex; rightptindexrow];
        rightcube = currentcube;
        rightcube(dim,:) = rightinterval;
        %cubes = [cubes,rightcube];
        if (rightcount > K*log(K)) & idimmin > idimbound
            rightstats(row,8) = size(stats,1) + 1;
            cubes = [cubes,rightcube];
            [ stats,pathindex, ptindex,cubes,M ] = IntDimPCtopdown5( data,SVDbound,idimbound,nextdim, maxscale,K, rightstats,rightpathindex, rightptindex,cubes, rightcube,varthresh,M );
        end
    end
    
  
end

