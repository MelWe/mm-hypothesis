function [ MMsummary,MM,u ] = MultiManifold(data,idim,SVDbound,maxscale,K, varthresh,eps)
%This function computes a linear multi-manifold structure for the data set
%data is a p xn  matrix with n columns, whose rows are the data points
%idim is a p x 1 whose values are the intrinsic dimension of the data points 
%The data is normalized to be in the unit cube. 
%The exploits the dyadic binary tree determined by recursively 
%splitting the unit cube in half along each dimension and then repeating
%for each resulting sub-cube. The nodes in the binary tree correspond to
%cubes. 
%For each value of idim, the subset dimdata of data with that intrinsic
%dimension is computed. The subset is analyzed top-down using the binary
%tree structure. For each node cube, local intrinsic dimension of the idim
%data in the cube is computed if there are not too many points (> SVDbound)
%in the cube and if there are not too few < cutoff = K*logK. The local
%intrinsic dimension is the dimension of the "best fitting
%subspace" (using the variance threshhold varthresh). The top down
%recursion is stopped when cubes are found whose local intrinsic dimension
%<= idim. For each node set the best fitting affine space is computed 
%(the EV of the points in the cube and a basis for the subpaces centere at
%EV) and the EV of the squared distances from the points to the linear spaces
%is computed. The affine spaces for the (disjoint) leaf noides form the multi-manifold for 
%the points with the specified value of idim. The defining information is
%recorded in the structure MM(i), for the ith value of idim.
% The linear multi-manifold for the whole data is the union of each of
% these multi-manifold structures. It is easy to understand approximately
%how the multi-manifolds for different values of idim fit together by
%using the binary tree structure. 
%Thus the Multi-Manifold Structure consists of tangent affine spaces to the
%to each of the idim data subsets
%EVsqdist is the expected value of the squared distance of points to
%their associated affine spaces.
%u is a list of the unique idim values for the data , in increaisng order.
ndata = normalizetounitcube( data,eps );
u = unique(idim);
%MMsummary = -2*ones(length(u),5);
dimdata3 = [];
for i = 1:length(u)
    d = u(i);
    ind1 = find(idim == d);
    dimdata = ndata(ind1,:);
    [m,n] = size(dimdata);
    
    dim = 1; stats = [0,0,m,-2,-2,-2,-2,-2,2*n,-2];pathindex = -1*ones(1,maxscale*n); 
    ptindex = ones(1,m); currentcube = [zeros(1,n);ones(1,n)]'; 
    cubes = currentcube; 
    M = 2*n*ones(0,n*(n+1));
    [ stats,pathindex, ptindex,cubes,M ] = IntDimPCtopdown5( dimdata,SVDbound,d, dim, maxscale,K, stats,pathindex, ptindex,cubes, currentcube,varthresh,M );
    
    %stats is a 10 column matrix:dim,scale,count,product coefficient,
    %idim,rowparent,rowleftchild,rowrightchild,minidim,TV
    ind2 = find(stats(:,7) == -2 & stats(:,8) == -2);
    EVsqdist = 0;
    %create the stucture defining the Multi-Manifold;; MM(i) defines
    %the multi-manifold for points with idim = i; 
    MM(i).idim = d;
    MM(i).stats = stats(ind2,:); 
    MM(i).pathindex = pathindex(ind2,:);
    MM(i).ptindex = ptindex(ind2,:);
    cubeindex = sort([2*ind2' - 1, 2*ind2']);
    MM(i).cubes = cubes(:,cubeindex);
    MM(i).M = M(ind2,:);
    MM(i).totalpts = length(ind1);%the number of data points with idim = i
    MM(i).cubepts = sum(MM(i).stats(:,3));  %number of points in MM(i)
    MM(i).components = size(MM(i).stats,1); %number of "linear subspaces" in MM(i) 
    TV = MM(i).stats(:,10);
    EVsqdist = EVsqdist + TV;
    MM(i).EVsqdist = (1-varthresh)*sum(TV)/MM(i).cubepts; %expected value of sum of squares of
    %distance of MM(i) points to the affine space in the same cube.
    %MMsummary(i,:) = [MM(i).idim, MM(i).totalpts,MM(i).cubepts,MM(i).components, MM(i).EVsqdist];
    MMsummary(u(i)).idim = u(i);
    MMsummary(u(i)).totalpts = MM(i).totalpts;
    MMsummary(u(i)).MMpoints = MM(i).cubepts;
    MMsummary(u(i)).components = MM(i).components;
    MMsummary(u(i)).EVsqdist  = MM(i).EVsqdist;
end




    
end

