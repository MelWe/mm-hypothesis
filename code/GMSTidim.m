function [kintdim, intdim, datadim] = GMSTidim(data, subset, neigh, param, varargin)
%========================================================================== 
% Syntax
%   [kintdim, intdim, datadim] = GMSTidim(data, subset, neigh, param)
%   [kintdim, intdim, datadim] = GMSTintdim(data, poi, neigh, param, aggr);
%
% GMSTintdimest -- GMST INTrinsic DIMension ESTimation 
%       This function estimates the intrinsic dimension at each point using 
%       a variant of the geodesic minimal spanning tree (GMST) to estimate 
%       the intrinsic dimension and entropy of the manifold on which the 
%       data lie.  For each of these 
%       sub-manifolds, the intrinsic dimension d is computed as follows. 
%       For each point x,
%       1. Form a neighborhood using either knn or balls
%       2. Compute the intrinsic dimension using the GMST algorithm
%       3. The intrinsic dimension at each point is obtained by taking the 
%       average ('mean') or the most common value ('mode') over all 
%       neighborhoods or scales. 
%==========================================================================
% Inputs: 
%       data    - m x n matrix: rows are the data points.
%       poi     - p x n matrix: rows points of interest, subset of data.
%       neigh   - 'kNN' or 'e-ball'
%       param   - k: vector that specifies number of nearest neighbors.
%                 For e.g., for the LIDAR data set, I used the set
%                 k = [70, 80, 100, 200, 250, 300, 500, 600, 750, 800].
%                 The SSV is normalized by the mean distance of the
%                 neighbors. 
%                 scales: vector that specifies scale range. Balls are 
%                 of radius 2^s, where s in scales
%       aggr    - aggregation method: 'mean' or 'mode' 
%                 Default is 'mean'.
% 
% Outputs:   
%       intdim  - intrinsic dimension at in each point based on aggregating 
%                 kintdim according to 'method'
%       kintdim - intrinsic dimension at in each point for each k
%       datadim - m x n+1 matrix [data,intdim]
%==========================================================================
% Reference     : GMST
% Author        : KYD
% Created       : Jul 27, 2017 by KYD
% Last revised	: Aug 11, 2017 by KYD
%==========================================================================

tstart = tic;

% Check if method is
c = length(varargin);
if c < 1
    aggr = 'mean';
else
    aggr = varargin{1};
end

% Check validity of neighborhood construction method
c1 = strcmp('knn', neigh);
c2 = strcmp('e-ball', neigh);
c = or(c1,c2);
if c == 0
    fprintf('\n Warning: Unknown method for constructing neighborhood at each data point! \n')
    fprintf('\n')
    return
end

% Check validity of method
c3 = strcmp('mode', aggr);
c4 = strcmp('mean', aggr);
c = or(c3,c4);
if c == 0
    fprintf('\n Warning: Unknown method for aggregating intrinsic dimension at each data point! \n')
    fprintf('\n')
    return
end

% Compute intrinsic dimension per using GMST
if c1 == 1 % Use kNN to construct manifolds
    % Declare sizes and initialize matrices and parameters
    [n, q] = size(subset);
    if q == 2
        X = subset(:,1); Y = subset(:,2);
    elseif q == 3
       X = subset(:,1); Y = subset(:,2); Z = subset(:,3);
    end
    k = param;
    kintdim = ones(n,length(k));
    no_neighs = length(k);
    % Compute intrinsic dimension for each k
    for s = 1:no_neighs
        [idx, ~] = knnsearch(data, subset,'dist','euclidean','k',k(s));
        for i = 1:n
            neighbors = data(idx(i,:),:);
            no_dims = intrinsic_dim(neighbors, 'GMST');
            if no_dims > 3
                kintdim(i,s) = 3;
            else
                kintdim(i,s) = no_dims;
            end
        end
        % Plot intrinsic dimension

        figure; scatter3(X,Y,Z,3,kintdim(:,s)); 
        title(strcat('Intrinsic Dimension for k =',{' '},num2str(k(s))));
    end
elseif c2 == 1 % Use e-balls to construct manifolds
    [n, q] = size(subset);
    scales = param;
    kintdim = ones(n,length(scales));
    no_neighs = length(scales);
    % Compute intrinsic dimension per scale
    for s = 1:no_neighs
        r = 2^(-scales(s));
        [idx, ~] = rangesearch(data, subset, r);  
        for i = 1:n
            neighbors = data(idx{i},:);     
            no_dims = intrinsic_dim(neighbors, 'GMST');
            if no_dims > 3
                kintdim(i,s) = 3;
            else
                kintdim(i,s) = no_dims;
            end
        end
        % Plot intrinsic dimension
        X = subset(:,1); Y = subset(:,2); Z = subset(:,3);
        figure; scatter3(X,Y,Z,3,kintdim(:,s));
        title(strcat('Intrinsic Dimension for scale =',{' '},num2str(2^(-scales(s)))));
    end
end 

% Intrinsic dimension for each point
if c3 == 1
    intdim = mode(kintdim,2);
elseif c4 == 1
    intdim = mean(kintdim,2);
end
    
% Format data for output
datadim  = [subset,intdim];

% Plot intrinsic dimension
if q == 2
    figure; scatter3(X,Y,intdim,3,intdim); 
    title('Intrinsic Dimension');
else
    figure; scatter3(X,Y,Z,3,intdim); 
    title('Intrinsic Dimension in 3D');
end

% Compute and print the time elapsed
tend = toc(tstart);
no_min = floor(tend/60);
no_sec = ceil(rem(tend,60));
disp(strcat('Time elapsed is',{' '},num2str(no_min), 'min',{' '}, ...
               num2str(no_sec),'sec'));

% Signal end of code    
load gong.mat; sound(y)
disp('Done!');

end