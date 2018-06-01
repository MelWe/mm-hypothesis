function [TD, TDstats] = TotalDisttoMM( data,stats,cubes,M)
%data is an mxn matrix - m data points in dimension n space
[m,n] = size(data);

%stats, cubes and M are outputs of IntDimPCtopdown5, for a different
%n-dimensional data set; cube, stats and M have the same number of rows.
%Each row of stats has one row corresponds to a dyadic cube.
%data set.It is assumed that the outputs determine a multi-manifold, i.e.
%none of the cubes have an ancestral relationship. 

c = size(stats,1);


%stats is a 10 column matrix:dim,scale,count,product coefficient,
%idim,rowparent,rowleftchild,rowrightchild,minidim,TV

%cubes is an n x 2*(rows in stats) matrix specifying the boundaries of
%each dyadic cube

% M is a matrix with n*(n+1) columns; EAch row defines an affine space, 
% the first n columns is the Expected Value of the points in the other data set 
% associated with the cube; the rests of the columns are ab orthogonal
% basis for a subspace centred at EV


%This function first computes  2 field structure TDstats with fields sqdist
%and count. TDstats(i).sqdist is the sum of the
%squared distances from the data points in the ith cube to the affine space
%defined by the ith row of M. TDstats(i).count is the number of data points
%in the ith cube. 

%TD is a 4-field structure with one set of 4 fields.
%TD(1).sqdist is the sum of TDstats(i).sqdist
%TD(1).MMcount is the sum of TDstats(i).count. TD(1).count = m. The
%difference of the last 2 is the number of data points not in a cube
%supporting the multi-manifold.TD(1).EVsqdist =
%TD(1).totalsqdist/TD(1).cubecount. 


TD(1).sqdist = 0;
TD(1).MMcount = 0;
TD(1).count = m;
TD(1).EVsqdist = NaN;
for i = 1:c
    currentcube = cubes(:,2*i-1:2*i);
    %find the data points in the current cube
    Pts = data;
    L = size(Pts,1);
    j= 1;
    while L > 0 & j <= n
        mn = currentcube(j,1);
        mx = currentcube(j,2);
        if mx < 1
            ind = find(Pts(:,j) >= mn & Pts(:,j) < mx);
        elseif mx == 1
            ind = find(Pts(:,j) >= mn & Pts(:,j) <= mx);
        end
        Pts = Pts(ind,:);
        L = size(Pts,1);
        j = j + 1;
    end
    TDstats(i).count = L;
    
    p = M(i,1:n);
    ind = find(~(M(i,n+1:n*(n+1)) == -2));
    l = length(ind);
    if rem(l,n) == 0
        idim = l/n;
        B = reshape(M(i,n + [1:n*idim]),n,idim);
        d2 = distancetosubspace2( data,p,B);
        TDstats(i).sqdist = d2;
    else
        TDstats(i).sqdist = NaN;
        %error 
    end
    TD(1).sqdist = TD(1).sqdist + TDstats(i).sqdist;
    TD(1).MMcount = TD(1).MMcount + TDstats(i).count;
end
if TD(1).MMcount > 0
    TD(1).EVsqdist = TD(1).sqdist/TD(1).MMcount;
end;
end

