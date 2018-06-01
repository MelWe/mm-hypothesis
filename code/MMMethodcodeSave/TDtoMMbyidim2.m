function [TDtoMM] = TDtoMMbyidim2( data,idim,MM,MMuidim)
%data is m xn matrix of m points with n coordinates
%idim is an m x 1 matrix with the intrinsic dimensions of the m points
%MM is a multi-manifold structure with fields: idim, stats,pathindex,ptindex,
%cubes,M,totalpts,cubepts,components, EVsqdist
%MMuidim is a list of the unique values of idim in MM; the order corresponds
%to the order of idim values in MM.

%Each value of idim determines a subset of the data, datai. The distance
%from data i to the multi-manifold in MM of with the same idim dimension. 

%TDtoMM is five field structure, with one set of fields for each possible unique value of
%idim. It is indexed by the value of idim. 
%TDtoMM(i).idim is value of idim, so it equals i. 
%The other four fields of TDtoMM are sqdist,MMcount,count, and EVsqdist
%These are the sum of the squared distances
%of the idim points in data to the idim multi-manifold, the number of idim
%points in the cubes supporting the multi-manifold, the number of idim
%data points, and the expected value of the squared distances of the supported 
%points. 



u = unique(idim);
%find the subset of u of positive dimensions
%idim == -1 if the neighborhoods were  too sparse (< cutoff # of points)
indu = find(u.*(u > 0));
u = u(indu,:);
L = length(u);
for i = 1:L
    TDtoMM(u(i)).idim = u(i);
    ind1 = find(MMuidim == u(i));
    ind2 = find(idim == u(i));
    datai = data(ind2,:);
    if length(ind1) > 0
        j = u(i);
        stats = MM(j).stats;
        cubes = MM(j).cubes;
        M = MM(j).M;
        [TD, TDstats] = TotalDisttoMM( datai,stats,cubes,M);
        TDtoMM(u(i)).sqdist = TD(1).sqdist;
        TDtoMM(u(i)).MMcount = TD(1).MMcount;
        TDtoMM(u(i)).count = TD(1).count;
        TDtoMM(u(i)).EVsqdist = TD(1).sqdist/TD(1).MMcount;
    else
        TDtoMM(u(i)).sqdist = 0;
        TDtoMM(u(i)).MMcount = 0;
        TDtoMM(u(i)).count = size(datai,1);
        TDtoMM(u(i)).EVsqdist = 0;
    end
end
end

