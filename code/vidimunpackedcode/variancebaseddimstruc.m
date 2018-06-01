function [ idimstats,SSVs] = variancebaseddimstruc( data,radii,K,varthreshhold,epszero,idimfigs,fname )
%This function has the functionality as variancebaseddim. 
%SSVs2 corrects an error in SSVs
%The idim related output variables are organized into a structure idimstats. 

%idim plots are shown if idimfigs == 1
%variance function pots are shown if TVfigs == 1
%radii is a decreasing sequence of radiuses determining a sequence of
%balls around each point. The index of each radius is referred to as the
%scale.
[m,n] = size(data); % m data points, n coordinates
maxdim = min(m,n);
q = length(radii);
%This program looks for points with intrinsic dimensions <= K (with respect
%to the family of point neighborhoods specified by radii
cutoff = K*log(K);
%This function computes the intrinsic dimension idim for each point by
%a concentration of variance analysis using varthreshhold (e.g. .95)
%SSVs is an m x maxdim*q matrix of q groups of squared singular values
%of matrices of points in each ball, the matrices are centered so EV = 0
%The squared singular values are then the variances of each centered coordinate.
%PeG is an m x maxdim*q matrix of q groups of discrete probability distriutions
%for each of the m points (relative variances); SG is an m x q matrix of
%total variances for each ball; CumG is an m x maxdim*q matrix of q groups
%of vectors (c1, ..., cn) defining the cumulative distribution for
%(p1,..,pn) so cj = p1 + ... pj 
% idim(p) = i, i = minimum over the scales of smallestindex j cj >= varthresshold
%e.g if varthreshhold = .95 ,idim(p) = i if there is a ball (scale) for which
%the sum of the variances of the first i centered coordinates >= .95
%scales is an m xq matrix of 0's and 1's. each row (pt) has 1's at the scales 
%i, such that ci >= varthresshold and i = idim(pt) 
%i.e. idim is computed by concentration of variance analysis  and the
%scales is the set of indices of balls where the variance accumulates to
%the fastest to varthreshhold

%Two summary matrices are computed: idimsummary and istatssummary
%rows of idimsummary (idim value, number of points with idim value)
%i.e. the discrete distribution of idim for the data set
%rows of istats summary (idim value, first scale index, number of points 
%with the specified idim value and first scale index.
%i.e. the discrete distribution of idim,firstscale for the data set

%If idimfigs == 1  scatter plots of the data points color-coded by idim
%are computed

%Different Total Variance functions are computed for each point, including
%the SSVEnergy function. 
%Summaries of idim statistics are computed; idim and the Total Varianc
%functions are plotted


[ SSVs,dense,ballcount ] = ComputeSSVs2( data,radii,cutoff );
M = SSVs; 
[ PrG,SG,CumG,indG ] = ProbDistGroups(M,q,maxdim,dense,epszero);
%PrG1599 = PrG(1599,:)
%SG199 = SG(1599,:)
%CumG1599 = CumG(1599,:)
groupdim = n;
t = varthreshhold;

[idim,firstscaleindex, scales,idimG,mxG,scaleprob,consecutive] = idimprob( CumG,q,groupdim,t,indG);
u = unique(idim);
idimsummary = zeros(length(u),2 );
for i = 1:length(u)
    ind = find(idim == u(i));
    lgth = length(ind);
    idimsummary(i,:) = [u(i),lgth];
end

[istats,irhs,ilhs] = unique([idim,firstscaleindex],'rows');

%istatssummary = zeros(length(istats),4);
istatssummary = zeros(length(istats),3);
for i = 1:size(istats,1)
    %x = abs(istats(i,1:3));
    x = istats(i,1:2);
    y = length(find(ilhs == i));
    istatssummary(i,:) = [x,y];
end

idimstats(1).idim = idim;
idimstats(1).firstscaleindex = firstscaleindex;
idimstats(1).indG = indG;
idimstats(1).scales = scales;
idimstats(1).idimG = idimG;
idimstats(1).scaleprob = scaleprob;
idimstats(1).scales = scales;
idimstats(1).consecutive = consecutive;
idimstats(1).dense = dense;
idimstats(1).SG = SG;
idimstats(1).idimsummary = idimsummary;
idimstats(1).istatssummary = istatssummary;
idimstats(1).ballcount = ballcount;

if idimfigs == 1
    ind = find(idim > 0);
    m2 = length(ind);
    
    h = figure;
    a = [1:m2];
    b = idim(ind,:);
    scatter(a,b,100);
    titlestr = sprintf('ordered points with y coordinate = idim');
    if ~length(fname) == 0
        fname1 = [fname,'orderedidim','.png'];
        saveas(h,fname1);
        fname1 = [fname,'orderedidim','.fig'];
        saveas(h,fname1);
    end
    
    if n == 3
        X = data(ind,1); Y = data(ind,2); Z = data(ind,3);
        idim2 = idim(ind,:);
        firstscaleindex2 = firstscaleindex(ind,:);

        h = figure;
        colormap(jet(m2));
        scatter3(X,Y,Z,10,idim2);
        titlestr = sprintf('points color-coded by idim');
        title(titlestr);
        if ~length(fname) == 0
            fname1 = [fname,'idim','.png'];
            saveas(h,fname1);
            fname1 = [fname,'idim','.fig'];
            saveas(h,fname1);
        end

        h = figure;
        colormap(jet(m));
        scatter3(X,Y,Z,10,ilhs(ind,:));
        %scatter3(X,Y,Z,100,firstscaleindex2);
        titlestr = sprintf('points color-coded by lexicographic ordering of idim,firstscaleindex');
        title(titlestr);
        if ~length(fname) == 0
            fname1 = [fname,'lexidim','.png'];
            saveas(h,fname1);
            fname1 = [fname,'lexidim','.fig'];
            saveas(h,fname1);
        end
    end
end
end
