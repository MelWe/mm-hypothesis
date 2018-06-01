function testvariancebaseddimLIDAR( testid )
%This file illustrates the variance based intrinsic dimension code and
%total variance functions on the Bridge_87K.txt file



if testid == 8
    fnamedata = '/Users/lindaness/Documents/MATLAB/LIDAR2/Bridge_87K.txt';
    data = dlmread(fnamedata);
    x = data(:,1);
    y = data(:,2);
    z = data(:,3);
    xdiam = max(x) - min(x);
    ydiam = max(y) - min(y);
    zdiam = max(z) - min(z);
    maxdiam = sqrt(xdiam^2 + ydiam^2 + zdiam^2)
   
    h1 = figure;
    scatter3(x,y,z);
    str = sprintf('scatterplot of data for testid %d',testid);
    title(str)
    dirname = '/Users/lindaness/Documents/MATLAB/MSVDLinda/testsLIDAR/';
    fnamescatter = [dirname,sprintf('test%d',testid)];
    saveas(h1,[fnamescatter,'.fig']);
    saveas(h1,[fnamescatter,'.png']);
    
    varthreshhold = .95; K = 3; epszero = 10^(-12); TVfigs = 1;
    radii = (maxdiam)*2.^(-1.*[4:7])  
    q = length(radii);
    fname = [dirname,sprintf('testvarbasedidim testid %d',testid)];
    idimfigs = 1;
    [ idimsummary,istatssummary,idim,firstscaleindex,scaleprob, scales,consecutive,dense,SG,indG,SSVs] = variancebaseddim( data,radii,K,varthreshhold,epszero,idimfigs,fname );
    idimsummary
    display('columheadings: idim,firstscaleindex, count');
    istatssummary
    save(fname,'varthreshhold', 'K', 'radii', 'epszero', 'TVfigs', 'fname');
    save(fname,'data','-append');
    save(fname, 'idimsummary','istatssummary','idim','firstscaleindex','scaleprob', 'scales','consecutive','dense','SG','indG','SSVs','-append');
    [SSVEnergy,idimSSVEnergy,EidimSSVEnergy,TV,ETV,idimTV,EidimTV ] = TotalVarianceFunctions(data,idim,SG,radii,q,scales,indG,TVfigs,fname);
    save(fname,'SSVEnergy','idimSSVEnergy','EidimSSVEnergy','TV','ETV','idimTV','EidimTV','-append');
elseif testid == 9
    fnamedata = '/Users/lindaness/Documents/MATLAB/wisdmreproducibility/datasets/Bridge_87K.txt';
    dirname = '/Users/lindaness/Documents/MATLAB/wisdmreproducibility/results/';
    fnamescatter = [dirname,sprintf('test%d-data',testid)];
    fname = [dirname,sprintf('test%d-results',testid)];
    
    data = dlmread(fnamedata);
    x = data(:,1);
    y = data(:,2);
    z = data(:,3);
    xdiam = max(x) - min(x);
    ydiam = max(y) - min(y);
    zdiam = max(z) - min(z);
    maxdiam = sqrt(xdiam^2 + ydiam^2 + zdiam^2)
   
    h1 = figure;
    scatter3(x,y,z);
    str = sprintf('scatterplot of data for testid %d',testid);
    title(str)
    saveas(h1,[fnamescatter,'.fig']);
    saveas(h1,[fnamescatter,'.png']);
    
    varthreshhold = .95; K = 3; epszero = 10^(-12); TVfigs = 1;
    radii = (maxdiam)*2.^(-1.*[4:7])  
    q = length(radii);
    
    idimfigs = 1;
    [ idimsummary,istatssummary,idim,firstscaleindex,scaleprob, scales,consecutive,dense,SG,indG,SSVs] = variancebaseddim( data,radii,K,varthreshhold,epszero,idimfigs,fname );
    idimsummary
    display('columheadings: idim,firstscaleindex, count');
    istatssummary
    save(fname,'varthreshhold', 'K', 'radii', 'epszero', 'TVfigs', 'fname');
    save(fname,'data','-append');
    save(fname, 'idimsummary','istatssummary','idim','firstscaleindex','scaleprob', 'scales','consecutive','dense','SG','indG','SSVs','-append');
    [SSVEnergy,idimSSVEnergy,EidimSSVEnergy,TV,ETV,idimTV,EidimTV ] = TotalVarianceFunctions(data,idim,SG,radii,q,scales,indG,TVfigs,fname);
    save(fname,'SSVEnergy','idimSSVEnergy','EidimSSVEnergy','TV','ETV','idimTV','EidimTV','-append');
elseif testid == 10
    fnamedata = '/Users/lindaness/Documents/MATLAB/wisdmreproducibility/datasets/Bridge_87K.txt';
    dirname = '/Users/lindaness/Documents/MATLAB/wisdmreproducibility/results/';
    fnamescatter = [dirname,sprintf('test%d-data',testid)];
    fname = [dirname,sprintf('test%d-vidim-',testid)];
    
    data = dlmread(fnamedata);
    x = data(:,1);
    y = data(:,2);
    z = data(:,3);
    xdiam = max(x) - min(x);
    ydiam = max(y) - min(y);
    zdiam = max(z) - min(z);
    maxdiam = sqrt(xdiam^2 + ydiam^2 + zdiam^2)
   
    h1 = figure;
    scatter3(x,y,z);
    str = sprintf('scatterplot of data for testid %d',testid);
    title(str)
    saveas(h1,[fnamescatter,'.fig']);
    saveas(h1,[fnamescatter,'.png']);
    
    varthreshhold = .95; K = 3; epszero = 10^(-12); TVfigs = 1;
    radii = (maxdiam)*2.^(-1.*[4:7])  
    q = length(radii)
    
    idimfigs = 1;
    [ idimsummary,istatssummary,idim,firstscaleindex,scaleprob, scales,consecutive,dense,SG,indG,SSVs] = variancebaseddim( data,radii,K,varthreshhold,epszero,idimfigs,fname );
    idimsummary
    display('columheadings: idim,firstscaleindex, count');
    istatssummary
    save(fname,'varthreshhold', 'K', 'radii', 'epszero', 'TVfigs', 'fname');
    save(fname,'data','-append');
    save(fname, 'idimsummary','istatssummary','idim','firstscaleindex','scaleprob', 'scales','consecutive','dense','SG','indG','SSVs','-append');
    [SSVEnergy,idimSSVEnergy,EidimSSVEnergy,TV,ETV,idimTV,EidimTV ] = TotalVarianceFunctions(data,idim,SG,radii,q,scales,indG,TVfigs,fname);
    save(fname,'SSVEnergy','idimSSVEnergy','EidimSSVEnergy','TV','ETV','idimTV','EidimTV','-append');
    
    
end
end


function [ idimsummary,istatssummary,idim,firstscaleindex,scaleprob, scales,consecutive,dense,SG,indG,SSVs] = variancebaseddim( data,radii,K,varthreshhold,epszero,idimfigs,fname )
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


[ SSVs,dense,ballcount ] = ComputeSSVs( data,radii,cutoff );
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

function [SSVEnergy,idimSSVEnergy,EidimSSVEnergy,TV,ETV,idimTV,EidimTV ] = TotalVarianceFunctions(data,idim, SG,radii,q,scales,indG,TVfigs,fname)
%data is an m xn matrix of m data points with n coordinate
%SG is an m x q matrix of total variances for data in each of the q balls
%with radii specified in radius (total variances of the centered data)
%SG entries for a ball are -1 if the intrinsic dimension was not computed for that ball 
%due to too few points etc
%scales is an m x q matrix of 0's and 1's specifying the scales at which idim was observed
%for each data point
%indG is an m x q matrix of 0's and 1's specifying the scales for which
%intrinsic dimension was computed 
%A variety of total variance functions are computed including the SSVEnergy
%function we computed at WiSDM
%if TVfig == 1, the total variance functions are plotted for the subset of
%points for which intrinsic dimension was computed.

[m,n] = size(data);
%SG is an m x q matrix containg the total variances for each ball 
%(or -1's for the ball where SG and hence the intrinsic dimension were
%not computed

%Compute different Versions of Total Variance Functions
%TV sum of total variances over the relevant balls
TV = sum(SG.*indG,2); 
ETV = TV./sum(indG,2); %m x 1 average of total variances of relevant balls

idimTV = sum(scales.*SG,2); %sum of scale variances for idim scales 
%total variance for each of these scale is concentrated near idim scales
EidimTV = idimTV./sum(scales,2); %Expected Value over the scales where 
%local idim = idim for the point

R = repmat(radii.^2,m,1);
RelSG = SG./R; %RelSG total variances for each ball are normalized by
%dividing by the square of the radius; functions analagous to above are 
%computed from the relative variances. 
SSVEnergy =sum(RelSG,2) ; %This is the SSVEnergy function computed at WiSDM
idimSSVEnergy = sum(scales.*RelSG,2);
EidimSSVEnergy = idimSSVEnergy./sum(scales,2);


if TVfigs == 1
     ind = find(idim > 0);
     SG = SG(ind,:);
     m = size(data(ind,:),1);
    
    %h = figure;
    %x = 1:m;
    %hold on
    %for i = 1:q
        %plot(x,SG(:,i));
    %end
    %titlestr = sprintf('Total Variance Curves for each scale');
    %title(titlestr)
    %hold off
    
    if n == 3
     X = data(ind,1); Y = data(ind,2); Z = data(ind,3);
     TV = TV(ind,:);
     idimTV = idimTV(ind,:);
     EidimTV = EidimTV(ind,:);
     SSVEnergy = SSVEnergy(ind,:);
     idimSSVEnergy = idimSSVEnergy(ind,:);
     EidimSSVEnergy = EidimSSVEnergy(ind,:);
        h = figure;
        colormap(jet(length(ind)));
        scatter3(X,Y,Z,100,TV);
        titlestr = sprintf('points color-coded by sum of Total Variances over scales where idim was computed');
        title(titlestr);
        if ~length(fname) == 0
            fname1 = [fname,'TV','.png'];
            saveas(h,fname1);
            fname1 = [fname,'TV','.fig'];
            saveas(h,fname1);
        end

        h = figure;
        colormap(jet(length(ind)));
        scatter3(X,Y,Z,100,idimTV);
        titlestr = sprintf('points color-coded by sum of Total Variance over idim scales for each point') ;
        title(titlestr);
        if ~length(fname) == 0
            fname1 = [fname,'idimTV','.png'];
            saveas(h,fname1);
            fname1 = [fname,'idimTV','.fig'];
            saveas(h,fname1);
        end

        h = figure;
        colormap(jet(length(ind)));
        scatter3(X,Y,Z,100,EidimTV);
        titlestr = sprintf('points color-coded by EV of Total Variance over idim scales for each point') ;
        title(titlestr);
        if ~length(fname) == 0
            fname1 = [fname,'EidimTV','.png'];
            saveas(h,fname1);
            fname1 = [fname,'EidimTV','.fig'];
            saveas(h,fname1);
        end

        h = figure;
        colormap(jet(length(ind)));
        scatter3(X,Y,Z,100,SSVEnergy);
        titlestr = sprintf('SSVEnergy - points color-coded by sum of normalized total variances over scales where idim was computed') ;
        title(titlestr);
        if ~length(fname) == 0
            fname1 = [fname,'SSVEnergy','.png'];
            saveas(h,fname1);
            fname1 = [fname,'SSVEnergy','.fig'];
            saveas(h,fname1);
        end

        h = figure;
        colormap(jet(length(ind)));
        scatter3(X,Y,Z,100,idimSSVEnergy);
        titlestr = sprintf('idmSSVEnergy - points color-coded by sum of normalized total variances over idim scales') ;
        title(titlestr);
        if ~length(fname) == 0
            fname1 = [fname,'idimSSVEnergy','.png'];
            saveas(h,fname1);
            fname1 = [fname,'idimSSVEnergy','.fig'];
            saveas(h,fname1);
        end

        h = figure;
        colormap(jet(length(ind)));
        scatter3(X,Y,Z,100,EidimSSVEnergy);
        titlestr = sprintf('EidmSSVEnergy - points color-coded by EV of normalized total variances over idim scales') ;
        title(titlestr); 
        if ~length(fname) == 0
            fname1 = [fname,'EidimSSVEnergy','.png'];
            saveas(h,fname1);
            fname1 = [fname,'EidimSSVEnergy','.fig'];
            saveas(h,fname1);
        end
    end
end

end

function [ SSVs,dense,ballcount ] = ComputeSSVs( data,radii,cutoff )
% data is an m x n matrix, whose rows are the data
%r is a 1 x q list of decreasing radii (e.g 2^0, 2^(-1), 2^(-2), ....);
%This function computes the squared singular values for each point p and each
%subset d(p,r) of data points at distance <= r from data point p 
%Results are returned in the m x q*maxdim matrix SSVs, maxdim = min(m,n);
%For the ith data point p = data(i,:) 
%SSVs for p and d(p,r(j)) are in row i,  columns [1:maxdim] + maxdim*(j-1)
%SSVs are not computed if |d(p,r)| < cutoff, instead -1's are stored in SSVs \
%dense is an m x q matrix of 0's and 1's
%dense(i,j)  = |d(p,r(j))| < cutoff

[m,n] = size(data);
q = length(radii);
%maxdim is the maximum number of non-zero singular values
maxdim = min(m,n);
%dense(i,j) = 1 if |d(p,r)| > cutoff else 0, p = data(i,:)
dense = ones(m,q);
SSVs = -1*ones(m,q*maxdim);
ballcount = zeros(m,q);

%for each point 
for i = 1:m
    if mod(i,1000) == 0
        i
    end
    p = data(i,:);
    ballmembership = zeros(m,q);
    for j = 1:q
        r = radii(j);
        [b,B] = ball(data,p,r);
        ballmembership(:,j) = B;
        %compute EV of b and translate by -EV so EV(c) = 0;
        c = center(b);
        cardc = size(c,1);
        if cardc > cutoff
            s = svd(c);
            sq = s.*s;
            inds2 = [1:maxdim] + maxdim*(j-1);
            SSVs(i,inds2) = sq; 
        else
            dense(i,j) = 0;
        end      
    end
    ballcount(i,:) = sum(ballmembership);
end
end

function [v,B] = ball( data,p,r) 
%% p=center; r=radius %%
%This functions finds the data points in the ball of radius r 
%centered at p. 
[m,n] = size(data);
x = sqrt(sum((data- ones(m,1)*p).^2,2));
B= x < r*ones(m,1);
v=data(B,:);
end

function [ centereddata ] = center( data )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
[m,n] = size(data);
EV = (1/m)*sum(data);
centereddata = data - ones(m,1)*EV;
end

function [ PrG,SG,CumG,indG ] = ProbDistGroups(M,q,maxdim,dense,epszero)
%M is an m x q*maxdim matrix
%dense is an m x q matrix
%epszero very small and positive and used to test for zero e.g 10^(-12)
%SG is m x q; PrG and SG are m x q*maxdim; indG is m x q
%For the jth group of maxdim columns, for each row i, if dense(i,j) == 1 and
%if the entries are >=0  with positive sum (> epszero) , indG(i,j) = 1 and
%the rowsum,probability distribution, and cumlative probability distribution
%of the row are computed and stored in the matrices SG, PrG and CumG.
%for the other rows in the group indG = 0 and the entries of 
%SG, PrG and CumG are -1's.
%if M does not have q*maxdim columns PrG,SG,CumG have -1 entries and inG = 0

[m,n] = size(M);
PrG = -1*ones(m,n);
CumG = -1*ones(m,n);
SG = -1*ones(m,q);
indG = zeros(m,q);
if n == q*maxdim
    for j = 1:q
        inddense = find(dense(:,j) == 1);
        indj = [1:maxdim] + maxdim*(j-1);
        [ Pr,S,Cum,ind ] = Probdist(M(inddense,indj),epszero);
        PrG(inddense,indj) = Pr;
        SG(inddense,j) = S;
        CumG(inddense,indj) = Cum;
        indG(inddense(ind),j) = 1;
    end
end
end


function [ Pr,S,Cum,ind ] = Probdist(M,epszero)
%M is a m x n matrix 
%S is an m x 1 matrix; Cum and Pr are m x n matrices. 
%S,Cum,Pr are initialized to -1;
%for the rows of M which are non-negative with positive sum (> epszero)
%S is the sum of the row
%Pr is the discrete probability distribution determined by the row
%Cum is the cumulative distribution determined by the row
%Example row = [1,2,3]; S = 6, Pr = [1/6,2/6,3/6], Cum = [1/6,1/2,1]
%For the other rows of M, the entries of S, Pr, and Cu are -1.
%ind is the index of the non-negative rows with positive sum


[m,n] = size(M);
S = -1*ones(m,1);
Pr = -1*ones(m,n);
Cum = -1*ones(m,n);
nonneg = M >= 0;
possum = sum(M,2);
N = sum(nonneg,2);
ind = find(N == n & possum > epszero);
S(ind,:)  = sum(M(ind,:),2);
k = length(ind);
if k > 0
    Tmp = zeros(k,1);
    for i = 1:n
        Pr(ind,i) = M(ind,i)./S(ind,1);
        Tmp = Tmp + Pr(ind,i);
        Cum(ind,i) = Tmp;
    end
end




end

function [idim,firstscaleindex, scales,idimG,mxG,scaleprob,consecutive] = idimprob( CumG,q,groupdim,t,indG)
%CumG is a matrix whose rows are q groups of cumulative distributions
%of size groupdim
%t is a threshhold 0 <= t <= 1
%The function returns column index of the first element in each row
%which is >= t; if no such index exists, the function returns -1.
[m,n] = size(CumG);
if ~(n == q*groupdim)
    idim = [];
    firstscaleindex = [];
    scales = [];
    idimG = [];
    mxG = [];
    scaleprob = [];
    consecutive = [];
    display('in idimprob -- CumG does not have the expected number of columns');
else
    idim = -1*ones(m,1);
    firstscaleindex = -1*ones(m,1);
    scales = zeros(m,q);
    idimG = -1*ones(m,q);
    mxG = -1*ones(m,q);
    scaleprob = -1*ones(m,1);
    consecutive = -1*ones(m,1);
    for j = 1:q
        ind = groupdim*(j-1) + [1:groupdim];
        C = CumG(:,ind);
        D = C >= t;
        %find first index in each row where C >= t
        [mx,ind] = max(D,[],2);
        for i = 1:m
            if mx(i) == 1
                idimG(i,j) = ind(i,1);
                mxG(i,j) = C(i,ind(i));
            elseif mx(i) == 0
                idimG(i,j) = -1;
                mxG(i,j) = -1;
            end
        end
        %col = mxG(:,j)
    end
end
%The intrinsic dimension for a point is the minimum of the intrinsic dimensions for
%each scale for which the intrinsic dimension is computed. 
%If it is not computed at any scale idim = -1.
for i = 1:m
    ind1 = find(indG(i,:) == 1);
    if length(ind1) > 0
        idim(i,1) = min(idimG(i,ind1));
        indscales = find(idimG(i,:) == idim(i,1));
        scales(i,indscales) = 1;
        scaleprob(i,1) = length(indscales)/length(ind1);
        if indscales == indscales(1) - 1 + [1:length(indscales)]
            consecutive(i,1) = 1;
        end
        firstscaleindex(i,1) = indscales(1);
    end
end
end

