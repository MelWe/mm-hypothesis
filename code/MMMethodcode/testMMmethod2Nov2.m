function [ MMbyidimstats, Runstats, MMsummary,MM,uidim,idimstats ] = testMMmethod2Nov2( testcaseid,dirname )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


if length(testcaseid) == length('LIDAR') & testcaseid == 'LIDAR'
     %run this in the directory '/Users/lindaness/Documents/MATLAB/wisdmreproducibility/code/MMMethodcode'
    fnamedata = '/Users/lindaness/Documents/MATLAB/wisdmreproducibility/datasets/Bridge_87K.txt';
    %dirname = '/Users/lindaness/Documents/MATLAB/wisdmreproducibility/results/';
    fnamescatter = [dirname,sprintf('test%s-data-',testcaseid)];
    fname1 = [dirname,sprintf('test%s-vidim-',testcaseid)];
    fname2 = [dirname,sprintf('test%s-MM-',testcaseid)];
    addpath('/Users/lindaness/Documents/MATLAB/wisdmreproducibility/code/vidim2code')
   

    
    %
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
    str = sprintf('scatterplot of data for testid %s',testcaseid);
    title(str)
    saveas(h1,[fnamescatter,'.fig']);
    saveas(h1,[fnamescatter,'.png']);
    
    varthreshhold = .95; K = 3; epszero = 10^(-12); TVfigs = 1;
    radii = (maxdiam)*2.^(-1.*[4:7])  
    idimfigs = 1;   
    [data,idim,idimstats] = computevidim2( data,radii,K,varthreshhold,epszero,idimfigs,TVfigs,fname1);
    
    %load(fname1,'data','idim');
    %only compute the Multi-Manifold for the points where idim > 0
    indpos = find(idim > 0);
    idim2 = idim(indpos,:);
    data2 = data(indpos,:);
    percentagenoidim = length(find(idim <= 0))/length(idim)
    
    testproportion = .4; runs = 20;seed = 5; eps = 10^(-14);SVDbound = 3000;
    maxscale = 3;K = 3; varthresh = .95; confidencebound = .95;figs = 0; 
    [MMbyidimstats, MMsummary,MM,uidim,Runstats ] = MMMethod(data2,idim2,testproportion,runs, seed,eps,SVDbound,maxscale,K,varthresh,confidencebound, figs,fname2 );
    
elseif length(testcaseid) == length('SphereLine') & testcaseid == 'SphereLine'
    %dirname = '/Users/lindaness/Documents/MATLAB/wisdmreproducibility/results/';
    fname1 = [dirname,sprintf('test%s-vidim-',testcaseid)];
    fname2 = [dirname,sprintf('test%s-MM-',testcaseid)];
    addpath('/Users/lindaness/Documents/MATLAB/wisdmreproducibility/code/vidim2code')
    addpath('/Users/lindaness/Documents/MATLAB/wisdmreproducibility/datasets')
    
    sampleparam = 800;rgparam = 5;
    [ sample,spheresample,linesample,samplesize ] = ExampleUpdated( sampleparam,rgparam);
    data = sample; radii = 2:-.1:.1;K = 3; varthreshhold = .95 ;epszero = 10^(-12);idimfigs = 1;TVfigs = 1;
    [data,idim,idimstats] = computevidim2( data,radii,K,varthreshhold,epszero,idimfigs,TVfigs,fname1);
    
    %load(fname1,'data','idim');
    %only compute the Multi-Manifold for the points where idim > 0
    indpos = find(idim > 0);
    idim2 = idim(indpos,:);
    data2 = data(indpos,:);
    percentagenoidim = length(find(idim <= 0))/length(idim)
    
    testproportion = .4; runs = 20;seed = 5; eps = 10^(-14);SVDbound = 3000;
    maxscale = 3;K = 3; varthresh = .95; confidencebound = .95;figs = 0; 
    
    [MMbyidimstats, MMsummary,MM,uidim,Runstats ] = MMMethod(data2,idim2,testproportion,runs, seed,eps,SVDbound,maxscale,K,varthresh,confidencebound, figs,fname2 );
    
end

end


