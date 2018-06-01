function [T] = idimtests( testid)
%This function tests the unpacked vidimcode on the common_prepreprocssed
%Word2Vec data set. The table T contains a list of the words, their 
%intrinsic dimension, the index of the largest neighborhood at which the
%idim was observed and the significance ratio. 
%No straightforward local correlation of the word semantics with the idim
%was observed. 


if testid == 10
    fnamedata = '/Users/lindaness/Documents/MATLAB/wisdmreproducibility/datasets/common_processed.xlsx';
    fnamedata2 = '/Users/lindaness/Documents/MATLAB/wisdmreproducibility/datasets/common_processed.mat';
    [num,txt,raw] = xlsread(fnamedata);
    [r,c] = size(num);
    numdata = num(:,1:(c-1));
    [p,q] = size(txt);
    words = txt(2:p,q);
    eps = 10^(-14);
    [ normalizeddata ] = normalizetounitcube( numdata,eps );
    size(normalizeddata)
    save(fnamedata2,'numdata','words','normalizeddata');
    display('numdata, words, and normalizeddata saved');
   
    addpath('/Users/lindaness/Documents/MATLAB/wisdmreproducibility/code/MMMethodcode');
    addpath('/Users/lindaness/Documents/MATLAB/wisdmreproducibility/code/vidimunpackedcode');
    dirname = '/Users/lindaness/Documents/MATLAB/wisdmreproducibility/results/';
    localfname = sprintf('W2V-test-%d',testid);
    fname = [dirname,localfname]
    tablename = [fname,'-wordsidim.txt'];
    
    load(fnamedata2,'words','normalizeddata');
    [m,n] = size(normalizeddata)
    eps = 10^(-14);
    diam = diameter(normalizeddata)
    diam = 3.9628
    radii = diam*[1,.9,.8,.7,.6,.5]
    K = 3;varthreshhold = .95; idimfigs = 0;
    [ idimstats,SSVs] = variancebaseddimstruc( normalizeddata,radii,K,varthreshhold,eps,idimfigs,fname );
    idimstats = moreidimstats( idimstats );
    J = idimstats(1).lexsortedindex;
    wordssorted = words(J,:);
    lexsorted = idimstats(1).lexsorted;
    save(fname,'idimstats','wordssorted','words');
    idim = idimstats(1).lexsorted(:,1);
    firstscaleindex = idimstats(1).lexsorted(:,2);
    sigratio = idimstats(1).lexsorted(:,3);
    T = table(wordssorted, idim,firstscaleindex,sigratio,'VariableNames',{'words', 'idim','index', 'sigratio'});
    writetable(T,tablename);
    
end
    
    
    


end

