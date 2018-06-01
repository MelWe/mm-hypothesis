function [MMbyidimstats,MMsummary,MM,uidim,Runstats ] = MMMethod2(data,idim,testproportion,runs, seed,eps,SVDbound,maxscale,K,varthresh,confidencebound, figs,fname )
%This function takes as input an m xn matrix data, whose rows are the
%data points, and an m x1 matrix idim whose values are the intrinsic 
%dimension of the data points (e.g. computed by variancebasedintrinsicdim
%uidim is a vector of the unique intrinsic dimensions idim.
%It normalizes the data points to be in the unit cube and then computes
%For each value of idim, linear multi-manifold MM(i) for the subset of 
%data points with idim = i.
%MM is a structure whose fields define the linear multi-manifold. 
%Each "component" of the multi-manifold is supported in a unique dyadic 
%cube. MMsummary is a structure summarizing MM
%A confidence interval is computed for the distribution of squared
%distances to each of this idim multi-manifolds,
%It is computed by randomly selecting runs number 
%of test subsets whose cardinality is ~testproportion*m and their
%complemenary training sets. For each pair, a multi-manifold is computed 
%using the training set and the Expected value of the sum of the squared 
%distances EVsquareddist from the testing set is computed, along with the Variance,
%Standard Deviation, and the proportion of the testing data set 
%in the dyadic cubes supporting the multi-manifold. 
%The confidence interval is computed for the set of  EVsquaredist values.
%The structure MMbyidimstats contains the results of the confidence interval
%calculation. Tables are generated for both MMsummary and MMybydimstats.
%The information generated is sufficient for a hypothesis test. 
%Given a sample from another multi-manifold MM',hypothesized to represent data,
%the intrinsic dimension for each of its points can be computed, and
%the Expected Value of the  Distance from the of MM'to MM can be computed using
%the function TDtoMMbyidim. The resulting value can be compared to see if
%if is outside the confidence interval. If it is, the hypothesis that it
%represents the data set well can be rejected. 

rng(seed);
[m,n] = size(data);
%m data points in n dimensions
[ ndata ] = normalizetounitcube(data,eps );
uidim = unique(idim);
%find the subset of uidim of positive dimensions
%idim == -1 if the neighborhoods were  too sparse (< cutoff # of points)
indu = find(uidim.*(uidim > 0));
uidim = uidim(indu,:);
L = length(uidim);
%Compute Multi-Manifold for the entire data set as a baseline
[ MMsummary,MM] = MultiManifold2(ndata,idim,SVDbound,maxscale,K, varthresh,eps);
Runstats = zeros(runs*L,6); %cols run number, idim, EVsqdist, MMproportion,mtrain,mtest
for i = 1:runs
    display(sprintf('runs = %d',i))
    mtest = floor(testproportion*m);
    testindex = randi(m,mtest,1);
    %the above samples with replacement, so repetitions must be removed
    testindex = unique(testindex);
    mtest = length(testindex);
    trainindex = setdiff(1:m,testindex);
    mtrain  = length(trainindex);
    traindata = ndata(trainindex,:);
    testdata = ndata(testindex,:);
    trainidim = idim(trainindex,:);
    testidim = idim(testindex,:);
    trainuidim = unique(trainidim);
    testuidim = unique(testidim);
    [ trainMMsummary,trainMM] = MultiManifold2(traindata,trainidim,SVDbound,maxscale,K, varthresh,eps);
    [TDtoMM] = TDtoMMbyidim2( testdata,testidim,trainMM,trainuidim);
    %TDtoMM is a seven field structure, with one set of fields for each possible value of
    %idim. TDtoMM(i).idim = i. The other 6 fields of TDtoMM
    %are sqdist,MMcount,count, EVsqdist -- the sum of the squared distances
    %of the idim points in the test data to the idim multi-manifold determined 
    %by the training data, the number of idim
    %points in the cubes supporting the multi-manifold, the number of idim
    %test data points, and the expected value of the squared distances of the supported 
    %test data points. TDtoMM(i) is empty if that value of idim does not occur.
    for j = 1:L
        row = (i-1)*runs + j;
        if TDtoMM(j).count == 0
            Runstats(row,:) = [i,j,0,0,mtrain,mtest];
            if ~(testm == 0)
                display('error testm in  MMMethod');
            end
        else
            Runstats(row,:) = [i,j,TDtoMM(j).EVsqdist, TDtoMM(j).MMcount/TDtoMM(j).count,mtrain,mtest];
        end
    end
end
%change name of table so it is saved in the intended directory
MMsummarytable = struct2table(MMsummary)
%write(MMsummarytable);
for j = 1:L
    ind = find(Runstats(:,2) == j);
    ind2 = find(Runstats(ind,4) > 0);
    MMbyidimstats(j).idim = j;
    MMbyidimstats(j).EVEVsqdist = sum(Runstats(ind,3))/length(ind2);
    MMbyidimstats(j).EVMMsupportproportion = sum(Runstats(ind,4))/length(ind);
    MMbyidimstats(j).EVtraincount = sum(Runstats(ind,5))/length(ind);
    MMbyidimstats(j).EVtestcount = sum(Runstats(ind,6))/length(ind);
    MMbyidimstats(j).runs = runs;
    
    x = Runstats(ind,3);
    [ zscorecutoff,EV,Var,SD ] = SDdist( x ,confidencebound);
    MMbyidimstats(j).SDEVsqdist = SD;
    MMbyidimstats(j).zscorecutoff = zscorecutoff;
    
end
MMbyidimstatstable = struct2table(MMbyidimstats)
%change name of table so it is save in the intended directory
%writetable(MMbyidimstatstable)

save(fname,'data','SVDbound','maxscale','K','varthresh','eps');
save(fname,'testproportion','runs','seed','confidencebound','-append');
save(fname,'figs','fname','-append');
save(fname,'m','n','ndata','uidim','L','MMsummary','MM','-append');
save(fname,'Runstats','MMbyidimstats','-append');
fnametable = [fname,'-MMsummarytable.txt'];
write(MMsummarytable,fnametable);
fnametable = [fname,'-MMbyidimstatstable.txt'];
write(MMbyidimstatstable,fnametable);
        
    
    


end

