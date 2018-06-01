function [SSVEnergy,idimSSVEnergy,EidimSSVEnergy,TV,ETV,idimTV,EidimTV ] = TotalVarianceFunctions(data,idim, SG,radii,q,scales,indG,TVfigs,fname)
%data is an m xn matrix of m data points with n coordinate SG is an m x q
%matrix of total variances for data in each of the q balls with radii
%specified in radius (total variances of the centered data) SG entries for
%a ball are -1 if the intrinsic dimension was not computed for that ball
%due to too few points etc scales is an m x q matrix of 0's and 1's
%specifying the scales at which idim was observed for each data point indG
%is an m x q matrix of 0's and 1's specifying the scales for which
%intrinsic dimension was computed A variety of total variance functions are
%computed including the SSVEnergy function we computed at WiSDM if TVfig ==
%1, the total variance functions are plotted for the subset of points for
%which intrinsic dimension was computed.

[m,n] = size(data);
%SG is an m x q matrix containg the total variances for each ball (or -1's
%for the ball where SG and hence the intrinsic dimension were not computed

%Compute different Versions of Total Variance Functions TV sum of total
%variances over the relevant balls
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
    
    %h = figure; x = 1:m; hold on for i = 1:q
        %plot(x,SG(:,i));
    %end titlestr = sprintf('Total Variance Curves for each scale');
    %title(titlestr) hold off
    
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
            display('line73 Total Variance Functions')
            fname1 = [fname,'idimTV','.png']
            saveas(h,fname1);
            fname1 = [fname,'idimTV','.fig']
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

