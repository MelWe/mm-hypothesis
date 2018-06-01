
function [idimsig, idimsigstats] = idimprob2( CumG,q,groupdim,t,indG,ballcount)
%CumG is a matrix whose rows are q groups of cumulative distributions
%of size groupdim
%t is a threshhold 0 <= t <= 1
%The function returns column index of the first element in each row
%which is >= t; if no such index exists, the function returns -1.
[m,n] = size(CumG);
if ~(n == q*groupdim)
    idimsig = [];
    firstscaleindex = [];
    scales = [];
    idimG = [];
    mxG = [];
    scaleprob = [];
    consecutive = [];
    display('in idimprob -- CumG does not have the expected number of columns');
else
    idimsig = -1*ones(m,1);
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
[ idimGsig,scalessig,totalsigscales ] = significantidim(idimG,ballcount );
%The intrinsic dimension for a point is the minimum of the intrinsic dimensions for
%each scale for which the intrinsic dimension is computed. 
%If it is not computed at any scale idim = -1.
for i = 1:m
    %ind1 = find(indG(i,:) == 1);
    ind1 = find(idimGsig(i,:)) > 0;
    if length(ind1) > 0
        idimsig(i,1) = min(idimGsig(i,ind1));
        indscales = find(idimGsig(i,:) == idimsig(i,1));
        scales(i,indscales) = 1;
        scaleprob(i,1) = length(indscales)/length(ind1);
        if indscales == indscales(1) - 1 + [1:length(indscales)]
            consecutive(i,1) = 1;
        end
        firstscaleindex(i,1) = indscales(1);
    end
end

totalidimscales = sum(scales,2);

idimsigstats(1).idimsig = idimsig;
idimsigstats(1).firstscaleindex = firstscaleindex;
idimsigstats(1).indG = indG;
idimsigstats(1).scales = scales;
idimsigstats(1).idimG = idimG;
idimsigstats(1).scaleprob = scaleprob;
idimsigstats(1).scales = scales;
idimsigstats(1).consecutive = consecutive;
idimsigstats(1).dense = dense;
idimsigstats(1).SG = SG;




end

