function [data,idim,idimstats] = computevidim2( data,radii,K,varthreshhold,epszero,idimfigs,TVfigs, fname)
%This computes variance-based intrinsic dimension for the 
%m x n data matrix. (The rows are the data points.)
%Two summaries idimsummary and istatssummary are displayed (and saved).
%The results are saved in the file [dirname,fname,'-vidim-']

[ idimsummary,istatssummary,idim,idimstats,dense,SG,SSVs] = variancebaseddim2( data,radii,K,varthreshhold,epszero,idimfigs,fname );
idimsummary
display('columheadings: idim,firstscaleindex, count');
istatssummary
save(fname,'varthreshhold', 'K', 'radii', 'epszero', 'TVfigs', 'fname');
save(fname,'data','-append');
save(fname, 'idimsummary','istatssummary','idim','idimstats','dense','SG','SSVs','-append');
q = length(radii);
scales = idimstats(1).scalessig;
indG = idimstats(1).idimGsig;
[SSVEnergy,idimSSVEnergy,EidimSSVEnergy,TV,ETV,idimTV,EidimTV ] = TotalVarianceFunctions(data,idim,SG,radii,q,scales,indG,TVfigs,fname);
save(fname,'SSVEnergy','idimSSVEnergy','EidimSSVEnergy','TV','ETV','idimTV','EidimTV','-append');   
end




