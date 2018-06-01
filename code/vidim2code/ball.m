function [v,B] = ball( data,p,r) 
%% p=center; r=radius %%
%This functions finds the data points in the ball of radius r 
%centered at p. 
[m,n] = size(data);
x = sqrt(sum((data- ones(m,1)*p).^2,2));
B= x < r*ones(m,1);
v=data(B,:);
end


