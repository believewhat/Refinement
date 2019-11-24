function [model_pos,model_neg ] = FindGuassianModel( x,y )
%FINDGUASSIANMODULE Summary of this function goes here
%   Detailed explanation goes here
x_pos = x(:,y==1);
model_pos.mu = mean(x_pos,2);
model_pos.var = cov(x_pos');
model_pos.prior = length(x_pos)/length(x);




x_neg = x(:,y~=1);
model_neg.mu = mean(x_neg,2);
model_neg.var = cov(x_neg');
model_neg.prior = length(x_neg)/length(x);

end

