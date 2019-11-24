function [err,h] = FindModelError(model_pos,model_neg, x,y )
%FINDGUASSIANMODULE Summary of this function goes here
%   Detailed explanation goes here
mu1 = model_pos.mu;
sigma1 = model_pos.var;
p1 = model_pos.prior;

mu2 = model_neg.mu;
sigma2 = model_neg.var;
p2 = model_neg.prior;

bias = 0.5*log(det(sigma2))-0.5*log(det(sigma1))+log(p1/p2);
err = 0;
h = zeros(size(y));
for i=1:length(y)
   c = bias + 0.5*(x(:,i)-mu2)'/sigma2*(x(:,i)-mu2) - 0.5*(x(:,i)-mu1)'/sigma1*(x(:,i)-mu1);
   if c > 0
       h(i) = 1;
   else
       h(i) = -1;
   end
   if h(i)~=y(i)
       err = err + 1;
   end   
    
end

end