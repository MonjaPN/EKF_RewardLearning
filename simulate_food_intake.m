function[cal] = simulate_food_intake(beta,ID)
% beta and ID have to have same length
% define parameters
%

intercept_mean = 1800;
intercept_sd = 1800/6;
b1_mean = 400;
b1_sd = 400/6;

switch nargin
    case 2
        IDs = unique(ID);
        intercept = zeros(length(beta),1);
        b1 = zeros(length(beta),1);
        
        for i_sub = 1:length(IDs)
            
            trials = sum(ID == IDs(i_sub));
            
            intercept(ID==IDs(i_sub)) = repmat(normrnd(intercept_mean,intercept_sd,[1,1]),trials,1);
            b1(ID==IDs(i_sub)) = repmat(normrnd(b1_mean,b1_sd,[1,1]),trials,1);
            
            
        end
    case 1
         intercept = normrnd(intercept_mean,intercept_sd,[length(beta),1]);
         b1 = normrnd(b1_mean,b1_sd,[length(beta),1]);
        
end

error_sd = 400;

eps = normrnd(0,error_sd,[length(beta),1]);



cal = intercept + b1.*beta + eps;



end 