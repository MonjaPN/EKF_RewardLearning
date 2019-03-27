function[cal] = simulate_food_intake(beta)

% define parameters
%
x1 = beta;

intercept_mean = 1800;
intercept_sd = 1800/6;
b1_mean = 400;
b1_sd = 400/6;

error_sd = 400;

eps = normrnd(0,error_sd,[length(x1),1]);
intercept = normrnd(intercept_mean,intercept_sd,[length(x1),1]);
b1 = normrnd(b1_mean,b1_sd,[length(x1),1]);



cal = intercept + b1.*x1 + eps;



end 