%% Model fitting for Influenca data / simulations
% Based on Reinforcement learning tutorial, 2017
%contact: nils.kroemer@uni-tuebingen.de @cornu_copiae
%v1.1 Dec-19 2017, Seminar Computational Psychiatry
%

function L = fit_model(x,D)

%Model Parameters, original units
alpha = x(1);
beta = x(2);

%Model Parameters, transformed for unconstrained search
%alpha = 1/(1+exp(-x(1)));
%beta  = exp(x(2)); 
gamma = 1;

%Parse Data
rew = D(:,2);
choice = D(:,1);
draw_blue = D(:,3);
rew_grid = D(:,4:5);
 
%Initialization of values for the model
L = 0;                          %initial value sum of squared error
RPE = zeros(length(choice),2);  %initial RPEs
RPE_alpha = zeros(length(choice),1);  %initial RPEs
Q = zeros(length(choice),2);    %initial Q values
est_prob = 0.5 * ones(length(choice),1); % initial estimates of probability 

%loop through trials
for t = 1:length(choice)

    %calculates RPEs
    %RPE(t,choice(t,1)) = rew(t,1) - Q(t,choice(t,1));
    
    % Compute reward prediction error RPE
    RPE(t,choice(t,1)) = rew(t,1) - Q(t,choice(t,1));

    % Compute outcome prediction error
    RPE_alpha(t,1) = draw_blue(t,1) - est_prob(t,1);
    
    
    %updates Q values
    if t < length(choice)
        %Q(t+1,1) = Q(t,1) + alpha * RPE(t,1);
        %Q(t+1,2) = Q(t,2) + alpha * RPE(t,2);
        est_prob(t+1,1) = est_prob(t,1) + alpha * RPE_alpha(t,1);
        
        Q(t+1,1) = max(min(gamma*(est_prob(t+1,1)-0.5)+0.5,1),0)*(rew_grid(t+1,1)/50);
        Q(t+1,2) = max(min(gamma*((1-est_prob(t+1,1))-0.5)+0.5,1),0)*(rew_grid(t+1,2)/50);
 
    end
    
    %calculates loglikelihood based on discrepancies between observed and
    %predicted choices
    %L = L + log( exp(Q(t,choice(t,1)) .* beta) / (exp(Q(t,1) .* beta) + exp(Q(t,2) .* beta)));
    L = L + log( exp(Q(t,choice(t,1)) .* beta) / (exp(Q(t,1) .* beta) + exp(Q(t,2) .* beta)));
    
end

L=-L; %negative due to function minimization
 
end