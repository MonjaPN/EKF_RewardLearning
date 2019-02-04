%% Reinforcement learning tutorial, WiSe2017
%contact: nils.kroemer@uni-tuebingen.de @cornu_copiae
%v1.2 Dec 18 2018, Seminar Computational Psychiatry

function L = fit_RWmodel(x,D)

%Model Parameters, transformed for unconstrained search
% alpha = 1/(1+exp(-x(1)));
% rs = exp(x(2));

alpha = x(1);
beta = x(2);

% Parse Data
reward = D(:,2);
choice = D(:,1);
%cond_stimulus = D(:,1); 
 
%Initialization of values for the model
L_t = zeros(length(choice),1);    %initial value sum of squared error
weights = zeros(length(choice),8); % initial action weights
RPE = zeros(length(choice),8);  %initial RPEs
Q = zeros(length(choice),8);    %initial Q values
%V = zeros(length(choice),4);





%loop through trials
for i = 1:length(choice)
    

    % calculate weights -> implement one action weight per condition and choice
%     weights(i,cond_stimulus(i)) = Q(i,cond_stimulus(i)); %nogo
%     weights(i,cond_stimulus(i)+4) = Q(i,cond_stimulus(i)+4) + a_bias; % go


    
    
    %compute reward prediction errors, RPEs, according to delta rule,
    %single update
    RPE(i,choice(i,1)) = reward(i,1) - Q(i,choice(i,1));

    
    
    %updates Q values

    if i < length(choice)
        
            Q(i+1,1) = Q(i,1) + alpha * RPE(i,1);
            Q(i+1,2) = Q(i,2) + alpha * RPE(i,2);
        
    end
    
    %???Calculate choic probability
    % ...
    
    %calculates loglikelihood based on discrepancies between observed and
    %predicted choices
    L_t(i) = log(action_probability(i));
   
 
end

L=sum(L_t);
L=-L; %negative due to function minimization
 
end