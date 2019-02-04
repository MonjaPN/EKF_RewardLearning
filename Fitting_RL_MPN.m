%% Model fittin RL
%contact: monja-pascale.neuser@uni-tuebingen.de
%v1.0 Feb 04 2019

%%Modification with Rescorla Wagner algorithm to model RL data

% 2 parameter model
% alpha (learning rate), reward sensitivity

clear

%% RL parameters, fit_par is a structure containing parameter starting values
fit_par.alpha       = 0.3;  % learning rate of the simulated agent
fit_par.beta_rew    = 2;    % reward sensitivity


% Data contains:

% data output from Simulation_RL_MPN
%
%(1) trial ID
%(2) Distribution ID
%(§) Q value Option A
%(4) Q value Option B
%(5) RPE Option A (updated if option A chosen)
%(6) RPE Option B (updated if option B chosen)
%(7) choice (1 = OptA, 2 = OptB)
%(8) reward (1 win, -1 lose)
%(9) correct_ground
%(10)good_opt


simulate_flag = 1;

if simulate_flag == 0
    
    n_subj = [1:1]; %number of app data sets
    Data = '...'; %add path of app data
    loat(Data)
    
else
    
    Simulation_RL_MPN;
    n_subj = 4; %number of beta distributions
    Data = '...'; %add path of app data
    load(Data)
    
end





