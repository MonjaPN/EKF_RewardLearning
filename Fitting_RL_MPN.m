%% Model fittin RL
%contact: monja-pascale.neuser@uni-tuebingen.de
%v1.0 Feb 04 2019

%%Modification with Rescorla Wagner algorithm to model RL data

% 2 parameter model
% alpha (learning rate), reward sensitivity



% =================
% TO DO
% 
%align simulation outoput and app data
%Conversion output-> table
%Export output
% =================
clear

%% RL parameters, fit_par is a structure containing parameter starting values
fit_par.alpha       = 0.3;  % learning rate of the simulated agent
fit_par.beta    = 2;    % reward sensitivity

count = 1;
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
    Data = ''; %add path of app data
    loat(Data)
    
else
    
    Simulation_RL_MPN;
    n_subj = 4; %number of beta distributions
    Data = 'C:\Users\Monja\PC_EKF\04_Review\RL\Output\simulation_RL_trials_1.csv'; %add path of app data
    load(Data)
    
end


for i_subj = 1 : n_subj
    
    %%Extract individual data for each subj i
    if simulate_flag == 1
        
        %Reduce to data of 1 distribution
        subj_data = simulation_RL_trials_1((simulation_RL_trials_1(:,2)==i_subj),:);
        %Reduce to relevant variables
        subj_data = subj_data(:,[7,8]);
    
    else
        
        subj_data = [];
        
    end
    
    
    
    
    %%Set initial values for fitting
    x0 = [fit_par.alpha, fit_par.beta]; %uses initial values as optimal starting values
    D =  subj_data;
    
    %Sets options for optimization. It might spam your console if you don't

    options = optimoptions('fmincon','Display','off','Algorithm','sqp');
    %options = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','off');
    mcon.A = [];
    mcon.b = [];
    mcon.lb = [0,0,0,0,0];
    mcon.ub = [1,Inf,Inf,1,1];
    
    
    [xout,fval,mcon.exitflag,mcon.out,mcon.lambda,mcon.grad,mcon.hessian] = fmincon(@(x)fit_RWmodel(x,D),x0,mcon.A,mcon.b,[],[],mcon.lb,mcon.ub,[], options);
    
    
    
%% Store output
% %original units, only to be used when fitting is run with original scale
% if simulate_flag == 0
% 
%     eval_fit(count,:) = [data_mat(i*150,1),data_mat(i*150,2),xout(:,1),xout(:,4),xout(:,2),xout(:,3),xout(:,5), fval];
%     
% else
%     
%     eval_fit(count,:) = [data_mat(i*150,1),xout(:,1),xout(:,4),xout(:,2),xout(:,3),xout(:,5), fval];
%     
% end
% 
% count = count + 1;
    

%transformed units
estp.alpha = xout(:,1); %1./(1+exp(xout(:,1)));
estp.beta = xout(:,2); %exp(xout(:,2));

eval_fit(i_subj,:) = [i_subj, sim_par.n_trials, sim_par.alpha, sim_par.beta, estp.alpha, estp.beta, estp.alpha - sim_par.alpha, estp.beta - sim_par.beta, fval, sum(reward)/sim_par.n_trials];

%out.part(:,:) = [ones(sim_par.n_trials,1) .* i_p, [1:sim_par.n_trials]', ones(sim_par.n_trials,1) .* sim_par.alpha, ones(sim_par.n_trials,1) .* sim_par.beta, choice, choice_p, reward_grid, reward, Q, rgw, draw_blue, est_prob, RPE_alpha, chosen_p, result];
out.part(:,:) = [ones(sim_par.n_trials,1) .* i_subj, [1:sim_par.n_trials]', ones(sim_par.n_trials,1) .* sim_par.alpha, ones(sim_par.n_trials,1) .* sim_par.beta, choice, choice_p, reward, Q, est_prob, RPE, chosen_p, result];

%out.mat = [out.mat; out.part];






end




% if simulate_flag == 0
% 
%     table_eval_fit = array2table(eval_fit,'VariableNames',{'Subject','Session','alpha','lapse','reward_sensitivity','punishment_sensitivity','action_bias','log_L'});
%     
% else
%     
%     table_eval_fit = array2table(eval_fit,'VariableNames',{'Subject','alpha','lapse','reward_sensitivity','punishment_sensitivity','action_bias','log_L'});
%     
% end

