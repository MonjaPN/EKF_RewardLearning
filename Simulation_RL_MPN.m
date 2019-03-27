%================================================================
% Simulation of RL
% based on Influenca performance with RW-Model and different
% distributions of temperature parameter
%
% Code by Anne K�hnel 01/2019
% Modified by Monja Neuser 01/2019 
% Modified Anne K�hnel 15.02.2019
%================================================================

%clear

% set =1, if plots and data .csv should be written
draw_output = 1;
simulate_data = 1;
fit_data = 1;

%data_out = pwd;
%data_out = 'C:\Users\Monja\PC_EKF\04_Review\RL\Output\';
data_out = 'C:\Users\Monja\PC_EKF\04_Review\App_Analysis';

%% Define simulation parameters
%  define parameters
sim_par.beta_mean = [3, 3, 1.5, 1.5];
sim_par.beta_std = [0, 0.9, 0, 0.9];
sim_par.alpha = 0.6;

sim_par.trials = 150; % trials per run
sim_par.simulations = 10; % agents per group
sim_par.runs = 30; % runs per agent
sim_par.gamma = 1;



%% Generate beta values
if simulate_data
dist = [];
data = [];
tmp = 1;
beta_values = zeros(sim_par.simulations*length(sim_par.beta_mean),sim_par.runs);


% Initialize beta values from distributions 2D Matrix of beta-values(rows =
% agents, columns = runs) + 1 vector (dist) with length (agents) tracking
% the group
for i_dist = 1:4
    if sim_par.beta_std(i_dist) > 0
    pd_start(1) = makedist('Normal','mu',sim_par.beta_mean(i_dist),'sigma',sim_par.beta_std(i_dist));
    %pd(i_dist) = truncate(pd_start,sim_par.beta_mean(i_dist)-2*sim_par.beta_std(i_dist),sim_par.beta_mean(i_dist)+2*sim_par.beta_std(i_dist));
    pd(tmp) = truncate(pd_start,sim_par.beta_mean(i_dist)-1*sim_par.beta_std(i_dist),sim_par.beta_mean(i_dist)+1*sim_par.beta_std(i_dist));
 
    %fill beta values 
    beta_values((i_dist-1)*sim_par.simulations+1:i_dist*sim_par.simulations,:) = random(pd(tmp),sim_par.simulations,sim_par.runs);
    tmp = tmp + 1;
    else
    beta_values((i_dist-1)*sim_par.simulations+1:i_dist*sim_par.simulations,:) = sim_par.beta_mean(i_dist);
   
    end
    %concatenated vector with i_dist
    dist = [dist; repmat(i_dist,sim_par.simulations,1)];
end

    
%% Generate random walk and input/correct sequence + reward grid
rdm_wlk = generate_random_walk(sim_par.runs, sim_par.trials);

sim_par.inputs = rdm_wlk.inputs;
sim_par.correct_option = rdm_wlk.correct_option;
sim_par.reward_grid = rdm_wlk.reward_grid;
sim_par.probs = rdm_wlk.probs;

%% Simulate Reward-Learning with separate Q values for left / right

% Loop through list of beta values
for i_p = 1:size(beta_values,1) % loop through agents
    
    for i_r = 1:size(beta_values,2) %loop through runs

    sim_par_run.alpha = sim_par.alpha;
    sim_par_run.beta = beta_values(i_p,i_r);
    sim_par_run.reward_grid = sim_par.reward_grid(:,:,i_r);
    sim_par_run.inputs = sim_par.inputs(:,i_r);
    sim_par_run.correct_option = sim_par.correct_option(:,i_r);
    sim_par_run.gamma = sim_par.gamma;
    sim_par_run.n_trials = sim_par.trials;
    sim_par_run.dist = dist(i_p);
    


out = Simulate_influenca(sim_par_run);
% 'out' contains 'data_sub' derived from Simulation
% (1)Trial ID / (2)Distribution ID / (3)Beta value [sim_par] / 
% (4)est_prob / (5)RPE_alpha /
% (6)choice / (7)reward/50 /
% (8)correct_ground / (9)double(inputs==1) / (10)inputs /
% (11&12)reward_grid optA and optB



% Get summary results of each agent/run
summary_simulation(i_p,:,i_r) = [i_p,i_r,dist(i_p),sim_par_run.alpha, sim_par_run.beta, sum(out(:,7))*50/sim_par.trials, sum(out(:,8))];

ID = ones(sim_par.trials,1)*i_p;
run = ones(sim_par.trials,1)*i_r;

data_sub = [ID,run,out];

%Store summary and trialwise simulation
data = [data;data_sub];



    end

end

data_fit = data(:,[8,9,11,13,14,1,2,3,4,5]); %Reduce simulated data to the coloumns relevant for fitting
data_table = array2table(data, 'VariableNames',{'ID','Run','Trial','Group','beta','Estimated_Prob','RPE_prob','choice','reward','correct_choice','draw_blue','rewarded_choice','reward_blue','reward_green'}); %Get a table to export for R / nicer plotting 
% !! correct choice = choice == real prob>.5, good_option is the rewarded option for each trial 
end
%% Fit model 
if fit_data
if simulate_data
 
    D = data_fit;
else
    % Load empirical data (make sure coloumns have same order/content as simulated data) 
    %data_path and IDs to analyze have to be defined
    IDs = [1,4];
    % IDs = []; for all IDs leave IDs empty
    data_path = 'C:\Users\Anne\Documents\GitHub\EKF_RewardLearning\app_data_trial-07-02-2019.xlsx'
    D = read_app_data(IDs,data_path);
end

runs = 1;
tmp = unique(D(:,6));
fitted_data_plot =[];
for i_p = 1:length(unique(D(:,6)))
   
    DID = D(D(:,6)==tmp(i_p),:);
    
    for i_r = 1:length(unique(DID(:,7)))
        
    cD = DID(DID(:,7)==i_r,:);
    ID = cD(1,6);
    Run = cD(1,7);
    x0 = [0.5, 2];
    %options = optimset('Display','off');
    %options = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','off');
    options = optimoptions('fmincon','Display','off','Algorithm','sqp');
    %[xout,fval,exitflag,output] = fminunc(@(x)fit_model(x,D),x0,options);
    
    mcon.A = [];
    mcon.b = [];
    mcon.lb = [0,0];
    mcon.ub = [1,Inf];
    
    [xout,fval,mcon.exitflag,mcon.out,mcon.lambda,mcon.grad,mcon.hessian] = fmincon(@(x)fit_model(x,cD),x0,mcon.A,mcon.b,[],[],mcon.lb,mcon.ub,[], options);
    estp.alpha = xout(:,1); %1./(1+exp(xout(:,1)));
    estp.beta = xout(:,2); %exp(xout(:,2));

  
 %   eval_fit(runs,:) = [ID, Run, sim_par.n_trials, sim_par.alpha, sim_par.beta, estp.alpha, estp.beta, estp.alpha - sim_par.alpha, estp.beta - sim_par.beta, fval];
    eval_fit(runs,:) = [ID, Run, size(cD,1), estp.alpha, estp.beta, fval, cD(1,10), cD(1,9)]; % This only works for simulated data
 
    fitted_data(runs) = get_fitted_values(cD, estp);
    fitted_data_plot = [fitted_data_plot;fitted_data(runs).run];
    runs = runs + 1;
    end
end



save_str = ['Fitted_results.mat']; %Differentiate simulated and real?
save(save_str,'fitted_data_plot')

end
%% Simulate food intake

simulated_cal = simulate_food_intake(eval_fit(:,5));




%% Plotting

% Some preliminary plotting
% barplots beta

for i_group = 1:4

    m_beta(i_group) = mean(beta_values(dist==i_group));
    std_beta(i_group) = std(beta_values(dist==i_group));
    performance(i_group) = mean(eval_fit(dist==i_group,4));
    std_performance(i_group) = std(eval_fit(dist==i_group,4));
    for i_trial = 1:150
        
        trialwise_choice (i_trial,i_group) = mean(data(data(:,2)==i_group&data(:,1)==i_trial,7));
        Q_1_trial (i_trial,i_group) = mean(data(data(:,2)==i_group&data(:,1)==i_trial,5));
    end    

end



 if draw_output == 1

%Plot simulations


dat = data(:,[1:5]);
simpar_betas = dat((dat(:,3)==1),:);


x = simpar_betas(:,5);
y = eval_fit(:,5);
group = simpar_betas(:,4);

gscatter(x,y,group)
axis([0 6 0 6])
    ylabel('estimated betas')
    xlabel('simulated betas')
    diag = refline([1 0]);
    %diag.Color = 'r';
    
%    refline(1,0)





% 
%     
%     out_fig1 = [data_out,'simulation_scatter_',num2str(i_inputs),'.png'];
%     %figure.Position=[0 0 [] []];
%     figure(1)
% 
%     subplot(3,2,1)
%     hold on
%     bar(1:4,m_beta)
%     errorbar(1:4,m_beta,std_beta,'.')
%     ylabel('beta')
% 
%     subplot(3,2,2)
%     hold on
%     bar(1:4,performance)
%     errorbar(1:4,performance,std_performance,'.')
%     ylabel('trials rewarded')
% 
%     subplot(3,2,3:6)
%     hold on
%     scatter(eval_fit(:,3),eval_fit(:,5),[],eval_fit(:,1),'filled')
%     legend({'high_b low_VAR','high_b high_VAR','low_b low_VAR','low_b highVAR'},'Location', 'NorthWest','FontSize',7,'Orientation','horizontal')
%     ylabel('correct choices')
%     xlabel('beta')
%      saveas(gcf,out_fig1)
% 
% 
% 
%  if i_dist <= 4
%  
%     out_fig2 = [data_out,'simulation_Qvalues',num2str(i_inputs),'.png'];
%     %figure.Position=[0 0 [] []];
%     figure(2)
%     
%     for i_plot=1:i_dist
%         
%         plot_data=data((data(:,2)==i_plot),:);
%         plot_data=plot_data(1:(sim_par.trials*2),:);
%         label_fig = ['beta_type_',num2str(i_plot)];
% 
%         subplot(i_dist,1,i_plot)
%         hold on
%         plot(1:sim_par.trials,plot_data(1:sim_par.trials,3))
%         plot(1:sim_par.trials,plot_data(sim_par.trials+1:sim_par.trials*2,3))
% 
%         ylabel('Q-value OptA')
%         xlabel(label_fig)
% 
%         % subplot(2,2,3:4)
%         % hold on
%         % plot(1:sim_par.trials,Q_1_trial)
%         % %legend({'high low','high high','low low','low high'},'Location', 'NorthWest','FontSize',7,'Orientation','horizontal')
%         % legend({'high_b low_VAR','high_b high_VAR','low_b low_VAR','low_b highVAR'},'Location', 'NorthWest','FontSize',7,'Orientation','horizontal')
%         % ylabel('Q_Value_trials')
% 
%         saveas(gcf,out_fig2)
% 
%     end
%     
%  end
%  
%  
 end




% Write output in csv format to designated diectory (set path at line 16/17)
% out_trials = [data_out,'simulation_RL_trials_',num2str(i_inputs),'.csv'];
% out = [data_out,'simulation_RL_',num2str(i_inputs),'.csv'];
out_trials = [data_out,'simulation_RL_trials.csv'];
out = [data_out,'simulation_RL_est.csv'];
csvwrite(out_trials, data)
csvwrite(out, eval_fit)

close all












