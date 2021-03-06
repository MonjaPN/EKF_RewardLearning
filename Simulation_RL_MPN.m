%================================================================
% Simulation of RL
% based on Influenca performance with RW-Model and different
% distributions of temperature parameter
%
% Code by Anne K�hnel 01/2019
% Modified by Monja Neuser 01/2019 
% Last correction 31.01.2019 correct vector
%================================================================

clear

% set =1, if plots and data .csv should be written
draw_output = 1;

%  define parameters
sim_par.beta_mean = [5, 5, 2, 2];
sim_par.beta_std = [0.1, 0.9, 0.1, 0.9];
sim_par.alpha = 0.6;

sim_par.trials = 150;
sim_par.simulations = 100;
% generate distributions

beta_values = [];
dist = [];
data = [];


% Initialize beta values from distributions
for i_dist = 1:4
    pd_start(1) = makedist('Normal','mu',sim_par.beta_mean(i_dist),'sigma',sim_par.beta_std(i_dist));
    %pd(i_dist) = truncate(pd_start,sim_par.beta_mean(i_dist)-2*sim_par.beta_std(i_dist),sim_par.beta_mean(i_dist)+2*sim_par.beta_std(i_dist));
    pd(i_dist) = truncate(pd_start,sim_par.beta_mean(i_dist)-2*sim_par.beta_std(i_dist),sim_par.beta_mean(i_dist)+2*sim_par.beta_std(i_dist));
    
    %concatenated vector of beta values of length 'simulations', drawn from disributions
    beta_values = [beta_values; random(pd(i_dist),sim_par.simulations,1)];
    %concatenated vector with i_dist
    dist = [dist; repmat(i_dist,sim_par.simulations,1)];
end


% Loop through number of ??
for i_inputs = 1:1
    
    
%% Generate random walk and input/correct sequence

%Initialize vector with probabilities for each trial
 probs = zeros(sim_par.trials,1);
 r = rand();
 probs(1) = 0.8 * (r > 0.5) + 0.2 * (r <= 0.5);

    for i_trial  = 1:sim_par.trials-1
        
    	step = randn*0.1; % random scalar drawn from the standard normal distribution, convert to probability value
        
        %check if step is ok
    	while (probs(i_trial) + step - 0.03 * (probs(i_trial) - 0.5)) < 0 || (probs(i_trial) + step - 0.03 * (probs(i_trial) - 0.5)) > 1
    		
            step = randn*0.1; %if value >1, draw again
            
        end
        
        %Concatenate probability values
    	probs(i_trial+1) = probs(i_trial) + step - 0.03 * (probs(i_trial) - 0.5);

    end

    
    
%% determine correct choices depending on probabilities

%rand returns a single uniformly distributed random number in the interval (0,1).
gpt = rand(sim_par.trials, 1) < probs; %1 left, 2 right

good_opt = ones(sim_par.trials,1);
good_opt(gpt) = 2;
%correct = probs > .5;
correct = double(probs > .5);
correct = correct + 1; 



%% Simulate Reward-Learning with separate Q values for left / right

% Loop through list of beta values
for i_p = 1:length(beta_values)

    %initialize variables
    RPE = zeros(sim_par.trials,2);
    Q = zeros(sim_par.trials,2);
    choice_p = zeros(sim_par.trials,2);
    reward = zeros(sim_par.trials,1);
    choice = zeros(sim_par.trials,1);
    good_opt = good_opt;
    beta = beta_values(i_p);

    
    %Loop through trials
    for i_t = 1:sim_par.trials

        %choice
        %1 = OptA, 2 = OptB
        if i_t == 1
            choice(i_t,1) = randi(2,1); %first choice is random
        else
            %[Q_max, best_opt] = max(Q(i_t,:));
            %choice(i_t,1) = best_opt;

            OptA = exp(Q(i_t,1)*beta);
            OptB = exp(Q(i_t,2)*beta);
            choice_p(i_t,1) = OptA / (OptA + OptB);
            choice_p(i_t,2) = OptB / (OptA + OptB);

            [Q_max, pref_opt] = max(choice_p(i_t,:));

            if rand < choice_p(i_t,1)
               choice(i_t,1) = 1;
            else
               choice(i_t,1) = 2;
            end

    %        choice(i_t,1) = pref_opt;

        end

  
        %calculate outcome "reward" based on p_win and p_lose
        if choice(i_t,1) == good_opt(i_t,1) %good option
            reward(i_t,1) = 1;
        else
            reward(i_t,1) = -1;
        end

        %compute reward prediction errors, RPEs, according to delta rule,
        %single update
        RPE(i_t,choice(i_t,1)) = reward(i_t,1) - Q(i_t,choice(i_t,1));

        %update Q values according to the learning rate, alpha, and the RPEs
        if i_t < sim_par.trials
            Q(i_t+1,1) = Q(i_t,1) + sim_par.alpha * RPE(i_t,1);
            Q(i_t+1,2) = Q(i_t,2) + sim_par.alpha * RPE(i_t,2);
        end

    end

 correct_ground = correct == choice;
    
%original units
eval_fit(i_p,:) = [dist(i_p),sim_par.alpha, beta, sum(reward)/sim_par.trials, sum(correct_ground)];
time = 1:sim_par.trials;
beta_vector = ones(sim_par.trials,1)*beta;


data_sub = [time', repmat(dist(i_p),sim_par.trials,1),Q,RPE,choice, reward, correct_ground,good_opt];
data = [data;data_sub];

end

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

    
out_fig1 = ['simulation_scatter_',num2str(i_inputs),'.png'];
%figure.Position=[0 0 [] []];
figure(1)


subplot(3,2,1)
hold on
bar(1:4,m_beta)
errorbar(1:4,m_beta,std_beta,'.')
ylabel('beta')

subplot(3,2,2)
hold on
bar(1:4,performance)
errorbar(1:4,performance,std_performance,'.')
ylabel('trials rewarded')

subplot(3,2,3:6)
hold on
scatter(eval_fit(:,3),eval_fit(:,5),[],eval_fit(:,1),'filled')
legend({'high_b low_VAR','high_b high_VAR','low_b low_VAR','low_b highVAR'},'Location', 'NorthWest','FontSize',7,'Orientation','horizontal')
ylabel('correct choices')
xlabel('beta')
 saveas(gcf,out_fig1)


out_fig2 = ['simulation_Qvalues',num2str(i_inputs),'.png'];
%figure.Position=[0 0 [] []];
figure(2)


if i_dist <= 4
    
    
for i_plot=1:i_dist
    
    
    
    plot_data=data((data(:,2)==i_plot),:);
    plot_data=plot_data(1:(sim_par.trials*2),:);
    label_fig = ['beta_type_',num2str(i_plot)];

subplot(i_dist,1,i_plot)
hold on
plot(1:sim_par.trials,plot_data(1:sim_par.trials,3))
plot(1:sim_par.trials,plot_data(sim_par.trials+1:sim_par.trials*2,3))

ylabel('Q-value OptA')
xlabel(label_fig)

% subplot(2,2,3:4)
% hold on
% plot(1:sim_par.trials,Q_1_trial)
% %legend({'high low','high high','low low','low high'},'Location', 'NorthWest','FontSize',7,'Orientation','horizontal')
% legend({'high_b low_VAR','high_b high_VAR','low_b low_VAR','low_b highVAR'},'Location', 'NorthWest','FontSize',7,'Orientation','horizontal')
% ylabel('Q_Value_trials')
 saveas(gcf,out_fig2)

end
end
 
 
end





out_trials = ['simulation_RL_trials_',num2str(i_inputs),'.csv'];
out = ['simulation_RL_',num2str(i_inputs),'.csv'];
csvwrite(out_trials, data)
csvwrite(out, eval_fit)

close all




end







