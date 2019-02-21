
function [data_sub] = Simulate_influenca(sim_par)

%function to simulate data in the influenca game given the parameters alpha
%and beta and inputs(rewarded choices and reward values for each trial)

    reward_grid = sim_par.reward_grid;
    inputs = sim_par.inputs;
    correct = sim_par.correct_option;
    %initialize variables
    RPE = zeros(sim_par.n_trials,2);
    Q = zeros(sim_par.n_trials,2);
    choice_p = zeros(sim_par.n_trials,2);
    est_prob = 0.5 * ones(sim_par.n_trials,1);
    reward = zeros(sim_par.n_trials,1);
    
    choice = zeros(sim_par.n_trials,1);
    
    choice_p(1,:) = [0.5, 0.5];

    
  
    for i_t = 1:sim_par.n_trials

   
        
    %% Make choice for blue or green
    
    reward_blue = reward_grid(i_t,1);
    reward_green = reward_grid(i_t,2);
        %choice
        %1 = blue, 2 = green
        if i_t == 1
            
            % First choice is based on higher reward
            if reward_blue > reward_green
            
                choice(i_t,1) = 1; 
                est_prob(i_t,1) = 0.5; % Initialize r with 0.5 in the first round  
            else
                
                choice(i_t,1) = 2;
                est_prob(i_t,1) = 0.5; % Initialize r with 0.5 in the first round    
            end
            
            chosen_p(i_t,1) = 0.5;
        else
            
            % Calculate predicted outcome Q depending on risk aversion,
            % outcome probability and current reward
            Q(i_t,1) = max(min(sim_par.gamma*(est_prob(i_t,1)-0.5)+0.5,1),0)*(reward_blue/50);
            Q(i_t,2) = max(min(sim_par.gamma*((1-est_prob(i_t,1))-0.5)+0.5,1),0)*(reward_green/50);
            
            % Calculate probability for each choice
            choice_p(i_t,1) = 1/(1+exp(-sim_par.beta*(Q(i_t,1)-Q(i_t,2))));
            choice_p(i_t,2) = 1/(1+exp(-sim_par.beta*(Q(i_t,2)-Q(i_t,1))));
            
            % Get option with highest probability
            [Q_max, pref_opt(i_t,1)] = max(choice_p(i_t,:));

            % Set choice of agent
            if rand < choice_p(i_t,1)
   
                choice(i_t,1) = 1;
            else 
                choice(i_t,1) = 2;
            end
            
            chosen_p(i_t,1) = choice_p(i_t,choice(i_t,1)); 
        end

%% Get outcome and update outcome probability

    
        % See if agent wins based on p_win and p_lose

        if choice(i_t,1) == inputs(i_t,1) %good option
            reward(i_t,1) = 1;
        else
            reward(i_t,1) = -1;
        end
        reward(i_t,1) = reward(i_t,1)*reward_grid(i_t,choice(i_t,1));
        % Compute reward prediction error RPE
        RPE(i_t,choice(i_t,1)) = reward(i_t,1) - Q(i_t,choice(i_t,1));
        
        % Compute outcome prediction error
        RPE_alpha(i_t,1) = double(inputs(i_t,1)==1) - est_prob(i_t,1);
        
        % Update outcome propability according to the former outcome probability, the learning rate alpha, and the RPE
        if i_t < sim_par.n_trials
    

            
                est_prob(i_t+1,1) = est_prob(i_t,1) + sim_par.alpha * RPE_alpha(i_t,1);
  
        end

    end
    correct_ground = choice == correct;
 
    time = 1:sim_par.n_trials;
    beta = ones(sim_par.n_trials,1)*sim_par.beta;
    data_sub = [time', repmat(sim_par.dist,sim_par.n_trials,1),beta,est_prob,RPE_alpha,choice, (reward/50), correct_ground,double(inputs==1), inputs, reward_grid];
 
 

end

