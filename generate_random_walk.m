function out = generate_random_walk(n_runs, n_trials)
% Generates the random walk as in the App together with rewarded choice and
% reward vectors for blue and green. Generates n_runs walks with n_trials
% trials


for i_run = 1:n_runs
%% Generate random walk
    %Generate as many random walks as runs
%Initialize vector with probabilities for each trial
 probs = zeros(n_trials,1);
 reward_grid = zeros(n_trials,2);
 r = rand();
 probs(1) = 0.8 * (r > 0.5) + 0.2 * (r <= 0.5);

    for i_trial  = 1:n_trials-1
        
    	step = randn*0.1; % random scalar drawn from the standard normal distribution, convert to probability value
        
        %check if step is ok
    	while (probs(i_trial) + step - 0.03 * (probs(i_trial) - 0.5)) < 0 || (probs(i_trial) + step - 0.03 * (probs(i_trial) - 0.5)) > 1
    		
            step = randn*0.1; %if value >1, draw again
            
        end
        
        %Concatenate probability values
    	probs(i_trial+1) = probs(i_trial) + step - 0.03 * (probs(i_trial) - 0.5);
        
        reward_blue = round(normrnd(50,16));
        if reward_blue < 1
            reward_blue = 1;
        elseif reward_blue > 99
            reward_blue = 99;
        end
        reward_green = 100-reward_blue;
        reward_grid(i_trial,1:2) = [reward_blue,reward_green];

    end

    
    
%% determine correct choices depending on probabilities

%rand returns a single uniformly distributed random number in the interval (0,1).
gpt = rand(n_trials, 1) < probs; %1 blue, 2 green

good_opt = ones(n_trials,1);
good_opt(gpt) = 2;

correct = double(probs > .5);
correct = correct + 1; 

out.inputs(:,i_run) = good_opt;
out.correct_option(:,i_run) = correct;
out.reward_grid(:,:,i_run)= reward_grid;
out.probs(:,i_run) = probs;
  
end    



end