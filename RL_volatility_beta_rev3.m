
%RL Tutorial 2017
%contact: nils.kroemer@tu-dresden.de

clear all
clc

%sim_par.alpha = 0.3; %learning rate of the simulated agent
%sim_par.beta = 0.5; %inverse temperature of the simulated agent
sim_par.gamma = 1; % risk sensitivity of the simulated agent
sim_par.n_trials = 150; %number of trials
% sim_par.p_win = .8; %probability of winning for the good option
% sim_par.p_lose = .2; %probability of losing for the good option
% sim_par.min_flip = 10; %minimum number of trials before another reversal can occur
sim_par.n_part = 50000; %number of simulated participants
out.mat = [];

for i_p = 1:sim_par.n_part
%% Generate probabilities of winning and losing for blue
rgw = nan(sim_par.n_trials,2);
rgw(1,1) = 0.80;
rgw(1,2) = 0.20;
sim_par.rgw.sigma = 0.10;        %0.075;
sim_par.rgw.reg =   0.03;        %0.05;

for t = 1:sim_par.n_trials-1

   step = normrnd(0,sim_par.rgw.sigma); 

   while rgw(t) + step - sim_par.rgw.reg * (rgw(t) - 0.5) < 0 || rgw(t) + step - sim_par.rgw.reg * (rgw(t) - 0.5) > 1

       step = normrnd(0,sim_par.rgw.sigma);

   end

   rgw(t+1,1) = rgw(t) + step - sim_par.rgw.reg * (rgw(t) - 0.5);
   rgw(t+1,2) = 1 - rgw(t+1,1);

end


%for i_p = 1:sim_par.n_part

    %initialize variables
    RPE = zeros(sim_par.n_trials,2);
    Q = zeros(sim_par.n_trials,2);
    choice_p = zeros(sim_par.n_trials,2);
    est_prob = 0.5 * ones(sim_par.n_trials,1);
    result = zeros(sim_par.n_trials,1);
    reward = zeros(sim_par.n_trials,1);
    reward_grid = zeros(sim_par.n_trials,2);
    choice = zeros(sim_par.n_trials,1);
    
    choice_p(1,:) = [0.5, 0.5];
%     good_opt = zeros(sim_par.n_trials,1);
%     count_flip = 0;
%    count_rev = 0;

    sim_par.beta =  (randi(10)/2-0.5);      %normrnd(3,0.75);
    sim_par.alpha = rand;                   %normrnd(0.4,0.15);
    
%     if sim_par.beta <  0.2
%         sim_par.beta = 0.2;
%     end
%     
%     if sim_par.alpha <  0.025
%         sim_par.alpha = 0.025;
%     end
   
    for i_t = 1:sim_par.n_trials

    %% Generate rewards for current trial
    reward_blue = round(normrnd(50,16));
    if reward_blue < 1
        reward_blue = 1;
    elseif reward_blue > 99
        reward_blue = 99;
    end
    reward_green = 100-reward_blue;
    reward_grid(i_t,1:2) = [reward_blue,reward_green];
    
        
    %% Make choice for blue or green
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
            %if reward_blue > reward_green    
                choice(i_t,1) = 1;
            else 
                choice(i_t,1) = 2;
            end
            
            chosen_p(i_t,1) = choice_p(i_t,choice(i_t,1)); 
        end

%% Get outcome and update outcome probability

    draw_blue(i_t,1) = double(rand < rgw(i_t,1));

        % See if agent wins based on p_win and p_lose
        if (choice(i_t,1) == 1)
            
            if (draw_blue(i_t,1) == 1)
                %result(i_t,1) = (double(rand < rgw(i_t,1)) - 0.5) * 2;
                %result(i_t,1) = double(rand < rgw(i_t,1));
                result(i_t,1) = 1;
            else
                result(i_t,1) = 0;
            end
            
        else
                
            if (draw_blue(i_t,1) == 1)
                result(i_t,1) = 0;
            else
                result(i_t,1) = 1;
            end
            %result(i_t,1) = (double(rand < rgw(i_t,2)) - 0.5) * 2;
            %result(i_t,1) = double(rand < rgw(i_t,2));

        end
        
        % Calculate reward based on outcome
        if result(i_t,1) == 1
            
            if choice(i_t,1) == 1
                
                reward(i_t,1) = reward_blue;
                
            else              
                reward(i_t,1) = reward_green;
                
            end
            
        elseif result(i_t,1) == 0
            if choice(i_t,1) == 1
                
                %reward(i_t,1) = 0;
                reward(i_t,1) = -reward_blue;
            else
                %reward(i_t,1) = 0;
                reward(i_t,1) = -reward_green;
            end
            
        end

        % Compute reward prediction error RPE
        RPE(i_t,choice(i_t,1)) = reward(i_t,1) - Q(i_t,choice(i_t,1));
        
        % Compute outcome prediction error
        RPE_alpha(i_t,1) = draw_blue(i_t,1) - est_prob(i_t,1);
        
        % Update outcome propability according to the former outcome probability, the learning rate alpha, and the RPE
        if i_t < sim_par.n_trials
            
        %    if choice(i_t,1) == 1
            
                est_prob(i_t+1,1) = est_prob(i_t,1) + sim_par.alpha * RPE_alpha(i_t,1);
                
        %    else
                
        %        est_prob(i_t+1,1) = (1 - est_prob(i_t,1)) + sim_par.alpha * RPE_alpha(i_t,1);
                %est_prob(i_t+1,1) = 1 - est_prob(i_t+1,1); % As outcome probability is defined in relation to the blue option,
                                                                                 % learning based on green choice needs to be reversed again
        %    end

        end

    end

    %set initial values for fitting
    x0 = [sim_par.alpha, sim_par.beta];
    D = [choice, (reward/50), draw_blue, reward_grid]; %concatenate data vector
    
    %options = optimset('Display','off');
    %options = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','off');
    options = optimoptions('fmincon','Display','off','Algorithm','sqp');
    %[xout,fval,exitflag,output] = fminunc(@(x)fit_model(x,D),x0,options);
    
    mcon.A = [];
    mcon.b = [];
    mcon.lb = [0,0];
    mcon.ub = [1,Inf];
    
    [xout,fval,mcon.exitflag,mcon.out,mcon.lambda,mcon.grad,mcon.hessian] = fmincon(@(x)fit_model(x,D),x0,mcon.A,mcon.b,[],[],mcon.lb,mcon.ub,[], options);

    %print output
    %fprintf('\n a=%.2f, b=%.2f, L=%.2f\n',xout(1),xout(2),fval);

%store output
%original units
%eval_fit(i_p,:) = [sim_par.alpha, sim_par.beta, xout(:,1), xout(:,2), fval, sum(reward)/sim_par.n_trials, count_rev*100/sim_par.n_trials];

%transformed units
estp.alpha = xout(:,1); %1./(1+exp(xout(:,1)));
estp.beta = xout(:,2); %exp(xout(:,2));

eval_fit(i_p,:) = [i_p, sim_par.n_trials, sim_par.alpha, sim_par.beta, estp.alpha, estp.beta, estp.alpha - sim_par.alpha, estp.beta - sim_par.beta, fval, sum(reward)/sim_par.n_trials];

out.part(:,:) = [ones(sim_par.n_trials,1) .* i_p, [1:sim_par.n_trials]', ones(sim_par.n_trials,1) .* sim_par.alpha, ones(sim_par.n_trials,1) .* sim_par.beta, choice, choice_p, reward_grid, reward, Q, rgw, draw_blue, est_prob, RPE_alpha, chosen_p, result];
%out.mat = [out.mat; out.part];
end

% fig_handle = figure('Position', [100, 100, 900, 700]);
% subplot(2,2,1)       % add first plot in 2 x 1 grid
% %histogram(eval_fit(:,1)-eval_fit(:,3))
% ecdf(eval_fit(:,1)-eval_fit(:,3),'bounds','on');
% title('ecdf of simulated vs. estimated alpha');
% 
% %figure
% %histogram(eval_fit(:,2)-eval_fit(:,4))
% subplot(2,2,2)       
% ecdf(eval_fit(:,2)-eval_fit(:,4),'bounds','on');
% title('ecdf of simulated vs. estimated beta');
% 
% %figure
% subplot(2,2,3)       
% ecdf(eval_fit(:,5),'bounds','on');
% title('ecdf of the obtained fit');
% 
% %figure
% subplot(2,2,4)  
% ecdf(eval_fit(:,6),'bounds','on');
% title('ecdf of win per trial');
% 
% figure
% ecdf(eval_fit(:,7),'bounds','on');
% title('ecdf of normalized reversal counts');

%  figure;histogram(out.mat(:,10),'normalization','pdf')
%  sum(out.mat(:,10))/30000
%  sum(out.mat(:,end))/30000
% sel_trials = find(out.mat(:,end)==1);
% sum(max(out.mat(sel_trials,8:9),[],2))/30000

save('example_NENA.mat','eval_fit');