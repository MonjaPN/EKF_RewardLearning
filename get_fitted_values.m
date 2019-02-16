function out = get_fitted_values(D, pars)


  

        cD = D;
        %p_winA = T.p_reward_a(T.subject_id_app==ID&T.run_ind==i,:);
        p_winA = cD(:,9);

        alpha = pars.alpha;
        beta = pars.beta;
        gamma = 1;

        rew = cD(:,2);
        choice = cD(:,1);
        draw_blue = cD(:,3);
        rew_grid = cD(:,4:5);

        RPE = zeros(length(choice),2);  %initial RPEs
        RPE_alpha = zeros(length(choice),1);  %initial RPEs
        Q = zeros(length(choice),2);    %initial Q values
        est_prob = 0.5 * ones(length(choice),1); % initial estimates of probability 

        for t = 1:length(choice)

            % Compute reward prediction error RPE
            RPE(t,choice(t,1)) = rew(t,1) - Q(t,choice(t,1));

            % Compute outcome prediction error
            RPE_alpha(t,1) = draw_blue(t,1) - est_prob(t,1);

            % Set expectations according to the initial values
            Q(t,1)=rew_grid(t,1)/50;
            Q(t,2)=rew_grid(t,2)/50;

            %updates Q values
            if t < length(choice)
                %Q(t+1,1) = Q(t,1) + alpha * RPE(t,1);
                %Q(t+1,2) = Q(t,2) + alpha * RPE(t,2);
                est_prob(t+1,1) = est_prob(t,1) + alpha * RPE_alpha(t,1);

                Q(t+1,1) = max(min(gamma*(est_prob(t+1,1)-0.5)+0.5,1),0)*(rew_grid(t+1,1)/50);
                Q(t+1,2) = max(min(gamma*((1-est_prob(t+1,1))-0.5)+0.5,1),0)*(rew_grid(t+1,2)/50);

            end

            p_choice(t,1) = exp(Q(t,choice(t,1)) .* beta) / (exp(Q(t,1) .* beta) + exp(Q(t,2) .* beta));
        end

        del_rew = rew_grid(:,1)/50-rew_grid(:,2)/50;
        out.run = [cD(:,6), cD(:,7), cD(:,8), draw_blue, rew, del_rew, abs(del_rew), choice, p_winA, est_prob, Q, Q(:,1)-Q(:,2), p_choice, choice-1];
     



end