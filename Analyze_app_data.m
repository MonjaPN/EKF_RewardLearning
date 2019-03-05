clear

str_path = 'C:\Users\Nils\Google Drive\TUE_general\Projects\Open Science Fellowship\App\App_Data\app_data_trial-07-02-2019.xlsx';
%str_path = 'C:\Users\Nils\Google Drive\TUE_general\Projects\Open Science Fellowship\App\App_Data\app_data_trial_Nils.xlsx';
%T = readtable(str_path,'ReadVariableNames',1,'Format','%u%u%u%f%f%f%s%u%u%s%u%f%u%f%s%s');
T = readtable(str_path,'ReadVariableNames',1);
%ID = T(:,1);
T = T(T.run_finished == 1,:);
tbl = tabulate(T.subject_id_app);
ID = [1];

for i_p = 1:length(ID)

    choice = ismember(T.choice,'B');
    draw_blue = ismember(T.drawn_outcome,'A');
    reward_grid = [T.reward_a, T.reward_b]; 
    for i=1:length(choice)
        reward(i,1) = reward_grid(i,choice(i)+1) * (2 * (T.win(i) - 0.5));
    end

    %set initial values for fitting
    x0 = [0.5, 2];
    D = [choice+1, (reward/50), draw_blue, reward_grid]; %concatenate data vector
    iDD = D(T.subject_id_app==ID(i_p),:);
    IDTab = T(T.subject_id_app==ID(i_p),:);

    for i=1:length(IDTab.trial_ind)

        if IDTab.trial_ind(i) == 1 && IDTab.run_finished(i)== 0

            if IDTab.run_finished(i+1) == 1
                IDTab.run_finished(i) = 1;
            end
        end   
    end

    IDTab = IDTab(IDTab.run_finished==1,:);
    sess_tbl = tabulate(IDTab.run_ind);

    for i=1:length(sess_tbl)
        %cD = D(T.subject_id_app==ID&T.run_ind==i,:);
        cD = iDD(IDTab.run_ind==i,:);
        %options = optimset('Display','off');
        %options = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','off');
        options = optimoptions('fmincon','Display','off','Algorithm','sqp');
        %[xout,fval,exitflag,output] = fminunc(@(x)fit_model(x,cD),x0,options);

        mcon.A = [];
        mcon.b = [];
        mcon.lb = [0,0];
        mcon.ub = [1,Inf];

        [xout,fval,mcon.exitflag,mcon.out,mcon.lambda,mcon.grad,mcon.hessian] = fmincon(@(x)fit_model(x,cD),x0,mcon.A,mcon.b,[],[],mcon.lb,mcon.ub,[], options);

        out.par(i,:) = [ID(i_p), i, xout, fval, mean(cD(:,2))*50]; 

    end

    save(['Fitted_results_ID' num2str(ID(i_p)) '.mat'],'out');
    out.plot = [];

    for i=1:length(sess_tbl)

        %cD = D(T.subject_id_app==ID&T.run_ind==i,:);
        cD = iDD(IDTab.run_ind==i,:);
        %p_winA = T.p_reward_a(T.subject_id_app==ID&T.run_ind==i,:);
        p_winA = IDTab.p_reward_a(IDTab.run_ind==i,:);

        alpha = out.par(i,3);
        beta = out.par(i,4);
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
        out.run = [ones(150,1)*ID(i_p), ones(150,1)*i, (1:150)', rew_grid/50, draw_blue, rew, del_rew, abs(del_rew), choice, p_winA, est_prob, Q, Q(:,1)-Q(:,2), p_choice, choice-1];
        out.plot = [out.plot; out.run];

    end

    save_str = ['Fitted_results_ID' num2str(ID(i_p)) '.mat'];
    save(save_str,'out')
    
end