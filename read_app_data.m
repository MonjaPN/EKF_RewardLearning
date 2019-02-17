function D = read_app_data(ID,data_path)


T = readtable(data_path,'ReadVariableNames',1);

T = T(T.run_finished == 1,:);
tbl = tabulate(T.subject_id_app);


%Define relevant columns

 choice = ismember(T.choice,'B');
 draw_blue = ismember(T.drawn_outcome,'A');
 reward_grid = [T.reward_a, T.reward_b]; 
 
 for i=1:length(choice)
        reward(i,1) = reward_grid(i,choice(i)+1) * (2 * (T.win(i) - 0.5));
 end

  Dcomplete = [choice+1, (reward/50), draw_blue, reward_grid,T.subject_id_app, T.run_ind, T.trial_ind,T.p_reward_a]; %concatenate data vector
  
  if ~isempty(ID)
  DIDs = Dcomplete(ismember(T.subject_id_app,ID),:);
  IDTab = T(ismember(T.subject_id_app,ID),:);
  else 
      DIDs = Dcomplete
      IDTab = T
      
  end

    for i=1:length(IDTab.trial_ind)

        if IDTab.trial_ind(i) == 1 && IDTab.run_finished(i)== 0

            if IDTab.run_finished(i+1) == 1
                IDTab.run_finished(i) = 1;
            end
        end   
    end
    
    D = DIDs(IDTab.run_finished==1,:);
      
   
end