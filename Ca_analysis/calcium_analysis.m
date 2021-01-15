function calcium_analysis(tablename,savename,fps,pacingfreq,n_paces,n_paces_needed,pk_threshold,...
    PassbandFrequency,StopbandFrequency,PassbandRipple,StopbandAttenuation,maxcells_per_plot,up_filter,down_filter,save_which,...
    cell_type,prepace_time)

clc;

d = designfilt('lowpassfir','PassbandFrequency',PassbandFrequency, ...
    'StopbandFrequency',StopbandFrequency,'PassbandRipple',PassbandRipple, ...
    'StopbandAttenuation',StopbandAttenuation);

paceable_index = 0;
paceable_cell = [];
unpaceable_index = 0;
unpaceable_cell = [];

Ca_table = readtable(tablename);
Ca_init = table2cell(Ca_table);
loop_length = size(Ca_init,2);

p = numSubplots(min([maxcells_per_plot, loop_length]));

for i = 1:loop_length
    
    try
        clearvars Ca_inter Ca_initial time_i Ca time
        
        Ca_inter = cell2mat(Ca_init(:,i));
        Ca_initial = Ca_inter(~isnan(Ca_inter));
        time_i = (1:1:length(Ca_initial))./fps;
        time_i = time_i - prepace_time;
        time = time_i(time_i>=0);
        Ca = Ca_initial(time_i>=0);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % STEP 1: Normalization: COMPLETE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Identify baseline points to nonlinear fit
        [fitting_time, fitting_mins] = find_mins(time,Ca,fps,pacingfreq,up_filter,down_filter);
        
        % Spline fit
        norm_spline = interp1(fitting_time,fitting_mins,time,'linear','extrap');
        
        % Plot Raw data
        JM_subplot_inloop(time,Ca,'Time (s)','Raw Ca',fignum(1,i,maxcells_per_plot),p,i,0,maxcells_per_plot,'Identifying Photobleaching Effect');
        
        % Plot minimums used in spline
        JM_subplot_inloop(fitting_time,fitting_mins,'Time (s)','Raw Ca',fignum(1,i,maxcells_per_plot),p,i,1,maxcells_per_plot,'Identifying Photobleaching Effect');
        JM_subplot_inloop(time,norm_spline,'Time (s)','Raw Ca',fignum(1,i,maxcells_per_plot),p,i,0,maxcells_per_plot,'Identifying Photobleaching Effect');
        
        % Normalization
        norm_Ca = Ca./norm_spline'-1;
        
        % Subtract by Baseline Noise Value
        noise = median(norm_Ca(norm_Ca<=pk_threshold.*max(norm_Ca)));
        norm_Ca = norm_Ca-noise;
        
        % Plot normalized data
        JM_subplot_inloop(time,norm_Ca,'Time (s)','Normalized Ca',fignum(2,i,maxcells_per_plot),p,i,0,maxcells_per_plot,'Normalized Data After Removing Photobleaching Effect');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % STEP 2: Filter: COMPLETE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        JM_subplot_inloop(time, norm_Ca, 'Time (s)','Paceable, Filtered Ca',fignum(3,i,maxcells_per_plot),p,i,0,maxcells_per_plot,'Paceable or Unpaceable?');
        hold on;
        filtered_Ca = filtfilt(d,norm_Ca);
        filtered_Ca(1) = 0;
        JM_subplot_inloop(time, filtered_Ca, 'Time (s)','Paceable, Filtered Ca',fignum(3,i,maxcells_per_plot),p,i,0,maxcells_per_plot,'Paceable or Unpaceable?');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % STEP 3: Check for paceable cells
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [binary_paceable(i), paced_pk, paced_pk_time] = check_paceable(time(time<=(n_paces-1/2)/pacingfreq)...
            ,filtered_Ca(time<=(n_paces-1/2)/pacingfreq),n_paces,n_paces_needed,fps,pacingfreq,cell_type);
        
        if binary_paceable(i) == 1
            paceable_index = paceable_index + 1;
            paceable_cell(paceable_index) = i;
            figure(fignum(3,i,maxcells_per_plot));
            subplot(p(1),p(2),1+mod(i-1,maxcells_per_plot));
            title('PACEABLE','Color','green')
        else
            unpaceable_index = unpaceable_index + 1;
            unpaceable_cell(unpaceable_index) = i;
            figure(fignum(3,i,maxcells_per_plot));
            subplot(p(1),p(2),1+mod(i-1,maxcells_per_plot));
            title('UNPACEABLE','Color','red')
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % STEP 4: Quantify Number and Duration of Post-Pacing Events: COMPLETE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if ismember(i,paceable_cell) == 1
            JM_subplot_inloop(time, filtered_Ca, 'Time (s)','Paceable, Filtered Ca',fignum(6,i,maxcells_per_plot),p,i,0,maxcells_per_plot,'Quantifying PPEs');
            [pk,pk_loc] = findpeaks(filtered_Ca,'MinPeakHeight',pk_threshold.*max(filtered_Ca),...
                'MinPeakProminence',pk_threshold*max(filtered_Ca));
            JM_subplot_inloop(time(pk_loc),pk,'Time (s)','Paceable, Filtered Ca',fignum(6,i,maxcells_per_plot),p,i,1,maxcells_per_plot,'Quantifying PPEs');
            [pk_type] = label_quantifypks(pk,pk_loc,time,paced_pk_time,n_paces,pacingfreq,fps,fignum(6,i,maxcells_per_plot),p,i,maxcells_per_plot);
            P = [pk,pk_loc, pk_type];
            normal_peak_locs = pk_loc(pk_type == 0);
            
            PPE_number_all(i) = length(pk_loc(pk_type ~= 0));
            PPE_number_post(i) = length(pk(time(pk_loc)>time(normal_peak_locs(end))));
            PPE_number_during(i) = PPE_number_all(i)-PPE_number_post(i);
            
            PPE_duration_all(i) = calculate_PPE_duration(time,filtered_Ca,pk,pk_loc,pk_type,pk_threshold,fps,p,i,maxcells_per_plot,1);
            PPE_duration_post(i) = calculate_PPE_duration(time,filtered_Ca,pk,pk_loc,pk_type,pk_threshold,fps,p,i,maxcells_per_plot,2);
            PPE_duration_during(i) = PPE_duration_all(i)-PPE_duration_post(i);

        
        else
            JM_subplot_inloop(NaN, NaN, 'Time (s)','Unpaceable',fignum(6,i,maxcells_per_plot),p,i,0,maxcells_per_plot,'Quantifying PPEs');
            P = [NaN,NaN, NaN];
            
            PPE_number_all(i) = NaN;
            PPE_number_post(i) = NaN;
            PPE_number_during(i) = NaN;
            
            PPE_duration_all(i) = NaN;
            PPE_duration_post(i) = NaN;
            PPE_duration_during(i) = NaN;

        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % STEP 5: Quantify Calcium Transient Parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Focus only on paced portion of data
        index_Ca = time<=n_paces/pacingfreq;
        t = time(index_Ca');
        Ca_paced = norm_Ca(index_Ca');
        
        if ismember(i,paceable_cell) == 1
            
            %%% amplitude
            JM_subplot_inloop(t, Ca_paced,'Time (s)','F/F0',fignum(4,i,maxcells_per_plot),p,i,0,maxcells_per_plot,'Calculating Calcium Transient Upstroke Parameters');
            [ca_pk,ca_pk_loc] = max_Ca(t,Ca_paced,pk,pk_loc,pk_type,pk_threshold,p,i,maxcells_per_plot);
            Amp_mean(i) = mean(ca_pk); % mean amplitude
            Amp_std(i) = std(ca_pk); % std amplitude
            
            %%% upstroke velocity
            [upslope_pk,upslope_pk_loc] = upslope(t,Ca_paced,fps,pacingfreq,ca_pk_loc,n_paces,p,i,maxcells_per_plot);
            UV_mean(i) = mean(upslope_pk); % mean upstroke velocity
            UV_std(i) = std(upslope_pk); % std upstroke velocity
            
            %%% T50, T90, C50
            JM_subplot_inloop(t, Ca_paced,'Time (s)','F/F0',fignum(5,i,maxcells_per_plot),p,i,0,maxcells_per_plot,'Calculating Calcium Transient Downstroke Parameters');
            [tau90_loc, preT90] = time_function(0.1,t,Ca_paced,ca_pk,ca_pk_loc,p,i);
            [tau50_loc, preT50] = time_function(0.5,t,Ca_paced,ca_pk,ca_pk_loc,p,i);
            [C50_loc, preC50] = C50_time_function(0.5,t,Ca_paced,ca_pk,ca_pk_loc,p,i);
            
            %%% Identify triggers during pacing
            nontriggered_index = check_triggered(time,ca_pk,ca_pk_loc,pk_loc, pk_type,n_paces,tau90_loc);
            percent_nontriggered(i) = length(nontriggered_index)/n_paces*100;
            
            %%% Only include T50, T90, C50 for nontriggered
            T90 = preT90(nontriggered_index);
            T50 = preT50(nontriggered_index);
            C50 = preC50(nontriggered_index) + T50;
            JM_subplot_inloop(tau90_loc(nontriggered_index), 0.1.*ca_pk(nontriggered_index),'Time (s)','F/F0',fignum(5,i,maxcells_per_plot),p,i,2,maxcells_per_plot,'Calculating Calcium Transient Downstroke Parameters');
            JM_subplot_inloop(tau50_loc(nontriggered_index), 0.5.*ca_pk(nontriggered_index),'Time (s)','F/F0',fignum(5,i,maxcells_per_plot),p,i,2,maxcells_per_plot,'Calculating Calcium Transient Downstroke Parameters');
            C50_x = [C50_loc(nontriggered_index)', tau50_loc(nontriggered_index)'];
            C50_y = [0.5.*ca_pk(nontriggered_index), 0.5.*ca_pk(nontriggered_index)];
            
            for C50_index = 1:length(nontriggered_index)
                JM_subplot_inloop(C50_x(C50_index,:), C50_y(C50_index,:),'Time (s)','F/F0',fignum(5,i,maxcells_per_plot),p,i,0,maxcells_per_plot,'Calculating Calcium Transient Downstroke Parameters');
            end
            
            %%% nontriggered downstroke velocity
            downslope_pk = downslope(t,-Ca_paced,fps,pacingfreq,ca_pk_loc,n_paces,p,i,nontriggered_index,maxcells_per_plot);
            
            T50_mean(i) = mean(T50); % mean T50
            T90_mean(i) = mean(T90); % mean T90
            C50_mean(i) = mean(C50); % mean C50
            T50_std(i) = std(T50); % std T50
            T90_std(i) = std(T90); % std T90
            C50_std(i) = std(C50); % std C50
            DV_mean(i) = mean(downslope_pk); % mean downstroke velocity
            DV_std(i) = std(downslope_pk); % std downstroke velocity
            
        else
            JM_subplot_inloop(NaN, NaN,'Time (s)','Unpaceable',fignum(4,i,maxcells_per_plot),p,i,0,maxcells_per_plot,'Calculating Calcium Transient Upstroke Parameters');
            JM_subplot_inloop(NaN, NaN,'Time (s)','Unpaceable',fignum(5,i,maxcells_per_plot),p,i,0,maxcells_per_plot,'Calculating Calcium Transient Downstroke Parameters');
            
            Amp_mean(i) = NaN;
            Amp_std(i) = NaN;
            UV_mean(i) = NaN;
            UV_std(i) = NaN;
            DV_mean(i) = NaN;
            DV_std(i) = NaN;
            T50_mean(i) = NaN;
            T90_mean(i) = NaN;
            C50_mean(i) = NaN;
            T50_std(i) = NaN;
            T90_std(i) = NaN;
            C50_std(i) = NaN;
            percent_nontriggered(i) = NaN;
            
        end
        
        
    catch
        warning(['Analysis of Cell',num2str(i),' skipped due to error']);
        
        JM_subplot_inloop(NaN, NaN,'Time (s)','Unpaceable',fignum(4,i,maxcells_per_plot),p,i,0,maxcells_per_plot,'Calculating Calcium Transient Upstroke Parameters');
        JM_subplot_inloop(NaN, NaN,'Time (s)','Unpaceable',fignum(5,i,maxcells_per_plot),p,i,0,maxcells_per_plot,'Calculating Calcium Transient Downstroke Parameters');
        
        Amp_mean(i) = NaN;
        Amp_std(i) = NaN;
        UV_mean(i) = NaN;
        UV_std(i) = NaN;
        DV_mean(i) = NaN;
        DV_std(i) = NaN;
        T50_mean(i) = NaN;
        T90_mean(i) = NaN;
        C50_mean(i) = NaN;
        T50_std(i) = NaN;
        T90_std(i) = NaN;
        C50_std(i) = NaN;
        percent_nontriggered(i) = NaN;
        binary_paceable(i) = NaN;
    end
    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 6: Output to Excel File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varNames = {'Cell_Number','Paceable_1yes_0no','Amp_mean_F_F0','Amp_std_F_F0','UV_mean_F_F0_s','UV_std_F_F0_s',...
    'DV_mean_F_F0_s','DV_std_F_F0_s','T50_mean_s','T50_std_s','T90_mean_s','T90_std_s','C50_mean_s','C50_std_s'...
    'PPE_number_during','PPE_number_post','PPE_number_all','PPE_duration_during_s','PPE_duration_post_s','PPE_duration_all_s',...
    'Percent_Paces_Nontriggered'};

relabel_table = readtable(tablename,'ReadVariableNames',0);

for names = 1:loop_length
    cell_name(names) = relabel_table{1,names};
end

T = table(cell_name',binary_paceable',Amp_mean',Amp_std',...
    UV_mean',UV_std',DV_mean',DV_std',...
    T50_mean',T50_std',T90_mean',T90_std',...
    C50_mean',C50_std',...
    PPE_number_during',PPE_number_post',PPE_number_all',...
    PPE_duration_during',PPE_duration_post',PPE_duration_all',percent_nontriggered',...
    'VariableNames',varNames);

if save_which == 1
    T2 =  T;
elseif save_which == 2
    T2 = rmmissing(T);
end

writetable(T2,[savename,'_Results.xlsx'],'WriteVariableNames',1)

varNames_real = {'Cell_Name','Paceable? (1=yes, 0=no)','Mean Amplitude (F/F0)','Std Amplitude (F/F0)',...
    'Mean Upstroke Velocity (F/F0/s)','Std Upstroke Velocity (F/F0/s)',...
    'Mean Downstroke Velocity (F/F0/s)','Std Downstroke Velocity (F/F0/s)',...
    'Mean T50 (s)','Std T50 (s)','Mean T90 (s)','Std T90 (s)','Mean C50 (s)','Std C50 (s)'...
    'PPE Number (during pacing)','PPE Number (post-pacing)','PPE Number (total)',...
    'PPE Duration (during pacing) (s)','PPE Duration (post-pacing) (s)','PPE Duration (total) (s)',...
    'Percent Paces Nontriggered (%)'};

T3 = table(varNames_real);
writetable(T3,[savename,'_Results.xlsx'],'WriteVariableNames',0)

FigList = findobj(allchild(0), 'flat', 'Type', 'figure'); savefig(FigList,'allfigures.fig');

end





