function during_pacing_duration_add = calculate_during_pacing_duration(trigger_loc_during_pacing,normal_loc,t,Ca,threshold,fps,p,i,maxcells_per_plot,ii)

pacing_PPE_loc = trigger_loc_during_pacing(trigger_loc_during_pacing>normal_loc(ii) & trigger_loc_during_pacing <normal_loc(ii+1));

if isempty(pacing_PPE_loc) == 0
    
    t_pacing_PPE_init = t(normal_loc(ii):normal_loc(ii+1));
    Ca_pacing_PPE_init = Ca(normal_loc(ii):normal_loc(ii+1));
    
    slope_Ca = diff(Ca_pacing_PPE_init);
    slope_Ca_flip = diff(flip(Ca_pacing_PPE_init));

    loc_of_interest_start = find(slope_Ca>0);
    loc_of_interest_end = find(flip(slope_Ca_flip)>0);
    
    t_pacing_PPE = t_pacing_PPE_init(loc_of_interest_start(1):loc_of_interest_end(end));
    Ca_pacing_PPE = Ca_pacing_PPE_init(loc_of_interest_start(1):loc_of_interest_end(end));
    
    t_pacing_threshold = find(Ca_pacing_PPE>=threshold.*max(Ca));
    during_pacing_duration_add = length(t_pacing_threshold)./fps;
    
    JM_subplot_inloop(t_pacing_PPE(t_pacing_threshold), zeros(length(t_pacing_threshold),1), 'Time (s)','Paceable, Filtered Ca',...
        fignum(6,i,maxcells_per_plot),p,i,3,maxcells_per_plot,'Quantifying PPEs');
else

during_pacing_duration_add = 0;

end

end
