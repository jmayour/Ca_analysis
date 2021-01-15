function after_pacing_duration = calculate_after_pacing_duration(normal_loc,t,Ca,threshold,fps,p,i,maxcells_per_plot)

t_PPE_init = t(normal_loc(end):end);
Ca_PPE_init = Ca(normal_loc(end):end);

slope_Ca = diff(Ca_PPE_init);

loc_start = find(slope_Ca>0);

t_PPE = t_PPE_init(loc_start(1):end);
Ca_PPE = Ca_PPE_init(loc_start(1):end);

t_threshold = find(Ca_PPE>=threshold.*max(Ca));
after_pacing_duration = length(t_threshold)./fps;

JM_subplot_inloop(t_PPE(t_threshold), zeros(length(t_threshold),1), 'Time (s)','Paceable, Filtered Ca',...
    fignum(6,i,maxcells_per_plot),p,i,3,maxcells_per_plot,'Quantifying PPEs');

end