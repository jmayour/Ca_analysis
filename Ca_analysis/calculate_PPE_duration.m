function PPE_duration = calculate_PPE_duration(t,Ca,pk,pk_loc,pk_type,threshold,fps,p,i,maxcells_per_plot,define_PPE)

normal_loc = pk_loc(pk_type == 0);
abnormal_loc = pk_loc(pk_type ~= 0);

PPE_duration = 0;
after_pacing_duration = 0;
during_pacing_duration = 0;
if isempty(abnormal_loc) == 0
    PPE_loc = pk_loc(pk_loc>normal_loc(end));
    if isempty(PPE_loc) == 0
        after_pacing_duration = calculate_after_pacing_duration(normal_loc,t,Ca,threshold,fps,p,i,maxcells_per_plot);
    end
    if define_PPE == 1
        trigger_loc_during_pacing = abnormal_loc(abnormal_loc<normal_loc(end));
        if isempty(trigger_loc_during_pacing) == 0
            for ii = 1:(length(normal_loc)-1)
                during_pacing_duration_add = calculate_during_pacing_duration(trigger_loc_during_pacing,normal_loc,t,Ca,threshold,fps,p,i,maxcells_per_plot,ii);
                during_pacing_duration = during_pacing_duration + during_pacing_duration_add;
            end
        end
    end
    
end

PPE_duration = during_pacing_duration + after_pacing_duration;

end
