function [paceable, paced_pk, paced_pk_time] = check_paceable(time,Ca,n_paces,n_paces_needed,fps,pacingfreq,cell_type)
Ca(1) = 0;
paceable = 0;

max_both = [max(Ca)/4,1e-2];
[pk,pk_time] = findpeaks(Ca,fps,'MinPeakHeight',max(max_both));

if cell_type == 1
    rise_filt = 1e-1;
elseif cell_type == 2
    rise_filt = 4e-1;
end

paced_pk_locs = mod(pk_time+2.5e-2,1/pacingfreq)<=rise_filt;
paced_pk = pk(paced_pk_locs);
paced_pk_time = pk_time(paced_pk_locs);

if isempty(paced_pk) == 0
    rounded_time = round(paced_pk_time);
    for ii = 1:length(rounded_time)
        same_peak(ii) = length(find(rounded_time(ii)==rounded_time));
    end
    paced_pk(same_peak>1) = [];
    paced_pk_time(same_peak>1) = [];
end

tau_achieved = 0;

if length(paced_pk) <= n_paces && length(paced_pk) >= n_paces_needed
    
    [min_val,min_time] = findpeaks(-Ca,fps);
    min_val = -min_val;
    
    for i = 2:length(paced_pk)
        
        check_min = min_time<paced_pk_time(i);
        min_val_time = min_time(check_min);
        min_val_val = min_val(check_min);
        
        if min_val_val(end)<=paced_pk(i-1)/4
            time_shift = time + 0.1;
            Ca_shift = Ca(time_shift>=min_val_time(end) & time<=min_val_time(end));
            
            if all(Ca_shift<=(paced_pk(i-1)/2)) == 1
                tau_achieved = tau_achieved + 1;
            end
            
        end
        
    end
    
    if tau_achieved >= (n_paces_needed - 1) && tau_achieved <= (n_paces - 1)
        paceable = 1;
    end
    
end

end
