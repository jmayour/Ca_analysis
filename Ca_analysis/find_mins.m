function [fitting_time fitting_mins] = find_mins(time,Ca,fps,pacingfreq,up_filter,down_filter)

interval = 1/pacingfreq;

n_intervals = time(end)/interval;

for i = 1:n_intervals
    [fitting_mins(i) fitting_loc(i)] = min(Ca(time>i-1 & time<i));
    fitting_time(i) = time(fitting_loc(i)+(i-1)*fps);
end

if median(diff(fitting_mins))>0
    [filt_min1 filt_time1] = filt_mins(fitting_time, fitting_mins, time, Ca, 5*up_filter);
    [filt_min2 filt_time2] = filt_mins(flip(max(time)-filt_time1), flip(filt_min1), time, Ca, 0.5*down_filter);
    fitting_mins = flip(filt_min2);
    fitting_time = flip(max(time)-filt_time2);
else
    [filt_min1 filt_time1] = filt_mins(fitting_time, fitting_mins, time, Ca, 0.5*down_filter);
    [filt_min2 filt_time2] = filt_mins(flip(max(time)-filt_time1), flip(filt_min1), time, Ca, 5*up_filter);
    fitting_mins = flip(filt_min2);
    fitting_time = flip(max(time)-filt_time2);
    
end

end

