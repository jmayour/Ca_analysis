function [filt_mins filt_time] = filt_mins(fitting_time, fitting_mins, time, Ca, percent)

q = 0;

for i = 1:length(fitting_mins)-1
    index = i-q;
    difference = (fitting_mins(index+1)-fitting_mins(index))./...
        (fitting_time(index+1)-fitting_time(index));
    if difference > abs((max(Ca)-min(Ca))/(max(time)-min(time))*percent)
        fitting_mins(index+1) = [];
        fitting_time(index+1) = [];
        q = q+1;
    end
end

filt_mins = fitting_mins;
filt_time = fitting_time;

end
