function [tau, timeconstant] = time_function(percent,t,Ca_upstroke,ca_pk,ca_pk_loc,p,i)

peak_t = t(ca_pk_loc);

for ii = 1:length(peak_t)
    timeconstant_init = find((t' >= peak_t(ii)) & (Ca_upstroke<=(percent.*ca_pk(ii))));
    if isempty(timeconstant_init) == 0
        timeconstant_index_1 = t(timeconstant_init(1));
        timeconstant_index_2 = t(timeconstant_init(1)-1);
        timeconstant_val_1 = Ca_upstroke(timeconstant_init(1));
        timeconstant_val_2 = Ca_upstroke(timeconstant_init(1)-1);
        timeconstant_val = [timeconstant_val_1; timeconstant_val_2];
        timeconstant_index = [timeconstant_index_1; timeconstant_index_2];
        
        tau(ii) = interp1(timeconstant_val,timeconstant_index,percent.*ca_pk(ii));

        timeconstant(ii) = tau(ii) - t(ca_pk_loc(ii));
    else
        timeconstant(ii) = 1;
        tau(ii) = ii.*1;
    end
end

end
