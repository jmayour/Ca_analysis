function [ca_pk,ca_pk_loc] = max_Ca(t,Ca_paced,pk,pk_loc,pk_type,pk_threshold,p,i,maxcells_per_plot)%(t, Ca_upstroke,fps,pacingfreq,p,i,maxcells_per_plot)

[allpk,allpk_loc] = findpeaks(Ca_paced,'MinPeakHeight',pk_threshold.*max(Ca_paced),...
    'MinPeakProminence',pk_threshold*max(Ca_paced));

normal_locs = pk_loc(pk_type==0);

for iii = 1:length(normal_locs)
    [M, I] = min(abs(normal_locs(iii)-allpk_loc));
    ca_pk_loc(iii) = allpk_loc(I);
end

ca_pk = Ca_paced(ca_pk_loc);
JM_subplot_inloop(t(ca_pk_loc), ca_pk,'Time (s)','F/F0',fignum(4,i,maxcells_per_plot),p,i,1,maxcells_per_plot,'Calculating Calcium Transient Upstroke Parameters');

end
