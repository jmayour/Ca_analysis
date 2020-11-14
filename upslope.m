function [slope_pk_final,slope_pk_loc_final] = upslope(t,Ca_upstroke,fps,pacingfreq,pk_loc,n_paces,p,i,maxcells_per_plot)

Ca_upstroke(1) = 0;

slope_init = diff(Ca_upstroke);
slope_Ca = zeros(length(Ca_upstroke),1);
slope_Ca(2:end) = slope_init;
[slope_pk,slope_pk_loc] = findpeaks(slope_Ca,'MinPeakHeight',max(slope_Ca)/8);

pk_loc_real = pk_loc;

for iii = 1:length(pk_loc_real)
    slope_init = slope_pk(t(slope_pk_loc)<=t(pk_loc_real(iii)));
    slope_pk_final(iii) = slope_init(end);
    slope_loc_init = slope_pk_loc(t(slope_pk_loc)<=t(pk_loc_real(iii)));
    slope_pk_loc_final(iii) = slope_loc_init(end);
end

slope_pk_final = slope_pk_final.*fps;

for ii = 1:length(slope_pk_loc_final)
    x = (t(slope_pk_loc_final(ii)-1)):1/fps:t(slope_pk_loc_final(ii)+1);
    y = slope_pk_final(ii).*(x-t(slope_pk_loc_final(ii)-1))+Ca_upstroke(slope_pk_loc_final(ii)-1);
    JM_subplot_inloop(x,y,'Time (s)','Paceable, F/F0',fignum(4,i,maxcells_per_plot),p,i,0,maxcells_per_plot,'Calculating Calcium Transient Upstroke Parameters');
    hold on;
    
end

end