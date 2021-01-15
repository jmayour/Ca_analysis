function slope_pk = downslope(t,Ca_upstroke,fps,pacingfreq,pk_loc,n_paces,p,i,nontriggered_index,maxcells_per_plot)

Ca_upstroke(1) = 0;

slope_init = diff(Ca_upstroke);
slope_Ca = zeros(length(Ca_upstroke),1);
slope_Ca(2:end) = slope_init;
[preslope_pk,preslope_pk_loc] = findpeaks(slope_Ca,'MinPeakHeight',max(slope_Ca)/8);

pk_loc_real = pk_loc;

for iii = 1:length(pk_loc_real)
    preslope_init = preslope_pk(t(preslope_pk_loc)>=t(pk_loc_real(iii)));
    preslope_pk(iii) = preslope_init(1);
    preslope_loc_init = preslope_pk_loc(t(preslope_pk_loc)>=t(pk_loc_real(iii)));
    preslope_pk_loc(iii) = preslope_loc_init(1);
end

preslope_pk = preslope_pk.*fps;
slope_pk = preslope_pk(nontriggered_index);
slope_pk_loc = preslope_pk_loc(nontriggered_index);

for ii = 1:length(slope_pk_loc)
    x = (t(slope_pk_loc(ii)-1)):1/fps:t(slope_pk_loc(ii)+1);
    y = -slope_pk(ii).*(x-t(slope_pk_loc(ii)-1))-Ca_upstroke(slope_pk_loc(ii)-1);
    JM_subplot_inloop(x,y,'Time (s)','Paceable, F/F0',fignum(5,i,maxcells_per_plot),p,i,0,maxcells_per_plot,'Calculating Calcium Transient Downstroke Parameters');
    hold on;

end

end