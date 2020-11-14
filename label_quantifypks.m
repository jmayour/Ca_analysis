function [pk_type] = label_quantifypks(pk,pk_loc,time,paced_pk_time,n_paces, pacingfreq,fps,figurenumber,p,i,maxcells_per_plot)

pk_type = ones(length(pk),1).*-1;
pk_type(1) = 0;
text_loc = pk;
JM_subplot_inloop_text(time(pk_loc(1)),text_loc(1),'N',figurenumber,p,i,maxcells_per_plot,'Quantifying PPEs');

for ii = 2:length(paced_pk_time)
    % 0 = normal
    % 1 = abnormal
    [M, I] = min(abs(time(pk_loc)-paced_pk_time(ii)));
    pk_type(I) = 0;
    JM_subplot_inloop_text(time(pk_loc(I)),text_loc(I),'N',figurenumber,p,i,maxcells_per_plot,'Quantifying PPEs');
end

pk_type(pk_type == -1) = 1;
JM_subplot_inloop_text(time(pk_loc(pk_type == 1)),text_loc(pk_type == 1),'A',figurenumber,p,i,maxcells_per_plot,'Quantifying PPEs');

end

