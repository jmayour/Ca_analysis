function nontriggered = check_triggered(t,pk,pk_loc,pk_all_loc,pk_all_type,n_paces,tau90)

normal_loc = pk_loc;
normal_times = t(normal_loc);
abnormal_loc = pk_all_loc(pk_all_type~=0);
abnormal_times = t(abnormal_loc);
empty = isempty(abnormal_loc);

q = 0;
nontriggered = [];
if empty == 0
    
    for i = 1:length(normal_loc)
        abnormal_pk_beat = abnormal_times(abnormal_times>=normal_times(i) ...
            & abnormal_times<= tau90(i));
        check_triggered = isempty(abnormal_pk_beat);
        if check_triggered == 1
            q = q+1;
            nontriggered(q) = i;
        end
    end
    
else
    nontriggered = 1:n_paces;
end


end