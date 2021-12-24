% This function takes theoverlap between the output sequences

function recover_rate=data_recovery_rate(m,L,d_sec,parity,breaks_points,X_seperated1,data_bits_total,frag_final)


compare_string=data_bits_total;

start_p=start_point_1(m,L,1,d_sec,parity);


frag_revover_data=[];
for L1=1:size(frag_final,1)
    data_combined=[];
    for jj=1:size(frag_final,2)-(L+2)
        data_combined=[data_combined  X_seperated1(X_seperated1(:,1)==frag_final(L1,jj),breaks_points+2:breaks_points+1+X_seperated1(X_seperated1(:,1)==frag_final(L1,jj),breaks_points+1))];
    end
    recov_data=[];
    for L2=1:sum(min(frag_final(:,end-1))+1>=start_p)-1
        recov_data=[recov_data  data_combined(start_p(L2):start_p(L2)+d_sec-1)];
    end
    L2=L2+1;
    recov_data=[recov_data  data_combined(start_p(L2):start_p(L2)+min(d_sec,min(frag_final(:,end-1))+1-start_p(L2))-1)];
    
    frag_revover_data(L1,:)=recov_data;
    
    if length(frag_revover_data(L1,:))<length(compare_string)
        
        data_bits_total=data_bits_total(1:length(frag_revover_data(L1,:)));
        compare_string=compare_string(1:length(frag_revover_data(L1,:)));
    end
   
    if sum(frag_revover_data(L1,:)==data_bits_total)<m^(L-1)*d_sec
        
        compare_string(frag_revover_data(L1,:)~=data_bits_total)=-1;   
    end
end

recover_rate=(length(compare_string)-sum(compare_string==-1))/(m^(L-1)*d_sec);

end

