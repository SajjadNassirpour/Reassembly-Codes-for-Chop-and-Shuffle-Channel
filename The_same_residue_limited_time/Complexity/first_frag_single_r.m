% This function finds candidate combinations including one fragment


function [frag,complex_i] = first_frag_single_r(m,L,d_sec,breaks_points,X_seperated,VT_section_length,parity,Main_r)
    complex_i=0;
    k=0;
    frag=[];
    for i=1:size(X_seperated,1)
        complex_i=complex_i+1;
        stop=0;
        La=L+1;
        r_mat=-1*ones(1,L);
        r=-1;
        sum_length=X_seperated(i,breaks_points+1);
        while La>1 && stop==0
            La=La-1;
            window_size=VT_section_length(La);
            start_p=start_point_1(m,L,La,d_sec,parity);
            if sum(X_seperated(i,breaks_points+1)+1>=start_p)>1
                valid_point_idx=sum(X_seperated(i,breaks_points+1)+1>=start_p);
                test_segment2_total=[X_seperated(i,breaks_points+2:breaks_points+1+X_seperated(i,breaks_points+1))];          
                if X_seperated(i,breaks_points+1)-start_p(valid_point_idx)+1>=window_size
                    valid_point=start_p(valid_point_idx)-1+window_size;
                    valid_point_idx=valid_point_idx+1;
                else

                    valid_point=start_p(valid_point_idx)-1;
                end
                test_segment2_total=test_segment2_total(1:valid_point);
                k0=0;
                while valid_point_idx>1 && stop==0
                    test_segment=test_segment2_total(1:window_size);
                    index=1:length(test_segment);
                    VT_sum=dot(test_segment,index);
                    r=mod(VT_sum,(length(test_segment)+1));

                    data_size=length(test_segment)-parity(La);
                    index1=1:data_size;
                    VT_sum1=dot(test_segment(1:data_size),index1);
                    r_prime=mod(VT_sum1,(length(test_segment)+1));

                    delta=dot(test_segment(data_size+1:end),parity(La):-1:1);


                    if r_prime>=r
                        delta_true=r_prime-r;
                    else
                        delta_true=length(test_segment)+1+r_prime-r;
                    end
                    if r==Main_r && delta_true==delta
                        k0=k0+1;
                        valid_point_idx=valid_point_idx-1;
                        if valid_point_idx>1

                            test_segment2_total(1:start_p(k0+1)-start_p(k0))=[];
                        end
                        r_mat(1,La)=r;
                    else
                        stop=1;
                    end
                end
                
            end
        end
        if stop==0
            k=k+1;
            if r==-1
               frag(k,:)=[X_seperated(i,1) r_mat sum_length   1];
            else
               frag(k,:)=[X_seperated(i,1) r_mat sum_length  0];
            end
        end
    end
end