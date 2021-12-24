% This function finds candidate combinations including more than one fragment

function frag = other_frag_single_r(m,L,d_sec,breaks_points,data_combined,X_seperated,VT_section_length,parity,fragment_combined,r_mat_orig,length_combined,rep_frag,Main_r,limited_time)

    time_lower=tic;
    k=0;
    frag=[];
    for i=1:size(X_seperated,1)
        time_upper=toc(time_lower);
        
        if time_upper<limited_time
            stop=0;
            La=L+1;

            r=-1;
            r_mat=r_mat_orig;
            while La>1 && stop==0 && time_upper<limited_time
                
                La=La-1;
                window_size=VT_section_length(La);
                start_p=start_point_1(m,L,La,d_sec,parity);

                residue_combined=r_mat(1,La);

                valid_point_idx=sum(length_combined+X_seperated(i,breaks_points+1)+1>=start_p);

                if length_combined+X_seperated(i,breaks_points+1)-start_p(valid_point_idx)+1>=window_size
                    valid_point=start_p(valid_point_idx)-1+window_size;
                    valid_point_idx=valid_point_idx+1;
                else

                    valid_point=start_p(valid_point_idx)-1;
                end

                if valid_point_idx>residue_combined+2

                    test_segment2_total=[data_combined   X_seperated(i,breaks_points+2:breaks_points+1+X_seperated(i,breaks_points+1))];
                    sum_length=length(test_segment2_total);

                    test_segment2_total=test_segment2_total(start_p(residue_combined+2):valid_point);
                    stop=0;
                    k0=residue_combined+1;

                    while valid_point_idx>residue_combined+2 && stop==0 && time_upper<limited_time
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
                            if valid_point_idx>residue_combined+2

                                test_segment2_total(1:start_p(k0+1)-start_p(k0))=[];
                            end
                            r_mat(1,La)=r;

                        else
                            stop=1;
                        end
                        time_upper=toc(time_lower);
                    end

                end
                time_upper=toc(time_lower);
            end
            if stop==0
                k=k+1;
                if r==-1
                   sum_length=length_combined+X_seperated(i,breaks_points+1);
                   frag(k,:)=[ fragment_combined   X_seperated(i,1) r_mat sum_length   rep_frag+1];
                else
                   frag(k,:)=[ fragment_combined   X_seperated(i,1) r_mat sum_length  0];
                end
            end
        end

    end
end