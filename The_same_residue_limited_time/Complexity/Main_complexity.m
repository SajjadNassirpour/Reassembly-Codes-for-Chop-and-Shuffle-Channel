clc
clear 
close all

% ------ Parameters ------- %

d_sec=50; % Length of a data section in layer 1
L=3; % Number of encoding layer
B_p=0.5; % Bernoulli Prob
m=2;  
alpha=0.05;
limited_time=20; % \Delta
Rt=0.65;  % Desired rate
Main_r=0; % The residue of all VT codewords
iters=1000;   % number of iterations


% Capacity of the channel
Capacity=exp(-alpha)

% Number of parity bits of a VT codeword in each layer
parity=zeros(1,L);
% Length of a VT codeword in each layer
VT_section_length=zeros(1,L);

iter=0;
rep_iter=0;
while iter<iters
    
    % If the decoder cannot find \hat{d}, we repeat the simulation with a new
    % breaking channel
    if rep_iter==0
        iter=iter+1;
        iter
    end
    
    % Generate input data bits
    [data_bits,data_bits_total] = rand_gen_nested_any_m(d_sec,m,L,B_p);
   
    % ----- Encoder ----- %
    for La=1:L
        parity(1,La)=ceil((sqrt(8*size(data_bits,2)+1)+1)/2);
        X_section=zeros(m^(L-La),size(data_bits,2)+parity(1,La));
        for kk=1:m^(L-La)
            [X_section(kk,:),p]=sloane_VT_encoder(data_bits(kk,:),parity(1,La),Main_r);
        end
        VT_section_length(La)=size(X_section,2);
        if La<L
        data_bits=reshape(X_section',m*VT_section_length(La),m^(L-La)/m)';
        end
    end
    X=X_section;
    
    % Codeword length
    n=length(X)
    
    Rate=length(data_bits_total)/n
    
    


    % ----- Chop-and-shuffle channel ------ %
   
    if alpha==0
        cut_rule=log2(n);
    else
        cut_rule=alpha*n/log2(n);
    end
    
    ran_vec=break_points_nested(cut_rule,n,5*ceil(cut_rule));



    if length(ran_vec)~=1 || ran_vec(1)~=0
        
    % Assign labels to different fragments
    [X_seperated,X_sep_orginal]=Data_frag_nested(X,ran_vec,cut_rule);
    
    
    % Save number of fragments in each iteration (M)
    row_cal(iter)=size(X_seperated,1);
    
    
    % ------ Decoder ------ %

    % Location of the first bit of each codeword in each layer 
    start_p=start_point_1(m,L,1,d_sec,parity);


    breaks_points=size(X_seperated,1);
    X_seperated1=[X_seperated(:,1).*ones(size(X_seperated,1),1)   zeros(size(X_seperated,1),size(X_seperated,1)-1)   X_seperated(:,2).*ones(size(X_seperated,1),1)   X_seperated(:,3:end).*ones(size(X_seperated,1),size(X_seperated,2)-2)];
    length_avg=n/cut_rule;


    tic
    done=0;
    
    % limited memory
    en_frag=1;
    
    % How much time our algorithm uses to find \hat{d}
    passed_time=0;
    
    
    % complexity
    complexity=0;
    
    while done==0 && passed_time<limited_time && en_frag<row_cal(iter)
        
    en_frag=en_frag+1;
    VT_num_check=1;
    choice_comb=1;
    
    % Find candidate combinations including one fragment
    [frag1,complex_i] = first_frag_single_r(m,L,d_sec,breaks_points,X_seperated1,VT_section_length,parity,Main_r);
    complexity=complexity+complex_i;
    frag_final=frag1;

    mul=1;
    for ii0=1:breaks_points-1
    kk=0;
    frag1=[];
    if ii0==mul*en_frag
        if ~isempty(frag_final)
            mul=mul+1;
            max_r=max(sum(frag_final(:,end-L-1:end-2),2));
            frag_final=frag_final(sum(frag_final(:,end-L-1:end-2),2)==max_r,:);
        end
    end
    for ii=1:size(frag_final,1)

        frag=[];
        X_seperated_test=X_seperated1;
        rep_frag=frag_final(ii,end);
        length_combined=frag_final(ii,end-1);
        r_mat=frag_final(ii,end-L-1:end-2);
        fragment_combined=frag_final(ii,1:end-L-2);
        data_combined=[];

        for jj=1:size(frag_final,2)-(L+2)
            data_combined=[data_combined  X_seperated_test(X_seperated_test(:,1)==frag_final(ii,jj),breaks_points+2:breaks_points+1+X_seperated_test(X_seperated_test(:,1)==frag_final(ii,jj),breaks_points+1))];
            X_seperated_test(X_seperated_test(:,1)==frag_final(ii,jj),:)=[];
        end
        % Find candidate combinations including more than one fragment
        [frag,complex_i] = other_frag_single_r(m,L,d_sec,breaks_points,data_combined,X_seperated_test,VT_section_length,parity,fragment_combined,r_mat,length_combined,rep_frag,Main_r);
        complexity=complexity+complex_i;
        if ~isempty(frag)
            valid_frag=frag_final;
            for ll=1:size(frag,1)
            kk=kk+1;
            frag1(kk,1:size(frag,2))=frag(ll,:);
            end
        end
    end
    frag_final=frag1;
    end
   

        if isempty(frag_final)

        else
            % Take the overlap between the reassembled sequences
            recover_rate(iter)=data_recovery_rate(m,L,d_sec,parity,breaks_points,X_seperated1,data_bits_total,frag_final);
            Number_failed_cases_plus_errors=sum(recover_rate<=(Rt/Rate))
            Details_fails_and_errors=sort(recover_rate(recover_rate<=(Rt/Rate)))
            done=1;
        end
        passed_time=toc;
    end
    else
        recover_rate(iter)=1;
        row_cal(iter)=0;
    end

if length(ran_vec)~=1 || ran_vec(1)~=0
    if passed_time>=20 && done==0
        recover_rate(iter)=0;
        rep_iter=1;
    else
        rep_iter=0;
    end
else
   complexity=0; 
end


total_complexity(iter)=complexity;

MM=row_cal(iter);
% Calculate the complexity of our approach / Brute-Force in each iteration
frac_comp(iter)=complexity/factorial(MM);

end

% Average value of M
mean_number_breaks=mean(row_cal)

% Average value of M(theoretically)
theory_break_number=cut_rule


% Average rate of recovery after taking the overlap
Average_recover_rate=mean(recover_rate(recover_rate>0))

% Heuristic Rt
Final_Rate=min(recover_rate(recover_rate>0))*Rate

% Error rate
Error_rate=sum(Details_fails_and_errors~=0)/iters

% Complexity of our approach / Brute-Force 
frac_comp=frac_comp(frac_comp~=0);
BF_our_ratio=mean(frac_comp)


