% This duction generates random breaking points based on Geometric(pn)

function ran_vec=break_points_nested(cut_rule,n,initial_random_length)
    N=floor(cut_rule)-1;
    sum_ran_vec=0;
    while sum_ran_vec<=n
        ran_vec=geornd(cut_rule/n,1,initial_random_length*(N+2));
        sum_ran_vec=sum(ran_vec);
    end
    if size(find(ran_vec==0),2)~=0
        ran_vec=ran_vec+ones(1,initial_random_length*(N+2));
    end
    
    pre_rran=0;
    ii=1;
    while ii<=length(ran_vec) && sum(ran_vec(1:ii))<n
        ii=ii+1;
    end
    if ii==1
        ran_vec=0;
    else
        rran_vec=zeros(1,(ii-1));
        for rran=1:length(rran_vec)
               rran_vec(1,rran)=ran_vec(1,rran)+pre_rran;
               pre_rran=rran_vec(1,rran);
        end
        ran_vec=rran_vec;
    end
end