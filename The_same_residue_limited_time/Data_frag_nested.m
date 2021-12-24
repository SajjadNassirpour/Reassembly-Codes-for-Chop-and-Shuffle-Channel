% This function assigns labels to different fragments

function [C_seperated,C_sep_org]=Data_frag_nested(C,ran_vec,cut_rule)
    n=length(C);
    if ran_vec(end)~=n
        section_number=length(ran_vec)+1;
    else
        section_number=length(ran_vec);
    end
        
    C_seperated=zeros(section_number,floor(cut_rule));
    po_0=0;
    for i=1:(section_number-1)
        po_1=ran_vec(1,i);
        C_seperated(i,1)=po_1-po_0;
        C_seperated(i,2:(po_1-po_0+1))=C(1,(po_0+1):po_1);
        po_0=po_1;
    end
    po_1=n;
    C_seperated(section_number,1)=po_1-po_0;
    C_seperated(section_number,2:(po_1-po_0+1))=C(1,(po_0+1):po_1);
    C_sep_org=C_seperated;
    %====== Sort C_seperated in descending order ======%
    abc=C_seperated(:,1);
    [aa,bb]=sort(abc,1,'descend');
    C_seperated_sorted=zeros(section_number,max(C_sep_org(:,1)+2));
    for i=1:section_number
        C_seperated_sorted(i,:)=[bb(i) C_seperated(bb(i),:)];
    end
    C_seperated=C_seperated_sorted;
end