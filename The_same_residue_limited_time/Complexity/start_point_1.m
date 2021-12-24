% This function determines the location of the first bi in each VT codeword

function start_p=start_point_1(m,L,Laa,d_sec,parity)

start_p=zeros(1,m^(L-Laa)-1);
coef=1:m^(L-Laa)-1;
for La=1:L-1
    
    if La==1
        start_p(coef)=start_p(coef)+m^(Laa-1)*(coef)*(d_sec+parity(La))+1;
    else
        
        start_p(coef)=start_p(coef)+(floor(m^(Laa-1)*coef/m^(La-1)))*(parity(La));
    end
    
end
start_p=[1 start_p];
end
    