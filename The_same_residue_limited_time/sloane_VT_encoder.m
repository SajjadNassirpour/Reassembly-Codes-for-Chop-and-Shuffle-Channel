% This function uses our VT approach to encode the data sting with respect to the desired residue 
function [X,p]=sloane_VT_encoder(d,p,r_desired)

    n=length(d)+p;
    index=1:length(d);
    VT_sum=dot(d,index);
    r_prime=mod(VT_sum,(n+1));
    
    if r_prime>=r_desired
        r=r_prime-r_desired;
    else
        r=(n+1)-(r_desired-r_prime);
    end
    j=1;
    answer=0;
    while j<p+1 && answer==0
        if(2*p+1)*j-j^2>=2*r
            answer=1;
            y=j-1;
        end
        j=j+1;
    end
    last_bit=r-(2*p*y-(y*(y-1)))/2;

    z=zeros(1,(p-y));
    if last_bit~=0
        z(end-last_bit+1)=1;
    end
    if y==0
        parity=z;
    else
        parity=[ones(1,y) z];
    end
    
    X=[d parity];
end
