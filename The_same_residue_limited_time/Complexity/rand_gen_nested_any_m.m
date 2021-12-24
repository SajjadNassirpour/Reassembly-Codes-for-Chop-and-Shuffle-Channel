% This function generates random data string based on Bernoulli (B_p)
% Length of data sting = d_sec*(m^(L-1))
function [data_bits,data_bits_total] = rand_gen_nested_any_m(d_sec,m,L,B_p)

data_bits_length=d_sec*(m^(L-1));
data_bits=double(rand(1,data_bits_length)<=B_p);

data_bits_total=data_bits;
data_bits=reshape(data_bits,d_sec,m^(L-1))';

end