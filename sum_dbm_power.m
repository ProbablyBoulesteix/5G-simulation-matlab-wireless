function total_dbm = sum_dbm_power(sig1, sig2)
%function takes two dbm-formated power measurements returns the power of
%the sum of those signals 
%  NOTE: A_db and B_db != A + B -> we SUM the power of the signals in
%  watts/mW and convert that to dbm

%both in mW
sig1_lin = 10^(sig1/10);
sig2_lin = 10^(sig2/10);

total_lin = sig1_lin + sig2_lin;
total_dbm = 10* log10(total_lin);

end