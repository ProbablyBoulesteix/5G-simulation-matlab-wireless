function [energy_per_symbole] = compute_energy_per_symbol(bit_energy,M)
%function finds the per-bit energy of an M-ary encoding scheme
%note: for most cases other than BPSK, symbol energy is log2(M) * bit
%energy

if M == 0 %disconnected/error case, no energy in signal
    energy = 0 ;
elseif (M == 2)
    energy = bit_energy ;
else
    energy = log2(M) * bit_energy; 
end
energy_per_symbole = energy;


end