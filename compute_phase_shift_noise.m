function [snr_of_distortion] = compute_phase_shift_noise(fq, subcarrier_frequencies, bw, M, phase_shift_mtx, freq_shift_mtx)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
sample_count = len(subcarrier_frequencies) * 10e3;
freq_div = floor(bw/sample_count);
freq_vector = 0:freq_div:bw
y = zeros(shape(freq_vector));
for k = 1:1:numel(subcarrier_frequencies)
    fc = subcarrier_frequencies()
    y = y + sinc(freq_vector - fc)
end
    

end