function [v, f, p] = find_spectral_peak(data, up)
freq_range = up.paramSet.rr_range/60;

cand_els = zeros(length(data.power),1);
for s = 2 : (length(data.power)-1)
    if data.power(s) > data.power(s-1) & data.power(s) > data.power(s+1) & data.freqs(s) > freq_range(1) & data.freqs(s) < freq_range(2)
        cand_els(s) = 1;
    end
end
clear s freq_range

temp = data.power - min(data.power);  % because the FFT and ACF have -ve values for the power spectra.
[~, r_el] = max(temp.*cand_els); clear cand_els
r_freq = data.freqs(r_el); clear r_el
v = 60*r_freq;

%% Store spectrum
f = data.freqs;
p = data.power;
end