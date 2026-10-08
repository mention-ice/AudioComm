function k = audioComm_prbs(n)
% Scrambling sequence of the audio communication demo:
% LFSR x^15 + x^14 + 1 with all-ones seed (same as acoustic_modem.html and audiocomm.py)
s = ones(1,15);              % s(1) = most recent bit
k = zeros(1,n);
for i = 1:n
    fb = xor(s(15),s(14));
    s = [fb s(1:14)];
    k(i) = fb;
end
end
