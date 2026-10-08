% audioComm_tx.m
% Transmitter of the two-device audio communication demo.
% Run audioComm_rx.m on the receiving computer FIRST, then run this script
% (or play the saved WAV file from a phone).
% Parameters are in audioComm_params.m (must be identical on both machines).

clear; close all; clc
audioComm_params

playSound = 1;           % play the signal on this computer's speaker
saveWav   = 1;           % save the signal as a WAV file (to be played from a phone)

%% modulation: one sinusoid per bit
L = round(T*Fe);                     % samples per symbol
t = (0:L-1)/Fe;
x0 = sin(2*pi*F0*t);                 % waveform for '0'
x1 = sin(2*pi*F1*t);                 % waveform for '1'
S = x0'*(1-s) + x1'*s;               % one column per symbol
x = 0.9*S(:);
x = [zeros(round(0.2*Fe),1); x; zeros(round(0.2*Fe),1)];  % short silences around the frame (as the web page)

if saveWav
    if istxt, tag = 'txt'; else, tag = 'img'; end
    wavName = sprintf('audioComm_%s_F%d-%d_%dbps.wav',tag,F0,F1,rate);
    audiowrite(wavName,x,Fe);
    disp(['saved ' wavName])
end

%% what the students see on the transmitter screen
figure('Name','Transmitter','Color','w')
if istxt
    subplot(2,1,1); axis off
    text(0.5,0.5,['sending: "' msg '"'],'FontSize',20,'HorizontalAlignment','center')
else
    subplot(2,2,1); imshow(img); title(sprintf('sending %dx%d image',size(img,1),size(img,2)))
    subplot(2,2,2)
end
nShow = 8;                           % first data symbols, in time
k0 = nh+np+NH;
plot((0:nShow*L-1)/Fe*1e3, reshape(S(:,k0+1:k0+nShow),1,[]))
xlabel('time (ms)'); title(['bits ' num2str(b(1:nShow))])
subplot(2,1,2)
nfft = 2^nextpow2(L);
specgram_img = abs(fft(reshape(x(1:floor(length(x)/L)*L),L,[]),nfft)).^2;
imagesc((0:size(specgram_img,2)-1)*T,(0:nfft/2-1)/nfft*Fe/1e3,10*log10(specgram_img(1:nfft/2,:)+1e-6))
axis xy; ylim([0 1.5*F1/1e3]); xlabel('time (s)'); ylabel('frequency (kHz)')
title(sprintf('%d bits, %g bit/s, %.1f s of sound',n0,1/T,length(x)/Fe))
drawnow

if playSound
    disp(sprintf('transmitting %d bits at %g bit/s (%.1f s)...',n0,1/T,length(x)/Fe))
    p = audioplayer(x,Fe);
    playblocking(p);
    disp('done.')
end
