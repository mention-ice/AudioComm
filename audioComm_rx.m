% audioComm_rx.m
% Receiver of the two-device audio communication demo.
% Start this script FIRST, then start the transmitter (audioComm_tx.m on
% another computer, or the WAV file played from a phone) while it listens.
% Parameters are in audioComm_params.m (must be identical on both machines).
%
% Receiver chain:
%  1. record (live spectrogram on screen)
%  2. non-coherent tone detection: energy at F0 and F1 over a sliding
%     window of one symbol, for every sample offset
%  3. frame + symbol synchronization: correlation with the known MLS pilot
%     gives the start of the preamble; the header (sent 3 times) gives the
%     message type and size, hence where to look for the postamble
%  4. clock offset: (measured pilot distance)/(expected distance) gives the
%     ratio between the two sound cards' sampling clocks (in ppm)
%  5. decision on each symbol at the corrected instants, descrambling

clear; close all; clc
audioComm_params

source  = 'mic';         % 'mic': record now, 'file': decode a recording (e.g. made with a phone)
recFile = 'recording.wav';
Fe_r    = 48000;         % receiver sampling frequency (does not need to equal Fe)
margin  = 8;             % extra listening time (s) to start the transmitter
                         % (listening time = local frame + margin: increase it if the
                         %  transmitter sends a longer message, e.g. from the web page)
saveRec = 1;             % save the recording to replay/analyse it later

%% 1. acquisition
if strcmp(source,'mic')
    listenTime = frameDuration + 0.4 + margin;
    recObj = audiorecorder(Fe_r,16,1);
    fig = figure('Name','Receiver - listening','Color','w');
    hIm = imagesc(0); axis xy; colormap(jet)
    xlabel('time (s)'); ylabel('frequency (kHz)')
    record(recObj,listenTime);
    disp(sprintf('listening for %.1f s: start the transmitter now',listenTime))
    Lw = round(T*Fe_r); nfft = 2^nextpow2(Lw); win = 3;   % show the last 3 s
    while isrecording(recObj)
        pause(0.2)
        try
            y = getaudiodata(recObj);
        catch
            continue
        end
        y = y(max(1,end-round(win*Fe_r)+1):end);
        nb = floor(length(y)/Lw);
        if nb < 2, continue, end
        P = abs(fft(reshape(y(end-nb*Lw+1:end),Lw,nb),nfft)).^2;
        set(hIm,'XData',(0:nb-1)*T,'YData',(0:nfft/2-1)/nfft*Fe_r/1e3,'CData',10*log10(P(1:nfft/2,:)+1e-9))
        axis tight; ylim([0 1.5*F1/1e3])
        title(sprintf('listening... %.1f / %.1f s',recObj.CurrentSample/Fe_r,listenTime))
        drawnow
    end
    y = getaudiodata(recObj);
    close(fig)
    if saveRec, audiowrite(recFile,y,Fe_r); end
else
    [y,Fe_r] = audioread(recFile);
end
y = y(:,1); y = y - mean(y);

%% 2. tone energies over a sliding window of one symbol
Lr = T*Fe_r;                         % samples per symbol at Rx (may be non-integer)
Lw = round(Lr);
n  = (0:length(y)-1)';
z0 = abs(filter(ones(Lw,1),1,y.*exp(-1i*2*pi*F0*n/Fe_r)));
z1 = abs(filter(ones(Lw,1),1,y.*exp(-1i*2*pi*F1*n/Fe_r)));
z0 = [z0(Lw:end); zeros(Lw-1,1)];    % z(m): window starting at sample m
z1 = [z1(Lw:end); zeros(Lw-1,1)];
d  = (z1-z0)./(z1+z0+eps);           % in [-1,1]: +1 for a '1', -1 for a '0'

%% 3. synchronization with the pilot (preamble), header, postamble
a   = 2*pilot-1;
off = round((0:np-1)*Lr);
M   = length(d) - off(end);
c   = zeros(M,1);
for j = 1:np
    c = c + a(j)*d(off(j)+(1:M));
end
c = c/np;                            % 1 = perfect match
thr = 0.3;                           % detection threshold
llr = @(m) log((z1(m)+eps)./(z0(m)+eps));
start = 1; hdr_rx = [];
while isempty(hdr_rx)
    cand = find(c(start:end) > thr, 1) + start - 1;
    if isempty(cand), error('no frame found: check rate/tones, volume and listening time'), end
    [cmax,i] = max(c(cand:min(M,cand+2*Lw)));
    p1 = cand + i - 1;               % start of the preamble
    % header: 24 bits sent 3 times each, soft combining of the repetitions
    hm = round(p1 + (np + (0:NH-1))*Lr);
    hb = sum(reshape(llr(hm),HREP,HBITS),1) > 0;
    typ = bin2dec(char(hb(1:2)+'0')); ha = bin2dec(char(hb(3:13)+'0')); hbb = bin2dec(char(hb(14:24)+'0'));
    if typ == 0 && ha >= 1 && hbb == 0
        hdr_rx = struct('type','text','n0',8*ha);
    elseif typ == 1 && ha >= 1 && hbb >= 1 && ha*hbb <= 20000
        hdr_rx = struct('type','image','n0',ha*hbb,'h',ha,'w',hbb);
    else
        start = p1 + Lw;             % false alarm: keep searching
    end
end
n0r = hdr_rx.n0;
Delta = (np+NH+n0r)*Lr;              % expected distance between the two pilots
tol   = round(0.005*Delta) + Lw;     % search +-0.5% around it
rng_  = max(1,round(p1+Delta)-tol):min(M,round(p1+Delta)+tol);
[best,i] = max(c(rng_)); p2 = rng_(i);
if isempty(best) || best < 0.5*cmax
    warning('postamble not found: the recording may not contain the whole frame')
    p2 = round(p1 + Delta); best = NaN;
end
rho = (p2-p1)/Delta;                 % ratio of Rx clock to Tx clock
ppm = (rho-1)*1e6;

%% 4. decision on each information symbol
mk = round(p1 + (np + NH + (0:n0r-1))*Lr*rho);
r  = llr(mk)';                       % soft output ('LLR')
b_hat = double(r > 0);
b0_hat = mod(b_hat + audioComm_prbs(n0r),2);   % descrambling
known = (n0r == n0);                 % same size as the local message: we can compute the BER
if known, ber = mean(b_hat ~= b); else, ber = NaN; end
pilot_err = sum((z1(round(p2+(0:np-1)*Lr*rho)) > z0(round(p2+(0:np-1)*Lr*rho)))' ~= pilot);

%% 5. display
fprintf('frame found at t = %.2f s, pilot match %.2f / %.2f\n',p1/Fe_r,cmax,best)
fprintf('clock offset Rx/Tx: %+.0f ppm (%.1f samples of drift over the frame)\n',ppm,p2-p1-Delta)
fprintf('bit rate %g bit/s, %d bits, postamble errors %d/%d, bit error rate %.4f\n',1/T,n0r,pilot_err,np,ber)

figure('Name','Receiver','Color','w')
subplot(2,2,1)
tt = (0:length(y)-1)/Fe_r;
plot(tt,y); hold on
yl = ylim; plot([1 1]*p1/Fe_r,yl,'g',[1 1]*(p2+np*Lr)/Fe_r,yl,'r'); hold off
xlabel('time (s)'); title('received signal (green: start, red: end of frame)')
subplot(2,2,2)
plot((1:M)/Fe_r,c); hold on; plot([p1 p2]/Fe_r,c([p1 p2]),'ro'); hold off
xlabel('time (s)'); title(sprintf('pilot correlation, clock offset %+.0f ppm',ppm))
subplot(2,2,3)
if known, s1 = b==1; else, s1 = r>0; end
histogram(r(s1),60,'FaceColor','r'); hold on; histogram(r(~s1),60,'FaceColor','b'); hold off
xlabel('log(E_{F1}/E_{F0})'); legend('sent 1','sent 0'); title('soft output')
subplot(2,2,4)
if strcmp(hdr_rx.type,'text')
    msg_hat = char(bin2dec(char(reshape(b0_hat,8,[])' + '0')))';
    axis off
    if istxt, text(0,0.7,['sent:     ' msg],'FontSize',14), end
    text(0,0.4,['received: ' msg_hat],'FontSize',14,'Color','b')
    text(0,0.1,sprintf('%g bit/s, BER = %.3f',1/T,ber),'FontSize',12)
    disp(['decoded message: ' msg_hat])
else
    img_hat = reshape(b0_hat,hdr_rx.h,hdr_rx.w);
    if ~istxt && isequal(size(img),size(img_hat))
        imshow([img ones(size(img,1),3) img_hat]); title(sprintf('sent | received   (%g bit/s, BER = %.3f)',1/T,ber))
    else
        imshow(img_hat); title(sprintf('received (%g bit/s)',1/T))
    end
end
