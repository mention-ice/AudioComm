% audioComm_params.m
% Shared parameters for the two-device version of the audio communication
% demo (audioComm_tx.m on the transmitter, audioComm_rx.m on the receiver).
% The frame is the same as in acoustic_modem.html (web page) and
% audiocomm.py (Python), so a WAV made or recorded by one can be decoded by
% the others:
%
%   0.1 s warm-up | MLS-127 preamble | 24-bit header x3 | scrambled data | MLS-127 postamble
%
%   header: type (2 bits: 0 text, 1 image), a (11 bits), b (11 bits), MSB
%           first, each bit sent 3 times. text: a = number of characters,
%           b = 0. image: a = rows, b = columns
%   data:   text = 8 bits per character; image = img(:) (column-major), 1 = white
%   scrambling: XOR with the LFSR x^15+x^14+1 sequence (all-ones seed)
%
% Differences with audioComm_v4.m:
%  - the receiver finds the frame in a recording that starts at any time
%    (preamble) and reads the message type and size from the header
%  - the distance between preamble and postamble gives the clock
%    (sample-rate) offset between the two sound cards, which is compensated
%  - default tones are lower (4/8 kHz) and symbols longer, to survive a
%    few meters of air, phone speakers and room echoes

%% message (transmitter side)
istxt   = 0;             % 1: send a text, 0: send an image
txt     = 'Hello CentraleSupelec!';
imgFile = 'bartS.png';
imgStep = 2;             % keep 1 pixel out of imgStep in each direction (1 = full image)

%% modulation (binary FSK, as in audioComm_v4.m) -- must match on both sides
Fe   = 48000;            % Tx sampling frequency (48 kHz is supported by all phones/laptops)
rate = 400;              % bit rate (200, 400, 800 or 1600 bit/s in the web page)
F0   = 4000;             % tone for bit '0'  (web page "8 / 16 kHz": F0 = 8000, F1 = 16000)
F1   = 8000;             % tone for bit '1'
T    = 1/rate;           % symbol duration (s)

%% frame constants
% MLS of length 127, as given by (GenerateMLS(7,1)+1)/2. (Without the flag,
% GenerateMLS seeds itself from the clock, so Tx and Rx pilots would differ.)
pilot = '0000000111111011111001111010111000011011101001100010101100000101111000111011011001001010010000100111001011010001000110011010101' - '0';
np    = length(pilot);
HBITS = 24; HREP = 3; NH = HBITS*HREP;

%% build the frame (used by the transmitter, and by the receiver to compute the BER)
if istxt
    msg = txt;
    b0 = reshape(dec2bin(double(msg),8)' - '0',1,[]);
    hdr = [dec2bin(0,2) dec2bin(length(msg),11) dec2bin(0,11)] - '0';
else
    I = imread(imgFile);
    if ndims(I) == 3, I = mean(double(I),3); else, I = double(I); end
    I = I(1:imgStep:end,1:imgStep:end);
    img = I > max(I(:))/2;           % black and white image
    b0 = double(img(:)');
    hdr = [dec2bin(1,2) dec2bin(size(img,1),11) dec2bin(size(img,2),11)] - '0';
end
n0 = length(b0);                     % number of information bits
key = audioComm_prbs(n0);
b = mod(b0+key,2);                   % scrambled data
sp_head = mod(1:round(0.1*rate),2);  % 0.1 s of alternating bits to wake up the sound card
s = [sp_head pilot kron(hdr,ones(1,HREP)) b pilot];   % transmitted bit sequence
nh = length(sp_head);
frameDuration = length(s)*T;         % duration of the transmitted sound (s)
