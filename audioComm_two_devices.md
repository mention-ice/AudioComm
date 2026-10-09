# Acoustic modem demo: web page, Python and MATLAB

Each version's page and Python script (`acoustic_modem_v0.html` and `audiocomm.py` for v0, `acoustic_modem_v1.html` and `audiocomm_v1.py` for v1, `acoustic_modem_v2.html` and `audiocomm_v2.py` for v2, `acoustic_modem_v3.html` and `audiocomm_v3.py` for v3) use the **same frame**: a WAV made or recorded by one is decoded by the other (the transmitted samples were checked to be identical). v1 also reads every v0 frame, v2 every v1 frame, and v3 every v2 frame. The MATLAB scripts still use the earlier frame (fixed settings, 127-bit preamble, no code) and do not read the current frames.

| File | Role |
|---|---|
| `acoustic_modem.html` | **Live demo**, the latest version (now v2, the same file as `acoustic_modem_v2.html`). Transmitter and receiver in one page, for laptops and phones. |
| `acoustic_modem_v0.html` | v0, the original page, frozen. |
| `audiocomm.py` | Python version for the labs (v0). Uses numpy only (`sounddevice` is needed only for live play and record). It does not read hopping or turbo frames: use `audiocomm_v1.py` for those. |
| `audioComm_tx.m`, `audioComm_rx.m`, `audioComm_params.m`, `audioComm_prbs.m` | MATLAB version (`audioComm_v4.m` is untouched) |
| `acoustic_modem_v1.html`, `audiocomm_v1.py` | **v1**: v0 plus room sounding, tone hopping, a turbo code, a fix for cut-off endings, fine timing with a fitted LLR scale, a hop band starting at 3.2 kHz, Reed–Solomon received but no longer sent, and turbo code rates 1/3 and 1/4 (the eight "v1" sections below). v0 = `acoustic_modem_v0.html` and `audiocomm.py`, frozen. |
| `acoustic_modem_v2.html`, `audiocomm_v2.py` | **v2**: v1 plus four tones per symbol (4-FSK) as a choice next to binary, 3200 and 6400 bit/s with four tones, and a narrower search for the postamble (see the two "v2" sections). v1 is unchanged. |
| `acoustic_modem_v3.html`, `audiocomm_v3.py` | **v3**: v2 plus OFDM as a third modulation, hundreds of subcarriers at once with a guard interval of 21 or 5 ms, 17 or 20 kbit/s before the code (see the last section). v2 is unchanged, and `acoustic_modem.html` still serves v2. |

## Running the live demo

**Laptop only:** open `acoustic_modem.html` in Chrome directly from the folder. The microphone works on local files.

**Phones:** the browser gives a page microphone access only over https, so the page must be hosted. The simplest way is GitHub Pages:

1. Create a public repository.
2. Upload the file as `index.html`.
3. Go to Settings → Pages → "Deploy from branch: main".

You get a link like `https://<user>.github.io/<repo>/`. Show it as a QR code in the slides.

This demo is published at `https://mention-ice.github.io/AudioComm/` (repository `mention-ice/AudioComm`): `acoustic_modem.html` is the version to use in class (now v2), and `acoustic_modem_v0.html` to `acoustic_modem_v3.html` keep each version (v3 once it is pushed). It was first published at `sheng-yang-cs.github.io/AudioComm`; GitHub does not redirect a Pages address when a repository moves, so links and QR codes printed before the move no longer work.

Receivers need no settings: the rate, tones, code and preamble length travel in the header of each frame, so the plain link (or the page's QR code) is enough. On a phone, add `?role=tx` to open the Transmit tab.

### In class

1. The students open the link and press **Start listening**. Their phones show a live spectrogram and the preamble match.
2. You press **Send** on your laptop: the demo image, a text typed live, a drawing, or a photo.
3. Every phone shows the message filling in bit by bit. At the end, each phone shows:
   - the clock offset of that phone relative to the laptop, in ppm
   - the rate, tones and preamble read from the header
   - the errors on the postamble
   - the BER (computed for the demo image)
   - the histogram of the soft output

What to vary:

- **Distance and position:** compare the BER between the front and the back of the room.
- **Rate:** at 200 or 400 bit/s it works across the room; at 1600 bit/s it fails because of the echoes.
- **Tones:** 8/16 kHz were the tones of `audioComm_v4.m`. They are weaker on phone speakers, and some adults no longer hear 16 kHz.
- **Noise:** have the students talk or clap during a transmission.
- **Swap roles:** a student's phone transmits a drawing and the projector laptop receives it.

Each phone can download the received recording as a WAV. The students then decode and analyse it in the lab with `audiocomm.py` or `audioComm_rx.m` (`source = 'file'`).

## Frame

Control part, always 4/8 kHz at 400 bit/s (the only thing receivers know in advance, like the SIGNAL field of Wi-Fi):
`0.1 s warm-up | MLS preamble (127, 255 or 511) | header (×1, ×2 or ×4)`
Message part, in the tones and rate the header announces:
`coded, scrambled data | MLS postamble (same sequence as the preamble)`

- **Bits:** bit 0 is a burst of F0 and bit 1 a burst of F1, as in v4.
- **Tones:** 4/8 kHz (default) or 8/16 kHz.
- **Rate:** 200, 400, 800 or 1600 bit/s.
- **Header:** type 2 | a 11 | b 11 | code 2 | rate 2 | tones 2 | 0 0, then a CRC-8 (x⁸+x²+x+1), all through the convolutional code (92 bits). The receiver needs no prior knowledge of the message or of the settings; the CRC rejects false detections and misread headers.
- **Scrambler:** the LFSR x¹⁵+x¹⁴+1, which replaces `randsrc` and is easy to reproduce in every language.
- **Error correction:** chosen on the transmitter ("Error correction" in the page, `--code` in `audiocomm.py`) and sent in the header, so receivers need no setting.
  - *Hamming (7,4)*, hard decisions: corrects 1 bit in 7. Rate 4/7.
  - *Convolutional (171,133), K = 7*, soft-decision Viterbi: the strongest against random errors. Rate 1/2.
  - *Reed–Solomon (255,223)* over GF(256): corrects 16 bytes in 255, good against bursts, poor against spread-out errors. Rate ≈ 7/8.
  - Hamming and convolutional bits are interleaved (rectangular, √n columns) to break up bursts; Reed–Solomon codewords are byte-interleaved.
  - The receiver shows the decoded message next to "same channel errors, without the code", and counts the channel errors the code corrected.

## Receiver

1. A sliding non-coherent energy detector at F0 and F1 (one symbol window, at every sample).
2. Correlation with the MLS pilot gives the frame start and the symbol timing.
3. The header gives the message length.
4. The postamble position gives the clock ratio between the two devices.
5. Decisions are taken at the corrected instants.
6. Soft output, chosen on the receiver ("Soft output" in the page, `--soft` in `audiocomm.py`):
   - *ln(E1/E0)*: the plain energy ratio; right sign when both tones arrive equally strong, but its size is not a reliability.
   - *Calibrated LLR* (default): the level A and noise σ² of each tone are measured on the 127 known preamble bits
     (E[E²] = A² + 2σ² when the tone is on, 2σ² when off), then λ = ln I0(A1E1/σ1²) − A1²/2σ1² − ln I0(A0E0/σ0²) + A0²/2σ0²,
     the LLR of non-coherent FSK (Rician / Rayleigh). Energies are first divided by the local signal level (±24 symbols),
     because the echo, the main disturbance in a room, follows the signal. It fixes tone imbalance and gives the Viterbi
     decoder correctly scaled reliabilities. In the page, switching after a reception re-decodes the same recording.

**Preamble length** (Transmit card, `--pre` in `audiocomm.py`): 127 (default), 255 or 511 bits, with the header sent 1, 2 or 4 times. The receiver runs one correlator per length, each with threshold 0.3·√(127/N), so all three have the same false-alarm rate while each doubling gives about 3 dB of detection margin. In simulation (conv code, 4 reflections + noise) 511 found every frame and read every header 2 to 4 dB below where 127 started missing them; for real long range combine 511 with 200 bit/s, since the message then becomes the weak link. The soft output is calibrated on the postamble, which is sent at the message's tones and rate.

The page also tracks timing on the fly with an early–late gate, so the image fills in correctly before the postamble arrives.

## Why audioComm_v4.m could not be used as-is across two devices

1. **Different pilots on each machine.** `GenerateMLS(7)` without the flag seeds from the clock, so two computers generate different pilots.
2. **Fixed recording window.** The recording starts with playback and lasts only 1.2× the message.
3. **No clock-drift handling.** Symbol timing is estimated once. About 100 ppm of drift over 10.5 s is roughly 1.7 symbols at 1600 bit/s.
4. **Tones hard for phones.** The 16 kHz tone at Fe = 80 kHz is badly suited to phones.

## Testing done

- **Simulated room channel:** ±150 ppm, Rx at 44.1 kHz versus Tx at 48 kHz, 3-path echo up to 30 ms, and noise. Text and image were received with 0–1 errors at 200–400 bit/s.
- **End to end in Chromium:** the page decoded the image with 1 error out of 4,218, using a simulated room recording as the microphone.
- **Cross-decoding:** web → Python, web → MATLAB (Octave), Python → web, Python → MATLAB and MATLAB → Python all decode.
- **Not yet tested:** real phones, real speakers and MATLAB itself (only Octave was used). Do one dry run in the classroom.

Echo test note: with echoes whose delay is an exact number of tone periods, one echo cancels the direct path and errors appear (multipath fading). That is a nice point to show by moving the phone a few centimetres.

## v1: room sounding (`acoustic_modem_v1.html`, `audiocomm_v1.py`)

v1 reads every v0 frame and adds a third frame type that measures the room. v1 link: `https://mention-ice.github.io/AudioComm/acoustic_modem_v1.html` (it was also served as `acoustic_modem.html` until v2). A v1 page's QR code points to its own online copy, also when the laptop runs it from a local file.

- **Sending:** "Sound the room" in the Transmit card, or `audiocomm_v1.py tx --sound [--order 14|15|16] [--periods 1..8] [--pre 127|255|511]`. After the usual control part (preamble, header with type 2, a = MLS order, b = periods) it plays an MLS of 2^order − 1 samples of ±0.35 at 48 kHz (period 0.34, 0.68 or 1.37 s) periods + 1 times, the first one to fill the room, then the postamble in the control mode. It sounds like a loud hiss of 3 to 13 s.
- **Receiving:** each period is resampled to the transmitter's clock (ratio from the preamble–postamble distance, refined by the drift of the direct-path peak from one period to the next) and circularly cross-correlated with the MLS by FFT. The MLS autocorrelation is N at lag 0 and −1 elsewhere, so h = (c + Σr)/(N + 1). The periods are averaged (3 dB less noise per doubling). If the audio skipped or repeated samples, each period is realigned on its own peak and the period cut by the jump is dropped.
- **Measures:** the impulse response h(t) of speaker + room + microphone; RT60 from the Schroeder decay curve (T20, or T10 when the decay range above the noise is short; the decay hidden under the noise is added back from the fitted slope, as in Lundeby's method); the direct-to-reverberant ratio (±2.5 ms around the direct sound against the rest); C50; the rms delay spread; the frequency response (1/6 octave) and the gains at 4, 8 and 16 kHz; the SNR per frequency against the noise heard just before the frame, turned into the Es/N0 a 400 bit/s tone would get; the clock offset.
- **Outputs:** the page plots the energy decay and the frequency response, lists the values and downloads the impulse response as a 48 kHz WAV. `audiocomm_v1.py rx` prints the values, saves `room_ir.wav`, and `--plot` draws the same plots (matplotlib).
- **Room simulator:** `audiocomm_v1.py sim tx.wav --ir room_ir.wav --snr 30 --ppm 50 -o heard.wav` plays any frame through a measured response. Without `--ir` it uses a model room: direct sound plus a Gaussian tail decaying by 60 dB in `--rt60` seconds, scaled for `--drr`.
- **Choosing the period:** it must outlast the reverberation, otherwise the tail folds back onto the start of the response. 0.68 s suits a classroom, 1.37 s an amphitheatre.

**Tests (simulation, page and Python):** rooms with RT60 0.4 to 1.2 s and DRR −3 to +8 dB, microphones at 44.1 and 48 kHz, clock offsets up to ±200 ppm. RT60, DRR, C50 and delay spread came within 3 % (or 0.2 dB) of the true response, the clock offset within 4 ppm, also with 10 ms of samples lost or repeated mid-sounding. The page in Chromium, fed a simulated room as its microphone, gives the same values as Python. All 96 v0 frame settings still decode. Not yet tried with real speakers, phones and rooms.

**What the room model says about the modem** (RT60 0.6 s, convolutional code, 8 trials each): with DRR ≥ +3 dB every frame is read. At 0 dB, the critical distance, the 127 preamble loses 2 frames in 8 and the 511 preamble none. At −3 dB at most half the frames are read, and at −6 dB none. Beyond the critical distance the limit is the room's own echo of earlier bits landing on the other tone, not the noise, so a longer preamble or a lower rate does not help. This motivates the next v1 step, tone hopping (next section). The exact boundary depends on the early reflections: a model with a diffuse tail only already fails at +3 dB. Since the sounding measures the DRR, the class can find where that boundary lies in their own room.

## v1: tone hopping

- **Sending:** "Hopping" in the Tones choice of the Transmit card, or `audiocomm_v1.py tx --tones hop`. At R bit/s the bits are spread over H = 6400 / R pairs of adjacent tones from 3.2 to 16 kHz (32 pairs at 200 bit/s, 16 at 400, 8 at 800, 4 at 1600): bit 0 on 3200 + 2pR Hz, bit 1 on 3200 + (2p + 1)R Hz. (The band started at 1.6 kHz until the change in the last section; the results in this section were measured with that band.) Bit k of a segment uses pair p = (H/2 − 1)k mod H (k mod 4 when H = 4), so a pair comes back only every H bits.
- **Why it helps:** in a room each bit rings on for a large part of a second, hundreds of times longer than the bit. With two tones, the echo of all the earlier bits sits on both tones; beyond the critical distance it is louder than the direct sound and buries the bit. Hopping spreads that echo over 2H tones, so the echo on the two tones in use is about H times weaker than with two tones: 12 dB with 16 pairs, about 3 dB per doubling of the pairs. A lower rate gives more pairs in the same band, so with hopping, slowing down helps again.
- **Frame:** the header's tones field takes a third value (2 = hop). A hopping message has a hopping control part (warm-up, preamble and header at 400 bit/s over 16 pairs), otherwise far receivers would not even find the frame; the postamble hops at the message's rate. A room sounding sent with Hopping selected uses the hopping control part and postamble, so it is found farther away too. v0 receivers ignore hopping frames.
- **Receiver:** next to the three preamble correlators, three more run on the contrast (E1 − E0)/(E1 + E0) of whichever pair each preamble bit uses (sliding energies of the 32 control tones). Header and message energies are computed symbol by symbol on the symbol's pair, with the same early–late gate; the postamble is searched first where the tracking expects it. Speaker, room and microphone give each tone its own level, so each tone's energies are divided by their rms over the symbols that used its pair before the calibrated LLR; the page shows the spread of these levels.
- **Cost:** the preamble search now steps by 1/12 of a bit, then sample by sample around a candidate. With the extra correlators the receiver still needs less computing than before (46 instead of 56 ms per second of audio on the test machine), with the same detection sensitivity.
- **Fix in `audiocomm_v1.py`:** a preamble candidate that runs past the end of a short recording no longer stops the search (v0's `audiocomm.py` could miss the frame in a short, noise-free WAV straight from `tx`).

**Tests (simulation):** all 48 hopping settings (rate × preamble × code) and the 96 v0 settings decode, with settings read from the header; hopping and plain frames back to back; no false detection in 120 s of noise, whistle and clicks; the page and `audiocomm_v1.py` give identical WAVs and decode each other's; the page in Chromium decodes a hopping frame heard at DRR −6 dB. Room soundings with a hopping control part are now measured in the cases where the 4 / 8 kHz one was lost (DRR −3 dB with RT60 0.8 s, DRR −8 dB with RT60 1.5 s).

**Plain FSK against hopping in a model room** (RT60 0.6 s with early reflections, convolutional code, the page's receiver, texts received right out of 8):

| DRR | 4 / 8 kHz, 400 bit/s, preamble 127 | Hopping, 400 bit/s, 127 | 4 / 8 kHz, 200 bit/s, 511 | Hopping, 200 bit/s, 511 |
|---|---|---|---|---|
| +3 dB | 8 | 8 | 8 | 8 |
| 0 dB | 6 | 8 | 8 | 8 |
| −3 dB | 3 | 8 | 1 | 8 |
| −6 dB | 0 | 8 | 0 | 8 |
| −9 dB | 0 | 8 | 0 | 8 |
| −12 dB | 0 | 2 | 0 | 8 |

`audiocomm_v1.py` on a diffuse-tail room gives the same picture. In a more reverberant room (RT60 1.2 s, 4 trials) hopping at 200 bit/s with the 511 preamble still reads 4 of 4 at −12 dB and 1 of 4 at −15 dB. Since the DRR falls by about 6 dB per doubling of the distance beyond the critical distance, this moves the limit from about 1.4 times the critical distance (plain) to about 3 times (hopping, 400 bit/s) or 4 times (hopping, 200 bit/s); this conversion is the textbook rule, not a measurement.

Try `audiocomm_v1.py sim hop.wav --ir room_ir.wav` with a measured response to see where the boundary lies in your own room. Not yet tried with real speakers, phones and rooms.

## v1: turbo code

- **Choosing it:** "Turbo" in the Error correction choice of the Transmit card, or `audiocomm_v1.py tx --code turbo`. Rate 1/2 and the same length as the convolutional code (2k + 12 channel bits for k message bits), so a frame lasts as long with either.
- **Code:** as in 4G, two 8-state recursive systematic convolutional encoders (feedback 1 + D² + D³, parity 1 + D + D³: 13 and 15 octal). The first encodes the message, the second the message shuffled by an interleaver; parity bits are taken alternately from each, and each encoder is brought back to state 0 by 3 tail steps (12 tail bits in all). The interleaver works for any length: a Fisher–Yates shuffle driven by xorshift32 seeded with k, then spread (S-random, S = min(16, ⌊√(k/2)⌋)); at 456 bits it does as well as 4G's QPP interleaver. The codeword then goes through the same rectangular interleaver as the convolutional code.
- **Decoder:** two log-MAP (BCJR) decoders take turns, each passing the other its extrinsic LLRs, up to 8 iterations, stopping when the decisions no longer change. The page shows the number of iterations.
- **Header:** the 2-bit code field was full, so the turbo code (code 4) is sent as code 0 with the CRC-8 inverted, the way 4G tells the number of transmit antennas by masking the CRC of its broadcast channel. Receivers that predate it (v0, and v1 before the turbo code) see a wrong CRC and ignore turbo frames; frames with the other codes are unchanged.
- **Gain:** on the textbook channel (BPSK, white noise, 456 message bits) turbo needs about 1.3 dB less Eb/N0 than the convolutional code for the same frame error rate; an LDPC code of the same length, tried for comparison (not in the modem), about 1.1 dB.

Through the model rooms (hopping, 400 bit/s, a 57-character text, `audiocomm_v1.py`'s receiver, 24 trials with the same rooms and noise for each code), texts received right:

| Phone far away, limited by | Convolutional | LDPC (not in the modem) | Turbo |
|---|---|---|---|
| echo, DRR −11 dB | 21 | 23 | 22 |
| echo, DRR −12 dB | 14 | 17 | 18 |
| echo, DRR −13 dB | 4 | 6 | 7 |
| noise, SNR −9 dB | 14 | 22 | 22 |
| noise, SNR −10 dB | 2 | 12 | 14 |

The page's own receiver, in a model room with early reflections (16 trials each): at DRR −11 dB the convolutional code reads 14 texts and turbo 16; at −12 dB, 8 and 13; at −13 dB, 3 and 1.

So turbo gains 0.5 to 1 dB in the room, roughly 10% more distance. The gain is small because the bit error rate climbs steeply with distance (4.5% at DRR −10 dB, 24% at −14 dB): the frames that fail have 15 to 19% wrong bits, where the soft outputs carry less than 0.5 bit per channel bit, beyond any rate-1/2 code; and about 1 dB farther the frame itself is lost. Lowering the rate to 200 bit/s with the 511 preamble gains about 3 dB. Turbo also lost one frame at DRR −10 dB that the convolutional code read, probably because log-MAP decoding relies more on well-calibrated soft values than the Viterbi algorithm does.

**Tests (simulation):** all 120 settings with 4 / 8 and 8 / 16 kHz and all 60 hopping settings (rate × preamble × code, now five codes) decode with the settings read from the header; the page and `audiocomm_v1.py` give byte-identical turbo frames and the same interleavers (12 lengths up to 20000 bits); receivers without the turbo code ignore turbo frames. Not yet tried with real speakers and phones.

## v1: cut-off endings and the calibrated LLR

- **Symptom:** at 1600 bit/s the calibrated LLR sometimes did worse than ln(E1/E0), while it did better at lower rates.
- **Cause, found on a real recording** (2026-10-08, hopping, 1600 bit/s, turbo, demo image): the sound stopped 0.2 s before the end of the frame. The page closed its AudioContext as soon as the buffer had been played, while the speaker still had about 0.4 s of audio to play (output latency, typical of Bluetooth); the 0.2 s of silence after the frame did not cover it. The postamble was lost, so the receiver calibrated the LLR on noise (SNR −62 / −8 dB) and the calibrated LLR got 38% of the channel bits wrong, against 9% for ln(E1/E0), which needs no calibration. The postamble lasts 255 symbols: 0.16 s at 1600 bit/s, 1.3 s at 200 bit/s. So the same cut removes all of it at 1600 bit/s, but only part of it at lower rates, where it is still found.
- **Fix, transmitter:** the page closes the AudioContext 2 s after the end of playback, plus whatever output latency the browser reports.
- **Fix, receiver (page and `audiocomm_v1.py`):** when the postamble is not found, A and σ² are estimated on the message itself by EM. Start from the hard decisions E1 > E0, compute each tone's level and noise with every symbol weighted by the probability that it is a 1 (from the current LLR), and repeat 8 times. The ±24-symbol level tracking is then relative to the message's average level. Frames whose postamble is found decode exactly as before (56 of 56 identical in a regression test). Next to the SNR, the page says "from the message: postamble not found"; the script says "postamble not found (soft output calibrated on the message)".
- **Results:** on that recording, the page's receiver now reads the image with 0 pixel errors using the calibrated LLR (22 when replayed through Chromium's fake microphone), against 1373 (1339 in Chromium) before; ln(E1/E0) gives 256 (199 in Chromium). In simulation (room with RT60 0.6 s, DRR −3 dB, last 0.4 s of the sound cut, 6 trials per setting), texts read right with the calibrated LLR, before / after, against ln(E1/E0): hopping at 1600 bit/s with turbo 1 / 6, against 0; with the convolutional code 1 / 4, against 5. At 800 bit/s and below, the same texts were read before and after.
- **Also seen in that recording:** the 1600 bit/s symbol timing, taken from the 400 bit/s preamble, was about 8 samples (0.17 ms, a quarter of a symbol) late. Reading the message 8 samples earlier gives 6% channel errors instead of 21%. The page's tracking corrects most of this when the postamble is missing (8%). `audiocomm_v1.py` does not track, so it does not read that recording either way (17% channel errors calibrated, 20% with the ratio).

## v1: fine timing on the message, and a fitted LLR scale

- **Fine timing (tone hopping only):** the instants of the message symbols come from the 400 bit/s preamble and the postamble. On the real 1600 bit/s recording they were about 8 samples (a quarter of a symbol) late for the message's short symbols, while the preamble's long symbols hardly cared (its correlation is flat within a few samples). Once the message is in, the receiver now cuts it into up to 8 blocks and, in each, tries shifts of up to half a symbol or 0.5 ms, keeping the one that makes the two tones most distinct (mean |E1 − E0| / (E1 + E0)). It fits a straight line through these shifts, which also follows a residual clock drift, and keeps the new instants only if they make the tones clearly more distinct, symbol by symbol, than the instants it had (paired t-test, t > 2). The page shows the shift at the start and end of the message; the script prints it.
- **Why only with hopping:** with fixed tones, the echo of the earlier bits lies on the same two tones. In simulated rooms the same measure then often picked a worse instant: at 1600 bit/s and DRR 0 dB, 4 of 6 texts were read instead of 6. With hopping the echo lands on other pairs.
- **LLR scale:** the noise model behind the calibrated LLR assumes Gaussian noise. On the real recording, even when it was calibrated on the true bits, its LLRs were about 3 times too large: the echoes of short symbols have heavier tails. When the postamble is found, the LLR is now multiplied by the factor α ≤ 1 (searched from 1 down to 0.05) that best predicts the postamble's own bits, P(bit) = 1 / (1 + e^(−αL)). Emulating a postamble on the recording (255 known message bits) gave α ≈ 0.3, and the information carried per channel bit at that scale rose from 0.28–0.59 to 0.73–0.79, which is the best any scale gives. In the model rooms α stays at 1. It matters for the turbo decoder; the Viterbi decoder's decisions do not depend on the scale. Without the postamble (calibration on the message by EM) there is no α, and the LLR stays too large.
- **Results:** the real recording now decodes with 0 pixel errors in both the page and `audiocomm_v1.py`, with either soft output. Its channel errors fell from 8.2% to 6.1% calibrated and from 9.3% to 6.6% with the ratio, and the turbo decoder needs 3 to 4 iterations instead of 8. Before, the script got 17–20% of the channel bits wrong and could not read it, and the page's ratio left 256 pixel errors. In model rooms (RT60 0.6 s, DRR 0 / −3 / −6 dB, 6 trials per setting, all rates, both codes, fixed tones and hopping), no text was lost. With the preamble's instant made 7 samples late, as in the real room, hopping at 1600 bit/s read 6 texts instead of 5 with the convolutional code, and 4 instead of 2 with turbo and the ratio. All 120 fixed-tone and 60 hopping settings still decode, and the page and the script read the same texts from the same simulated recordings.

## v1: tone hopping from 3.2 kHz

- **Change:** the hopping pairs now start at 3.2 kHz instead of 1.6 kHz and run up to 16 kHz. The band keeps its width (12.8 kHz), so each rate keeps its number of pairs. At 1600 bit/s the four pairs are 3.2 / 4.8, 6.4 / 8.0, 9.6 / 11.2 and 12.8 / 14.4 kHz; the hopping control part (400 bit/s, 16 pairs) runs from 3.2 to 15.6 kHz.
- **Why, from the real 1600 bit/s recording:** checked against the true bits, nearly all the channel errors came from the lowest pair. The 1.6 kHz tone was read wrong 43% of the time and its partner at 3.2 kHz 8%, against 0 to 0.2% for every other pair. In the symbols that used it, the 1.6 kHz tone arrived about 10 dB weaker than the tones from 4.8 kHz up, compared with what leaks onto unrelated tones: it is low in a small speaker's range. The room's noise in that recording is concentrated below 1 kHz, and between 1 and 2 kHz it is still about 8 dB above its level between 2 and 5 kHz. A 1600 bit/s symbol lasts 0.6 ms, so each tone's detector is about 1.6 kHz wide, and the one at 1.6 kHz also collects that low-frequency noise. The preamble said the same at 400 bit/s: an SNR of 7 dB at 1.6 kHz against 13 to 27 dB for the tones above 2 kHz.
- **A coincidence the new band also avoids:** at 1600 bit/s the second harmonic of the 1.6 kHz tone fell exactly on its partner at 3.2 kHz, so a distorting speaker would push bit 0 towards bit 1. On this recording the 3.2 kHz energy during the 1.6 kHz symbols was no higher than on unrelated tones, so it was not the cause here. A harmonic of the bit-0 tone can only land on its partner if the bit-0 tone is at most the bit rate, which can no longer happen once the band starts at 3.2 kHz.
- **Results in simulation** (room model with RT60 0.6 s, the background noise of the real recording added at the given SNR over 1 to 14 kHz, 57-character text, 6 trials per setting, the page's receiver, old and new band on the same rooms and noise): over 88 settings (200 to 1600 bit/s, convolutional and turbo, DRR 0 / −3 / −6 dB, SNR 0 to 22 dB, with a flat speaker or one rolled off below 2.8 kHz like the recording's), 500 of 528 texts were read instead of 472, and only one setting lost a text (1600 bit/s, turbo, rolled-off speaker, DRR 0 dB, SNR 3 dB: 5 instead of 6). The channel bit errors fell in nearly every setting, e.g. from 4.1% to 0.7% at 1600 bit/s, convolutional, DRR 0 dB, SNR 9 dB. Texts read (old band → new band, out of 6):

| Setting | Convolutional | Turbo |
|---|---|---|
| flat speaker, DRR 0 dB, 1600 bit/s, SNR 3 dB | 2 → 6 | 5 → 6 |
| flat speaker, DRR 0 dB, 1600 bit/s, SNR 0 dB | 1 → 3 | 1 → 6 |
| rolled-off speaker, DRR −3 dB, 1600 bit/s, SNR 16 / 12 / 9 dB | 4 / 3 / 2 → 6 / 4 / 5 | 5 / 5 / 4 → 6 / 6 / 5 |
| flat speaker, DRR −6 dB, 1600 bit/s, SNR 22 / 16 / 12 dB | 3 / 1 / 0 → 3 / 1 / 0 | 4 / 3 / 1 → 5 / 4 / 3 |

At 800 bit/s and below every text was read with either band, with about 2 to 10 times fewer channel errors. Far from the speaker at 1600 bit/s (DRR −6 dB) the limit is the room's echo, which the band's position does not change.
- **Compatibility:** the hopping frames changed, control part included, so earlier v1 pages and scripts do not find the new hopping frames and the new ones do not find the old. Every phone needs to reload the page once this version is on GitHub Pages. Frames with fixed tones are unchanged, and v0 is untouched.
- **Tests (simulation):** all 120 fixed-tone and 60 hopping settings decode with the settings read from the header; the page and `audiocomm_v1.py` give byte-identical hopping WAVs at every rate and decode each other's; the page in Chromium reads a 1600 bit/s hopping frame. The new top tones (14.4 kHz at 1600 bit/s, up to 15.8 kHz at 200 bit/s) are not yet tried with real speakers and phones; in the recording, the control part's tones at 13.2 to 14 kHz still had an SNR of 25 to 27 dB.

## v1: Reed–Solomon no longer sent

- **Change:** the v1 page and `audiocomm_v1.py tx` no longer offer Reed–Solomon (255,223). Their receivers still decode it, so Reed–Solomon frames from the v0 page or `audiocomm.py` are read as before. The header is unchanged: code 3 stays Reed–Solomon, and the turbo code keeps its masked CRC. Giving code 3 to the turbo code would make v0 receivers decode turbo frames as Reed–Solomon and show garbage instead of ignoring them.
- **Why:** in a room Reed–Solomon fails as early as Hamming. Both decide each bit before decoding, so they fail once about 1–2% of the channel bits are wrong, while the convolutional and turbo decoders work on the soft values and still read texts with 5–9% wrong. Texts read right out of 6 (57-character text, room model with RT60 0.6 s, the real recording's noise, the page's receiver):

| Setting | None | Hamming | Reed–Solomon | Convolutional | Turbo |
|---|---|---|---|---|---|
| hopping 400 bit/s, DRR −6 dB, SNR 20 dB | 2 | 6 | 6 | 6 | 6 |
| hopping 400 bit/s, DRR −9 dB, SNR 20 dB | 0 | 1 | 1 | 6 | 6 |
| hopping 400 bit/s, DRR −12 dB, SNR 20 dB | 0 | 0 | 0 | 2 | 3 |
| hopping 400 bit/s, DRR 0 dB, SNR 3 dB | 0 | 3 | 3 | 6 | 6 |
| hopping 400 bit/s, DRR 0 dB, SNR 0 dB | 0 | 0 | 0 | 6 | 6 |
| hopping 1600 bit/s, DRR 0 dB, SNR 6 dB | 0 | 1 | 4 | 6 | 6 |
| 4 / 8 kHz 400 bit/s, DRR 0 dB, SNR 20 dB (out of the 5 frames found) | 0 | 1 | 3 | 5 | 5 |

- **What is lost:** its high rate (about 7/8). The demo image at 1600 bit/s lasts 4.9 s with Reed–Solomon, 7.0 s with the convolutional or turbo code, 6.4 s with Hamming and 4.4 s uncoded.
- **Tests (simulation):** Reed–Solomon frames built by the v0 page's own code (4 / 8 and 8 / 16 kHz, 200 to 1600 bit/s) are read right by the v1 receiver; the page shows four codes and falls back to None when its address asks for `code=rs`; all 120 fixed-tone and 60 hopping settings still decode.

## v1: turbo code rates 1/3 and 1/4

- **Choosing it:** with Turbo selected, a "Turbo code rate" row offers 1/2 (as before), 1/3 and 1/4, and the page's address keeps the choice (`code=turbo3`, `code=turbo4`); `audiocomm_v1.py tx --code turbo3` or `--code turbo4`. Receivers read the rate from the header.
- **Code:** the same two 8-state encoders and interleaver as at rate 1/2. Rate 1/3 sends every parity bit of both encoders, as cdma2000 does at rate 1/3 and as 4G's mother code before rate matching: k message bits give 3k + 12 channel bits. Rate 1/4 adds cdma2000's second parity output to each encoder, 1 + D + D² + D³ (17 octal), taken from each encoder in turn: 4k + 18 channel bits, the extra 6 being the tails' second parity bits. The three rates are nested (each one's bits include those of the higher rates); cdma2000's own rate-1/4 puncturing pattern may differ, it was not checked. The decoder is the same, with the extra parity LLRs added to each log-MAP decoder's branch metrics.
- **Header:** codes 5 and 6 are sent as code fields 1 and 2 with the CRC inverted, the same trick as for code 4. Receivers that predate them (v0, and v1 until now) see a wrong CRC or a code they do not know, and ignore the frame. Frames with the other codes, rate-1/2 turbo included, are unchanged bit for bit.
- **Textbook channel** (BPSK, white noise, 456 message bits, 400 frames per point), Eb/N0 for 10% frame errors: convolutional 2.7 dB, turbo 1/2 1.4 dB, 1/3 0.7 dB, 1/4 0.4 dB. At the same energy per message bit, rates 1/3 and 1/4 gain 0.7 and 1 dB over rate 1/2.

**In the model rooms, what counts is the airtime.** Tone hopping, 511-bit preamble, 57-character text, the page's receiver, room with RT60 0.6 s. "Message bit/s" is the bit rate times the code rate, and "Frame" is the length of the sound.

Limited by noise (DRR +6 dB, white noise, SNR measured over 3.2–16 kHz, 12 trials per point), SNR at which half the texts are read:

| Bit rate | Code rate | Message bit/s | Frame | SNR for half the texts |
|---|---|---|---|---|
| 800 | 1/4 | 200 | 5.6 s | −8.7 dB |
| 400 | 1/2 | 200 | 6.3 s | −9.0 dB |
| 400 | 1/3 | 133 | 7.4 s | −10.6 dB |
| 400 | 1/4 | 100 | 8.6 s | −11.7 dB |
| 200 | 1/2 | 100 | 9.9 s | −11.3 dB |
| 200 | 1/3 | 67 | 12.2 s | −12.4 dB (frame detection fails first) |
| 200 | 1/4 | 50 | 14.5 s | −12.5 dB (frame detection fails first) |

Limited by echo (little noise, 24 rooms per point), texts read out of 24:

| Bit rate, code rate | Message bit/s | Frame | DRR −12 | −14 | −16 | −18 | −20 | −22 | −24 dB |
|---|---|---|---|---|---|---|---|---|---|
| 400, 1/2 | 200 | 6.3 s | 19 | 2 | 0 | 0 | 0 | 0 | 0 |
| 400, 1/3 | 133 | 7.4 s | 22 | 16 | 2 | 1 | 3 | 1 | 1 |
| 400, 1/4 | 100 | 8.6 s | 23 | 22 | 11 | 7 | 9 | 8 | 5 |
| 200, 1/2 | 100 | 9.9 s | 24 | 22 | 13 | 16 | 16 | 16 | 17 |
| 200, 1/3 | 67 | 12.2 s | 24 | 24 | 21 | 23 | 24 | 24 | 24 |
| 200, 1/4 | 50 | 14.5 s | 24 | 24 | 24 | 24 | 24 | 24 | 24 |

- **Noise:** at the same bit rate, rates 1/3 and 1/4 tolerate 1.6 and 2.7 dB more noise than rate 1/2, but the frame lasts longer. Spending the same airtime on a lower bit rate does as well: 400 bit/s at rate 1/4 and 200 bit/s at rate 1/2 both carry 100 message bits per second and are 0.4 dB apart; 800 bit/s at 1/4 and 400 bit/s at 1/2 both carry 200 and are 0.3 dB apart. The 1 dB that rate 1/4 gains on the textbook channel is lost here, probably because a non-coherent detector loses more when each channel bit carries less energy: at low SNR the information it extracts per bit falls with the square of the SNR, where a coherent detector's falls in proportion.
- **Echo:** far beyond the critical distance the bit rate counts more than the code rate. In the same airtime, from DRR −16 to −24 dB, 200 bit/s at rate 1/2 read 78 of 120 texts and 400 bit/s at rate 1/4 read 40. 800 bit/s at 1/4 did worse than 400 bit/s at 1/2 (1 of 12 texts at DRR −12 dB, against 8 of 12, in a separate run). Below about −16 dB the counts stop falling. The phone then hears mostly reverberation, whose first milliseconds still carry the current symbol, and the interference is the reverberation of the earlier symbols sent on the same tone pair. Probably, then, what matters is how long a pair waits before it is used again: H symbols, 40 ms at 400 bit/s (16 pairs) and 160 ms at 200 bit/s (32 pairs), by which time the reverberation has decayed by 4 and 16 dB (60 dB per RT60). Halving the bit rate doubles that wait; lowering the code rate does not change it.
- **Both:** 200 bit/s with rate 1/3 or 1/4 read 116 and 120 of these 120 far texts, at 12.2 and 14.5 s for 57 characters.
- **Preamble:** with the 127-bit preamble the frame itself is lost at about the noise (SNR −9 to −10 dB) and echo (DRR −12 to −14 dB) where rate 1/2 already fails, so rates 1/3 and 1/4 only help with the 511-bit preamble.
- **For the farthest phones**, in a reverberant room: hopping, 200 bit/s, 511 preamble, turbo 1/3, or 1/4 for the last bit of margin. When airtime is short, halve the bit rate before lowering the code rate; in a quiet, dry room the two are worth the same.
- **Tests (simulation):** all 108 turbo settings (fixed tones and hopping, 4 bit rates, 3 preambles, 3 code rates) decode with the settings read from the header; the page and `audiocomm_v1.py` give byte-identical WAVs at rates 1/3 and 1/4, and the script decodes the page's; rate-1/2 codewords, decisions and iteration counts are identical to before; the previous page and script ignore the new frames (in Chromium the previous page stays on "Listening"); the page shows the rate row only for Turbo, keeps the rate in its address, and reads a rate-1/4 frame in Chromium. Not yet tried with real speakers and phones.

## v2: four tones per symbol (4-FSK)

v2 is `acoustic_modem_v2.html` and `audiocomm_v2.py`, made from v1; the v1 files are unchanged. Since 2026-10-09 `acoustic_modem.html` is a copy of `acoustic_modem_v2.html`, so the main GitHub Pages link serves v2; v1 stays at `acoustic_modem_v1.html`. Phones that had the page open need to reload it to read 4-tone frames. The study that led to it is `mfsk_study.md` (an ideal receiver); the results below come from the real v2 receiver.

- **Choosing it:** the Transmit card has a new "Modulation" row: "2 tones, 1 bit per symbol" (as before) or "4 tones, 2 bits per symbol". With 4 tones the rates are 400, 800, 1600 and 3200 bit/s: 200 bit/s is binary only and 3200 bit/s is 4 tones only, and switching moves 200 to 400 and 3200 to 1600. The page's address keeps the choice (`mod=4`); `audiocomm_v2.py tx --mod 4 --rate 3200 --tones hop`. Receivers read it from the header. Room soundings have no modulation choice.
- **Modulation:** a symbol carries 2 bits and lasts 2/R seconds, twice as long as a binary one at the same bit rate. Bit pairs 00, 01, 11, 10 go to tones 0, 1, 2, 3 (Gray code), so that neighbouring tones differ in one bit. The tones are spaced by the symbol rate R/2, the smallest orthogonal spacing, so 4 tones take the band of a binary pair.
  - Fixed tones: low 4 / 5.6 / 7.2 / 8.8 kHz, high 8 / 9.6 / 11.2 / 12.8 kHz. Being 1.6 kHz apart, they are orthogonal at every rate (the binary 4 / 8 kHz pair is not, at 1600 bit/s).
  - Hopping: the 6400/R pairs become 6400/R groups of 4 adjacent tones (3200 + 2pR + mR/2 Hz), with the same hopping pattern. The top tone is 15.8, 15.6, 15.2 and 14.4 kHz at 400, 800, 1600 and 3200 bit/s.
- **Frame and header:** the control part is unchanged (binary, 400 bit/s), so every receiver finds every frame and reads the header. A 4-tone header has tones = 3, unused until now, the tone set (0 low, 1 high, 2 hopping) in the two reserved bits, and its rate field read as 400 / 800 / 1600 / 3200. v0 and v1 pages and scripts know no tone set 3 and ignore the frame. Binary frames are unchanged bit for bit, and a v2 receiver rejects a binary header whose reserved bits are not 0. The message and the postamble are sent with 4 tones, each padded with one 0 to an even number of bits (the 127-bit postamble becomes 64 symbols). A frame lasts as long as a binary one at the same bit rate; at 3200 bit/s the 57-character text with the convolutional code takes 0.98 s and the demo image 3.3 s (127-bit preamble).
- **Receiver:** the same steps as v1, with 4 energies per symbol instead of 2. Each one reduces exactly to v1's with 2 tones, and the v2 page and script decode binary frames to the same bits as v1.
  - Calibrated LLR: level and noise of each of the 4 tones from the postamble, then each tone's log-likelihood λ_m = ln I0(A_m E_m/σ_m²) − A_m²/2σ_m². The LLR of each bit is ln Σ e^λ over the 2 tones whose bit is 1, minus the same over the 2 tones whose bit is 0. The scale α, the blind estimate (postamble missing) and the level tracking work the same way on 4 tones. The plain ratio becomes the strongest tone whose bit is 1 against the strongest whose bit is 0.
  - Tone levelling with hopping: a tone is sent in only a quarter of its group's symbols, too few in a short frame for v1's rms. A tone never sent in a group would have been raised to the level of the others, and in simulation this made errors at 400 bit/s even with no noise. Each tone's level is now the rms over the symbols where it is the strongest tone (its group's level if it never is). Binary keeps v1's levelling.
  - Timing and postamble: the contrast is (strongest − second strongest) / sum of the 4 energies; the postamble correlator measures how much the known tone stands out, (4 E_s − ΣE) / (3 ΣE). With 2 tones both are v1's (E1 − E0)/(E1 + E0).
- **Postamble search narrowed to 300 ppm (binary frames too):** the postamble is now looked for within ±(300 ppm of the frame + a quarter of a symbol) of where it should be, instead of ±(0.5% + 1 symbol). In strong echo the correlation can peak a symbol or more late, and the wide search took such peaks for clock offsets of 1000 to 4300 ppm, after which the message was re-read at the wrong instants. The clocks of real devices are usually within ±50 ppm of nominal. With 4 tones at 400 bit/s and DRR −14.6 dB, this raised the texts read from 7 to 12 of 12 (page) and from 5 to 11 (script). Binary frames decoded as with v1 in every other test. Room soundings keep the wide search.

**Results (simulation).** Tone hopping, turbo code rate 1/2, 511-bit preamble, 57-character text, room model with RT60 0.6 s, the script's receiver, 16 trials per point on a 1 dB grid. The page's receiver read the same number of texts within 1 or 2 of 12 on the same recordings.

Limited by noise (DRR +6 dB, white noise, SNR over 3.2–16 kHz), SNR at which half the texts are read (lower is better):

| Bit rate | 2 tones | 4 tones | Gain |
|---|---|---|---|
| 400 | −8.4 dB | −10.3 dB | 1.9 dB |
| 800 | −5.4 dB | −7.8 dB | 2.4 dB |
| 1600 | −2.5 dB | −4.5 dB | 2.0 dB |
| 3200 | (binary not offered) | −1.7 dB | |

Limited by echo (SNR 40 dB), DRR at which half / 9 in 10 of the texts are read (lower is better):

| Bit rate | 2 tones | 4 tones | Gain |
|---|---|---|---|
| 400 | −12.6 / −11.6 dB | −15.5 / −14.6 dB | 2.9 / 3.0 dB |
| 800 | −9.1 / −7.8 dB | −10.5 / −9.7 dB | 1.4 / 1.9 dB |
| 1600 | −7.0 / above −4.8 dB | −7.3 / −6.5 dB | 0.3 / over 1.7 dB |
| 3200 | (binary not offered) | −5.2 / −3.5 dB | |

- **At the same bit rate, 4 tones are better:** about 2 dB in noise and, for 9 texts in 10, 2 to 3 dB in echo. At 1600 bit/s the half-way echo gain is small (0.3 dB) because binary fails gradually: 13 of 16 texts were still read 2 dB above its half-way point. By the rule of 6 dB per doubling of the distance (in the direct field for noise, beyond the critical distance for the DRR), 2 to 3 dB is 25 to 40% more distance. That is a rule of thumb, not a measurement.
- **Why:** each decision gets a symbol twice as long, hence twice the energy, and a non-coherent detector gains more from that than a coherent one. With hopping, a tone group is used again only after twice the time, so its echo has decayed more. The timing tolerance also doubles.
- **4 tones at twice the bit rate** (twice the throughput) need only 0.6 to 0.9 dB more SNR than 2 tones at half that rate, but about 2 dB more DRR. 3200 bit/s with 4 tones needs 0.8 dB more SNR and 1.8 dB more DRR than binary at 1600 bit/s.
- **For class:** use 4 tones whenever the bit rate is fixed; to double the throughput of a binary setting, 4 tones at twice the rate cost about 1 dB in a quiet room but 2 dB in an echoing one.
- In deep echo at 1600 and 3200 bit/s the 4-tone postamble is often not found; the receiver then uses the blind calibration and still reads the texts down to the thresholds above.
- **Tests (simulation):**
  - 168 frames (2 and 4 tones, every rate, tone set, preamble and code) decode with the settings read from the header.
  - v2's binary frames are bit-identical to v1's (84 of 84), and v1 decodes none of the 4-tone frames.
  - The page and `audiocomm_v2.py` give byte-identical WAVs (13 settings) and decode each other's (52 of 52 each way, with both soft outputs).
  - `audiocomm_v2.py` decodes v1's binary WAVs to the same bits as `audiocomm_v1.py` (32 of 32). `audiocomm_v1.py` and the v0 page ignore all 36 four-tone WAVs.
  - In Chromium, the page:
    - switches its rates and tone labels with the modulation, keeps `mod=4` in its address, and hides the row for room soundings;
    - receives a 4-tone 1600 bit/s hopping image (turbo code) and a 4-tone 400 bit/s text on 4 / 5.6 / 7.2 / 8.8 kHz (convolutional code), with no errors.
  - **Not yet tried with real speakers and phones.** A first classroom test could send the same text with 2 and 4 tones at 1600 bit/s with hopping, to the same phones.

## v2: 6400 bit/s

Added to v2 (both `acoustic_modem.html` and `acoustic_modem_v2.html`, and `audiocomm_v2.py`) after Sheng asked how to raise the bit rate further. Phones that had the page open need to reload it to read 6400 bit/s frames.

- **The signal:** 4 tones at 3.2, 6.4, 9.6 and 12.8 kHz, spaced by the symbol rate of 3200 symbols/s, so a symbol lasts 0.31 ms (15 samples at 48 kHz). One group of 4 tones fills the whole 3.2–16 kHz band, so there is nothing left to hop over: it is the hopping rule with H = 6400/R = 1 group. This is as fast as FSK goes in this band, since sending one tone at a time carries at most half a bit per second per hertz. Going further needs many tones at once, each with its own phase (OFDM).
- **Choosing it:** the Rate row shows 6400 with 4 tones. Choosing it selects Hopping, the only tone set that spans the band, and choosing fixed tones moves the rate back to 3200 (`rate=6400` in the page's address). `audiocomm_v2.py tx --mod 4 --tones hop --rate 6400`.
- **Frame and header:** tone set 3, unused until now, with the rate field at 3 (any other rate value with tone set 3 is rejected). v2 receivers made before this change, and v0 and v1, ignore these frames (checked on the page and the script). All other frames are unchanged bit for bit (13 settings compared). The control part still runs at 400 bit/s, so it becomes a large share of a fast frame: the 57-character text with the convolutional code takes 0.81 s (0.98 s at 3200 bit/s) and the demo image 1.99 s (3.33 s), of which 0.65 s is preamble and header (127-bit preamble).
- **Receiver:** the hopping receiver with one group. One change: the postamble counts as found if its correlation reaches a fifth of the preamble's, instead of a half. With every symbol on the same 4 tones, the early echoes (a few milliseconds, many symbols at this rate) cap the postamble's contrast at about 0.39 in the model room even with little echo (0.69 at 3200 bit/s). With the old threshold the postamble was never found at 6400 bit/s. Its side lobes, a symbol or more away, stay near 0.02, so a fifth is still safe. It is now found in 14 to 16 of 16 frames from a DRR of −3 dB up, which gives the clock offset and the postamble calibration; the texts are read about as often as with the blind calibration (within two in 16).

**Results (simulation)**, same conditions as the 4-tone tables above (tone hopping, turbo code rate 1/2, 511-bit preamble, 57-character text, RT60 0.6 s, the script's receiver, 16 trials per point), SNR or DRR at which half / 9 in 10 of the texts are read:

| Bit rate, 4 tones | Noise (DRR +6 dB) | Echo (SNR 40 dB) |
|---|---|---|
| 3200 | −1.7 / −0.7 dB | −5.2 / −3.5 dB |
| 6400 | +3.0 / +3.8 dB | −1.6 / −0.5 dB |
| Cost of 6400 | 4.7 / 4.5 dB | 3.6 / 3.0 dB |

- Each earlier doubling of the 4-tone bit rate cost about 3 dB. In echo 6400 bit/s costs about the same, but in noise it costs 4.5 to 4.7 dB, because with a single group the echo of the earlier symbols always lands on the same 4 tones and adds to the noise. By the 6 dB per doubling of distance rule of thumb, 6400 bit/s reaches 60 to 70% of the distance of 3200 bit/s.
- A clock offset of ±50 ppm and a 44.1 kHz microphone changed these counts by at most two texts in 16.
- **Beyond 6400 bit/s:** a wider band (up to 20 kHz, if phones and laptops reproduce it) would add about 30%. After that comes OFDM, many tones at once with a phase on each. A simulation study in the same room model (`ofdm_study.md`: DQPSK on 546 subcarriers from 3.2 to 16 kHz, 64 ms symbols with a 21 ms cyclic prefix, peaks clipped 6 dB above the rms) found:
  - in echo, about 2.7 times FSK's information rate for 1 dB less range (in that model room, whose early reflections fall inside the guard; in a room whose echo is diffuse from the start the gain is smaller, see v3): 4.3 kbit/s (turbo 1/4) read down to DRR −4.3 dB, against −5.2 dB for FSK at 3200 bit/s (1.6 kbit/s); 8.5 kbit/s (turbo 1/2) down to −0.5 dB, against −1.6 dB for FSK at 6400 bit/s (3.2 kbit/s);
  - in noise, worse: 4.3 kbit/s needs an SNR of +6.0 dB, 3 dB more than FSK at 6400 bit/s, because OFDM plays 3.4 dB quieter for the same peaks and differential detection costs about 3 dB;
  - a phone moving at 0.1 m/s scales time by 300 ppm, which breaks OFDM unless the receiver estimates the scale and resamples the recording.
- **Tests (simulation):**
  - The page and `audiocomm_v2.py` give byte-identical WAVs for 18 settings (4 of them at 6400 bit/s, with no code, the convolutional code and turbo codes 1/2 and 1/3), and each decodes all 72 WAVs (clean at 48 kHz, and at 44.1 kHz with ±50–70 ppm and an SNR of 12 dB).
  - The v2 page from before this change ignores all 16 WAVs at 6400 bit/s and decodes the other 56.
  - In Chromium, the page shows 6400 only with 4 tones, selects Hopping with it, goes back to 3200 with fixed tones and to 1600 with 2 tones, and receives a 6400 bit/s demo image (turbo code) and text (convolutional code) with no errors, postamble found.
  - **Not yet tried with real speakers and phones.** A first test could send the same text at 3200 and 6400 bit/s to phones at several distances.

## v3: OFDM

v3 is `acoustic_modem_v3.html` and `audiocomm_v3.py`, made from v2 (with 6400 bit/s); the v2 files are unchanged and `acoustic_modem.html` still serves v2. FSK frames and room sounding are unchanged bit for bit (28 FSK frames and a sounding frame compared with `audiocomm_v2.py`). To read OFDM frames, phones need the v3 page.

- **The signal:** FSK sends one tone at a time; OFDM sends hundreds at once, each with its own phase. An inverse FFT of N samples at 48 kHz turns one complex number per subcarrier k (at k × 48000 / N Hz, from 3.2 to 16 kHz) into a symbol. Before each symbol its last CP samples are sent again: the cyclic prefix, or **guard interval**. An echo that arrives within the guard only multiplies each subcarrier by a complex gain (the room's frequency response) instead of spilling into the next symbol; echo that arrives later acts as noise.

  | Guard | N | Guard interval | Symbol + guard | Subcarriers | Before the code | Turbo 1/2 | Turbo 1/4 |
  |---|---|---|---|---|---|---|---|
  | long | 2048 | 1024 samples, 21.3 ms | 64 ms | 546, 23.4 Hz apart | 17.1 kbit/s | 8.5 kbit/s | 4.3 kbit/s |
  | short | 1024 | 256 samples, 5.3 ms | 26.7 ms | 273, 46.9 Hz apart | 20.5 kbit/s | 10.2 kbit/s | 5.1 kbit/s |

- **Modulation:** DQPSK. Each subcarrier carries 2 bits per symbol as a phase step from the symbol before: bits b0 b1 give a step of π/4 + m π/2 with m = (b0 xor b1) + 2 b1 (Gray: 00, 10, 11, 01), so the sign of the real part of the step gives b0 and that of its imaginary part b1. Bits 2j and 2j + 1 of each symbol go on subcarrier j. The receiver compares each symbol with the one before on the same subcarrier, so it needs no estimate of the room's response.
- **Level:** the rms is 0.45, half the 0.9 peak of FSK (each subcarrier a cosine of amplitude 0.45 √(2/K)); peaks above 0.9 are clipped, which touches 4.5% of the samples and leaves the clipping noise 19 dB below the signal. With the same peaks, OFDM plays about 3 dB quieter than FSK.
- **Choosing it:** the Modulation row has a third button, "OFDM". It hides the Tones and Rate rows and shows a "Guard interval" row: "21 ms · 17 kbit/s" (long) or "5 ms · 20 kbit/s" (short). In the page's address: `mod=ofdm` and `guard=short`. Script: `audiocomm_v3.py tx --mod ofdm --guard long --code turbo4 --text "..."`. Every code can be used; without one, a single wrong bit spoils a text, so use a turbo code.
- **Frame and header:** a hopping control part as in every hopping frame (warm-up, preamble and header at 400 bit/s), then a known reference symbol, the ns = ⌈(coded bits) / (2 K)⌉ data symbols, and the reference symbol again, each with its guard. The reference symbol has the phases π j² / K on its K subcarriers (Newman phases: its peaks stay 5.4 dB above its rms, below the clipping level). The header uses tone set 3 (as 6400 bit/s) with the rate field 0 for the long guard, 1 for the short one (2 is rejected, 3 stays 6400 bit/s). Older receivers (v0, v1, v2 page and script) ignore these frames. The 57-character text with the convolutional code lasts 0.84 s, of which 0.65 s is the control part (127-bit preamble) and 0.19 s the three OFDM symbols. The demo image lasts 1.29 s with turbo 1/2 (8 data symbols), 1.80 s with turbo 1/4 and 1.13 s with the short guard and turbo 1/2, against 1.99 s with FSK at 6400 bit/s and 3.33 s at 3200 bit/s.
- **Receiver** (`ofdm_receive` in the script, `ofdmReceive` in the page; it starts once the header has been read and the whole OFDM part is in):
  1. It reads the OFDM part on a 48 kHz grid, whatever the microphone's rate, with a windowed-sinc interpolation (Hann window, 16 samples each side).
  2. Timing: the first reference symbol, divided by what was sent, gives the room's impulse response as seen in the band. The FFT windows are placed so that the guard covers the stretch with the most energy: the latest stretch of CP + 1 samples holding at least 99% of the most energy any stretch holds, within ±N/4 of the preamble's timing.
  3. Time scale: a clock offset, or a phone moving at 0.1 m/s (about 300 ppm), stretches the recording by ε. Between the two reference symbols, (ns + 1)(N + CP) samples apart, the phase on subcarrier f then turns by 2π f (ns + 1)(N + CP) ε / 48000. The receiver tries ε from −500 to +500 ppm in steps of 1 ppm, keeps the one that best lines up the subcarriers, refines it with a parabola, and reads the part again on a grid stretched by 1 + ε. The timing is then found again from both reference symbols.
  4. The FFT windows start 8 samples before the guard's end, so that the direct sound never falls outside the guard.
  5. Differential detection: z = Y<sub>i</sub> Y<sub>i−1</sub>* on each subcarrier, LLR(b0) = −√2 Re z / N0 and LLR(b1) = −√2 Im z / N0, with N0 measured from the decisions on each group of 32 subcarriers: the noise plus the echo that arrives after the guard. These LLRs go to the decoder whatever soft output is chosen.
- **What the page shows:** the impulse response as the receiver's FFT windows see it, with the guard shaded (echo inside it is harmless, echo after it is noise); the share of the energy that arrives after the guard; the SNR per subcarrier (lowest, highest and median over the groups of 32); the clock offset; and the timing offset from the preamble's. The WAVs are named `acoustic_text_ofdm_long.wav` and the like.

**Results (simulation)**, 57-character text, 511-bit preamble, the script's receiver, RT60 0.6 s. DRR at which half of the texts are read (SNR 40 dB), in two model rooms:

- room A: the page tests' room (used for v2's tables): direct sound, 4 early reflections between 3 and 14 ms, a diffuse tail, through a speaker filter; its DRR counts the early reflections as echo (16 trials, 1 dB steps);
- room B: the script's `sim` room: direct sound, then a diffuse tail from 1 ms on, no distinct early reflections (8 trials, 2 dB steps).

| | Message bit/s | Room A | Room B |
|---|---|---|---|
| FSK, 4 tones, 3200 bit/s, turbo 1/2 | 1.6 k | −5.2 dB | −5.0 dB |
| FSK, 4 tones, 6400 bit/s, turbo 1/2 | 3.2 k | −1.6 dB | −1.6 dB |
| OFDM, long guard, turbo 1/4 | 4.3 k | −4.4 dB | −1.3 dB |
| OFDM, long guard, turbo 1/2 | 8.5 k | −0.5 dB | +3.1 dB |
| OFDM, short guard, turbo 1/2 | 10.2 k | +1.3 dB | +5.3 dB |

- FSK does not care how the echo is spread in time, so both rooms give it the same numbers. OFDM does: in room A much of the "echo" is the early reflections, which fall inside the guard and do no harm, so at the same DRR OFDM sees 3 to 4 dB less harmful echo than in room B. The gain of OFDM in a reverberant room therefore depends on how much of the room's echo arrives within 21 ms: in room A it carries 2.7 times FSK's information rate for about 1 dB less range (4.3 kbit/s against 1.6, 8.5 against 3.2), in room B 1.3 times FSK 6400's rate at the same range (4.3 kbit/s against 3.2) and 2.7 times FSK 3200's for 3.7 dB less.
- In both rooms the short guard needs about 2 dB more than the long one for 20% more bit rate: in a reverberant room much of the echo arrives after 5.3 ms.
- **Noise** (room A at DRR +6 dB, so no diffuse tail; the noise level is the same for every frame, set relative to a 4-tone 3200 bit/s frame), SNR at which half of the texts are read: FSK 3200 bit/s −1.7 dB, FSK 6400 +3.0 dB, OFDM long guard with turbo 1/4 +5.5 dB, with turbo 1/2 +9.4 dB, short guard with turbo 1/2 +10.7 dB. Where noise limits the link, OFDM with turbo 1/4 needs 2.5 dB more than FSK at 6400 bit/s for 1.3 times its rate, because it plays 3 dB quieter and differential detection costs about 3 dB.
- A time scale of ±300 ppm (a phone moving at 0.1 m/s) and a 44.1 kHz microphone with ±50 ppm gave the same counts as none (long guard, turbo 1/4 at DRR −4 and −3 dB in room A: 13 and 16 of 16; turbo 1/2 and the short guard at +1 dB: 16 and 4 of 16 in each case).
- These match the earlier study (`ofdm_study.md`, which assumed perfect timing and time scale and used room A) within 0.5 dB.

**Sheng's recording of 2026-10-08** (hopping, 1600 bit/s, turbo, demo image, a phone at a distance; old hop band from 1.6 kHz). From the decoded image the frame that was sent can be rebuilt, and compared with what the phone heard:

- **Echo:** in 1.6–14 kHz, 81% of the received energy arrives within 1 ms of the direct sound, 16% between 1 and 5.3 ms, 1% between 5.3 and 21.3 ms and about 1% later. At that spot even the short guard would have held nearly all the echo.
- **Background noise:** 36 dB below the signal in 3.2–16 kHz (it is strong only below 1 kHz).
- **Motion:** the delay changed by less than 0.1 sample per 64 ms (20 ppm) most of the time, and by up to 0.6 sample per 64 ms (about 200 ppm, a few cm/s) for about a second, which the time-scale estimate handles when it lasts the whole OFDM part.
- **What limited the link:** only 55 to 80% of the received energy in 3.2–16 kHz is a linear echo of what was sent (55% in the second where the phone moved). The rest, 1 to 6 dB below the signal, is neither echo nor background noise. Probably distortion in the speaker or the microphone (inferred, not measured): all the tones of that frame were multiples of 1.6 kHz, so harmonics land on other tones, and the second harmonic of the 1.6 kHz tone is its pair, 3.2 kHz, which fits the 43% of 0s read as 1s on that pair. For OFDM this acts as noise, and in noise OFDM needs more signal than FSK: at that spot FSK at 3200 bit/s should reach farther than OFDM. With the signal 36 dB above the background noise, a lower volume may reduce the distortion more than it costs.

**Tests (simulation):**

- The page and `audiocomm_v3.py` give WAVs equal within one LSB for 13 settings (8 with OFDM: both guards, no code, Hamming, convolutional, turbo 1/2, 1/3 and 1/4; 5 with FSK, including 6400 bit/s), and each decodes all 52 WAVs (made by the page and by the script, clean at 48 kHz and at 44.1 kHz with −70 ppm and an SNR of 12 dB) with the same clock offsets and the same number of bits corrected.
- The v2 page and `audiocomm_v2.py` ignore all OFDM frames and decode the FSK ones.
- In Chromium (fake microphone), the page received an OFDM demo image (long guard, turbo 1/2, 773 channel bits corrected, +21 ppm, no pixel errors) and a short-guard text; the Modulation row hides the tones and rates with OFDM and shows the guard.
- The OFDM receiver takes about 0.2 s per frame on the test computer, the turbo decoder about as long again; a phone may take a few times longer.
- **Not yet tried with real speakers and phones.** A first test could send the same text with FSK 3200, FSK 6400 and OFDM (long guard, turbo 1/4 and 1/2) to phones at several distances, at full and at reduced volume.
