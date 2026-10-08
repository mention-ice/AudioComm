# Acoustic modem demo: web page, Python and MATLAB

Each version's page and Python script (`acoustic_modem_v0.html` and `audiocomm.py` for v0, `acoustic_modem_v1.html` and `audiocomm_v1.py` for v1) use the **same frame**: a WAV made or recorded by one is decoded by the other (the transmitted samples were checked to be identical). v1 also reads every v0 frame. The MATLAB scripts still use the earlier frame (fixed settings, 127-bit preamble, no code) and do not read the current frames.

| File | Role |
|---|---|
| `acoustic_modem.html` | **Live demo**, the latest version (now v1, the same file as `acoustic_modem_v1.html`). Transmitter and receiver in one page, for laptops and phones. |
| `acoustic_modem_v0.html` | v0, the original page, frozen. |
| `audiocomm.py` | Python version for the labs (v0). Uses numpy only (`sounddevice` is needed only for live play and record). It does not read hopping or turbo frames: use `audiocomm_v1.py` for those. |
| `audioComm_tx.m`, `audioComm_rx.m`, `audioComm_params.m`, `audioComm_prbs.m` | MATLAB version (`audioComm_v4.m` is untouched) |
| `acoustic_modem_v1.html`, `audiocomm_v1.py` | **v1**: v0 plus room sounding, tone hopping, a turbo code, a fix for cut-off endings, fine timing with a fitted LLR scale, a hop band starting at 3.2 kHz, and Reed–Solomon received but no longer sent (see the last seven sections). v0 = `acoustic_modem_v0.html` and `audiocomm.py`, frozen. |

## Running the live demo

**Laptop only:** open `acoustic_modem.html` in Chrome directly from the folder. The microphone works on local files.

**Phones:** the browser gives a page microphone access only over https, so the page must be hosted. The simplest way is GitHub Pages:

1. Create a public repository.
2. Upload the file as `index.html`.
3. Go to Settings → Pages → "Deploy from branch: main".

You get a link like `https://<user>.github.io/<repo>/`. Show it as a QR code in the slides.

This demo is published at `https://mention-ice.github.io/AudioComm/` (repository `mention-ice/AudioComm`): `acoustic_modem.html` is always the latest version, and `acoustic_modem_v0.html` and `acoustic_modem_v1.html` keep each version. It was first published at `sheng-yang-cs.github.io/AudioComm`; GitHub does not redirect a Pages address when a repository moves, so links and QR codes printed before the move no longer work.

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

v1 reads every v0 frame and adds a third frame type that measures the room. v1 link: `https://mention-ice.github.io/AudioComm/acoustic_modem_v1.html`, also served as `acoustic_modem.html`, the latest version. A v1 page's QR code points to its own online copy, also when the laptop runs it from a local file.

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
