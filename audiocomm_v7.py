"""audiocomm_v7.py - acoustic modem (FSK with 2 or 4 tones, OFDM with DQPSK, QPSK, 8PSK, 16, 64 or 256QAM, or v5's single
carrier with QPSK to 256QAM and turbo equalization) and room sounding, Python version for the labs. v6: a short
control part (0.19 s instead of 0.65 s) for the OFDM and single-carrier frames. v7: LDPC codes (5G NR) at rates 1/4 to
3/4; the Hamming code is gone. v7.1: pictures as JPEG files, with a CRC-32. v7.2: the single carrier tries a
frame a second time when it fails, with a longer room (2048 taps, 213 ms) learned from the frame itself. v7.4: the
single carrier's wide band, 13 714 symbols/s from 2.3 to 20.5 kHz (43% more than 3.2 to 16 kHz). v7.5: 1024QAM on
the single carrier (10 bits per symbol), for phones close together in a quiet room. v7.6: LDPC rate 5/6.

© 2026 Sheng Yang. All rights reserved.

Same frames as acoustic_modem_v7.html, so a WAV produced or recorded by either can be
decoded by the other. The receiver needs no settings: they travel in the header.

  control part, at 400 bit/s, 4 / 8 kHz or tone hopping (two modes known to every receiver in advance):
    0.1 s warm-up | MLS preamble (127, 255 or 511 bits) | header (sent 1, 2 or 4 times)
  message part, in the tones and rate announced by the header:
    coded, scrambled data | MLS postamble (same sequence as the preamble)
  or, with OFDM (v3; its control part always hops):
    reference symbol | data symbols (coded, scrambled data) | reference symbol, each preceded by its guard interval
  or, with a single carrier (v5; its control part always hops):
    0.1 s silence | TA (1280 QPSK symbols) TC | data | TC | data | ... | TC, the training TC of 640, 1280 or 2560
    symbols (short, medium or long: the longest echo it lets the receiver learn, 20, 47 or 100 ms after the first strong arrival;
    v7.4, the wide band: TA 1792, TC 896, 1792 or 3584, about the same durations)
  v6, the short control part of the OFDM and single-carrier frames (their default; the 400 bit/s one remains, for range
    and for v5 receivers, which do not read the short one):
    P P P | H, 0.192 s: three periods of P, 2048 samples at 48 kHz (OFDM's reference symbol for N = 2048 with its
    phases negated, 3.2 to 16 kHz), then the header symbol H (a guard of 1024 samples, then 2048) carrying the same
    92 coded header bits by BPSK on P's subcarriers: subcarrier j is P's value times 1 - 2 c, c = bit (37 j) mod 92

  header : type 2 | a 11 | b 11 | code 2 | rate 2 | tones 2 | 0 0  (MSB first), then CRC-8 (x^8 + x^2 + x + 1),
           all through the convolutional code;  type 0 text: a = number of characters, b = 0;
           type 1 image: a = rows, b = columns;  type 2 room sounding: a = MLS order (14, 15, 16), b = periods (1..8);
           rate: 200 400 800 1600;  tones: low (4/8 kHz), high (8/16 kHz), hop (tone hopping, 3.2 to 16 kHz)
           v2: tones 3 = 4 tones per symbol (4-FSK), the tone set (0 low, 1 high, 2 hop) in the last two bits, and the
           rate field read in 400 800 1600 3200; receivers older than v2 find no tone set 3 and ignore the frame;
           tone set 3 with rate field 3 = 6400 bit/s, one group of 4 tones over the whole band (3.2 / 6.4 / 9.6 /
           12.8 kHz, tone hopping with a single group); receivers older than 6400 bit/s ignore it
           v3: tone set 3 with rate field 0 or 1 = OFDM with the long or the short guard (see OFDM_MODES);
           receivers older than v3 ignore it
           v4: an OFDM frame with a coherent mapping xors its CRC with 0x33 (QPSK), 0x55 (8PSK), 0x66 (16QAM), 0x3C (64QAM)
           or 0x0F (256QAM), on top of the 0xFF of codes 4 to 6; receivers older than v4 ignore it (64QAM and 256QAM:
           receivers older than their addition)
           v5: tone set 3 with rate field 2 = single carrier (see SC_FIELD), its mapping (QPSK, 16QAM, 64QAM or
           256QAM) in the CRC mask as for OFDM, turbo codes only; receivers older than v5 ignore it; the training
           length in the top bits of a and b (unused: at most 1023 characters, images at most 512 x 512): 0 0 medium
           (1280), 0 1 short (640), 1 0 long (2560)
           v7: type 3 = a message coded with an LDPC code (5G NR, LDPC_CODES): a text if b = 0 (a = number of
           characters), else an image (a rows, b columns); the code field (+ 4 with the CRC inverted, as for codes 4
           to 6) gives the rate: 0 1/4, 1 1/3, 2 1/2, 3 2/3, 4 3/4; receivers older than v7 ignore it
           v7.1: an image with a > 512 is a JPEG file of n = 512 (a - 513) + b bytes (b from 1 to 512, n at most
           JPEG_MAX), with any code (type 1, or type 3 with LDPC); receivers older than v7.1 reject images over 512 rows
           v7.3: type 3 with the code field 5, 6 or 7 (1, 2, 3 with the CRC inverted) = a polar code (5G NR, POLAR_CODES)
           of rate 1/3, 1/2, 2/3; receivers older than v7.3 ignore it
           v7.4: type 2 on the single carrier = the wide band (SC_BANDS), read as type 3 (an LDPC or a polar code, a
           text if b = 0, else an image); receivers older than v7.4 ignore it
           v7.5: on the single carrier, 8PSK's mask (0x55) = 1024QAM (SC_SLOT); receivers older than v7.5 ignore it
           v7.6: type 0 or 1 with the code field 7 (3 with the CRC inverted) = LDPC 5/6 (LDPC56), on the single carrier's
           wide band type 0 or 1 with the code field 2 (convolutional, which the single carrier never takes); receivers
           older than v7.6 ignore it
  code   : 0 none, 1 Hamming (7,4) + interleaver (until v6; v7 neither sends nor reads it), 2 convolutional (171,133) K = 7 + interleaver,
           3 Reed-Solomon (255,223) over bytes (received only: v1 no longer sends it), 4 turbo (13, 15 octal,
           rate 1/2, as in 4G) + interleaver, 5 and 6 the same turbo code at rates 1/3 and 1/4 (1/4 adds cdma2000's second parity);
           codes 4, 5, 6 are sent as 0, 1, 2 with the CRC inverted (as 4G masks a CRC), so receivers that predate
           them ignore them
  data   : text = 8 bits per character (Latin-1, MSB first);
           image = pixels in column-major order (MATLAB img(:)), 1 = white;
           JPEG (v7.1) = the file's bytes, then its CRC-32 (as zlib.crc32, 4 bytes), each byte MSB first
  scrambling: XOR with the LFSR x^15 + x^14 + 1 sequence, all-ones seed
  bit 0 -> sin(2 pi F0 t), bit 1 -> sin(2 pi F1 t), one symbol = 1/rate seconds
  4 tones (v2): 2 bits per symbol, one symbol = 2/rate seconds; bits (b0, b1) -> the tone m whose Gray code
           m xor (m >> 1) is b0 b1, so that neighbouring tones differ in one bit; tones low 4 / 5.6 / 7.2 / 8.8 kHz,
           high 8 / 9.6 / 11.2 / 12.8 kHz, or with hopping H = 6400 / R groups of 4 tones spaced R / 2
           (3200 + 2pR + mR/2 Hz). The message and the postamble are each padded with a 0 to a whole number of symbols.
  tone hopping: at R bit/s, H = 6400 / R pairs of adjacent tones (3200 + 2pR, 3200 + (2p + 1)R Hz); symbol k of a
           segment uses pair (H/2 - 1) k mod H (k mod 4 when H = 4), counting from the first preamble bit, from the
           first message bit through the postamble, or from the first bit of a sounding's postamble. A hopping
           message has a hopping control part (at 400 bit/s), so that far receivers find it in a reverberant room.
  room sounding (type 2): after the header, an MLS of 2^order - 1 samples of +-0.35 at 48 kHz, played periods + 1
           times (the first one fills the room), then the postamble in the control mode. The receiver correlates
           each period with the MLS to get the impulse response h(t) of speaker + room + microphone, then derives
           RT60 (Schroeder, T20 or T10), the direct-to-reverberant ratio, C50, the rms delay spread, the frequency
           response and the SNR per frequency (against the noise heard before the frame).
  OFDM (v3): an inverse FFT of N points at 48 kHz, subcarriers k 48000 / N from 3.2 to 16 kHz, each preceded by a
           cyclic prefix of CP samples (N 2048 and CP 1024, or N 1024 and CP 256); a reference symbol with phases
           pi j^2 / K, the data symbols, then the reference symbol again; DQPSK: bits b0 b1 of subcarrier j (in increasing frequency, 2K bits per symbol, the last
           symbol padded with 0) turn its phase by pi/4 + m pi/2, m = (b0 xor b1) + 2 b1; rms 6 dB below the 0.9 peak,
           peaks clipped at 0.9. The receiver times its windows and measures the time scale (clock offset, motion)
           on the two reference symbols.
  OFDM, coherent (v4): subcarrier j of data symbol i (from 1) is a pilot if (j - 3 i) mod 12 = 0: it carries the
           reference symbol's value; the other subcarriers, in increasing frequency, symbol after symbol, carry m bits
           each (MSB first, Gray: neighbouring points differ in 1 bit): QPSK (1 - 2 b0 + j (1 - 2 b1)) / sqrt 2;
           8PSK exp(j 2 pi p / 8) with bits p xor (p >> 1); 16QAM (A[u] + j A[v]) / sqrt 10, A = (-3, -1, 1, 3), with
           bits g(u) g(v), g = (00, 01, 11, 10); 64QAM and 256QAM the same with A = (-7, -5, ..., 7) or (-15, ..., 15),
           3 or 4 Gray bits per axis, divided by sqrt 42 or sqrt 170; v7.5, the single carrier's 1024QAM: A = (-31, ..., 31),
           5 Gray bits per axis, divided by sqrt 682. The last data symbol is filled up with scrambled zeros (the scrambling
           sequence goes on), then the reference symbol ends the part. The receiver estimates each subcarrier's gain
           and phase from the reference symbols and the pilots (see ofdm_receive).
  single carrier (v5): 9600 symbols/s at 9.6 kHz, root-raised-cosine pulses of roll-off 1/3 (3.2 to 16 kHz), rms 0.45;
           v7.4, the wide band: 96000 / 7 = 13 714 symbols/s at 11.4 kHz (2.3 to 20.5 kHz), 7 samples per symbol at 96 kHz,
           sent as every other sample (3.5 samples per symbol at 48 kHz); training 10/7 longer in symbols;
           points as OFDM's coherent mappings; the training blocks are QPSK exp(j pi/2 (q + 1/2)), q the top 2 bits of
           xorshift32 (seeds 0x5C0A for TA, 0x5C0C for TC); the first block is TA TC, every later one TC, after at most
           4800 data symbols; channel bits: the code's output padded with 0, scrambled. The receiver
           learns the room (512 taps at the symbol rate) on the training blocks and decodes by turbo equalization:
           soft SMAP (Ramos, Senecal, Yang) or linear MMSE, and the turbo (v7: or LDPC) decoder, exchanging beliefs (see sc_receive).
  The MATLAB scripts audioComm_tx.m / audioComm_rx.m use the earlier frame and do not read these frames.

Usage
  python audiocomm_v7.py tx --text "Hello" -o hello.wav     # write a WAV (add --play to play it)
  python audiocomm_v7.py tx --image greatwave_80x54.png -o img.wav   # the page's demo image (Hokusai's Great Wave); greatwave_160x108.png: the large one
  python audiocomm_v7.py tx --image greatwave_320x216.png --mod ofdm --map 256qam --code turbo -o big.wav   # the 4x one
  python audiocomm_v7.py tx --text "Hello" --mod sc --map 64qam --code ldpc --wide -o wide.wav   # v7.4: 2.3 to 20.5 kHz
  python audiocomm_v7.py tx --image greatwave_160x108.png --mod sc --map 1024qam --code ldpc23 -o q1024.wav   # v7.5: 1024QAM
  python audiocomm_v7.py tx --image greatwave_160x108.png --mod sc --map 256qam --code ldpc56 -o r56.wav   # v7.6: LDPC 5/6
  (the earlier demo image: --image ../../Lab/bartS.png --step 2, without --step, or with --up 2 for the large and 4x ones)
  python audiocomm_v7.py tx --jpeg greatwave_160x108_q50.jpg --mod sc --code ldpc -o jpg.wav   # a JPEG file as it is (v7.1)
  python audiocomm_v7.py tx --image photo.jpg --quality 50 --width 160 --mod sc --code ldpc -o photo.wav   # a picture made into a JPEG (v7.1)
  python audiocomm_v7.py tx --text "Hello" --code conv -o coded.wav   # codes: none conv turbo turbo3 turbo4 ldpc4 ldpc3 ldpc ldpc23 ldpc34 ldpc56 polar3 polar polar23
  python audiocomm_v7.py rx recording.wav                    # decode a recording
  python audiocomm_v7.py rx recording.wav --soft ratio       # soft output ln(E1/E0) instead of the calibrated LLR
  python audiocomm_v7.py tx --text "Hello" --code conv --rate 200 --pre 511 -o far.wav   # long range
  python audiocomm_v7.py tx --text "Hello" --tones high --rate 800 -o fast.wav
  python audiocomm_v7.py tx --text "Hello" --tones hop --code conv -o hop.wav   # tone hopping, for reverberant rooms
  python audiocomm_v7.py tx --text "Hello" --tones hop --code turbo -o turbo.wav   # turbo code (turbo3, turbo4: rates 1/3, 1/4)
  python audiocomm_v7.py tx --text "Hello" --tones hop --code ldpc -o ldpc.wav   # LDPC code (v7; ldpc4 ldpc3 ldpc23 ldpc34 ldpc56: 1/4 1/3 2/3 3/4 5/6)
  python audiocomm_v7.py tx --text "Hello" --mod sc --code polar -o polar.wav   # polar code (v7.3; polar3 polar23: 1/3 2/3)
  python audiocomm_v7.py tx --text "Hello" --mod 4 --tones hop --rate 3200 --code turbo -o fast4.wav   # 4 tones, 2 bits per symbol
  python audiocomm_v7.py tx --text "Hello" --mod 4 --tones hop --rate 6400 --code turbo -o fast6.wav   # the whole band, no hopping
  python audiocomm_v7.py tx --image greatwave_160x108.png --mod ofdm --code turbo -o ofdm.wav   # OFDM, 17 kbit/s before the code
  python audiocomm_v7.py tx --text "Hello" --mod ofdm --guard short --code turbo -o ofdm5.wav   # OFDM, 5 ms guard
  python audiocomm_v7.py tx --image greatwave_160x108.png --mod ofdm --map 16qam --code turbo -o qam.wav   # 16QAM with pilots, 31 kbit/s
  python audiocomm_v7.py tx --image greatwave_160x108.png --mod ofdm --map 256qam --code turbo -o qam256.wav   # 256QAM, 63 kbit/s
  python audiocomm_v7.py tx --image greatwave_160x108.png --mod sc --map 256qam --code turbo -o sc256.wav   # single carrier (v5)
  python audiocomm_v7.py tx --text "Hello" --mod ofdm --map qpsk --code turbo --control long -o v5ctrl.wav   # v5's control part
  python audiocomm_v7.py rx sc256.wav --det lmmse            # the same with the linear MMSE detector instead of SMAP
  python audiocomm_v7.py rx sc256.wav --est mmse             # the room estimated by MMSE instead of least squares
  python audiocomm_v7.py rx --listen 20                      # record 20 s from the mic, then decode
  python audiocomm_v7.py tx --sound -o sound.wav             # room sounding: MLS of 0.68 s, 4 periods averaged
  python audiocomm_v7.py tx --sound --order 16 --periods 8 --pre 511 -o sound.wav   # large or reverberant room
  python audiocomm_v7.py rx recorded_sound.wav --ir room_ir.wav --plot   # room measures; saves the impulse response
  python audiocomm_v7.py sim hello.wav --ir room_ir.wav --snr 30 -o heard.wav   # play a frame through a measured room
  python audiocomm_v7.py sim hello.wav --rt60 0.6 --drr 0 --snr 30 --ppm 50 -o heard.wav   # or through a model room
  python audiocomm_v7.py sim hop.wav --rt60 0.6 --drr -6 -o far.wav && python audiocomm_v7.py rx far.wav   # beyond the critical distance
--play / --listen need the `sounddevice` package, --plot needs matplotlib.
"""
import argparse
import math
import os
import wave
import zlib
import numpy as np

PILOT = np.array([int(c) for c in
    "0000000111111011111001111010111000011011101001100010101100000101111000111011011001001010010000100111001011010001000110011010101"])
NP, HBITS = 127, 40          # HBITS: header, 32 bits + CRC-8
NHC = 2 * (HBITS + 6)        # header after the convolutional code


def mls(deg, taps):
    """Maximum-length sequence: Fibonacci LFSR, all-ones seed; taps = exponents of the feedback polynomial."""
    s, out = (1 << deg) - 1, []
    for _ in range((1 << deg) - 1):
        out.append(s & 1)
        fb = 0
        for t in taps:
            fb ^= (s >> (deg - t)) & 1
        s = (s >> 1) | (fb << (deg - 1))
    return np.array(out)


# Longer preambles for longer distances; the header is repeated more so that it is read as far away
PILOTS = {127: PILOT, 255: mls(8, (8, 6, 5, 4)), 511: mls(9, (9, 5))}
HREPS = {127: 1, 255: 2, 511: 4}
CTRL = dict(tones="low", rate=400)               # mode of the control part (preamble and header)
RATES, TONE_SETS = [200, 400, 800, 1600], ["low", "high", "hop"]
RATES4 = [400, 800, 1600, 3200]                  # v2, 4 tones per symbol
# 6400 bit/s, 4 tones only: one group of 4 tones spaced 3.2 kHz fills the band 3.2-12.8 kHz, so there is nothing left to
# hop over (H = 1). It is as fast as FSK goes in this band: one tone at a time carries at most 1/2 bit/s per Hz.
RATE_WIDE = 6400
TONES = {"low": (4000, 8000), "high": (8000, 16000)}
# v2: 4-FSK. With 4 tones a symbol carries 2 bits, so at the same bit rate it lasts twice as long; tones spaced by a
# multiple of the symbol rate stay orthogonal (here 1.6 kHz, at every rate up to 3200 bit/s).
TONES4 = {"low": (4000, 5600, 7200, 8800), "high": (8000, 9600, 11200, 12800)}
HOP_F0, HOP_BAND = 3200, 12800                  # tone hopping: pairs from 3.2 kHz, over 12.8 kHz
# The postamble is looked for only where the clocks of two devices can put it: they differ by far less than 300 ppm
# (0.03 %), plus a quarter of a symbol. In a strong echo the contrast can peak a symbol or more after the true
# instant, and with a wider search such a peak was taken for a clock offset of several thousand ppm.
MAX_PPM = 300
FS_TX = 48000
SOUND_TAPS = {14: (14, 13, 12, 2), 15: (15, 14), 16: (16, 15, 13, 4)}   # room sounding: MLS orders, primitive polynomials
SOUND_AMP = 0.35


def prbs(n):
    """Scrambling sequence: LFSR x^15 + x^14 + 1, all-ones seed."""
    s, k = 0x7FFF, np.zeros(n, dtype=int)
    for i in range(n):
        fb = ((s >> 14) ^ (s >> 13)) & 1
        s = ((s << 1) | fb) & 0x7FFF
        k[i] = fb
    return k


def msb(v, nb):
    return [(v >> i) & 1 for i in range(nb - 1, -1, -1)]


def tone_set(tones, M=2):
    return (TONES4 if M == 4 else TONES)[tones]


def symbol_bits(m, K):
    """Gray mapping: tone m carries the K bits of the Gray code of m, m xor (m >> 1), most significant first, so that
    neighbouring tones, the likeliest mix-up, differ in one bit only. With 2 tones, tone m is bit m."""
    return msb(m ^ (m >> 1), K)


def to_symbols(bits, M):
    """bits -> tones, K = log2 M bits per symbol, the last symbol padded with 0."""
    K = int(math.log2(M))
    bits = [int(b) for b in bits] + [0] * (-len(bits) % K)
    out = []
    for i in range(0, len(bits), K):
        v = int("".join(map(str, bits[i:i + K])), 2)
        m, sh = v, v >> 1
        while sh:                                  # the tone whose Gray code is v
            m ^= sh
            sh >>= 1
        out.append(m)
    return np.array(out, dtype=int)


def bit_sets(M):
    """For each bit i of a symbol: (the tones whose bit i is 0, the tones whose bit i is 1)."""
    K = int(math.log2(M))
    lab = [symbol_bits(m, K) for m in range(M)]
    return [tuple([m for m in range(M) if lab[m][i] == b] for b in (0, 1)) for i in range(K)]


# ---------------------------------------------------------------- channel codes
# Same codes as acoustic_modem.html; the code is announced in the header, so the receiver finds it there.
CODES = ["none", None, "conv", "rs", "turbo", "turbo3", "turbo4"]   # 1: Hamming (7,4) until v6
TURBO_R = {"turbo": 2, "turbo3": 3, "turbo4": 4}   # turbo codes: rate 1/r
RS_K, RS_P = 223, 32


def coded_length(code, n0):
    if code in POLAR_R:
        return polar_params(n0, POLAR_R[code])["G"]
    if code in LDPC_R:
        num, den = LDPC_R[code]
        return -(-n0 * den // num)
    if code in ("conv", "turbo"):
        return 2 * (n0 + 6)
    if code == "turbo3":
        return 3 * n0 + 12
    if code == "turbo4":
        return 4 * n0 + 18
    if code == "rs":
        K = math.ceil(n0 / 8)
        return 8 * (K + RS_P * math.ceil(K / RS_K))
    return n0


def ilv_perm(n):
    """Rectangular interleaver: write row by row into ceil(sqrt(n)) columns, read column by column."""
    C = max(1, math.ceil(math.sqrt(n)))
    R = math.ceil(n / C)
    idx = np.arange(R * C).reshape(R, C).T.ravel()
    return idx[idx < n]


def interleave(x):
    return np.asarray(x)[ilv_perm(len(x))]


def deinterleave(y):
    y = np.asarray(y)
    x = np.empty_like(y)
    x[ilv_perm(len(y))] = y
    return x


# Convolutional code, rate 1/2, K = 7, generators 171 and 133 (octal), 6 zero tail bits
CC_G = (0o171, 0o133)


def _parity(v):
    return bin(v).count("1") & 1


def conv_encode(d):
    s, out = 0, []
    for b in list(d) + [0] * 6:
        reg = (int(b) << 6) | s
        out += [_parity(reg & CC_G[0]), _parity(reg & CC_G[1])]
        s = reg >> 1
    return np.array(out)


def conv_decode(llr, n0, lim=4):
    """Soft-decision Viterbi, 64 states. llr > 0 means 1, clipped to +-lim."""
    T = len(llr) // 2
    L = np.clip(np.asarray(llr, dtype=float), -lim, lim).reshape(T, 2)
    o = np.array([[2 * _parity(r & g) - 1 for g in CC_G] for r in range(128)])
    ns = np.arange(64)
    sa = (ns & 31) << 1                       # the two predecessors of state ns are sa and sa + 1
    ra = ((ns >> 5) << 6) | sa
    rb = ra | 1
    pm = np.full(64, -1e18)
    pm[0] = 0
    dec = np.zeros((T, 64), dtype=np.uint8)
    for t in range(T):
        ma, mb = pm[sa] + o[ra] @ L[t], pm[sa | 1] + o[rb] @ L[t]
        dec[t] = mb > ma
        pm = np.maximum(ma, mb)
    out, s = np.zeros(T, dtype=int), 0
    for t in range(T - 1, -1, -1):
        out[t] = s >> 5
        s = ((s & 31) << 1) | int(dec[t, s])
    return out[:n0]


# Turbo code (as in 4G): two 8-state recursive systematic encoders, feedback 1 + D^2 + D^3 and parity 1 + D + D^3
# (13 and 15 octal); the second one encodes the message in the order of a spread pseudo-random interleaver. Each encoder
# returns to state 0 with 3 tail steps. Rate 1/2: parity bits taken alternately from each encoder,
# k bits -> [k message bits | k parity bits | 12 tail bits] = 2k + 12 bits, as long as with the convolutional code.
# Rate 1/3 sends every parity bit of both encoders (as cdma2000, and 4G before rate matching), [k | k | k | 12 tail bits];
# 1/4 also sends cdma2000's second parity, 1 + D + D^2 + D^3 (17 octal), taken alternately from each encoder,
# [k | k | k | k | 18 tail bits]. The three rates are nested: each one's bits include those of the higher rates.
# The decoder runs two log-MAP (BCJR) decoders in turn, each passing the other what it learnt (extrinsic LLRs).
TB_NEXT, TB_PAR, TB_PAR2 = (np.zeros((8, 2), dtype=int) for _ in range(3))   # next state, parity bits 15, 17 (state, input)
for _s in range(8):
    for _u in (0, 1):
        _a = _u ^ ((_s >> 1) & 1) ^ (_s & 1)
        TB_PAR[_s, _u], TB_NEXT[_s, _u] = _a ^ (_s >> 2) ^ (_s & 1), (_a << 2) | (_s >> 1)
        TB_PAR2[_s, _u] = _a ^ (_s >> 2) ^ ((_s >> 1) & 1) ^ (_s & 1)
_TB_PERMS = {}


def turbo_perm(k):
    """Fisher-Yates shuffle driven by xorshift32 seeded with k, then spread (S-random): positions closer than S stay
    more than S apart where possible, which makes low-weight codewords rare. The page computes the same."""
    if k not in _TB_PERMS:
        M, p = 0xFFFFFFFF, list(range(k))
        x = (0x9E3779B9 ^ k) & M or 1
        for i in range(k - 1, 0, -1):
            x = (x ^ (x << 13)) & M
            x ^= x >> 17
            x = (x ^ (x << 5)) & M
            j = x % (i + 1)
            p[i], p[j] = p[j], p[i]
        S = min(16, int(math.sqrt(k / 2)))
        for i in range(k):
            for m in range(i, k):
                if all(abs(p[m] - p[j]) > S for j in range(max(0, i - S), i)):
                    p[i], p[m] = p[m], p[i]
                    break
        _TB_PERMS[k] = np.array(p, dtype=int)
    return _TB_PERMS[k]


def _rsc(u):
    """Parity bits (15 and 17 octal), then the 3 tail inputs that bring the register back to 0 and their parity bits."""
    s, p, q, tu, tp, tq = 0, [], [], [], [], []
    for b in u:
        p.append(TB_PAR[s, b])
        q.append(TB_PAR2[s, b])
        s = TB_NEXT[s, b]
    for _ in range(3):
        b = ((s >> 1) & 1) ^ (s & 1)
        tu.append(b)
        tp.append(TB_PAR[s, b])
        tq.append(TB_PAR2[s, b])
        s = TB_NEXT[s, b]
    return np.array(p, dtype=int), np.array(q, dtype=int), tu, tp, tq


def turbo_encode(d, r=2):
    """Rate 1/r, r = 2, 3 or 4."""
    d = np.asarray(d, dtype=int)
    p1, q1, tu1, tp1, tq1 = _rsc(d)
    p2, q2, tu2, tp2, tq2 = _rsc(d[turbo_perm(len(d))])
    even = np.arange(len(d)) % 2 == 0
    if r == 2:
        parts = [d, np.where(even, p1, p2), tu1, tp1, tu2, tp2]
    elif r == 3:
        parts = [d, p1, p2, tu1, tp1, tu2, tp2]
    else:
        parts = [d, p1, p2, np.where(even, q1, q2), tu1, tp1, tq1, tu2, tp2, tq2]
    return np.concatenate(parts).astype(int)


def _bcjr(ls, la, lp, lq=None):
    """Log-MAP over the steps of ls, from state 0 to state 0; LLRs are ln P(0)/P(1) here. A-posteriori LLR of each input.
    lq: LLRs of the second parity bits (17 octal), at rate 1/4 only."""
    T = len(ls)
    g = 0.5 * ((ls + la)[:, None, None] * np.array([1, -1]) + lp[:, None, None] * (1 - 2 * TB_PAR))   # (T, 8, 2)
    if lq is not None:
        g += 0.5 * lq[:, None, None] * (1 - 2 * TB_PAR2)
    pred = [[(s, u) for s in range(8) for u in (0, 1) if TB_NEXT[s, u] == n] for n in range(8)]
    ps, pu = np.array([[q[0][0], q[1][0]] for q in pred]), np.array([[q[0][1], q[1][1]] for q in pred])
    A = np.full((T + 1, 8), -1e9)
    A[0, 0] = 0
    for t in range(T):
        a = A[t][ps] + g[t][ps, pu]
        A[t + 1] = np.logaddexp(a[:, 0], a[:, 1])
        A[t + 1] -= A[t + 1].max()
    B, out = np.full(8, -1e9), np.empty(T)
    B[0] = 0
    for t in range(T - 1, -1, -1):
        gb = g[t] + B[TB_NEXT]
        m = A[t][:, None] + gb
        out[t] = np.logaddexp.reduce(m[:, 0]) - np.logaddexp.reduce(m[:, 1])
        B = np.logaddexp(gb[:, 0], gb[:, 1])
        B -= B.max()
    return out


def turbo_decode(llr, n0, lim=np.inf, iters=8, r=2):
    """llr > 0 means 1; rate 1/r. Returns the n0 message bits and the number of iterations run (it stops when the
    decisions no longer change)."""
    k, pi = n0, turbo_perm(n0)
    L = -np.clip(np.clip(np.asarray(llr, dtype=float), -lim, lim), -30, 30)
    part = lambda j: L[j * k:(j + 1) * k]
    ls, tl, z3, t2 = L[:k], L[r * k:], np.zeros(3), 9 if r == 4 else 6   # t2: tail of encoder 2
    odd = np.arange(k) % 2 == 1
    if r == 2:
        lp1, lp2 = np.where(odd, 0, part(1)), np.where(odd, part(1), 0)
    else:
        lp1, lp2 = part(1), part(2)
    lq1 = np.r_[np.where(odd, 0, part(3)), tl[6:9]] if r == 4 else None
    lq2 = np.r_[np.where(odd, part(3), 0), tl[15:18]] if r == 4 else None
    ls2 = ls[pi]
    le2, x = np.zeros(k), None
    for it in range(1, iters + 1):
        a1 = _bcjr(np.r_[ls, tl[0:3]], np.r_[le2, z3], np.r_[lp1, tl[3:6]], lq1)[:k]
        le1 = a1 - ls - le2
        a2 = _bcjr(np.r_[ls2, tl[t2:t2 + 3]], np.r_[le1[pi], z3], np.r_[lp2, tl[t2 + 3:t2 + 6]], lq2)[:k]
        le2, xn = np.empty(k), np.empty(k, dtype=int)
        le2[pi], xn[pi] = a2 - ls2 - le1[pi], a2 < 0
        same = x is not None and np.array_equal(xn, x)
        x = xn
        if same:
            break
    return x, it


# Reed-Solomon over GF(256) (primitive polynomial 0x11d), 32 parity bytes: up to 16 wrong bytes per codeword
GF_EXP, GF_LOG = [0] * 512, [0] * 256
_x = 1
for _i in range(255):
    GF_EXP[_i], GF_LOG[_x] = _x, _i
    _x <<= 1
    if _x & 256:
        _x ^= 0x11D
for _i in range(255, 512):
    GF_EXP[_i] = GF_EXP[_i - 255]


def gmul(a, b):
    return GF_EXP[GF_LOG[a] + GF_LOG[b]] if a and b else 0


def gdiv(a, b):
    return GF_EXP[(GF_LOG[a] + 255 - GF_LOG[b]) % 255] if a else 0


def gpow(e):
    return GF_EXP[e % 255]


RS_GEN = [1]
for _i in range(RS_P):
    _h = RS_GEN + [0]
    for _j, _c in enumerate(RS_GEN):
        _h[_j + 1] ^= gmul(_c, gpow(_i))
    RS_GEN = _h


def rs_encode_block(msg):
    """Systematic: the message bytes, then 32 parity bytes."""
    r = [0] * RS_P
    for m in msg:
        f = m ^ r[0]
        r = r[1:] + [0]
        if f:
            for j in range(RS_P):
                r[j] ^= gmul(RS_GEN[j + 1], f)
    return list(msg) + r


def _poly_eval(p, x):
    acc = 0
    for c in p:
        acc = gmul(acc, x) ^ c
    return acc


def rs_decode_block(cw):
    """Berlekamp-Massey + Chien + Forney. Returns (message, number of corrected bytes or None if it failed)."""
    n, k = len(cw), len(cw) - RS_P
    S = [_poly_eval(cw, gpow(j)) for j in range(RS_P)]
    if not any(S):
        return list(cw[:k]), 0
    C, B, L, m, b = [1], [1], 0, 1, 1
    for i in range(RS_P):
        d = S[i]
        for j in range(1, L + 1):
            d ^= gmul(C[j] if j < len(C) else 0, S[i - j])
        if not d:
            m += 1
            continue
        T, coef = C[:], gdiv(d, b)
        C += [0] * max(0, len(B) + m - len(C))
        for j, v in enumerate(B):
            C[j + m] ^= gmul(coef, v)
        if 2 * L <= i:
            L, B, b, m = i + 1 - L, T, d, 1
        else:
            m += 1
    lam = C[:L + 1]
    omega = [0] * RS_P
    for i in range(RS_P):
        for j in range(min(L, i) + 1):
            omega[i] ^= gmul(lam[j], S[i - j])
    out, found = list(cw), 0
    for p in range(n):
        val = 0
        for j, c in enumerate(lam):
            val ^= gmul(c, gpow(-p * j))
        if val:
            continue
        om = 0
        for i in range(RS_P):
            om ^= gmul(omega[i], gpow(-p * i))
        dl = 0
        for j in range(1, L + 1, 2):
            dl ^= gmul(lam[j], gpow(-p * (j - 1)))
        if not dl:
            return list(cw[:k]), None
        out[n - 1 - p] ^= gmul(gpow(p), gdiv(om, dl))
        found += 1
    if found != L or any(_poly_eval(out, gpow(j)) for j in range(RS_P)):
        return list(cw[:k]), None
    return out[:k], L


def _rs_layout(K):
    B = math.ceil(K / RS_K)
    return [K // B + (1 if i < K % B else 0) for i in range(B)]


def _rs_order(lens):
    """Byte interleaving between codewords: byte j of every codeword, then byte j + 1, ..."""
    return [(i, j) for j in range(max(lens)) for i, n in enumerate(lens) if j < n]


def rs_encode(d):
    by = list(np.packbits(np.asarray(d, dtype=np.uint8)))
    s, cws = 0, []
    for k in _rs_layout(len(by)):
        cws.append(rs_encode_block([int(v) for v in by[s:s + k]]))
        s += k
    return np.unpackbits(np.array([cws[i][j] for i, j in _rs_order([len(c) for c in cws])], dtype=np.uint8)).astype(int)


def rs_decode(c, n0):
    lens = [k + RS_P for k in _rs_layout(math.ceil(n0 / 8))]
    by = np.packbits(np.asarray(c, dtype=np.uint8))
    cws = [[0] * n for n in lens]
    for q, (i, j) in enumerate(_rs_order(lens)):
        cws[i][j] = int(by[q])
    res = [rs_decode_block(cw) for cw in cws]
    msg = np.array([v for r in res for v in r[0]], dtype=np.uint8)
    return np.unpackbits(msg)[:n0].astype(int), sum(r[1] is None for r in res), len(res)


# ---- v7: LDPC codes of 5G NR (3GPP TS 38.212: 5.2.2 segmentation, 5.3.2 encoding, 5.4.2.1 rate matching), rates 1/4
# to 3/4 (v7.6: and 5/6). Differences from NR: no transport-block or code-block CRC (the parity checks say whether a block decoded),
# no bit interleaving (the modem's channel interleaver follows), the coded bits split evenly between the blocks,
# redundancy version 0 only. Base graphs: Tables 5.3.2-2 and 5.3.2-3, entries "row col V(iLS=0..7)".
LDPC_R = {"ldpc4": (1, 4), "ldpc3": (1, 3), "ldpc": (1, 2), "ldpc23": (2, 3), "ldpc34": (3, 4), "ldpc56": (5, 6)}   # code -> rate num/den
LDPC_ITERS = 25          # layered belief propagation: at most this many iterations (about 50 of the flooding schedule)
LDPC_CLIP, LDPC_CLIP1 = 30.0, 20.0   # input LLRs clipped (as the turbo decoder); the extension parity bits, which one check
# alone protects, lower than what that check can send (about 28), or one read wrong with near certainty would stay wrong
NR_BG1 = ('0 0 250 307 73 223 211 294 0 135;0 1 69 19 15 16 198 118 0 227;0 2 226 50 103 94 188 167 0 126;'
          '0 3 159 369 49 91 186 330 0 134;0 5 100 181 240 74 219 207 0 84;0 6 10 216 39 10 4 165 0 83;'
          '0 9 59 317 15 0 29 243 0 53;0 10 229 288 162 205 144 250 0 225;0 11 110 109 215 216 116 1 0 205;'
          '0 12 191 17 164 21 216 339 0 128;0 13 9 357 133 215 115 201 0 75;0 15 195 215 298 14 233 53 0 135;'
          '0 16 23 106 110 70 144 347 0 217;0 18 190 242 113 141 95 304 0 220;0 19 35 180 16 198 216 167 0 90;'
          '0 20 239 330 189 104 73 47 0 105;0 21 31 346 32 81 261 188 0 137;0 22 1 1 1 1 1 1 0 1;0 23 0 0 0 0 0 0 0 0;'
          '1 0 2 76 303 141 179 77 22 96;1 2 239 76 294 45 162 225 11 236;1 3 117 73 27 151 223 96 124 136;'
          '1 4 124 288 261 46 256 338 0 221;1 5 71 144 161 119 160 268 10 128;1 7 222 331 133 157 76 112 0 92;'
          '1 8 104 331 4 133 202 302 0 172;1 9 173 178 80 87 117 50 2 56;1 11 220 295 129 206 109 167 16 11;'
          '1 12 102 342 300 93 15 253 60 189;1 14 109 217 76 79 72 334 0 95;1 15 132 99 266 9 152 242 6 85;'
          '1 16 142 354 72 118 158 257 30 153;1 17 155 114 83 194 147 133 0 87;1 19 255 331 260 31 156 9 168 163;'
          '1 21 28 112 301 187 119 302 31 216;1 22 0 0 0 0 0 0 105 0;1 23 0 0 0 0 0 0 0 0;1 24 0 0 0 0 0 0 0 0;'
          '2 0 106 205 68 207 258 226 132 189;2 1 111 250 7 203 167 35 37 4;2 2 185 328 80 31 220 213 21 225;'
          '2 4 63 332 280 176 133 302 180 151;2 5 117 256 38 180 243 111 4 236;2 6 93 161 227 186 202 265 149 117;'
          '2 7 229 267 202 95 218 128 48 179;2 8 177 160 200 153 63 237 38 92;2 9 95 63 71 177 0 294 122 24;'
          '2 10 39 129 106 70 3 127 195 68;2 13 142 200 295 77 74 110 155 6;2 14 225 88 283 214 229 286 28 101;'
          '2 15 225 53 301 77 0 125 85 33;2 17 245 131 184 198 216 131 47 96;2 18 205 240 246 117 269 163 179 125;'
          '2 19 251 205 230 223 200 210 42 67;2 20 117 13 276 90 234 7 66 230;2 24 0 0 0 0 0 0 0 0;2 25 0 0 0 0 0 0 0 0;'
          '3 0 121 276 220 201 187 97 4 128;3 1 89 87 208 18 145 94 6 23;3 3 84 0 30 165 166 49 33 162;'
          '3 4 20 275 197 5 108 279 113 220;3 6 150 199 61 45 82 139 49 43;3 7 131 153 175 142 132 166 21 186;'
          '3 8 243 56 79 16 197 91 6 96;3 10 136 132 281 34 41 106 151 1;3 11 86 305 303 155 162 246 83 216;'
          '3 12 246 231 253 213 57 345 154 22;3 13 219 341 164 147 36 269 87 24;3 14 211 212 53 69 115 185 5 167;'
          '3 16 240 304 44 96 242 249 92 200;3 17 76 300 28 74 165 215 173 32;3 18 244 271 77 99 0 143 120 235;'
          '3 20 144 39 319 30 113 121 2 172;3 21 12 357 68 158 108 121 142 219;3 22 1 1 1 1 1 1 0 1;3 25 0 0 0 0 0 0 0 0;'
          '4 0 157 332 233 170 246 42 24 64;4 1 102 181 205 10 235 256 204 211;4 26 0 0 0 0 0 0 0 0;'
          '5 0 205 195 83 164 261 219 185 2;5 1 236 14 292 59 181 130 100 171;5 3 194 115 50 86 72 251 24 47;'
          '5 12 231 166 318 80 283 322 65 143;5 16 28 241 201 182 254 295 207 210;5 21 123 51 267 130 79 258 161 180;'
          '5 22 115 157 279 153 144 283 72 180;5 27 0 0 0 0 0 0 0 0;6 0 183 278 289 158 80 294 6 199;'
          '6 6 22 257 21 119 144 73 27 22;6 10 28 1 293 113 169 330 163 23;6 11 67 351 13 21 90 99 50 100;'
          '6 13 244 92 232 63 59 172 48 92;6 17 11 253 302 51 177 150 24 207;6 18 157 18 138 136 151 284 38 52;'
          '6 20 211 225 235 116 108 305 91 13;6 28 0 0 0 0 0 0 0 0;7 0 220 9 12 17 169 3 145 77;7 1 44 62 88 76 189 103 88 146;'
          '7 4 159 316 207 104 154 224 112 209;7 7 31 333 50 100 184 297 153 32;7 8 167 290 25 150 104 215 159 166;'
          '7 14 104 114 76 158 164 39 76 18;7 29 0 0 0 0 0 0 0 0;8 0 112 307 295 33 54 348 172 181;8 1 4 179 133 95 0 75 2 105;'
          '8 3 7 165 130 4 252 22 131 141;8 12 211 18 231 217 41 312 141 223;8 16 102 39 296 204 98 224 96 177;'
          '8 19 164 224 110 39 46 17 99 145;8 21 109 368 269 58 15 59 101 199;8 22 241 67 245 44 230 314 35 153;'
          '8 24 90 170 154 201 54 244 116 38;8 30 0 0 0 0 0 0 0 0;9 0 103 366 189 9 162 156 6 169;'
          '9 1 182 232 244 37 159 88 10 12;9 10 109 321 36 213 93 293 145 206;9 11 21 133 286 105 134 111 53 221;'
          '9 13 142 57 151 89 45 92 201 17;9 17 14 303 267 185 132 152 4 212;9 18 61 63 135 109 76 23 164 92;'
          '9 20 216 82 209 218 209 337 173 205;9 31 0 0 0 0 0 0 0 0;10 1 98 101 14 82 178 175 126 116;'
          '10 2 149 339 80 165 1 253 77 151;10 4 167 274 211 174 28 27 156 70;10 7 160 111 75 19 267 231 16 230;'
          '10 8 49 383 161 194 234 49 12 115;10 14 58 354 311 103 201 267 70 84;10 32 0 0 0 0 0 0 0 0;'
          '11 0 77 48 16 52 55 25 184 45;11 1 41 102 147 11 23 322 194 115;11 12 83 8 290 2 274 200 123 134;'
          '11 16 182 47 289 35 181 351 16 1;11 21 78 188 177 32 273 166 104 152;11 22 252 334 43 84 39 338 109 165;'
          '11 23 22 115 280 201 26 192 124 107;11 33 0 0 0 0 0 0 0 0;12 0 160 77 229 142 225 123 6 186;'
          '12 1 42 186 235 175 162 217 20 215;12 10 21 174 169 136 244 142 203 124;12 11 32 232 48 3 151 110 153 180;'
          '12 13 234 50 105 28 238 176 104 98;12 18 7 74 52 182 243 76 207 80;12 34 0 0 0 0 0 0 0 0;'
          '13 0 177 313 39 81 231 311 52 220;13 3 248 177 302 56 0 251 147 185;13 7 151 266 303 72 216 265 1 154;'
          '13 20 185 115 160 217 47 94 16 178;13 23 62 370 37 78 36 81 46 150;13 35 0 0 0 0 0 0 0 0;'
          '14 0 206 142 78 14 0 22 1 124;14 12 55 248 299 175 186 322 202 144;14 15 206 137 54 211 253 277 118 182;'
          '14 16 127 89 61 191 16 156 130 95;14 17 16 347 179 51 0 66 1 72;14 21 229 12 258 43 79 78 2 76;14 36 0 0 0 0 0 0 0 0;'
          '15 0 40 241 229 90 170 176 173 39;15 1 96 2 290 120 0 348 6 138;15 10 65 210 60 131 183 15 81 220;'
          '15 13 63 318 130 209 108 81 182 173;15 18 75 55 184 209 68 176 53 142;15 25 179 269 51 81 64 113 46 49;'
          '15 37 0 0 0 0 0 0 0 0;16 1 64 13 69 154 270 190 88 78;16 3 49 338 140 164 13 293 198 152;'
          '16 11 49 57 45 43 99 332 160 84;16 20 51 289 115 189 54 331 122 5;16 22 154 57 300 101 0 114 182 205;'
          '16 38 0 0 0 0 0 0 0 0;17 0 7 260 257 56 153 110 91 183;17 14 164 303 147 110 137 228 184 112;'
          '17 16 59 81 128 200 0 247 30 106;17 17 1 358 51 63 0 116 3 219;17 21 144 375 228 4 162 190 155 129;'
          '17 39 0 0 0 0 0 0 0 0;18 1 42 130 260 199 161 47 1 183;18 12 233 163 294 110 151 286 41 215;'
          '18 13 8 280 291 200 0 246 167 180;18 18 155 132 141 143 241 181 68 143;18 19 147 4 295 186 144 73 148 14;'
          '18 40 0 0 0 0 0 0 0 0;19 0 60 145 64 8 0 87 12 179;19 1 73 213 181 6 0 110 6 108;19 7 72 344 101 103 118 147 166 159;'
          '19 8 127 242 270 198 144 258 184 138;19 10 224 197 41 8 0 204 191 196;19 41 0 0 0 0 0 0 0 0;'
          '20 0 151 187 301 105 265 89 6 77;20 3 186 206 162 210 81 65 12 187;20 9 217 264 40 121 90 155 15 203;'
          '20 11 47 341 130 214 144 244 5 167;20 22 160 59 10 183 228 30 30 130;20 42 0 0 0 0 0 0 0 0;'
          '21 1 249 205 79 192 64 162 6 197;21 5 121 102 175 131 46 264 86 122;21 16 109 328 132 220 266 346 96 215;'
          '21 20 131 213 283 50 9 143 42 65;21 21 171 97 103 106 18 109 199 216;21 43 0 0 0 0 0 0 0 0;'
          '22 0 64 30 177 53 72 280 44 25;22 12 142 11 20 0 189 157 58 47;22 13 188 233 55 3 72 236 130 126;'
          '22 17 158 22 316 148 257 113 131 178;22 44 0 0 0 0 0 0 0 0;23 1 156 24 249 88 180 18 45 185;'
          '23 2 147 89 50 203 0 6 18 127;23 10 170 61 133 168 0 181 132 117;23 18 152 27 105 122 165 304 100 199;'
          '23 45 0 0 0 0 0 0 0 0;24 0 112 298 289 49 236 38 9 32;24 3 86 158 280 157 199 170 125 178;'
          '24 4 236 235 110 64 0 249 191 2;24 11 116 339 187 193 266 288 28 156;24 22 222 234 281 124 0 194 6 58;'
          '24 46 0 0 0 0 0 0 0 0;25 1 23 72 172 1 205 279 4 27;25 6 136 17 295 166 0 255 74 141;25 7 116 383 96 65 0 111 16 11;'
          '25 14 182 312 46 81 183 54 28 181;25 47 0 0 0 0 0 0 0 0;26 0 195 71 270 107 0 325 21 163;'
          '26 2 243 81 110 176 0 326 142 131;26 4 215 76 318 212 0 226 192 169;26 15 61 136 67 127 277 99 197 98;'
          '26 48 0 0 0 0 0 0 0 0;27 1 25 194 210 208 45 91 98 165;27 6 104 194 29 141 36 326 140 232;'
          '27 8 194 101 304 174 72 268 22 9;27 49 0 0 0 0 0 0 0 0;28 0 128 222 11 146 275 102 4 32;28 4 165 19 293 153 0 1 1 43;'
          '28 19 181 244 50 217 155 40 40 200;28 21 63 274 234 114 62 167 93 205;28 50 0 0 0 0 0 0 0 0;'
          '29 1 86 252 27 150 0 273 92 232;29 14 236 5 308 11 180 104 136 32;29 18 84 147 117 53 0 243 106 118;'
          '29 25 6 78 29 68 42 107 6 103;29 51 0 0 0 0 0 0 0 0;30 0 216 159 91 34 0 171 2 170;30 10 73 229 23 130 90 16 88 199;'
          '30 13 120 260 105 210 252 95 112 26;30 24 9 90 135 123 173 212 20 105;30 52 0 0 0 0 0 0 0 0;'
          '31 1 95 100 222 175 144 101 4 73;31 7 177 215 308 49 144 297 49 149;31 22 172 258 66 177 166 279 125 175;'
          '31 25 61 256 162 128 19 222 194 108;31 53 0 0 0 0 0 0 0 0;32 0 221 102 210 192 0 351 6 103;'
          '32 12 112 201 22 209 211 265 126 110;32 14 199 175 271 58 36 338 63 151;32 24 121 287 217 30 162 83 20 211;'
          '32 54 0 0 0 0 0 0 0 0;33 1 2 323 170 114 0 56 10 199;33 2 187 8 20 49 0 304 30 132;33 11 41 361 140 161 76 141 6 172;'
          '33 21 211 105 33 137 18 101 92 65;33 55 0 0 0 0 0 0 0 0;34 0 127 230 187 82 197 60 4 161;'
          '34 7 167 148 296 186 0 320 153 237;34 15 164 202 5 68 108 112 197 142;34 17 159 312 44 150 0 54 155 180;'
          '34 56 0 0 0 0 0 0 0 0;35 1 161 320 207 192 199 100 4 231;35 6 197 335 158 173 278 210 45 174;'
          '35 12 207 2 55 26 0 195 168 145;35 22 103 266 285 187 205 268 185 100;35 57 0 0 0 0 0 0 0 0;'
          '36 0 37 210 259 222 216 135 6 11;36 14 105 313 179 157 16 15 200 207;36 15 51 297 178 0 0 35 177 42;'
          '36 18 120 21 160 6 0 188 43 100;36 58 0 0 0 0 0 0 0 0;37 1 198 269 298 81 72 319 82 59;'
          '37 13 220 82 15 195 144 236 2 204;37 23 122 115 115 138 0 85 135 161;37 59 0 0 0 0 0 0 0 0;'
          '38 0 167 185 151 123 190 164 91 121;38 9 151 177 179 90 0 196 64 90;38 10 157 289 64 73 0 209 198 26;'
          '38 12 163 214 181 10 0 246 100 140;38 60 0 0 0 0 0 0 0 0;39 1 173 258 102 12 153 236 4 115;'
          '39 3 139 93 77 77 0 264 28 188;39 7 149 346 192 49 165 37 109 168;39 19 0 297 208 114 117 272 188 52;'
          '39 61 0 0 0 0 0 0 0 0;40 0 157 175 32 67 216 304 10 4;40 8 137 37 80 45 144 237 84 103;'
          '40 17 149 312 197 96 2 135 12 30;40 62 0 0 0 0 0 0 0 0;41 1 167 52 154 23 0 123 2 53;41 3 173 314 47 215 0 77 75 189;'
          '41 9 139 139 124 60 0 25 142 215;41 18 151 288 207 167 183 272 128 24;41 63 0 0 0 0 0 0 0 0;'
          '42 0 149 113 226 114 27 288 163 222;42 4 157 14 65 91 0 83 10 170;42 24 137 218 126 78 35 17 162 71;'
          '42 64 0 0 0 0 0 0 0 0;43 1 151 113 228 206 52 210 1 22;43 16 163 132 69 22 243 3 163 127;'
          '43 18 173 114 176 134 0 53 99 49;43 25 139 168 102 161 270 167 98 125;43 65 0 0 0 0 0 0 0 0;'
          '44 0 139 80 234 84 18 79 4 191;44 7 157 78 227 4 0 244 6 211;44 9 163 163 259 9 0 293 142 187;'
          '44 22 173 274 260 12 57 272 3 148;44 66 0 0 0 0 0 0 0 0;45 1 149 135 101 184 168 82 181 177;'
          '45 6 151 149 228 121 0 67 45 114;45 10 167 15 126 29 144 235 153 93;45 67 0 0 0 0 0 0 0 0;')
NR_BG2 = ('0 0 9 174 0 72 3 156 143 145;0 1 117 97 0 110 26 143 19 131;0 2 204 166 0 23 53 14 176 71;0 3 26 66 0 181 35 3 165 21;'
          '0 6 189 71 0 95 115 40 196 23;0 9 205 172 0 8 127 123 13 112;0 10 0 0 0 1 0 0 0 1;0 11 0 0 0 0 0 0 0 0;'
          '1 0 167 27 137 53 19 17 18 142;1 3 166 36 124 156 94 65 27 174;1 4 253 48 0 115 104 63 3 183;'
          '1 5 125 92 0 156 66 1 102 27;1 6 226 31 88 115 84 55 185 96;1 7 156 187 0 200 98 37 17 23;'
          '1 8 224 185 0 29 69 171 14 9;1 9 252 3 55 31 50 133 180 167;1 11 0 0 0 0 0 0 0 0;1 12 0 0 0 0 0 0 0 0;'
          '2 0 81 25 20 152 95 98 126 74;2 1 114 114 94 131 106 168 163 31;2 3 44 117 99 46 92 107 47 3;'
          '2 4 52 110 9 191 110 82 183 53;2 8 240 114 108 91 111 142 132 155;2 10 1 1 1 0 1 1 1 0;2 12 0 0 0 0 0 0 0 0;'
          '2 13 0 0 0 0 0 0 0 0;3 1 8 136 38 185 120 53 36 239;3 2 58 175 15 6 121 174 48 171;3 4 158 113 102 36 22 174 18 95;'
          '3 5 104 72 146 124 4 127 111 110;3 6 209 123 12 124 73 17 203 159;3 7 54 118 57 110 49 89 3 199;'
          '3 8 18 28 53 156 128 17 191 43;3 9 128 186 46 133 79 105 160 75;3 10 0 0 0 1 0 0 0 1;3 13 0 0 0 0 0 0 0 0;'
          '4 0 179 72 0 200 42 86 43 29;4 1 214 74 136 16 24 67 27 140;4 11 71 29 157 101 51 83 117 180;4 14 0 0 0 0 0 0 0 0;'
          '5 0 231 10 0 185 40 79 136 121;5 1 41 44 131 138 140 84 49 41;5 5 194 121 142 170 84 35 36 169;'
          '5 7 159 80 141 219 137 103 132 88;5 11 103 48 64 193 71 60 62 207;5 15 0 0 0 0 0 0 0 0;6 0 155 129 0 123 109 47 7 137;'
          '6 5 228 92 124 55 87 154 34 72;6 7 45 100 99 31 107 10 198 172;6 9 28 49 45 222 133 155 168 124;'
          '6 11 158 184 148 209 139 29 12 56;6 16 0 0 0 0 0 0 0 0;7 1 129 80 0 103 97 48 163 86;7 5 147 186 45 13 135 125 78 186;'
          '7 7 140 16 148 105 35 24 143 87;7 11 3 102 96 150 108 47 107 172;7 13 116 143 78 181 65 55 58 154;'
          '7 17 0 0 0 0 0 0 0 0;8 0 142 118 0 147 70 53 101 176;8 1 94 70 65 43 69 31 177 169;8 12 230 152 87 152 88 161 22 225;'
          '8 18 0 0 0 0 0 0 0 0;9 1 203 28 0 2 97 104 186 167;9 8 205 132 97 30 40 142 27 238;9 10 61 185 51 184 24 99 205 48;'
          '9 11 247 178 85 83 49 64 81 68;9 19 0 0 0 0 0 0 0 0;10 0 11 59 0 174 46 111 125 38;10 1 185 104 17 150 41 25 60 217;'
          '10 6 0 22 156 8 101 174 177 208;10 7 117 52 20 56 96 23 51 232;10 20 0 0 0 0 0 0 0 0;11 0 11 32 0 99 28 91 39 178;'
          '11 7 236 92 7 138 30 175 29 214;11 9 210 174 4 110 116 24 35 168;11 13 56 154 2 99 64 141 8 51;11 21 0 0 0 0 0 0 0 0;'
          '12 1 63 39 0 46 33 122 18 124;12 3 111 93 113 217 122 11 155 122;12 11 14 11 48 109 131 4 49 72;12 22 0 0 0 0 0 0 0 0;'
          '13 0 83 49 0 37 76 29 32 48;13 1 2 125 112 113 37 91 53 57;13 8 38 35 102 143 62 27 95 167;'
          '13 13 222 166 26 140 47 127 186 219;13 23 0 0 0 0 0 0 0 0;14 1 115 19 0 36 143 11 91 82;'
          '14 6 145 118 138 95 51 145 20 232;14 11 3 21 57 40 130 8 52 204;14 13 232 163 27 116 97 166 109 162;'
          '14 24 0 0 0 0 0 0 0 0;15 0 51 68 0 116 139 137 174 38;15 10 175 63 73 200 96 103 108 217;'
          '15 11 213 81 99 110 128 40 102 157;15 25 0 0 0 0 0 0 0 0;16 1 203 87 0 75 48 78 125 170;'
          '16 9 142 177 79 158 9 158 31 23;16 11 8 135 111 134 28 17 54 175;16 12 242 64 143 97 8 165 176 202;'
          '16 26 0 0 0 0 0 0 0 0;17 1 254 158 0 48 120 134 57 196;17 5 124 23 24 132 43 23 201 173;'
          '17 11 114 9 109 206 65 62 142 195;17 12 64 6 18 2 42 163 35 218;17 27 0 0 0 0 0 0 0 0;'
          '18 0 220 186 0 68 17 173 129 128;18 6 194 6 18 16 106 31 203 211;18 7 50 46 86 156 142 22 140 210;'
          '18 28 0 0 0 0 0 0 0 0;19 0 87 58 0 35 79 13 110 39;19 1 20 42 158 138 28 135 124 84;19 10 185 156 154 86 41 145 52 88;'
          '19 29 0 0 0 0 0 0 0 0;20 1 26 76 0 6 2 128 196 117;20 4 105 61 148 20 103 52 35 227;20 11 29 153 104 141 78 173 114 6;'
          '20 30 0 0 0 0 0 0 0 0;21 0 76 157 0 80 91 156 10 238;21 8 42 175 17 43 75 166 122 13;21 13 210 67 33 81 81 40 23 11;'
          '21 31 0 0 0 0 0 0 0 0;22 1 222 20 0 49 54 18 202 195;22 2 63 52 4 1 132 163 126 44;22 32 0 0 0 0 0 0 0 0;'
          '23 0 23 106 0 156 68 110 52 5;23 3 235 86 75 54 115 132 170 94;23 5 238 95 158 134 56 150 13 111;'
          '23 33 0 0 0 0 0 0 0 0;24 1 46 182 0 153 30 113 113 81;24 2 139 153 69 88 42 108 161 19;24 9 8 64 87 63 101 61 88 130;'
          '24 34 0 0 0 0 0 0 0 0;25 0 228 45 0 211 128 72 197 66;25 5 156 21 65 94 63 136 194 95;25 35 0 0 0 0 0 0 0 0;'
          '26 2 29 67 0 90 142 36 164 146;26 7 143 137 100 6 28 38 172 66;26 12 160 55 13 221 100 53 49 190;'
          '26 13 122 85 7 6 133 145 161 86;26 36 0 0 0 0 0 0 0 0;27 0 8 103 0 27 13 42 168 64;27 6 151 50 32 118 10 104 193 181;'
          '27 37 0 0 0 0 0 0 0 0;28 1 98 70 0 216 106 64 14 7;28 2 101 111 126 212 77 24 186 144;'
          '28 5 135 168 110 193 43 149 46 16;28 38 0 0 0 0 0 0 0 0;29 0 18 110 0 108 133 139 50 25;'
          '29 4 28 17 154 61 25 161 27 57;29 39 0 0 0 0 0 0 0 0;30 2 71 120 0 106 87 84 70 37;30 5 240 154 35 44 56 173 17 139;'
          '30 7 9 52 51 185 104 93 50 221;30 9 84 56 134 176 70 29 6 17;30 40 0 0 0 0 0 0 0 0;31 1 106 3 0 147 80 117 115 201;'
          '31 13 1 170 20 182 139 148 189 46;31 41 0 0 0 0 0 0 0 0;32 0 242 84 0 108 32 116 110 179;32 5 44 8 20 21 89 73 0 14;'
          '32 12 166 17 122 110 71 142 163 116;32 42 0 0 0 0 0 0 0 0;33 2 132 165 0 71 135 105 163 46;'
          '33 7 164 179 88 12 6 137 173 2;33 10 235 124 13 109 2 29 179 106;33 43 0 0 0 0 0 0 0 0;'
          '34 0 147 173 0 29 37 11 197 184;34 12 85 177 19 201 25 41 191 135;34 13 36 12 78 69 114 162 193 141;'
          '34 44 0 0 0 0 0 0 0 0;35 1 57 77 0 91 60 126 157 85;35 5 40 184 157 165 137 152 167 225;'
          '35 11 63 18 6 55 93 172 181 175;35 45 0 0 0 0 0 0 0 0;36 0 140 25 0 1 121 73 197 178;'
          '36 2 38 151 63 175 129 154 167 112;36 7 154 170 82 83 26 129 179 106;36 46 0 0 0 0 0 0 0 0;'
          '37 10 219 37 0 40 97 167 181 154;37 13 151 31 144 12 56 38 193 114;37 47 0 0 0 0 0 0 0 0;38 1 31 84 0 37 1 112 157 42;'
          '38 5 66 151 93 97 70 7 173 41;38 11 38 190 19 46 1 19 191 105;38 48 0 0 0 0 0 0 0 0;39 0 239 93 0 106 119 109 181 167;'
          '39 7 172 132 24 181 32 6 157 45;39 12 34 57 138 154 142 105 173 189;39 49 0 0 0 0 0 0 0 0;'
          '40 2 0 103 0 98 6 160 193 78;40 10 75 107 36 35 73 156 163 67;40 13 120 163 143 36 102 82 179 180;'
          '40 50 0 0 0 0 0 0 0 0;41 1 129 147 0 120 48 132 191 53;41 5 229 7 2 101 47 6 197 215;41 11 118 60 55 81 19 8 167 230;'
          '41 51 0 0 0 0 0 0 0 0;')
NR_GRAPHS = {g: [tuple(int(v) for v in e.split()) for e in s.strip(";").split(";")] for g, s in ((1, NR_BG1), (2, NR_BG2))}
NR_ZSETS = [[a << j for j in range(8) if a << j <= 384] for a in (2, 3, 5, 7, 9, 11, 13, 15)]   # Table 5.3.2-1: Z = a 2^j
NR_GEOM = {1: (22, 68, 8448), 2: (10, 52, 3840)}   # information columns, columns, largest block (bits)


def ldpc_params(nbits, rate):
    """How nbits message bits are coded at rate (num, den): base graph (7.2.2: graph 2 for short messages or low
    rates), C blocks of kp message bits, lifting Z (set ils), K = kb Z bits per block (message, then filler zeros),
    E coded bits sent per block, the columns and rows of the graph those bits reach."""
    num, den = rate
    bg = 2 if nbits <= 292 or (nbits <= 3824 and num * 100 <= 67 * den) or num * 4 <= den else 1
    kbase, ncol, kcb = NR_GEOM[bg]
    C = max(1, -(-nbits // kcb))
    kp = [nbits // C + (i < nbits % C) for i in range(C)]
    kb = 22 if bg == 1 else 10 if nbits > 640 else 9 if nbits > 560 else 8 if nbits > 192 else 6
    Z = min(z for s in NR_ZSETS for z in s if kb * z >= max(kp))
    ils = next(i for i, s in enumerate(NR_ZSETS) if Z in s)
    K, G = kbase * Z, -(-nbits * den // num)
    E = [G // C + (i < G % C) for i in range(C)]
    cols = min(ncol, max(kbase + 4, -(-(max(E) + K - min(kp)) // Z) + 2))
    return dict(bg=bg, C=C, kp=kp, Z=Z, ils=ils, K=K, E=E, cols=cols, rows=cols - kbase, kb=kbase)


_LDPC_H = {}
def ldpc_graph(p):
    """The lifted graph's rows: [(column, shift)] per row (shift = V mod Z). Block P^s: check t of a row sees bit
    (t + s) mod Z of the column, so gathering is roll(x, -s) and scattering roll(y, s)."""
    key = (p["bg"], p["ils"], p["Z"], p["rows"])
    if key not in _LDPC_H:
        H = [[] for _ in range(p["rows"])]
        for e in NR_GRAPHS[p["bg"]]:
            if e[0] < p["rows"]:
                H[e[0]].append((e[1], e[2 + p["ils"]] % p["Z"]))
        _LDPC_H[key] = H
    return _LDPC_H[key]


def ldpc_encode_block(info, p):
    """K bits (message, filler zeros) -> the codeword as (columns, Z). The four core rows summed leave only the first
    parity column (its shifts cancel in pairs but one: the dual diagonal); then each row with a single unknown
    column gives that column."""
    Z, kb, H = p["Z"], p["kb"], ldpc_graph(p)
    x = np.zeros((p["cols"], Z), dtype=np.uint8)
    x[:kb] = np.asarray(info, dtype=np.uint8).reshape(kb, Z)
    tot = np.zeros(Z, dtype=np.uint8)
    for r in range(4):
        for c, s in H[r]:
            if c < kb:
                tot ^= np.roll(x[c], -s)
    sh = [s for r in range(4) for c, s in H[r] if c == kb]
    x[kb] = np.roll(tot, [s for s in sh if sh.count(s) % 2][0])
    known, todo = set(range(kb + 1)), list(range(p["rows"]))
    while todo:
        for r in todo:
            unk = [(c, s) for c, s in H[r] if c not in known]
            if len(unk) <= 1:
                if unk:
                    acc = np.zeros(Z, dtype=np.uint8)
                    for c, s in H[r]:
                        if c != unk[0][0]:
                            acc ^= np.roll(x[c], -s)
                    x[unk[0][0]] = np.roll(acc, unk[0][1])
                    known.add(unk[0][0])
                todo.remove(r)
                break
    return x


def ldpc_sent(p, i):
    """Positions (in the codeword's columns x Z bits) of block i's sent bits, in order: from column 2 on (the first
    2Z message bits are never sent), filler bits skipped."""
    idx = np.arange(2 * p["Z"], p["cols"] * p["Z"])
    return idx[(idx < p["kp"][i]) | (idx >= p["K"])][:p["E"][i]]


def ldpc_encode(bits, rate):
    bits, p, out, q = np.asarray(bits, dtype=np.uint8), ldpc_params(len(bits), rate), [], 0
    for i in range(p["C"]):
        info = np.zeros(p["K"], dtype=np.uint8)
        info[:p["kp"][i]] = bits[q:q + p["kp"][i]]
        q += p["kp"][i]
        out.append(ldpc_encode_block(info, p).ravel()[ldpc_sent(p, i)])
    return np.concatenate(out).astype(int)


def _ldpc_phi(x):
    return -np.log(np.tanh(np.clip(x, 1e-12, 50.0) / 2))


def ldpc_decode(llr, nbits, rate, iters=LDPC_ITERS):
    """Layered belief propagation, every block at once. llr: the sent bits' LLRs (> 0 means 1), clipped to +-30 (+-20). Returns
    the message bits, the a posteriori LLRs of the sent bits (with the clipped part added back), the iterations run,
    whether every block's checks hold, and the decoder's own a posteriori LLRs (the clipped part not added back)."""
    p = ldpc_params(nbits, rate)
    Z, C, cols, H = p["Z"], p["C"], p["cols"], ldpc_graph(p)
    lam, pos, q, clip = np.zeros((C, cols * Z)), [], 0, []      # ln P(0)/P(1); 0 for the bits never sent
    for i in range(C):
        pos.append(ldpc_sent(p, i))
        clip.append(np.where(pos[i] >= p["K"] + 4 * Z, LDPC_CLIP1, LDPC_CLIP))   # as the turbo decoder: a few wrong bits
        # read with near certainty (a quiet channel) must stay within what the checks can overrule
        lam[i, pos[i]] = -np.clip(np.asarray(llr[q:q + p["E"][i]], dtype=float), -clip[i], clip[i])
        lam[i, p["kp"][i]:p["K"]] = 60.0                         # filler bits: known zeros
        q += p["E"][i]
    L = lam.reshape(C, cols, Z)
    R = [[np.zeros((C, Z)) for _ in row] for row in H]
    ok, it = False, 0
    for it in range(1, iters + 1):
        for r, row in enumerate(H):
            Q = np.stack([np.roll(L[:, c], -s, axis=1) - R[r][k] for k, (c, s) in enumerate(row)])
            sg = np.where(Q < 0, -1.0, 1.0)
            ph = _ldpc_phi(np.abs(Q))
            tot, sgt = ph.sum(axis=0), np.prod(sg, axis=0)
            for k, (c, s) in enumerate(row):
                R[r][k] = sgt * sg[k] * _ldpc_phi(tot - ph[k])
                L[:, c] = np.roll(Q[k] + R[r][k], s, axis=1)
        hard = (L < 0).astype(np.uint8)
        ok = all(not np.bitwise_xor.reduce(np.stack([np.roll(hard[:, c], -s, axis=1) for c, s in row]), axis=0).any() for row in H)
        if ok:
            break
    flat = L.reshape(C, cols * Z)
    bits = np.concatenate([(flat[i, :p["kp"][i]] < 0).astype(int) for i in range(C)])
    post = np.concatenate([-flat[i, pos[i]] for i in range(C)])
    llr, cl = np.asarray(llr, dtype=float)[:len(post)], np.concatenate(clip)
    # app - llr: what the code adds, the clipping undone; post: the decoder's own a posteriori LLRs (its decisions)
    return bits, post + llr - np.clip(llr, -cl, cl), it, ok, post


# ---- v7.3: polar codes of 5G NR (3GPP TS 38.212: 5.3.1 code construction and encoding, 5.4.1 rate matching, the
# CRC of the uplink control channel, 6.3.1.2.1), rates 1/3, 1/2, 2/3, decoded by CRC-aided list decoding. Differences
# from NR: the message is split into as many blocks as needed for at most 1024 sent bits each (NR: at most two blocks),
# every block carries an 11-bit CRC (NR: 6 bits and 3 parity-check bits under 20 message bits), no bit interleaving of
# NR's own (the modem's channel interleaver follows). The rate counts the CRC as part of the block: E = (a + 11) / R.
POLAR_R = {"polar3": (1, 3), "polar": (1, 2), "polar23": (2, 3)}   # code -> rate num/den
POLAR_Q = [int(v) for v in (
    '0 1 2 4 8 16 32 3 5 64 9 6 17 10 18 128 12 33 65 20 256 34 24 36 7 129 66 512 11 40 68 130 19 13 48 14 72 257 '
    '21 132 35 258 26 513 80 37 25 22 136 260 264 38 514 96 67 41 144 28 69 42 516 49 74 272 160 520 288 528 192 544 '
    '70 44 131 81 50 73 15 320 133 52 23 134 384 76 137 82 56 27 97 39 259 84 138 145 261 29 43 98 515 88 140 30 146 '
    '71 262 265 161 576 45 100 640 51 148 46 75 266 273 517 104 162 53 193 152 77 164 768 268 274 518 54 83 57 521 '
    '112 135 78 289 194 85 276 522 58 168 139 99 86 60 280 89 290 529 524 196 141 101 147 176 142 530 321 31 200 90 '
    '545 292 322 532 263 149 102 105 304 296 163 92 47 267 385 546 324 208 386 150 153 165 106 55 328 536 577 548 '
    '113 154 79 269 108 578 224 166 519 552 195 270 641 523 275 580 291 59 169 560 114 277 156 87 197 116 170 61 531 '
    '525 642 281 278 526 177 293 388 91 584 769 198 172 120 201 336 62 282 143 103 178 294 93 644 202 592 323 392 '
    '297 770 107 180 151 209 284 648 94 204 298 400 608 352 325 533 155 210 305 547 300 109 184 534 537 115 167 225 '
    '326 306 772 157 656 329 110 117 212 171 776 330 226 549 538 387 308 216 416 271 279 158 337 550 672 118 332 579 '
    '540 389 173 121 553 199 784 179 228 338 312 704 390 174 554 581 393 283 122 448 353 561 203 63 340 394 527 582 '
    '556 181 295 285 232 124 205 182 643 562 286 585 299 354 211 401 185 396 344 586 645 593 535 240 206 95 327 564 '
    '800 402 356 307 301 417 213 568 832 588 186 646 404 227 896 594 418 302 649 771 360 539 111 331 214 309 188 449 '
    '217 408 609 596 551 650 229 159 420 310 541 773 610 657 333 119 600 339 218 368 652 230 391 313 450 542 334 233 '
    '555 774 175 123 658 612 341 777 220 314 424 395 673 583 355 287 183 234 125 557 660 616 342 316 241 778 563 345 '
    '452 397 403 207 674 558 785 432 357 187 236 664 624 587 780 705 126 242 565 398 346 456 358 405 303 569 244 595 '
    '189 566 676 361 706 589 215 786 647 348 419 406 464 680 801 362 590 409 570 788 597 572 219 311 708 598 601 651 '
    '421 792 802 611 602 410 231 688 653 248 369 190 364 654 659 335 480 315 221 370 613 422 425 451 614 543 235 412 '
    '343 372 775 317 222 426 453 237 559 833 804 712 834 661 808 779 617 604 433 720 816 836 347 897 243 662 454 318 '
    '675 618 898 781 376 428 665 736 567 840 625 238 359 457 399 787 591 678 434 677 349 245 458 666 620 363 127 191 '
    '782 407 436 626 571 465 681 246 707 350 599 668 790 460 249 682 573 411 803 789 709 365 440 628 689 374 423 466 '
    '793 250 371 481 574 413 603 366 468 655 900 805 615 684 710 429 794 252 373 605 848 690 713 632 482 806 427 904 '
    '414 223 663 692 835 619 472 455 796 809 714 721 837 716 864 810 606 912 722 696 377 435 817 319 621 812 484 430 '
    '838 667 488 239 378 459 622 627 437 380 818 461 496 669 679 724 841 629 351 467 438 737 251 462 442 441 469 247 '
    '683 842 738 899 670 783 849 820 728 928 791 367 901 630 685 844 633 711 253 691 824 902 686 740 850 375 444 470 '
    '483 415 485 905 795 473 634 744 852 960 865 693 797 906 715 807 474 636 694 254 717 575 913 798 811 379 697 431 '
    '607 489 866 723 486 908 718 813 476 856 839 725 698 914 752 868 819 814 439 929 490 623 671 739 916 463 843 381 '
    '497 930 821 726 961 872 492 631 729 700 443 741 845 920 382 822 851 730 498 880 742 445 471 635 932 687 903 825 '
    '500 846 745 826 732 446 962 936 475 853 867 637 907 487 695 746 828 753 854 857 504 799 255 964 909 719 477 915 '
    '638 748 944 869 491 699 754 858 478 968 383 910 815 976 870 917 727 493 873 701 931 756 860 499 731 823 922 874 '
    '918 502 933 743 760 881 494 702 921 501 876 847 992 447 733 827 934 882 937 963 747 505 855 924 734 829 965 938 '
    '884 506 749 945 966 755 859 940 830 911 871 639 888 479 946 750 969 508 861 757 970 919 875 862 758 948 977 923 '
    '972 761 877 952 495 703 935 978 883 762 503 925 878 735 993 885 939 994 980 926 764 941 967 886 831 947 507 889 '
    '984 751 942 996 971 890 509 949 973 1000 892 950 863 759 1008 510 979 953 763 974 954 879 981 982 927 995 765 '
    '956 887 985 997 986 943 891 998 766 511 988 1001 951 1002 893 975 894 1009 955 1004 1010 957 983 958 987 1012 '
    '999 1016 767 989 1003 990 1005 959 1011 1013 895 1006 1014 1017 1018 991 1020 1007 1015 1019 1021 1022 1023').split()]   # Table 5.3.1.2-1: Q_0 .. Q_1023, least reliable first
POLAR_P = [0, 1, 2, 4, 3, 5, 6, 7, 8, 16, 9, 17, 10, 18, 11, 19, 12, 20, 13, 21, 14, 22, 15, 23, 24, 25, 26, 28, 27, 29, 30, 31]
# Table 5.4.1.1-1, the sub-block interleaver: 32 sub-blocks
POLAR_CRC = [1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1]   # CRC11: D^11 + D^10 + D^9 + D^5 + 1
POLAR_NCRC, POLAR_EMAX = 11, 1024
POLAR_LIST = 8           # list decoding: the 8 likeliest paths are kept; the first whose CRC holds wins
POLAR_SCAN = 2           # soft cancellation iterations: the soft output that the single carrier's passes need
POLAR_CLIP, POLAR_BIG, POLAR_SURE = 30.0, 1000.0, 100.0   # input LLRs clipped (as LDPC); known bits; a block whose CRC holds


def polar_crc(bits):
    """The 11 CRC bits of bits (shift register from zero, MSB first)."""
    r = [0] * POLAR_NCRC
    for b in bits:
        fb = r[0] ^ int(b)
        r = r[1:] + [0]
        if fb:
            r = [x ^ g for x, g in zip(r, POLAR_CRC[1:])]
    return r


_POLAR_M = {}
def polar_crc_matrix(a):
    """The CRC is linear: its bits are M m (mod 2) for a message m of a bits; M is 11 x a."""
    if a not in _POLAR_M:
        _POLAR_M[a] = np.array([polar_crc(np.eye(a, dtype=int)[i]) for i in range(a)], dtype=np.int64).reshape(a, POLAR_NCRC).T
    return _POLAR_M[a]


_POLAR_B = {}
def polar_block(K, E):
    """A block of K bits (message + CRC) sent as E bits: mother code length N = 2^n (5.3.1), the K information
    positions (the most reliable ones that rate matching leaves free, 5.4.1.1), and each sent bit's codeword position:
    the sub-block interleaver (y_k = d_J(k)), then repetition (E >= N), puncturing (the first N - E bits of y,
    K/E <= 7/16) or shortening (the last N - E bits of y, which the frozen bits make zeros)."""
    if (K, E) not in _POLAR_B:
        lg = (E - 1).bit_length()
        n1 = lg - 1 if 8 * E <= 9 * (1 << (lg - 1)) and 16 * K < 9 * E else lg
        n = max(min(n1, (8 * K - 1).bit_length(), 10), 5)
        N = 1 << n
        J = [POLAR_P[32 * k // N] * (N // 32) + k % (N // 32) for k in range(N)]
        mode = "repeat" if E >= N else "puncture" if 16 * K <= 7 * E else "shorten"
        frozen = set()
        if mode == "puncture":
            frozen |= set(J[:N - E]) | set(range(-(-(3 * N - 2 * E) // 4) if 4 * E >= 3 * N else -(-(9 * N - 4 * E) // 16)))
        elif mode == "shorten":
            frozen |= set(J[E:])
        free = [q for q in POLAR_Q if q < N and q not in frozen]
        isinfo = np.zeros(N, dtype=bool)
        isinfo[free[len(free) - K:]] = True
        sent = [J[k % N] for k in range(E)] if mode != "puncture" else J[N - E:]
        _POLAR_B[(K, E)] = dict(K=K, E=E, N=N, n=n, mode=mode, isinfo=isinfo, info=np.flatnonzero(isinfo),
                                sent=np.array(sent), known=np.array(J[E:] if mode == "shorten" else [], dtype=int))
    return _POLAR_B[(K, E)]


def polar_params(nbits, rate):
    """C blocks, the fewest that send at most 1024 bits each; a message bits and E sent bits per block (split evenly)."""
    num, den = rate
    C = 1
    while -(-(nbits + POLAR_NCRC * C) * den // num) > POLAR_EMAX * C:
        C += 1
    G = -(-(nbits + POLAR_NCRC * C) * den // num)
    a = [nbits // C + (i < nbits % C) for i in range(C)]
    E = [G // C + (i < G % C) for i in range(C)]
    return dict(C=C, a=a, E=E, G=G, blocks=[polar_block(a[i] + POLAR_NCRC, E[i]) for i in range(C)])


def polar_transform(u):
    """d = u G_N along the last axis, G_N = F (x) ... (x) F, F = [1 0; 1 1] (NR has no bit reversal): in stage s,
    bit i (bit s of i clear) adds bit i + s."""
    x = np.array(u, dtype=np.uint8)
    N, s = x.shape[-1], 1
    while s < N:
        v = x.reshape(x.shape[:-1] + (N // (2 * s), 2, s))
        v[..., 0, :] ^= v[..., 1, :]
        s *= 2
    return x


def polar_codeword(m, b):
    u = np.zeros(b["N"], dtype=np.uint8)
    u[b["info"]] = np.r_[m, polar_crc(m)]
    return polar_transform(u)[b["sent"]]


def polar_encode(bits, rate):
    bits, p, out, q = np.asarray(bits, dtype=np.uint8), polar_params(len(bits), rate), [], 0
    for i, b in enumerate(p["blocks"]):
        out.append(polar_codeword(bits[q:q + p["a"][i]], b))
        q += p["a"][i]
    return np.concatenate(out).astype(int)


def _bplus(a, b):
    """The check rule on LLRs: the LLR of x xor y from those of x and y (2 atanh(tanh(a/2) tanh(b/2)))."""
    return (np.sign(a) * np.sign(b) * np.minimum(np.abs(a), np.abs(b))
            + np.log1p(np.exp(-np.abs(a + b))) - np.log1p(np.exp(-np.abs(a - b))))


# Decoding works on the tree of the code: a node of size s holds bits (v1 xor v2, v2), v1 and v2 being the codewords of
# its two halves; from the LLRs a = (a1, a2) of its bits (ln P(0)/P(1)), v1's are a1 [+] a2, then, v1 known, v2's are
# a2 + (1 - 2 v1) a1 (successive cancellation, SC). The leaves are the bits u_i, decided in order.
def polar_scl(lam, b, L=POLAR_LIST):
    """SC list decoding, G blocks at once (lam: G x N codeword LLRs ln P(0)/P(1)): at each information bit every path
    tries 0 and 1 and the L likeliest stay (path metric: the sum of ln(1 + exp(-(1 - 2 u) LLR)) over the leaves). Each
    depth keeps one array of LLRs (and of left-half decisions) for the L paths; a path copied only moves its index
    (perm), and an array is read through it, so nothing is copied (Tal and Vardy). Returns the K bits (message + CRC)
    of every block's L paths, likeliest first (G x L x K)."""
    G, N, n = lam.shape[0], b["N"], b["n"]
    ninfo = np.concatenate([[0], np.cumsum(b["isinfo"])])
    alpha, beta = [None] * (n + 1), [None] * (n + 1)
    perm = np.tile(np.arange(L), (2 * (n + 1), G, 1))   # perm[d]: each path's row of alpha[d]; perm[n + 1 + d]: of beta[d]
    alpha[0] = np.repeat(lam[:, None, :], L, axis=1)
    PM = np.full((G, L), np.inf)
    PM[:, 0] = 0.0
    hist = []                                           # per information bit: each path's parent and bit

    def get(x, d):
        return np.take_along_axis(x, perm[d][:, :, None], axis=1)

    def fork(lr):
        nonlocal PM, perm
        cand = np.concatenate([PM + np.logaddexp(0, -lr), PM + np.logaddexp(0, lr)], axis=1)
        order = np.argsort(cand, axis=1, kind="stable")[:, :L]
        par, bit = order % L, (order // L).astype(np.uint8)
        PM = np.take_along_axis(cand, order, axis=1)
        perm = np.take_along_axis(perm, np.broadcast_to(par, perm.shape), axis=2)
        hist.append((par, bit))
        return bit[:, :, None]

    def frozen(a):                                      # the leaves' LLRs under a node with no information bit
        G_, L_, s = a.shape
        x = a.reshape(G_, L_, 1, s)
        while x.shape[-1] > 1:
            h = x.shape[-1] // 2
            a1, a2 = x[..., :h], x[..., h:]
            x = np.stack([_bplus(a1, a2), a1 + a2], axis=3).reshape(G_, L_, -1, h)
        return x.reshape(G_, L_, s)

    def node(d, off):
        nonlocal PM
        s = N >> d
        a = get(alpha[d], d)
        if ninfo[off + s] == ninfo[off]:                # no information bit below: all zeros; the paths pay for their LLRs
            PM = PM + np.logaddexp(0, -frozen(a)).sum(axis=2)
            return np.zeros((G, L, s), dtype=np.uint8)
        if s == 1:
            return fork(a[:, :, 0])
        h = s // 2
        alpha[d + 1], perm[d + 1] = _bplus(a[..., :h], a[..., h:]), np.arange(L)
        b1 = node(d + 1, off)
        if ninfo[off + h] > ninfo[off]:
            a = get(alpha[d], d)                        # the paths have moved
        alpha[d + 1], perm[d + 1] = a[..., h:] + (1.0 - 2.0 * b1) * a[..., :h], np.arange(L)
        beta[d], perm[n + 1 + d] = b1, np.arange(L)
        b2 = node(d + 1, off + h)
        return np.concatenate([get(beta[d], n + 1 + d) ^ b2, b2], axis=2)

    node(0, 0)
    idx = np.argsort(PM, axis=1, kind="stable")
    u = np.zeros((G, L, len(hist)), dtype=np.uint8)
    for k in range(len(hist) - 1, -1, -1):              # back through the parents
        par, bit = hist[k]
        u[:, :, k] = np.take_along_axis(bit, idx, axis=1)
        idx = np.take_along_axis(par, idx, axis=1)
    return u


def polar_scan(lam, b, iters=POLAR_SCAN):
    """Soft cancellation (SCAN, Fayyaz and Barry 2014), G blocks at once: SC with soft decisions, a node returning the
    LLRs that the code (its frozen bits) gives its bits, the right halves' from the last iteration reused. Returns the
    codeword's extrinsic LLRs (ln P(0)/P(1): what the code says about each bit, the bit's own LLR left out)."""
    G, N = lam.shape[0], b["N"]
    ninfo, right = np.concatenate([[0], np.cumsum(b["isinfo"])]), {}

    def node(a, off, s):
        k = ninfo[off + s] - ninfo[off]
        if k == 0:
            return np.full((G, s), POLAR_BIG)           # all frozen: known zeros
        if k == s:
            return np.zeros((G, s))                     # all free: the code says nothing
        h = s // 2
        a1, a2 = a[:, :h], a[:, h:]
        b1 = node(_bplus(a1, a2 + right.get((off, s), 0.0)), off, h)
        b2 = right[(off, s)] = node(a2 + _bplus(a1, b1), off + h, h)
        return np.concatenate([_bplus(b1, a2 + b2), b2 + _bplus(a1, b1)], axis=1)

    for _ in range(iters):
        ext = node(lam, 0, N)
    return ext


def polar_decode(llr, nbits, rate, soft=False):
    """CRC-aided list decoding; llr: the sent bits' LLRs (> 0 means 1), clipped to +-30. Returns the message bits and
    every block's CRC verdict; with soft (turbo equalization) also the decoder's a-posteriori LLRs of the sent bits
    (> 0 means 1): +-100 on the re-encoded message where the CRC holds, else soft cancellation's (repeated bits share
    their total)."""
    p, llr = polar_params(nbits, rate), np.asarray(llr, dtype=float)
    q0 = np.concatenate([[0], np.cumsum(p["E"])])
    bits, oks, post, groups = [None] * p["C"], np.zeros(p["C"], dtype=bool), np.zeros(p["G"]), {}
    for i, b in enumerate(p["blocks"]):
        groups.setdefault((p["a"][i], b["K"], b["E"]), []).append(i)
    for (a, K, E), ii in groups.items():                # blocks of the same size, all at once
        b = p["blocks"][ii[0]]
        lam = np.zeros((len(ii), b["N"]))               # the codeword's LLRs: repeated bits add up, punctured ones 0
        np.add.at(lam, (slice(None), b["sent"]), -np.clip(np.stack([llr[q0[i]:q0[i] + E] for i in ii]), -POLAR_CLIP, POLAR_CLIP))
        lam[:, b["known"]] = POLAR_BIG                  # shortened bits: known zeros
        u = polar_scl(lam, b)
        ok = ((u[..., :a].astype(np.int64) @ polar_crc_matrix(a).T) % 2 == u[..., a:]).all(axis=2)
        best = u[np.arange(len(ii)), np.where(ok.any(axis=1), ok.argmax(axis=1), 0), :a]   # the first whose CRC holds
        ext = polar_scan(lam, b) if soft and not ok.any(axis=1).all() else None
        for j, i in enumerate(ii):
            bits[i], oks[i] = best[j], ok[j].any()
            if soft:
                post[q0[i]:q0[i] + E] = (np.where(polar_codeword(best[j], b) == 1, POLAR_SURE, -POLAR_SURE) if oks[i]
                                         else -(lam[j] + ext[j])[b["sent"]])
    return np.concatenate(bits).astype(int), oks, post


LDPC_CODES = ["ldpc4", "ldpc3", "ldpc", "ldpc23", "ldpc34"]   # v7, header type 3: the code field (+ 4 with the CRC inverted)
POLAR_CODES = ["polar3", "polar", "polar23"]   # v7.3: header type 3, the code field 5 to 7
TYPE3_CODES = LDPC_CODES + POLAR_CODES
LDPC56 = "ldpc56"   # v7.6: LDPC 5/6, header type 0 or 1 with the code field 7 (single carrier, wide band: 2), see header_bits
SC_CODES = list(TURBO_R) + LDPC_CODES + [LDPC56] + POLAR_CODES   # the codes with soft output, which turbo equalization needs
WIDE_CODES = TYPE3_CODES + [LDPC56]   # v7.4: the codes of the single carrier's wide band (v7.6: and LDPC 5/6)
OFDM_NO56 = ["64qam", "256qam"]   # v7.6: OFDM's mappings without LDPC 5/6 (see build_frame)


def encode_fec(code, d):
    if code in POLAR_R:
        return interleave(polar_encode(d, POLAR_R[code]))
    if code in LDPC_R:
        return interleave(ldpc_encode(d, LDPC_R[code]))
    if code == "conv":
        return interleave(conv_encode(d))
    if code in TURBO_R:
        return interleave(turbo_encode(d, TURBO_R[code]))
    if code == "rs":
        return rs_encode(d)
    return np.asarray(d, dtype=int)


def decode_fec(code, n0, llr, lim=4):
    """llr: channel LLRs after descrambling (> 0 means 1). Returns (n0 message bits, failed RS blocks)."""
    llr = np.asarray(llr, dtype=float)
    if code in POLAR_R:
        return polar_decode(deinterleave(llr), n0, POLAR_R[code])[0], None
    if code in LDPC_R:
        return ldpc_decode(deinterleave(llr), n0, LDPC_R[code])[0], None
    if code == "conv":
        return conv_decode(deinterleave(llr), n0, lim), None
    if code in TURBO_R:
        return turbo_decode(deinterleave(llr), n0, lim, 8, TURBO_R[code])[0], None
    if code == "rs":
        bits, fails, _ = rs_decode((llr > 0).astype(int), n0)
        return bits, fails
    return (llr > 0).astype(int), None


# ---------------------------------------------------------------- soft output
def ln_i0(x):
    """ln I0(x), modified Bessel function (Abramowitz & Stegun 9.8.1-2)."""
    x = np.asarray(x, dtype=float)
    t = (np.minimum(x, 3.75) / 3.75) ** 2
    small = np.log(1 + t * (3.5156229 + t * (3.0899424 + t * (1.2067492 + t * (0.2659732 + t * (0.0360768 + t * 0.0045813))))))
    u = 3.75 / np.maximum(x, 3.75)
    large = np.maximum(x, 3.75) - 0.5 * np.log(np.maximum(x, 3.75)) + np.log(
        0.39894228 + u * (0.01328592 + u * (0.00225319 + u * (-0.00157565 + u * (0.00916281 + u * (-0.02057706
        + u * (0.02635537 + u * (-0.01647633 + u * 0.00392377))))))))
    return np.where(x < 3.75, small, large)


def calibrate(pE, S):
    """Level A (mu) and noise variance per real dimension (s2) of each tone, from its energies on the postamble:
    E[E^2] = A^2 + 2 s2 when the tone is sent (Rician), 2 s2 when it is not (Rayleigh).
    pE[m][j]: energy of tone m in postamble symbol j; S: the postamble's symbols (with 2 tones, its bits)."""
    pE, S = np.asarray(pE, dtype=float), np.asarray(S)

    def tone(e, on):
        p_on, p_off = np.mean(e[on] ** 2), np.mean(e[~on] ** 2)
        a2, s2 = max(p_on - p_off, 1e-6 * p_on), max(p_off / 2, 1e-12 * p_on)
        return dict(mu=np.sqrt(a2), s2=s2, snr_db=10 * np.log10(a2 / (2 * s2)))
    cal = tuple(tone(pE[m], S == m) for m in range(len(pE)))
    # The model assumes Gaussian noise. The echoes of short symbols have heavier tails, which makes its LLR too sure of
    # itself: shrink it by the factor alpha (at most 1) that best predicts the postamble's own bits.
    B = np.array([symbol_bits(int(v), int(math.log2(len(pE)))) for v in S])
    a = fit_scale((soft_llr("cal", cal, pE) * (2 * B - 1)).ravel())
    for t in cal:
        t["alpha"] = a
    return cal


def fit_scale(r):
    """Factor alpha <= 1 maximizing the likelihood of the known bits, P(bit) = 1 / (1 + e^(-alpha L)), on a grid from 1
    down to 0.05 (r: LLRs times the sign of the known bits)."""
    best, fb = 1.0, -np.inf
    for k in range(40):
        a = 0.05 ** (k / 39)
        f = -np.sum(np.logaddexp(0, -a * np.asarray(r)))
        if f > fb:
            best, fb = a, f
    return best


def calibrate_blind(E, rounds=8):
    """Without the postamble (not found, e.g. when the end of the sound was cut off): the same levels and noise,
    estimated on the message. Its symbols are unknown, so each tone counts as sent with the probability given by the
    current model, and the estimate is refined a few times (EM), starting from the hard decisions (strongest tone)."""
    E = np.asarray(E, dtype=float)

    def fit(e, w):                             # w[j]: probability that the tone is on in symbol j
        p_on = np.sum(w * e ** 2) / max(np.sum(w), 1e-9)
        p_off = np.sum((1 - w) * e ** 2) / max(np.sum(1 - w), 1e-9)
        a2, s2 = max(p_on - p_off, 1e-6 * p_on), max(p_off / 2, 1e-12 * p_on)
        return dict(mu=np.sqrt(a2), s2=s2, snr_db=10 * np.log10(a2 / (2 * s2)))
    W = np.zeros_like(E)
    W[np.argmax(E, axis=0), np.arange(E.shape[1])] = 1
    for r in range(rounds + 1):
        cal = tuple(fit(E[m], W[m]) for m in range(len(E)))
        if r < rounds:                         # probability of each tone: softmax of the symbol metrics
            lam = symbol_metrics(cal, E)
            W = np.exp(lam - lam.max(axis=0))
            W /= W.sum(axis=0)
    return cal


def fine_timing(read, nc, R):
    """Fine timing at the message's rate, with tone hopping. The symbol instants come from the 400 bit/s preamble and
    the postamble. The preamble's long symbols tolerate an error of several samples, the short symbols of a fast
    message do not: in a real room their best instant was a quarter of a symbol away at 1600 bit/s. So, in up to 8
    blocks of the message, find the shift (at most R samples) that makes the tones most distinct, mean contrast, and
    fit a straight line through these shifts. read(s): energies of the tones of all message symbols, symbol j read
    s[j] samples after its instant. Returns the shift of each symbol (all 0 unless it makes the tones clearly more
    distinct). Used with tone hopping only: with fixed tones the echo of the earlier bits lies on the same tones, and
    in simulated rooms this measure then often favoured a worse instant."""
    con = contrast
    K, step, cache = max(1, min(8, nc // 128)), max(1, round(R / 6)), {}

    def c_at(s):
        if s not in cache:
            cache[s] = con(read(np.full(nc, s)))
        return cache[s]
    pts = []
    for k in range(K):
        j0, j1 = k * nc // K, (k + 1) * nc // K
        score = lambda s: c_at(s)[j0:j1].mean()
        s0 = max(range(-R, R + 1, step), key=score)          # coarse grid, then sample by sample
        s0 = max([s0] + [s for s in range(max(-R, s0 - step + 1), min(R, s0 + step - 1) + 1) if s != s0], key=score)
        pts.append(((j0 + j1 - 1) / 2, s0))
    x, sh = np.array(pts, dtype=float).T
    b = np.sum((x - x.mean()) * (sh - sh.mean())) / np.sum((x - x.mean()) ** 2) if len(x) > 1 else 0.0
    shift = np.clip(np.floor(sh.mean() + b * (np.arange(nc) - x.mean()) + 0.5), -2 * R, 2 * R).astype(int)
    d = con(read(shift)) - c_at(0)                          # kept if the tones are clearly more distinct (paired t-test)
    return shift if d.mean() > 2 * d.std(ddof=1) / np.sqrt(nc) else np.zeros(nc, dtype=int)


def level_track(cal, R, E, W=24):
    """Signal level (whichever tone is on) averaged over +-W symbols, relative to the energies R (the postamble)."""
    mu2 = np.array([t["mu"] ** 2 for t in cal])[:, None]
    q = lambda X: np.max(np.asarray(X, dtype=float) ** 2 / mu2, axis=0)
    qref = np.mean(q(R))
    cs = np.concatenate([[0], np.cumsum(q(E))])
    n = np.shape(E)[1]
    j = np.arange(n)
    a, b = np.maximum(0, j - W), np.minimum(n, j + W + 1)
    return np.sqrt(np.clip((cs[b] - cs[a]) / (b - a) / qref, 0.01, 4))


def symbol_metrics(cal, E, g=1.0):
    """Log-likelihood of each tone being the one sent, up to a term common to all: the tone sent is Rician, the
    others Rayleigh, so lambda_m = ln I0(A_m E_m / s2_m) - A_m^2 / 2 s2_m. In a room most of the disturbance is the
    echo of the signal, so it follows the signal level g: the energies are divided by g before using the
    postamble's levels and noise."""
    E = np.asarray(E, dtype=float) / g
    return np.array([ln_i0(t["mu"] * e / t["s2"]) - t["mu"] ** 2 / (2 * t["s2"]) for t, e in zip(cal, E)])


def soft_llr(soft, cal, E, g=1.0):
    """Soft output of the K bits of each symbol, shape (symbols, K); E[m]: energies of tone m.
    'cal': LLR, ln of the sum of e^lambda over the tones whose bit is 1 minus the same over the tones whose bit is 0
    (with 2 tones, lambda_1 - lambda_0, the LLR of non-coherent FSK).  'ratio': the strongest tone whose bit is 1
    against the strongest whose bit is 0 (with 2 tones, ln(E1/E0))."""
    E = np.asarray(E, dtype=float)
    sets = bit_sets(len(E))
    if soft == "ratio":
        return np.stack([np.log((E[s1].max(axis=0) + 1e-9) / (E[s0].max(axis=0) + 1e-9)) for s0, s1 in sets], axis=1)
    lam = symbol_metrics(cal, E, g)
    return cal[0].get("alpha", 1.0) * np.stack([np.logaddexp.reduce(lam[s1], axis=0) - np.logaddexp.reduce(lam[s0], axis=0)
                                                for s0, s1 in sets], axis=1)


def post_errors(pE, S):
    """Bit errors on the known postamble, from the hard decisions (strongest tone)."""
    K = int(math.log2(len(pE)))
    return int(sum(np.sum(np.array(symbol_bits(int(h), K)) != np.array(symbol_bits(int(v), K))) for h, v in zip(np.argmax(pE, axis=0), S)))


def contrast(E):
    """How distinct the tones of each symbol are: (strongest - second) / sum; |E1 - E0| / (E1 + E0) with 2 tones."""
    E = np.asarray(E, dtype=float)
    top = np.sort(E, axis=0)
    return (top[-1] - top[-2]) / (E.sum(axis=0) + 1e-12)


def known_contrast(E, s):
    """How much the known tone s stands out: 1 if only it is heard, about 0 in noise; +-(E1 - E0) / (E1 + E0) with 2 tones."""
    E = np.asarray(E, dtype=float)
    if len(E) == 2:
        return (2 * s - 1) * ((E[1] - E[0]) / (E[1] + E[0] + 1e-12))
    t = E.sum(axis=0)
    return (len(E) * E[s] - t) / ((len(E) - 1) * t + 1e-12)


# ---------------------------------------------------------------- transmitter
def crc8(bits):
    """CRC-8, generator x^8 + x^2 + x + 1."""
    r = 0
    for b in bits:
        fb = ((r >> 7) & 1) ^ int(b)
        r = (r << 1) & 0xFF
        if fb:
            r ^= 0x07
    return msb(r, 8)


def header_bits(typ, a, b, code, rate, tones, M=2, ofdm=None, omap="dqpsk", sc=False, nt=None, wide=False):
    if code in TYPE3_CODES:                      # v7: type 3, the LDPC (v7.3: or polar) code's index in the code field (b = 0: a text)
        typ, c = 3, TYPE3_CODES.index(code)
    elif code == LDPC56:                         # v7.6: type 0 or 1, the code field 7 (the single carrier's wide band: 2)
        c = 2 if sc and wide else 7
    else:
        c = CODES.index(code)                    # codes 4 to 7: field c mod 4 and the CRC inverted
    if sc:                                       # v5, single carrier: tone set 3, rate field 2, the training in a's and b's top bits
        t, four, r, T = 3, True, SC_FIELD, SC_TCODE[nt or SC_NC]
        a, b = a | (T >> 1) << 10, b | (T & 1) << 10
        if wide and code != LDPC56:              # v7.4: type 2, read as type 3
            if typ != 3:
                raise ValueError(f"the wide band takes an LDPC or a polar code: one of {WIDE_CODES}")
            typ = 2
    elif ofdm is not None:                         # OFDM: tones 3, tone set 3, rate field = the mode (0 long guard, 1 short)
        t, four, r = 3, True, OFDM_MODES.index(ofdm)
    else:
        t, four = TONE_SETS.index(tones), M == 4
        r = 3 if rate == RATE_WIDE else (RATES4 if four else RATES).index(rate)
    if rate == RATE_WIDE:
        t = 3                                    # 6400 bit/s: tone set 3, rate field 3
    v = (msb(typ, 2) + msb(a, 11) + msb(b, 11) + msb(c & 3, 2) + msb(r, 2)
         + msb(3 if four else t, 2) + msb(t if four else 0, 2))   # 4 tones: tones 3, the tone set in the last two bits
    slot = SC_SLOT.get(omap, omap) if sc else omap   # v7.5: 1024QAM in 8PSK's place
    mask = msb((0xFF if c >> 2 else 0) ^ (OFDM_MASKS[OFDM_MAPS.index(slot)] if ofdm is not None or sc else 0), 8)
    return encode_fec("conv", np.array(v + [x ^ k for x, k in zip(crc8(v), mask)]))


# ---------------------------------------------------------------- v7.1: JPEG pictures
# A JPEG travels as its file, then a CRC-32 of the file (as in Ethernet or Wi-Fi), so that every receiver can tell an
# intact file from a corrupted one whatever the code; the file is shown as it arrives. The header marks it as an image
# taller than 512 rows, which receivers before v7.1 reject: a = 513 + (n - 1) // 512, b = 1 + (n - 1) % 512 for a file
# of n bytes (at most JPEG_MAX). jpeg_encode makes the same files as the page: baseline, with the example tables of the
# standard scaled by the quality as libjpeg does, and a restart marker after every row of blocks. Huffman codes have
# variable lengths and the DC values are sent as differences, so a wrong bit garbles the picture up to the next restart
# marker; with a marker on every row, a wrong bit costs one stripe of 16 (grey: 8) rows instead of the rest of the picture.
JPEG_MAX = 20480


def jpeg_fields(n):
    """The header's a and b for a JPEG file of n bytes."""
    return 513 + (n - 1) // 512, 1 + (n - 1) % 512


def jpeg_payload(f):
    """The message bits: the file, then its CRC-32, each byte most significant bit first."""
    f = bytes(f)
    return np.unpackbits(np.frombuffer(f + zlib.crc32(f).to_bytes(4, "big"), dtype=np.uint8)).astype(int)


def jpeg_unpack(bits, n):
    """The received file (as it is, errors included) and whether its CRC-32 matches."""
    u = np.packbits(np.asarray(bits[:8 * (n + 4)], dtype=np.uint8)).tobytes()
    return u[:n], zlib.crc32(u[:n]) == int.from_bytes(u[n:n + 4], "big")



def jpeg_resync(f):
    """The file with every marker that cannot stand in the coded data (0xFF followed by anything but 0 or a restart
    marker) turned into data. Wrong bits make such markers; a decoder (libjpeg, so Chrome, Firefox and Pillow) stops at
    the first one and leaves the rest of the picture grey. Without them it reads the corrupted stripes as garbage, then
    starts again at the next restart marker: only the stripes that held wrong bits are lost."""
    f, i = bytearray(f), 2
    while i + 4 <= len(f) and f[i] == 0xFF and f[i + 1] not in (0xDA, 0xD9):   # the header's segments, up to the scan
        i += 2 + (f[i + 2] << 8 | f[i + 3])
    if i + 4 > len(f) or f[i] != 0xFF or f[i + 1] != 0xDA:
        return bytes(f)
    for j in range(i + 2 + (f[i + 2] << 8 | f[i + 3]), len(f) - 3):   # the coded data, up to the end of image marker
        if f[j] == 0xFF and f[j + 1] and not 0xD0 <= f[j + 1] <= 0xD7:
            f[j + 1] = 0
    return bytes(f)

# The example tables of the JPEG standard (ITU-T T.81, Annex K): quantization (row by row) and Huffman codes (number of
# codes of each length 1 to 16, then the symbols): DC luminance, AC luminance, DC chrominance, AC chrominance
JPEG_QL = [16, 11, 10, 16, 24, 40, 51, 61, 12, 12, 14, 19, 26, 58, 60, 55, 14, 13, 16, 24, 40, 57, 69, 56, 14, 17, 22, 29, 51, 87, 80, 62,
           18, 22, 37, 56, 68, 109, 103, 77, 24, 35, 55, 64, 81, 104, 113, 92, 49, 64, 78, 87, 103, 121, 120, 101, 72, 92, 95, 98, 112, 100, 103, 99]
JPEG_QC = [17, 18, 24, 47, 99, 99, 99, 99, 18, 21, 26, 66, 99, 99, 99, 99, 24, 26, 56, 99, 99, 99, 99, 99, 47, 66, 99, 99, 99, 99, 99, 99] + [99] * 32
JPEG_HUFF = [bytes.fromhex(s) for s in (
    "00010501010101010100000000000000000102030405060708090a0b",
    "0002010303020403050504040000017d01020300041105122131410613516107227114328191a1082342b1c11552d1f02433627282090a161718191a25262728292a3435363738393a434445464748494a535455565758595a636465666768696a737475767778797a838485868788898a92939495969798999aa2a3a4a5a6a7a8a9aab2b3b4b5b6b7b8b9bac2c3c4c5c6c7c8c9cad2d3d4d5d6d7d8d9dae1e2e3e4e5e6e7e8e9eaf1f2f3f4f5f6f7f8f9fa",
    "00030101010101010101010000000000000102030405060708090a0b",
    "00020102040403040705040400010277000102031104052131061241510761711322328108144291a1b1c109233352f0156272d10a162434e125f11718191a262728292a35363738393a434445464748494a535455565758595a636465666768696a737475767778797a82838485868788898a92939495969798999aa2a3a4a5a6a7a8a9aab2b3b4b5b6b7b8b9bac2c3c4c5c6c7c8c9cad2d3d4d5d6d7d8d9dae2e3e4e5e6e7e8e9eaf2f3f4f5f6f7f8f9fa")]
JPEG_ZZ = [0, 1, 8, 16, 9, 2, 3, 10, 17, 24, 32, 25, 18, 11, 4, 5, 12, 19, 26, 33, 40, 48, 41, 34, 27, 20, 13, 6, 7, 14, 21, 28,
           35, 42, 49, 56, 57, 50, 43, 36, 29, 22, 15, 23, 30, 37, 44, 51, 58, 59, 52, 45, 38, 31, 39, 46, 53, 60, 61, 54, 47, 55, 62, 63]
JPEG_C16 = [1, 0.9807852804032304, 0.9238795325112867, 0.8314696123025452, 0.7071067811865476, 0.5555702330196022, 0.3826834323650898, 0.19509032201612825, 0]   # cos(k pi / 16)
JPEG_COS = np.array([[(0.5 if u else 0.7071067811865476 / 2) * (lambda a: JPEG_C16[a] if a <= 8 else -JPEG_C16[16 - a] if a <= 16 else -JPEG_C16[a - 16] if a <= 24
                      else JPEG_C16[32 - a])((2 * x + 1) * u % 32) for x in range(8)] for u in range(8)])   # C(u) / 2 cos((2x + 1) u pi / 16): the page's numbers


def jpeg_codes(t):
    """The code of each symbol, as a string of bits (T.81, Annex C)."""
    code, c, k = {}, 0, 16
    for n in range(1, 17):
        for _ in range(t[n - 1]):
            code[t[k]] = format(c, f"0{n}b")
            c, k = c + 1, k + 1
        c <<= 1
    return code


def jpeg_encode(rgb, q=50, grey=False):
    """A baseline JPEG (bytes) of the picture rgb (rows x columns x 3, 0 to 255) at quality q (1 to 100): colour with
    the chrominance at half the width and height (4:2:0), or grey (luminance only)."""
    q = int(min(100, max(1, round(q))))
    s = 5000 // q if q < 50 else 200 - 2 * q
    QT = [np.clip((np.array(T) * s + 50) // 100, 1, 255) for T in ([JPEG_QL] if grey else [JPEG_QL, JPEG_QC])]
    h, w = rgb.shape[:2]
    m = 8 if grey else 16
    mx, my = -(-w // m), -(-h // m)
    P = np.asarray(rgb, dtype=float)[np.minimum(np.arange(my * m), h - 1)][:, np.minimum(np.arange(mx * m), w - 1)]   # the last row and column repeated
    r, g, b = P[..., 0], P[..., 1], P[..., 2]
    planes = [0.299 * r + 0.587 * g + 0.114 * b - 128]
    if not grey:
        planes += [-0.168736 * r - 0.331264 * g + 0.5 * b, 0.5 * r - 0.418688 * g - 0.081312 * b]
        planes[1:] = [0.25 * (C[0::2, 0::2] + C[0::2, 1::2] + C[1::2, 0::2] + C[1::2, 1::2]) for C in planes[1:]]

    def coefs(X, T):                                               # the quantized DCT of each 8 x 8 block, in zigzag order
        B = X.reshape(X.shape[0] // 8, 8, X.shape[1] // 8, 8).transpose(0, 2, 1, 3)
        R = sum((B[..., x:x + 1] * JPEG_COS[:, x] for x in range(1, 8)), B[..., 0:1] * JPEG_COS[:, 0])   # rows, then columns, the
        F = sum((JPEG_COS[:, y:y + 1] * R[..., y:y + 1, :] for y in range(1, 8)), JPEG_COS[:, 0:1] * R[..., 0:1, :])   # terms added in the page's order
        return np.floor(F / T.reshape(8, 8) + 0.5).astype(int).reshape(B.shape[0], B.shape[1], 64)[:, :, JPEG_ZZ]
    Z = [coefs(X, QT[min(i, 1)]) for i, X in enumerate(planes)]
    seg = lambda mk, body: bytes([0xFF, mk, (len(body) + 2) >> 8, (len(body) + 2) & 255]) + bytes(body)
    out = b"\xff\xd8" + seg(0xE0, [0x4A, 0x46, 0x49, 0x46, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0])
    out += seg(0xDB, [v for k, T in enumerate(QT) for v in [k] + [int(T[z]) for z in JPEG_ZZ]])
    out += seg(0xC0, [8, h >> 8, h & 255, w >> 8, w & 255] + ([1, 1, 0x11, 0] if grey else [3, 1, 0x22, 0, 2, 0x11, 1, 3, 0x11, 1]))
    out += seg(0xC4, [v for k in ((0, 1) if grey else (0, 1, 2, 3)) for v in [((k & 1) << 4) | (k >> 1)] + list(JPEG_HUFF[k])])
    out += seg(0xDD, [mx >> 8, mx & 255])                        # a restart interval of one row of blocks
    out += seg(0xDA, [1, 1, 0, 0, 63, 0] if grey else [3, 1, 0, 2, 0x11, 3, 0x11, 0, 63, 0])
    DC, AC = [jpeg_codes(JPEG_HUFF[0]), jpeg_codes(JPEG_HUFF[2])], [jpeg_codes(JPEG_HUFF[1]), jpeg_codes(JPEG_HUFF[3])]
    val = lambda v: format((v - 1 if v < 0 else v) & ((1 << abs(v).bit_length()) - 1), f"0{abs(v).bit_length()}b") if v else ""
    for row in range(my):
        bits, pred = [], [0, 0, 0]
        for k in range(mx):
            blocks = [(0, Z[0][row, k])] if grey else [(0, Z[0][2 * row + i, 2 * k + j]) for i in (0, 1) for j in (0, 1)] + [(1, Z[1][row, k]), (2, Z[2][row, k])]
            for c, z in blocks:
                t, d = min(c, 1), int(z[0]) - pred[c]
                pred[c] = int(z[0])
                bits += [DC[t][abs(d).bit_length()], val(d)]
                run = 0
                for v in z[1:].tolist():
                    if not v:
                        run += 1
                        continue
                    while run > 15:
                        bits.append(AC[t][0xF0])
                        run -= 16
                    bits += [AC[t][(run << 4) | abs(v).bit_length()], val(v)]
                    run = 0
                if run:
                    bits.append(AC[t][0])
        s = "".join(bits)
        s += "1" * (-len(s) % 8)                                  # the row's last byte filled with 1s
        out += int(s, 2).to_bytes(len(s) // 8, "big").replace(b"\xff", b"\xff\x00") if s else b""
        if row < my - 1:
            out += bytes([0xFF, 0xD0 + (row & 7)])                # restart marker RSTr mod 8
    return out + b"\xff\xd9"


def load_rgb(path, width=None):
    """The picture's colours (rows x columns x 3), scaled down (pixels averaged) to the given width if it is wider."""
    from PIL import Image
    im = Image.open(path).convert("RGB")
    if width and im.width > width:
        im = im.resize((width, max(1, round(im.height * width / im.width))), Image.BOX)
    return np.asarray(im)


def build_frame(text=None, img=None, code="none", rate=400, tones="low", pre=127, sound=None, mod=2, guard="long", omap="dqpsk",
                train="medium", control="short", jpeg=None, wide=False):
    """Return (segments, data): a list of (bits, tones, rate[, index of the first hop[, M]]), then the message bits;
    with 4 tones (mod=4) the message segment holds tones 0..3 instead of bits, and M = 4. With mod="ofdm" the message
    segment is the OFDM part (samples at 48 kHz, None, None), guard "long" or "short" (see OFDM_MODES), omap one of OFDM_MAPS.
    With mod="sc" (v5) it is the single-carrier part, omap one of SC_MAPS, code one of SC_CODES, train one of SC_TRAIN;
    wide (v7.4): in the wide band, 2.3 to 20.5 kHz (an LDPC or a polar code).
    img: 2-D boolean array (True = white). jpeg (v7.1): the bytes of a JPEG file, sent as they are. code: none, conv, turbo, turbo3, turbo4, an LDPC code (LDPC_R) or a polar code (POLAR_R).
    sound = (order, periods): room sounding frame; its MLS segment is (samples at 48 kHz, None, None).
    control (v6), OFDM and single carrier only: "short" (the short control part, c6_wave; receivers older than v6 do
    not read it) or "long" (the 400 bit/s one, with the preamble pre; it reaches farther)."""
    if sound is not None:
        order, periods = sound
        warm = [(i + 1) % 2 for i in range(round(0.1 * CTRL["rate"]))]
        ctl = "hop" if tones == "hop" else CTRL["tones"]   # control part and postamble: hopping or 4 / 8 kHz
        hdr = header_bits(2, order, periods, "none", CTRL["rate"], ctl)
        ctrl = np.concatenate([warm, PILOTS[pre], np.tile(hdr, HREPS[pre])]).astype(int)
        wave = SOUND_AMP * (2.0 * np.tile(mls(order, SOUND_TAPS[order]), periods + 1) - 1)
        return [(ctrl, ctl, CTRL["rate"], -len(warm)), (wave, None, None), (PILOTS[pre], ctl, CTRL["rate"])], np.zeros(0, dtype=int)
    if text is not None:
        codes = [ord(c) if ord(c) < 256 else 63 for c in text[:1023]]
        typ, a, b = 0, len(codes), 0
        data = np.array([x for c in codes for x in msb(c, 8)], dtype=int)
    elif jpeg is not None:                                 # v7.1: a JPEG file, then its CRC-32
        if not 1 <= len(jpeg) <= JPEG_MAX:
            raise ValueError(f"a JPEG file of {len(jpeg)} bytes: at most {JPEG_MAX}")
        typ, (a, b) = 1, jpeg_fields(len(jpeg))
        data = jpeg_payload(jpeg)
    else:
        typ, (a, b) = 1, img.shape
        data = img.astype(int).flatten(order="F")
    warm = [(i + 1) % 2 for i in range(round(0.1 * CTRL["rate"]))]
    if mod == "sc":                                        # v5: a hopping control part, then the single carrier
        if omap not in SC_MAPS or code not in SC_CODES:
            raise ValueError(f"single carrier: mapping one of {SC_MAPS}, code one of {SC_CODES} (turbo equalization)")
        sc_set_band(wide)                                  # v7.4: the usual band or the wide one
        nt = SC_TRAIN[train]
        hdr = header_bits(typ, a, b, code, None, None, omap=omap, sc=True, nt=nt, wide=wide)
        part = (sc_part(sc_channel_bits(code, data, omap)[0], omap, nt), None, None)
        if control == "short":                             # v6
            return [(c6_wave(hdr), None, None), part], data
        ctrl = np.concatenate([warm, PILOTS[pre], np.tile(hdr, HREPS[pre])]).astype(int)
        return [(ctrl, "hop", CTRL["rate"], -len(warm)), part], data
    if mod == "ofdm":                                      # OFDM: a hopping control part, then the OFDM symbols
        if omap not in OFDM_MAPS:
            raise ValueError(f"OFDM: mapping one of {OFDM_MAPS}")
        if omap in OFDM_NO56 and code == LDPC56:         # v7.6: OFDM's own noise (its clipping, about 21 dB below the signal) is too loud
            raise ValueError("OFDM with 64QAM or 256QAM: LDPC 3/4 at most (OFDM's clipping leaves 64QAM 3 to 5 % of wrong bits on a short "
                             "message, more than LDPC 5/6 corrects, and 256QAM with LDPC 5/6 needs 22.5 dB)")
        mode = ofdm_mode(guard)
        hdr = header_bits(typ, a, b, code, None, None, ofdm=mode, omap=omap)
        ch = encode_fec(code, data)
        ctrl = np.concatenate([warm, PILOTS[pre], np.tile(hdr, HREPS[pre])]).astype(int)
        n = len(ch) if omap == "dqpsk" else ofdm_layout(mode, omap, len(ch))[1]   # coherent: fill the last symbol up
        b = np.zeros(n, dtype=int)                         # with scrambled zeros (plain zeros, the same point on many
        b[:len(ch)] = ch                                   # subcarriers, would add up to a peak)
        part = (ofdm_symbols(b ^ prbs(n), mode, omap), None, None)
        if control == "short":                             # v6
            return [(c6_wave(hdr), None, None), part], data
        return [(ctrl, "hop", CTRL["rate"], -len(warm)), part], data
    M = 4 if mod == 4 else 2
    if rate not in (RATES4 + [RATE_WIDE] if M == 4 else RATES):
        raise ValueError(f"rate {rate} bit/s not available with {M} tones: {RATES4 + [RATE_WIDE] if M == 4 else RATES}")
    if rate == RATE_WIDE and tones != "hop":
        raise ValueError(f"{RATE_WIDE} bit/s uses one group of 4 tones over the whole band: choose tones hop")
    hdr = header_bits(typ, a, b, code, rate, tones, M)
    ch = encode_fec(code, data)
    ctrl = np.concatenate([warm, PILOTS[pre], np.tile(hdr, HREPS[pre])]).astype(int)
    body = np.concatenate([ch ^ prbs(len(ch)), PILOTS[pre]]).astype(int)
    ctl = "hop" if tones == "hop" else CTRL["tones"]       # a hopping message has a hopping control part
    if M == 4:                                             # 2 bits per symbol, message and postamble each padded
        syms = np.concatenate([to_symbols(ch ^ prbs(len(ch)), M), to_symbols(PILOTS[pre], M)])
        return [(ctrl, ctl, CTRL["rate"], -len(warm)), (syms, tones, rate, 0, M)], data
    return [(ctrl, ctl, CTRL["rate"], -len(warm)), (body, tones, rate)], data


def modulate(segs, fs=FS_TX):
    out = [np.zeros(round(0.2 * fs))]
    for seg in segs:
        bits, tones, rate = seg[:3]                           # bits: the tone of each symbol (with 2 tones, its bit)
        M = seg[4] if len(seg) > 4 else 2
        if tones is None:                                     # a waveform at 48 kHz (room sounding)
            out.append(bits if fs == FS_TX else np.interp(np.arange(round(len(bits) * fs / FS_TX)) * FS_TX / fs, np.arange(len(bits)), bits))
            continue
        t = np.arange(round(fs * math.log2(M) / rate)) / fs
        if tones == "hop":                                    # tone hopping: symbol j on the pair (group) of k0 + j (k0 < 0: warm-up)
            F = np.array(hop_tones(rate, (seg[3] if len(seg) > 3 else 0) + np.arange(len(bits)), M))
            out.append((0.9 * np.sin(2 * np.pi * F[np.asarray(bits), np.arange(len(bits))][:, None] * t)).ravel())
            continue
        w = 0.9 * np.stack([np.sin(2 * np.pi * f * t) for f in tone_set(tones, M)])
        out.append(w[bits].ravel())
    out.append(np.zeros(round(0.2 * fs)))
    return np.concatenate(out)


def load_image(path, step=1, up=1):
    """Black (False) or white (True) pixels around half the brightest grey: every step-th pixel, or up = 2 for twice
    the width and height (the grey levels interpolated between pixels, then the threshold: smooth edges). The
    page's demo images are greatwave_80x54.png, greatwave_160x108.png and greatwave_320x216.png (Hokusai's Great Wave,
    public domain); before 2026-10-10 they were bartS.png with step 2, step 1 and up 2."""
    from PIL import Image
    I = np.asarray(Image.open(path).convert("RGB"), dtype=float).mean(axis=2)[::step, ::step]
    while up > 1:
        I, up = upsample2(upsample2(I).T).T, up // 2
    return I > I.max() / 2


def upsample2(I):
    """Twice as many rows, interpolated linearly between rows (row r of the result sits at r / 2 - 1/4)."""
    n = I.shape[0]
    t = np.arange(2 * n) / 2 - 0.25
    r0 = np.clip(np.floor(t).astype(int), 0, n - 1)
    f = np.where(t < 0, 0, t - np.floor(t))[:, None]
    return I[r0] * (1 - f) + I[np.minimum(r0 + 1, n - 1)] * f


def write_wav(path, x, fs):
    with wave.open(path, "wb") as w:
        w.setnchannels(1); w.setsampwidth(2); w.setframerate(fs)
        w.writeframes((np.clip(x, -1, 1) * 32767).astype("<i2").tobytes())


def read_wav(path):
    """16-, 24- or 32-bit PCM WAV, or 32- or 64-bit float WAV (v7: the page's telemetry uploads) -> (first channel as
    float, fs). Reads the RIFF chunks itself: the wave module knows integer PCM only."""
    with open(path, "rb") as f:
        b = f.read()
    if b[:4] != b"RIFF" or b[8:12] != b"WAVE":
        raise ValueError(f"{path}: not a WAV file")
    u = lambda a, n: int.from_bytes(b[a:a + n], "little")
    p, fmt, data = 12, None, None
    while p + 8 <= len(b):
        cid, n = b[p:p + 4], u(p + 4, 4)
        if cid == b"fmt ":
            tag, ch, fs, bits = u(p + 8, 2), u(p + 10, 2), u(p + 12, 4), u(p + 22, 2)
            if tag == 0xFFFE:                                 # WAVE_FORMAT_EXTENSIBLE: the sub-format's first two bytes
                tag = u(p + 32, 2)
            fmt = (tag, ch, fs, bits)
        elif cid == b"data":
            data = b[p + 8:p + 8 + n]
            break
        p += 8 + n + (n & 1)
    if fmt is None or data is None:
        raise ValueError(f"{path}: no fmt or data chunk")
    tag, ch, fs, bits = fmt
    sw = bits // 8
    data = data[:len(data) // (sw * ch) * sw * ch]
    if tag == 3 and sw in (4, 8):                             # IEEE float
        return np.frombuffer(data, dtype="<f4" if sw == 4 else "<f8").reshape(-1, ch)[:, 0].astype(float), fs
    if tag == 1 and sw == 3:                                  # 24-bit PCM: three bytes, little-endian, signed
        a = np.frombuffer(data, dtype=np.uint8).reshape(-1, ch, 3)[:, 0, :].astype(np.int32)
        v = a[:, 0] | (a[:, 1] << 8) | (a[:, 2] << 16)
        return np.where(v >= 1 << 23, v - (1 << 24), v) / float(2 ** 23), fs
    if tag == 1 and sw in (2, 4):
        y = np.frombuffer(data, dtype={2: "<i2", 4: "<i4"}[sw]).reshape(-1, ch)[:, 0]
        return y / float(2 ** (8 * sw - 1)), fs
    raise ValueError(f"{path}: format {tag}, {bits} bits not supported")


# ---------------------------------------------------------------- tone hopping
# At R bit/s there are H = 6400 / R pairs of adjacent tones, bit 0 -> 3200 + 2pR Hz and bit 1 -> 3200 + (2p + 1)R Hz,
# p = 0 .. H - 1 (3.2 to 16 kHz). The pairs started at 1.6 kHz at first, but in a real room that tone, low in a small
# speaker's range and close to the room's noise, caused nearly all the errors. Symbol k of a segment (counted from the first preamble bit, from the first message bit
# through the postamble, or from the first bit of a sounding's postamble) uses pair p = s k mod H, s = H/2 - 1 (1 if H = 4):
# a pair comes back only every H symbols, so the echo it hears is the one left after H bits of decay.
# v2, 4 tones: the same H groups, each now of 4 tones spaced by the symbol rate R / 2, so they fill the same band.
def hop_set(rate):
    H = HOP_BAND // (2 * rate)
    return H, (H // 2 - 1 if H > 4 else 1)


def hop_pair(rate, p, M=2):
    """The M tones of pair (group) p; with 2 tones, (bit-0 tone, bit-1 tone)."""
    p, K = np.asarray(p), int(math.log2(M))
    return tuple(HOP_F0 + 2 * p * rate + m * rate // K for m in range(M))


def hop_tones(rate, k, M=2):
    """The M tones of the symbols k of a hopping segment (tone m of each)."""
    H, s = hop_set(rate)
    return hop_pair(rate, (s * np.asarray(k)) % H, M)


def tone_desc(tones, rate, M=2):
    if tones != "hop":
        return "tones " + " / ".join(f"{f / 1000:g}" for f in tone_set(tones, M)) + " kHz"
    H, _ = hop_set(rate)
    if H == 1:
        return "tones " + " / ".join(f"{f / 1000:g}" for f in hop_pair(rate, 0, M)) + " kHz (the whole band, no hopping)"
    top = hop_pair(rate, H - 1, M)[-1]
    return f"tone hopping, {H} {'pairs' if M == 2 else f'groups of {M} tones'} from {HOP_F0 / 1000:g} to {top / 1000:g} kHz"


def sliding_energy(y, fs, F, Lw):
    """|sum of y(n) exp(-j 2 pi F n / fs)| over a window of Lw samples starting at each sample."""
    n = np.arange(len(y))
    c = np.concatenate([[0], np.cumsum(y * np.exp(-2j * np.pi * F * n / fs))])
    return np.abs(c[Lw:] - c[:-Lw])


def energies_at(y, fs, f, pos, Lw):
    """Same energy for one window per start m in pos, the tone f[i] for window i."""
    n, f, pos = np.arange(Lw), np.asarray(f, dtype=float), np.asarray(pos, dtype=int)
    out = np.empty(len(pos))
    for a in range(0, len(pos), 2048):
        b = slice(a, a + 2048)
        out[b] = np.abs(np.sum(y[pos[b, None] + n] * np.exp(-2j * np.pi * f[b, None] * n / fs), axis=1))
    return out


def hop_energies(y, fs, rate, k0, pos, Lw, M=2):
    """Energies of the M tones of the symbols k0, k0 + 1, ... of a hopping segment, at the window starts pos."""
    return tuple(energies_at(y, fs, f, pos, Lw) for f in hop_tones(rate, k0 + np.arange(len(pos)), M))


def hop_normalize(rate, E):
    """Speaker, room and microphone give each tone its own level: divide each tone's energies (E[m]: symbols 0, 1,
    2 ... of a hopping segment) by its level, so that one model fits all the tones. With 2 tones, that level is the rms
    over the symbols that used its pair (half of them carry the tone, half only the echo). With 4 tones a tone is sent
    in only a quarter of the symbols of its group, too few in a short frame for that rms to be fair (a tone never sent
    in a group would be raised to the level of the others): its level is the rms over the symbols where it is the
    strongest, or its group's if it never is. Also returns the spread of the tone levels in dB."""
    E = np.asarray(E, dtype=float)
    M, n = E.shape
    H, s = hop_set(rate)
    p = (s * np.arange(n)) % H
    cnt = np.bincount(p, minlength=H)
    if M == 2:
        G = np.array([np.sqrt(np.bincount(p, e ** 2, H) / np.maximum(cnt, 1)) for e in E])
        G[G == 0] = 1.0
        used = np.tile(cnt > 0, (M, 1))
    else:
        w, c, S = np.argmax(E, axis=0), np.zeros((M, H)), np.zeros((M, H))
        np.add.at(c, (w, p), 1)
        np.add.at(S, (w, p), E[w, np.arange(n)] ** 2)
        grp = np.sqrt(S.sum(axis=0) / np.maximum(cnt, 1))
        grp[grp == 0] = 1.0
        G = np.sqrt(S / np.maximum(c, 1))
        G[G == 0] = 1.0
        G = np.where(c > 0, G, grp[None, :])
        used = c > 0
    lv = G[used]
    return E / G[:, p], 20 * np.log10(lv.max() / lv.min())


def hop_correlation(y, fs, Lr):
    """Preamble correlation of a hopping control part at every sample: for each pair, the contrast (E1 - E0) / (E1 + E0)
    of its two tones, added at the offsets of the preamble bits that use that pair."""
    H, s = hop_set(CTRL["rate"])
    Lw = round(Lr)
    off = {n: np.round(np.arange(n) * Lr).astype(int) for n in PILOTS}
    out = {n: np.zeros(max(0, len(y) - Lw + 1 - off[n][-1])) for n in PILOTS}
    for p in range(H):
        f0, f1 = hop_pair(CTRL["rate"], p)
        e0, e1 = sliding_energy(y, fs, f0, Lw), sliding_energy(y, fs, f1, Lw)
        d = (e1 - e0) / (e1 + e0 + 1e-12)
        for n, P in PILOTS.items():
            M = len(out[n])
            for j in np.nonzero((s * np.arange(n)) % H == p)[0]:
                out[n] += (2 * P[j] - 1) * d[off[n][j]:off[n][j] + M]
    return {n: c / n for n, c in out.items()}


def hop_postamble(y, fs, rate, k0, P, lo, hi, Ld, M=2):
    """Correlation with a hopping postamble (symbols P) whose first symbol is symbol k0, for every start lo .. hi."""
    Lw = round(Ld)
    off = np.round(np.arange(len(P)) * Ld).astype(int)
    seg = y[lo:hi + off[-1] + Lw]
    F = hop_tones(rate, k0 + np.arange(len(P)), M)
    E, c = {}, np.zeros(hi - lo + 1)
    for j in range(len(P)):
        fj = [f[j] for f in F]
        for f in fj:
            if f not in E:
                E[f] = sliding_energy(seg, fs, f, Lw)
        c += known_contrast([E[f][off[j]:off[j] + hi - lo + 1] for f in fj], P[j])
    return c / len(P)


# ---------------------------------------------------------------- OFDM (v3)
# Many tones at once, each with its own phase. An inverse FFT of N points at 48 kHz puts one complex number on each
# subcarrier k (at k 48000 / N Hz) from 3.2 to 16 kHz, and the last CP samples of each symbol are sent again just before
# it: the cyclic prefix, or guard interval. An echo that arrives within the guard then only multiplies each subcarrier
# by a complex gain (the room's frequency response) instead of spilling into the next symbol. Each subcarrier carries
# 2 bits per symbol as a phase step from one symbol to the next (DQPSK), so the receiver needs no estimate of that gain:
# it compares each symbol with the one before on the same subcarrier. A known reference symbol comes first and again
# last: the first starts the phase steps, and the pair gives the receiver its timing and the time scale.
#   long guard:  N 2048, CP 1024 (21.3 ms), 546 subcarriers 23.4 Hz apart, a symbol every 64 ms: 17.1 kbit/s before the code
#   short guard: N 1024, CP 256 (5.3 ms), 273 subcarriers 46.9 Hz apart, a symbol every 26.7 ms: 20.5 kbit/s
OFDM_MODES = [dict(name="long", N=2048, CP=1024), dict(name="short", N=1024, CP=256)]   # header rate field 0, 1
OFDM_F0, OFDM_F1 = 3200, 16000
OFDM_RMS = 0.45          # half the 0.9 peak of FSK (6 dB below it); higher peaks (4.5% of the samples) are clipped at 0.9
OFDM_GROUP = 32          # subcarriers per noise estimate
OFDM_BACK = 8            # FFT windows this many samples early: the direct sound must not fall outside the guard
OFDM_MAX_PPM = 500       # time-scale search: clock offset plus the motion of the phone (0.1 m/s is 290 ppm)


def ofdm_mode(name):
    return next(m for m in OFDM_MODES if m["name"] == name)


def ofdm_carriers(N):
    k = np.arange(N // 2 + 1)
    return k[(k * FS_TX >= OFDM_F0 * N) & (k * FS_TX < OFDM_F1 * N)]


def ofdm_rate(mode, omap="dqpsk"):
    """Channel bits per second: m bits per subcarrier (but the pilots) per symbol of N + CP samples."""
    K = len(ofdm_carriers(mode["N"]))
    return OFDM_BITS[omap] * (K if omap == "dqpsk" else K * (1 - 1 / OFDM_PSTEP)) * FS_TX / (mode["N"] + mode["CP"])


# ---------------------------------------------------------------- coherent OFDM (v4)
# DQPSK needs no knowledge of the channel, but it pays for it: comparing two noisy symbols doubles the noise, and a
# phase step carries 2 bits at most in practice. With pilots, known values on some subcarriers of every data symbol,
# the receiver estimates each subcarrier's gain and phase (the room's frequency response) and decides on the symbol
# itself: QPSK (2 bits, 2 dB less SNR than DQPSK), 8PSK (3 bits) or 16QAM (4 bits: the amplitude carries bits too),
# 64QAM (6 bits) or 256QAM (8 bits; in simulation about 5 and 11 dB more SNR than 16QAM).
#   pilots: subcarrier j of data symbol i (1 .. ns) carries the reference symbol's value if (j - 3 i) mod 12 = 0,
#   1 subcarrier in 12, shifted by 3 from one symbol to the next (8% of the subcarriers)
#   header: the CRC xor OFDM_MASKS[mapping] (and xor 0xFF for codes 4 to 6), so that older receivers ignore the frame
OFDM_MAPS = ["dqpsk", "qpsk", "8psk", "16qam", "64qam", "256qam"]
OFDM_NAMES = dict(dqpsk="DQPSK", qpsk="QPSK", **{"8psk": "8PSK", "16qam": "16QAM", "64qam": "64QAM", "256qam": "256QAM", "1024qam": "1024QAM"})
OFDM_MASKS = [0x00, 0x33, 0x55, 0x66, 0x3C, 0x0F]   # with xor 0xFF: 12 distinct masks, any two at least 4 bits apart
OFDM_BITS = {"dqpsk": 2, "qpsk": 2, "8psk": 3, "16qam": 4, "64qam": 6, "256qam": 8, "1024qam": 10}   # 1024QAM: single carrier only (v7.5)
OFDM_PSTEP, OFDM_PSHIFT = 12, 3
OFDM_DMAX = 8            # delay change of a symbol (moving phone) searched within +-8 samples
OFDM_ITERS = 3           # rounds of channel estimate <-> per-symbol gain, phase and delay
OFDM_ALPHA = 0.01        # channel smoothing A (A + alpha)^-1, see ofdm_smoother
OFDM_SLIP = 16           # v4.1: a symbol whose timing jumps by more than CP / 16 from the symbols before it: samples lost or repeated
OFDM_GAPS = 4            # v4.1: at most 4 candidate gaps per frame are weighed against the references (more: none kept)
_SMOOTH = {}


def ofdm_smoother(mode):
    """The matrix that smooths a channel estimate across the band: an echo that arrives within the guard (-16 .. CP
    samples after the FFT window starts) makes a response that varies slowly across the subcarriers. A: to the delay
    domain (inverse FFT), keep -16 .. CP, back (FFT); A[j, l] depends on j - l only. A alone would bend the estimate
    near the band's edges; A (A + alpha I)^-1 keeps what the guard holds well (eigenvalues of A near 1) and drops what
    it does not (near 0): about half the noise goes, with the long guard."""
    if mode["name"] not in _SMOOTH:
        N, CP = mode["N"], mode["CP"]
        K = len(ofdm_carriers(N))
        keep = np.zeros(N)
        keep[np.arange(-16, CP + 1) % N] = 1
        a = np.fft.fft(keep)[:K] / N
        d = np.subtract.outer(np.arange(K), np.arange(K))
        A = np.where(d >= 0, a[np.abs(d)], np.conj(a[np.abs(d)]))
        _SMOOTH[mode["name"]] = A @ np.linalg.inv(A + OFDM_ALPHA * np.eye(K))
    return _SMOOTH[mode["name"]]


def ofdm_points(omap):
    """A coherent mapping's points (mean energy 1) and the bits of each (MSB first, Gray)."""
    if omap == "qpsk":
        lab = np.arange(4)
        return ((1 - 2 * (lab >> 1)) + 1j * (1 - 2 * (lab & 1))) / math.sqrt(2), lab
    if omap == "8psk":
        p = np.arange(8)
        return np.exp(2j * np.pi * p / 8), p ^ (p >> 1)
    h = OFDM_BITS[omap] // 2                       # square QAM: h Gray bits for the real part, h for the imaginary
    L = 1 << h
    A, g = np.arange(-(L - 1), L, 2), np.arange(L) ^ (np.arange(L) >> 1)
    u, v = np.repeat(np.arange(L), L), np.tile(np.arange(L), L)
    return (A[u] + 1j * A[v]) / math.sqrt(2 * (L * L - 1) / 3), g[u] << h | g[v]


def ofdm_pilots(K, i):
    """True on the pilot subcarriers of data symbol i (1 .. ns)."""
    return (np.arange(K) - OFDM_PSHIFT * i) % OFDM_PSTEP == 0


def ofdm_layout(mode, omap, nc):
    """Data symbols for nc channel bits, and the channel bits they hold."""
    K, m = len(ofdm_carriers(mode["N"])), OFDM_BITS[omap]
    if omap == "dqpsk":
        ns = -(-nc // (2 * K))
        return ns, 2 * K * ns
    ns, cap = 0, 0
    while cap < nc:
        ns += 1
        cap += m * int(np.sum(~ofdm_pilots(K, ns)))
    return ns, cap


def ofdm_reference(K):
    """Phases of the reference symbol, pi j^2 / K (Newman): its peaks stay low, 5.4 dB above its rms, so that it is never clipped."""
    j = np.arange(K)
    return np.pi * ((j * j) % (2 * K)) / K


def ofdm_symbols(bits, mode, omap="dqpsk"):
    """Channel bits (coded, scrambled) -> the OFDM part at 48 kHz: the reference symbol, ceil(len / 2K) data symbols
    and the reference symbol again, each preceded by its cyclic prefix. Bits b0 b1 of a subcarrier -> phase step
    pi/4 + m pi/2, m = (b0 xor b1) + 2 b1 (Gray: 00, 10, 11, 01), so that the sign of the real part of the step gives
    b0 and that of its imaginary part b1. Coherent mappings (v4): pilots and points, see ofdm_points and ofdm_pilots."""
    N, CP = mode["N"], mode["CP"]
    kk = ofdm_carriers(N)
    K = len(kk)
    a = OFDM_RMS * math.sqrt(2 / K)
    if omap != "dqpsk":
        pts, lab = ofdm_points(omap)
        m, inv = OFDM_BITS[omap], np.argsort(lab)              # inv[bits] = the point
        ns, cap = ofdm_layout(mode, omap, len(bits))
        b = np.zeros(cap, dtype=int)
        b[:len(bits)] = bits
        ref = np.exp(1j * ofdm_reference(K))
        X = np.zeros((ns + 2, N // 2 + 1), dtype=complex)
        X[0, kk] = X[ns + 1, kk] = ref
        q = 0
        for i in range(1, ns + 1):
            row, dsc = ref.copy(), np.nonzero(~ofdm_pilots(K, i))[0]
            v = np.zeros(len(dsc), dtype=int)
            for k in range(m):
                v = 2 * v + b[q + k:q + m * len(dsc):m]
            row[dsc] = pts[inv[v]]
            X[i, kk] = row
            q += m * len(dsc)
        x = np.fft.irfft(X, N, axis=1) * (N / 2 * a)
        return np.clip(np.hstack([x[:, N - CP:], x]).ravel(), -0.9, 0.9)
    ns = -(-len(bits) // (2 * K))
    b = np.zeros(ns * 2 * K, dtype=int)
    b[:len(bits)] = bits
    b0, b1 = b[0::2].reshape(ns, K), b[1::2].reshape(ns, K)
    step = 1 + 2 * ((b0 ^ b1) + 2 * b1)                     # in units of pi/4
    P = np.vstack([np.zeros(K, dtype=int), np.cumsum(step, axis=0) % 8, np.zeros(K, dtype=int)])
    ph = ofdm_reference(K)[None, :] + P * np.pi / 4
    X = np.zeros((ns + 2, N // 2 + 1), dtype=complex)
    X[:, kk] = np.cos(ph) + 1j * np.sin(ph)
    # each subcarrier a cosine of amplitude a: the rms of a symbol is a sqrt(K / 2), set to OFDM_RMS
    x = np.fft.irfft(X, N, axis=1) * (N / 2 * a)
    return np.clip(np.hstack([x[:, N - CP:], x]).ravel(), -0.9, 0.9)


def frac_read(y, t, A=16):
    """y at the fractional sample positions t: band-limited interpolation, a sinc with a Hann window, A samples each side."""
    t = np.asarray(t, dtype=float)
    k0 = np.floor(t).astype(int)
    out = np.zeros(len(t))
    for k in range(-A + 1, A + 1):
        i = k0 + k
        u = t - i
        w = np.sinc(u) * (0.5 + 0.5 * np.cos(np.pi * u / A))
        ok = (i >= 0) & (i < len(y))
        out[ok] += w[ok] * y[i[ok]]
    return out


def ofdm_receive(y, fs, t0, mode, ns, nc, omap="dqpsk"):
    """The OFDM part of a frame whose first reference symbol (its guard first) starts near sample t0 of the recording y.
    1. Read it on a 48 kHz grid, whatever the microphone's rate.
    2. Timing: a reference symbol gives the room's impulse response as seen in the band (the inverse FFT of the
       received subcarriers divided by the sent ones); the FFT windows are placed so that the guard covers the stretch
       of that response with the most energy, from the direct sound on.
    3. Time scale: a clock offset, or a moving phone, stretches the recording by eps. Between the two reference
       symbols, (ns + 1) (N + CP) samples apart, the phase on subcarrier f then turns by 2 pi f (ns + 1) (N + CP) eps /
       48000: eps is the slope that best lines up the subcarriers. The part is then read again on a grid stretched by
       1 + eps.
    4. Differential detection: z = Y_i conj(Y_i-1) on each subcarrier, LLR(b0) = -sqrt(2) Re z / N0 and LLR(b1) =
       -sqrt(2) Im z / N0 (> 0 means 1), N0 measured from the decisions per group of 32 subcarriers: the noise and the
       echo that arrives after the guard. v4: a phone that moves changes its delay d from one symbol to the next,
       which turns subcarrier f by 2 pi f d / 48000 (90 degrees at 16 kHz for 0.75 sample, a 2 cm swing); z^4 has
       no data left (every step is pi/4 + m pi/2), so d is the delay that best lines up -z^4 across the band, and z
       is turned back by it before the decisions.
    4'. Coherent detection (v4): the channel of symbol i is H(f) c_i exp(-j 2 pi f d_i / 48000), the room's response
       H times a gain and phase c_i and a delay change d_i (a moving phone), fitted on the symbol's pilots (weighted
       least squares, d_i searched within +-8 samples). H is the average of the two reference symbols and of all the
       pilots once their c_i and d_i are removed, smoothed across the band (ofdm_smoother: keep what an echo within
       the guard can make). Three rounds. N0 per group of 32 subcarriers from what the
       references and pilots leave (x 1.25: the estimate absorbs some of the noise). LLR of each bit:
       ln sum over the points s with that bit 1 of exp(-|Y - H_i s|^2 / N0), minus the same over bit 0.
    2'. Samples lost or repeated (v4.1): the audio path of a phone may drop a block of samples (128, one block of the
       browser's audio) in the middle of a frame. Every symbol after the gap then arrives that much early: its
       subcarriers turn by 2 pi f 128 / 48000, far beyond what the delay tracking follows, and the reference symbols
       no longer give the time scale. Each symbol's coarse timing is the shift (within +-N/4) that best lines its guard
       up with its end; a jump of more than CP / OFDM_SLIP that the next symbol confirms is a gap. Its size is refined
       on the symbols on either side of it (pilots or the reference: their delay against the first reference; DQPSK:
       -z^4 across the gap), every symbol after it is read that much later (or earlier), and the symbol that holds
       the gap (its guard no longer matches its end) is erased: its LLRs are 0, the code fills them in. A faint frame
       can make the guards of a few symbols look shifted, so of the gaps found, ofdm_receive keeps the set with which
       the two reference symbols agree best at their best time scale (the fewest gaps when several sets agree as well)."""
    N, CP = mode["N"], mode["CP"]
    kk = ofdm_carriers(N)
    K, L = len(kk), N + CP
    ref = np.exp(1j * ofdm_reference(K))
    sc = fs / FS_TX
    pad = N // 2
    n = np.arange(-pad, L + pad)                              # each symbol read on its own: grid position i L + m is r[i, m + pad]

    def read(eps, S=None):                                    # S: each symbol's shift, the samples lost (< 0) or repeated before it
        S = np.zeros(ns + 2) if S is None else S
        return np.array([frac_read(y, t0 + sc * (S[i] + (i * L + n) * (1 + eps * 1e-6))) for i in range(ns + 2)])

    def spec(r, i, start):                                    # the subcarriers of symbol i, window from i L + start
        return np.fft.rfft(r[i, pad + start:pad + start + N])[kk]

    b0 = CP // 2                                              # a reference window in the middle of its guard
    taper = np.hanning(K + 2)[1:-1]                           # across the band: lower side lobes in time

    def response(r, refs, start):                             # |h|^2 at (CP - start + delay) mod N, from the reference symbols
        e = np.zeros(N)
        for i in refs:
            G = np.zeros(N, dtype=complex)
            G[kk] = spec(r, i, start) * np.conj(ref) * taper
            e += np.abs(np.fft.ifft(G)) ** 2
        return e

    def timing(e):                                            # the latest stretch of CP + 1 samples with 99% of the most energy
        c = np.concatenate([[0], np.cumsum(np.concatenate([e, e]))])
        s = np.arange(b0 - N // 4, b0 + N // 4 + 1)
        E = c[(s % N) + CP + 1] - c[s % N]
        best = int(np.nonzero(E >= 0.99 * E.max())[0][-1])
        return int(s[best] - b0), float(E[best] / max(e.sum(), 1e-30))

    r = read(0.0)
    tau, _ = timing(response(r, [0], CP - b0))
    f = kk * FS_TX / N
    cand = gaps(r, tau, f, ref, mode, omap, ns, N, CP, pad, spec, read)   # v4.1: samples lost or repeated, candidates
    S, slips, erased = np.zeros(ns + 2), [], set()
    if 0 < len(cand) <= OFDM_GAPS:                            # keep the gaps with which the two references agree best
        R = np.exp(1j * 2 * np.pi * (ns + 1) * L / FS_TX * 1e-6 * np.outer(np.arange(-OFDM_MAX_PPM, OFDM_MAX_PPM + 1, 1.0), f))
        Y0 = spec(r, 0, CP // 2 + tau)

        def agree(tot):                                       # 1 if the last reference matches the first one exactly
            x = frac_read(y, t0 + sc * (tot + (ns + 1) * L + n))
            z = spec(x[None], 0, CP // 2 + tau) * np.conj(Y0)
            return float(np.abs(R @ z).max() / max(np.abs(z).sum(), 1e-30))
        tries = []
        for mk in range(1 << len(cand)):
            sub = [c for k, c in enumerate(cand) if mk >> k & 1]
            tries.append((agree(sum(c[1] for c in sub)), len(sub), sub))
        top = max(t[0] for t in tries)
        _, _, sub = min((t for t in tries if t[0] >= 0.9 * top), key=lambda t: (t[1], -t[0]))
        for i, d, e in sub:
            S[i:] += d
            slips.append((i, d))
            erased |= set(e)
    erased = sorted(erased)
    if slips:
        r = read(0.0, S)
    # time scale, from the two reference symbols (windows in the middle of the guard: a drift of up to CP/2 fits)
    z = spec(r, ns + 1, CP // 2 + tau) * np.conj(spec(r, 0, CP // 2 + tau))
    grid = np.arange(-OFDM_MAX_PPM, OFDM_MAX_PPM + 1, 1.0)
    score = np.abs(np.exp(1j * 2 * np.pi * (ns + 1) * L / FS_TX * 1e-6 * np.outer(grid, f)) @ z)
    i = int(np.argmax(score))
    eps = float(grid[i])
    if 0 < i < len(grid) - 1:                                 # parabola through the peak and its neighbours
        a0, a1, a2 = score[i - 1], score[i], score[i + 1]
        eps += 0.5 * (a0 - a2) / (a0 - 2 * a1 + a2 - 1e-30)
    if abs(eps) > 0.05 or slips:
        r = read(eps, S)
    tau, _ = timing(response(r, [0, ns + 1], CP - b0))        # both references, now aligned
    w = CP + tau - OFDM_BACK                                  # start of the FFT windows
    # the response as these windows see it: what arrives within the guard (index 0 .. CP) only scales each
    # subcarrier, the rest (echo after the guard, noise) spreads over all the indices
    h = response(r, [0, ns + 1], w)
    inside = float(h[:CP + 1].sum() / max(h.sum(), 1e-30))
    Y = np.array([spec(r, i, w) for i in range(ns + 2)])     # the two references and the data
    gap = dict(slips=[(i, d * sc) for i, d in slips], erased=erased)   # the gaps in samples of the recording
    if omap != "dqpsk":
        return dict(ofdm_coherent(Y, ns, nc, mode, omap, erased), eps=eps, tau=tau, inside=inside, K=K, h=h, **gap)
    z = Y[1:ns + 1] * np.conj(Y[:ns])
    dgrid = np.arange(-4, 4.001, 0.01)                        # delay change between two symbols (moving phone)
    E = np.exp(8j * np.pi * np.outer(dgrid, f) / FS_TX)
    for i in range(ns):
        w4 = -z[i] ** 4 / np.maximum(np.abs(z[i]) ** 2, 1e-30)
        z[i] *= np.exp(2j * np.pi * f * dgrid[int(np.argmax((E @ w4).real))] / FS_TX)
    d = (np.sign(z.real) + 1j * np.sign(z.imag)) / math.sqrt(2)
    q = np.imag(z * np.conj(d)) ** 2                          # noise across the decision
    pw = 0.5 * (np.abs(Y[1:ns + 1]) ** 2 + np.abs(Y[:ns]) ** 2)
    ok = np.array([i not in erased and i - 1 not in erased for i in range(1, ns + 1)])   # v4.1: z of an erased symbol is 0
    if not ok.any():
        ok[:] = True
    N0, P, snr = np.empty(K), np.empty(K), []
    for g in range(0, K, OFDM_GROUP):
        sl = slice(g, g + OFDM_GROUP)
        P[sl] = max(pw[ok, sl].mean(), 1e-30)
        N0[sl] = max(1e-30, q[ok, sl].mean() / P[g])
        snr.append(10 * math.log10(max(1e-3, P[g] / N0[g] - 1)))
    llr = np.empty((ns, 2 * K))
    llr[:, 0::2] = -math.sqrt(2) * z.real / N0
    llr[:, 1::2] = -math.sqrt(2) * z.imag / N0
    llr[~ok] = 0.0
    return dict(llr=llr.ravel()[:nc], eps=eps, tau=tau, inside=inside, snr=np.array(snr), K=K, h=h,
                z=(z / np.maximum(pw.mean(axis=0), 1e-30)).ravel(),   # z over each subcarrier's power, to plot
                slips=[(i, d * sc) for i, d in slips], erased=erased)


def gaps(r, tau, f, ref, mode, omap, ns, N, CP, pad, spec, read):
    """Step 2' of ofdm_receive (v4.1): samples lost or repeated in the recording. r holds the symbols read on the
    nominal grid. Returns the candidate gaps [(first symbol after it, size in samples on the 48 kHz grid, symbols
    to erase)]; ofdm_receive keeps those with which the two references agree best."""
    ks = np.arange(-N // 4, N // 4 + 1)                       # gaps of up to N / 4 samples (512 with the long guard, 256 with the short one)
    a, o, q = pad + ks, [], []
    for x in r:                                               # coarse timing: the shift that best lines the guard up with the end
        c = [np.concatenate([[0.0], np.cumsum(u)]) for u in (x[:-N] * x[N:], x[:-N] ** 2, x[N:] ** 2)]
        v = (c[0][a + CP] - c[0][a]) / np.sqrt(np.maximum((c[1][a + CP] - c[1][a]) * (c[2][a + CP] - c[2][a]), 1e-30))
        o.append(float(ks[np.argmax(v)]))
        q.append(float(v.max()))
    o, q = np.array(o), np.array(q)
    med = float(np.median(q))
    good = q >= 0.6 * med                                     # a symbol that holds a gap: its guard no longer matches its end
    clean = q >= 0.9 * med                                    # its guard matches its end as well as most symbols' do
    jumps, seg, last = [], [o[0]], 0                          # last: the last clean symbol that agrees with the ones before it
    for i in range(1, ns + 2):
        if not good[i]:
            continue
        nxt = next((j for j in range(i + 1, ns + 2) if good[j]), None)
        m = float(np.median(seg))
        off = abs(o[i] - m) > CP / OFDM_SLIP
        if off and (nxt is None or abs(o[nxt] - o[i]) <= CP / OFDM_SLIP):   # a jump that the next symbol confirms
            a = i if nxt is None else nxt                     # sure to be after the gap: i may still hold its end
            jumps.append((i, o[a] - m, last, a))
            seg, last = [o[i]], i if clean[i] else a
        elif not off:
            seg.append(o[i])
            if clean[i]:
                last = i
    S = np.zeros(ns + 2)
    if not jumps:
        return []
    # erased: the (at most two) symbols between the last clean one before the gap and the first one after it (one
    # of them holds the gap, the guard of the next one may reach back across it), and the first one after it unless
    # its guard is clean
    ers = [[k for k in range(max(b + 1, i - 2), i + 1) if 0 < k <= ns and (k < i or not clean[k])] for i, _, b, a in jumps]
    for i, d, _, _ in jumps:
        S[i:] += d
    rr, st = read(0.0, S), CP // 2                            # windows in the middle of the guard (tau from one reference may be off)
    W = 0.4 * N / OFDM_PSTEP                                  # within half the pilots' period in time (N / 12 samples)
    dgrid = np.arange(-W, W + 0.001, 0.05)
    E = np.exp(2j * np.pi * np.outer(dgrid, f) / FS_TX)
    K = len(f)
    if omap != "dqpsk":
        H0 = ofdm_smoother(mode) @ (spec(rr, 0, st) / ref)

        def delay(k):                                         # of symbol k against the first reference: its pilots, or all of a reference
            p = np.arange(K) if k in (0, ns + 1) else np.nonzero(ofdm_pilots(K, k))[0]
            return dgrid[int(np.argmax(np.abs(E[:, p] @ (spec(rr, k, st)[p] * np.conj(ref[p] * H0[p])))))]
    out = []
    for (i, d, b, a), e in zip(jumps, ers):                   # refined on b, before the gap, and a, after it
        if omap != "dqpsk":
            fine = delay(a) - delay(b)
        elif a <= ns:                                         # DQPSK: z^4 of the two symbols has no data left (sign (-1)^(a - b))
            z = spec(rr, a, st) * np.conj(spec(rr, b, st))
            w4 = (-1) ** (a - b) * z ** 4 / np.maximum(np.abs(z) ** 2, 1e-30)
            fine = dgrid[int(np.argmax((E ** 4 @ w4).real))]
        else:                                                 # DQPSK, a gap just before the last reference: the coarse size
            fine = 0.0
        S[i:] += fine
        out.append((i, d + fine, e))
    return out


def ofdm_coherent(Y, ns, nc, mode, omap, erased=()):
    """Step 4' of ofdm_receive: Y holds the subcarriers of the reference, the ns data symbols and the reference.
    erased (v4.1): data symbols that hold a gap; they are left out and their LLRs are 0."""
    N, CP = mode["N"], mode["CP"]
    kk = ofdm_carriers(N)
    K, f = len(kk), kk * FS_TX / N
    ref = np.exp(1j * ofdm_reference(K))
    pts, lab = ofdm_points(omap)
    m = OFDM_BITS[omap]
    W = ofdm_smoother(mode)
    smooth = lambda G: W @ G
    pil = [np.arange(K)] + [np.nonzero(ofdm_pilots(K, i))[0] for i in range(1, ns + 1)] + [np.arange(K)]
    for i in erased:
        pil[i] = np.zeros(0, dtype=int)
    dgrid = np.arange(-OFDM_DMAX, OFDM_DMAX + 0.001, 0.02)

    def fit(i, H):                                            # c_i exp(-j 2 pi f d_i / 48000) on the pilots of symbol i
        p = pil[i]
        if len(p) == 0:
            return np.ones(K, dtype=complex)
        e = ref[p] * H[p]
        q = Y[i, p] * np.conj(e)
        dl = dgrid[int(np.argmax(np.abs(np.exp(2j * np.pi * np.outer(dgrid, f[p]) / FS_TX) @ q)))]
        c = np.sum(q * np.exp(2j * np.pi * f[p] * dl / FS_TX)) / max(np.sum(np.abs(e) ** 2), 1e-30)
        return c * np.exp(-2j * np.pi * f * dl / FS_TX)

    H = smooth(Y[0] / ref)
    for _ in range(OFDM_ITERS):
        C = [np.ones(K, dtype=complex)] + [fit(i, H) for i in range(1, ns + 2)]
        acc, cnt = np.zeros(K, dtype=complex), np.zeros(K)
        for i in range(ns + 2):
            p = pil[i]
            acc[p] += Y[i, p] / (ref[p] * C[i][p])
            cnt[p] += 1
        H = smooth(acc / np.maximum(cnt, 1))
    num, den = np.zeros(K), np.zeros(K)
    for i in range(ns + 2):
        p = pil[i]
        num[p] += np.abs(Y[i, p] - ref[p] * H[p] * C[i][p]) ** 2
        den[p] += 1
    N0, snr = np.empty(K), []
    for g in range(0, K, OFDM_GROUP):
        sl = slice(g, g + OFDM_GROUP)
        N0[sl] = max(1e-30, 1.25 * num[sl].sum() / max(den[sl].sum(), 1))
        snr.append(10 * math.log10(max(1e-3, np.mean(np.abs(H[sl]) ** 2) / N0[g])))
    llr, zs = [], []
    for i in range(1, ns + 1):
        dsc = np.nonzero(~ofdm_pilots(K, i))[0]
        Hd = H[dsc] * C[i][dsc]
        d2 = -np.abs(Y[i, dsc][:, None] - Hd[:, None] * pts[None, :]) ** 2 / N0[dsc][:, None]
        out = np.empty((len(dsc), m))
        for k in range(m):
            bit = (lab >> (m - 1 - k)) & 1
            out[:, k] = np.logaddexp.reduce(d2[:, bit == 1], axis=1) - np.logaddexp.reduce(d2[:, bit == 0], axis=1)
        if i in erased:
            out[:] = 0.0
        llr.append(out.ravel())
        zs.append(Y[i, dsc] / Hd)
    return dict(llr=np.concatenate(llr)[:nc], snr=np.array(snr), z=np.concatenate(zs))


def ofdm_desc(mode, omap="dqpsk"):
    N, CP = mode["N"], mode["CP"]
    K = len(ofdm_carriers(N))
    how = OFDM_NAMES[omap] + ("" if omap == "dqpsk" else f" with pilots (1 subcarrier in {OFDM_PSTEP})")
    return (f"OFDM, {K} subcarriers from {OFDM_F0 / 1000:g} to {OFDM_F1 / 1000:g} kHz {FS_TX / N:.1f} Hz apart, {how}, "
            f"symbol {1000 * N / FS_TX:.1f} ms + guard {1000 * CP / FS_TX:.1f} ms, {ofdm_rate(mode, omap):.0f} bit/s before the code")


# ---------------------------------------------------------------- single carrier (v5)
# One carrier at 9.6 kHz, 9600 symbols per second, each shaped by a root-raised-cosine pulse of roll-off 1/3: the
# signal fills 3.2 to 16 kHz as OFDM does, but every symbol uses the whole band, and the room's echoes make each
# received sample a mix of the symbols sent over the last tens of ms (intersymbol interference). The receiver
# writes this mix as a convolution at the symbol rate, y = c * s + noise, with c the room seen through the two
# pulse filters (512 taps, 53 ms, of which 64 before the first strong arrival), learns c from known training symbols,
# and then undoes it:
#   detector: soft SMAP (Ramos, Senecal and Yang): the most likely symbols under a Gaussian-mixture prior around the
#     constellation's levels, found by gradient steps (FFTs only), with a soft output per bit; or a linear MMSE
#     equalizer (with soft interference cancellation once the decoder speaks)
#   turbo equalization: the detector and the turbo decoder exchange beliefs on the coded bits, a few passes.
# Part (after the control part): 0.1 s of silence, then
#   TA (1280 QPSK symbols) TC (1280) | data | TC | data | TC ... | data | TC
# with at most 4800 data symbols (0.5 s) between training blocks: the room changes as people move (the channel
# drifts by about -24 dB in 0.25 s and -18 dB in 1 s in our recordings), and each block gives the receiver the
# room again, blended linearly over the data between two blocks. Every block ends with the same 1280 symbols TC,
# so their arrival times give the clock (and a moving phone) along the whole part.
# Channel bits: the code's output (turbo codes only, with their rectangular interleaver), padded with 0 to whole
# symbols, then scrambled; points as OFDM's (ofdm_points). A random bit interleaver on top, which helps OFDM's
# 256QAM by 0.5 to 2 dB (one pass, no feedback), made turbo equalization worse here: in simulation 0 of 4 frames at
# the level where 4 of 4 decode without it (256QAM, turbo 1/2), so there is none.
SC_RS, SC_SPS, SC_FC, SC_BETA, SC_SPAN = 9600, 5, 9600, 1 / 3, 12   # baud, samples per symbol, carrier, roll-off, pulse half length (symbols)
SC_FS = FS_TX            # v7.4: the single carrier's own sampling rate (its pulses, carrier and time map)
SC_BANDS = {False: (9600, 5, 9600, 1), True: (96000 / 7, 7, 11400, 2)}   # v7.4: (baud, samples per symbol, carrier,
                         # SC_FS / FS_TX) of the usual band (3.2 to 16 kHz) and of the wide one (2.3 to 20.5 kHz), the
                         # latter on a 96 kHz grid whose every other sample is sent: room soundings in Impasse Carnot heard
                         # 16 to 20.5 kHz as clearly as the usual band (the phone's speaker 6 to 10 dB weaker there, the
                         # noise as much lower), its echo die out fastest, and 2.3 to 3.2 kHz 3 to 5 dB less clearly; the
                         # iPhone's microphone stops at 21 kHz
SC_RMS = 0.45            # rms at the speaker, as OFDM's; the peaks above 0.9 are clipped (about 1 % of the samples
                         # at 256QAM, the distortion about 35 dB below the symbols)
SC_NA = 1280             # TA: the first block's own training, before its TC
SC_NC = 1280             # TC, the training block, by default (medium)
SC_TRAIN = {"short": 640, "medium": 1280, "long": 2560}   # TC: the sender's choice; like OFDM's guard, it sets the
SC_TCODE = {1280: 0, 640: 1, 2560: 2}                     # longest room the receiver can learn, for some rate (header code)
SC_TAPS = {640: (256,), 1280: (512,), 2560: (512, 1024)}  # rooms fitted (taps, 64 before the first strong arrival): each
                         # later TC gives TC - taps + 1 rows that see it alone, at least 1.5 taps' worth; with the long
                         # TC the receiver keeps 1024 taps when they leave less (their own error counted)
SC_D = 4800              # at most this many data symbols between two training blocks
SC_GAP = 0.1             # silence after the control part: its echoes fade before the training
SC_L, SC_PRE = 512, 64   # channel taps, of which SC_PRE before the first strong arrival
SC_FIELD = 2             # header: tone set 3, rate field 2 (OFDM uses 0 and 1, 6400 bit/s 3)
SC_MAPS = ["qpsk", "16qam", "64qam", "256qam", "1024qam"]
SC_WIDTH = {"qpsk": 0.25, "16qam": 0.25, "64qam": 0.4, "256qam": 1.0, "1024qam": 1.0}   # SMAP: prior width tau^2 / d^2 (d: spacing of the levels)
SC_SLOT = {"1024qam": "8psk"}   # v7.5: the header mask a single-carrier mapping borrows (8PSK is never sent on the single
                         # carrier, and every older receiver rejects it there, so they ignore 1024QAM frames)
SC_PASSES = 8            # turbo equalization: at most this many detector <-> decoder passes (full effort)
SC_EFFORT = {"fast": ((30, 20), 4), "full": ((150, 100), SC_PASSES)}   # the receiver's choice: SMAP's steps (first pass,
                         # later ones) and the passes; fast costs at most 1 dB at the threshold on our recordings, and
                         # takes a third (a frame that fails) to two thirds of the time
SC_RELEARN_TAPS = 2048   # v7.2, the second try: when the passes end without a codeword, one room of this many taps
                         # (213 ms; the training blocks allow 1024 at most) is learned from the whole frame, its data
                         # as the decoder's soft symbols; in a living room at 2 to 4 m the echo after 107 ms was only
                         # 15 to 19 dB down, so the training's room capped the SNR there (14 -> 29 of 34 frames)
SC_RELEARN_PASSES = 4    # v7.2: at most this many more passes, the room learned again after each (0: no second try)
SC_REL = 1.0             # ridge prior of a training block's room: the one before, give or take its own size (a weak
                         # pull; a tight one, -20 dB, cost 2 to 7 dB on our recordings when people moved)
SC_MAX_PPM = 500         # clock offset plus a moving phone
SC_EST = ("ls", "mmse")  # the receiver's estimate of the room on the first block: least squares (default) or MMSE with
                         # a prior power per tap learned from a least-squares pilot; on Sheng's recordings with noise
                         # added, MMSE decoded at the same SNR as least squares or 1 dB lower (v5/taps/tap_estimation.md)
SC_SMOOTH = 33           # MMSE: taps over which the pilot's power is averaged into the prior
SC_GAP_MAX = 1024        # v5.1: the time map also looks this many samples farther, for samples the microphone lost or
SC_GAP_COH = 1.5         # repeated between two training blocks; such a jump is a gap if the blocks agree this much
                         # better there than in the usual window (coherence of their cross-spectrum)
SC_GAP_AR = 32           # order of the linear predictor that finds where the samples went missing (sc_click)
SC_GAP_W = 16            # sc_click: samples on each side of the join over which the prediction errors are summed
SC_GAP_SPAN = SC_SPAN * SC_SPS + 8   # samples on each side of a jump that the matched filter and the interpolation smear it over
SC_GAP_EDGE = 600        # a jump inside the first block is looked for only this far from its end (closer, the next block
                         # sees it) and past its first L symbols (before, the timing takes the part after it)
SC_GAP_MIN = 4           # inside the first block: a jump of more samples than this


def sc_rrc(beta=SC_BETA, sps=SC_SPS, span=SC_SPAN):
    t = np.arange(-span * sps, span * sps + 1) / sps
    h = np.zeros_like(t)
    for i, tt in enumerate(t):
        if abs(tt) < 1e-9:
            h[i] = 1 - beta + 4 * beta / np.pi
        elif abs(abs(4 * beta * tt) - 1) < 1e-9:
            h[i] = beta / np.sqrt(2) * ((1 + 2 / np.pi) * np.sin(np.pi / (4 * beta)) + (1 - 2 / np.pi) * np.cos(np.pi / (4 * beta)))
        else:
            h[i] = (np.sin(np.pi * tt * (1 - beta)) + 4 * beta * tt * np.cos(np.pi * tt * (1 + beta))) / (np.pi * tt * (1 - (4 * beta * tt) ** 2))
    return h / np.sqrt((h ** 2).sum())


SC_G = sc_rrc()
SC_DLY = SC_SPAN * SC_SPS          # 60 samples: the pulse's peak after its start
SC_GS = {w: sc_rrc(sps=v[1]) for w, v in SC_BANDS.items()}
# v7.4: TA, the training blocks, the rooms fitted, the data between blocks and the second try's room, in symbols, then
# the gap search's reach and margin in samples of SC_FS, for each band; the wide band's last about as long as the usual
# ones (10/7 the symbols, rooms in multiples of 64 taps): with 3.2 to 19.2 kHz in Impasse Carnot's echo (rooms measured
# by sounding), the same symbols lost 1.7 to 3.3 dB of SNR, the same durations gained 0.9 to 2.9 dB over the usual band
SC_SIZES = {False: (SC_NA, SC_NC, SC_TRAIN, SC_TCODE, SC_TAPS, SC_D, SC_RELEARN_TAPS, SC_GAP_MAX, SC_GAP_EDGE),
            True: (1792, 1792, {"short": 896, "medium": 1792, "long": 3584}, {1792: 0, 896: 1, 3584: 2},
                   {896: (352,), 1792: (704,), 3584: (704, 1408)}, 6720, 2944, 2048, 1200)}


def sc_set_band(wide):
    """v7.4: the band of the single-carrier frame about to be built or read (wide: 13 714 symbols/s at 11.4 kHz, on
    a 96 kHz grid), and its sizes (SC_SIZES)."""
    global SC_RS, SC_SPS, SC_FC, SC_FS, SC_G, SC_DLY, SC_GAP_SPAN, SC_NA, SC_NC, SC_TRAIN, SC_TCODE, SC_TAPS, SC_D, \
        SC_RELEARN_TAPS, SC_GAP_MAX, SC_GAP_EDGE
    wide = bool(wide)
    SC_RS, SC_SPS, SC_FC, u = SC_BANDS[wide]
    SC_FS, SC_G, SC_DLY, SC_GAP_SPAN = FS_TX * u, SC_GS[wide], SC_SPAN * SC_SPS, SC_SPAN * SC_SPS + 8
    SC_NA, SC_NC, SC_TRAIN, SC_TCODE, SC_TAPS, SC_D, SC_RELEARN_TAPS, SC_GAP_MAX, SC_GAP_EDGE = SC_SIZES[wide]


def xorshift32(seed, n):
    M, x, out = 0xFFFFFFFF, (seed & 0xFFFFFFFF) or 1, []
    for _ in range(n):
        x = (x ^ (x << 13)) & M
        x ^= x >> 17
        x = (x ^ (x << 5)) & M
        out.append(x)
    return np.array(out, dtype=np.int64)


def sc_training(nt=SC_NC):
    """TA and TC (nt symbols): QPSK symbols exp(j pi/2 (q + 1/2)), q = the top two bits of xorshift32 seeded 0x5C0A
    and 0x5C0C (a shorter TC is the start of a longer one)."""
    return [np.exp(1j * np.pi / 2 * ((xorshift32(s, n) >> 30) + 0.5)) for s, n in ((0x5C0A, SC_NA), (0x5C0C, nt))]


def sc_axis(omap):
    """A square mapping axis by axis: its levels (increasing), the h bits of each level (MSB first), h."""
    pts, lab = ofdm_points(omap)
    h = OFDM_BITS[omap] // 2
    amp = np.unique(np.round(pts.real, 12))
    g = np.array([lab[np.argmin(np.abs(pts.real - a))] >> h for a in amp])
    return amp, g, h


def sc_layout(ns, nt=SC_NC):
    """Training blocks and data segments as (first, end) symbol indices, and the symbols in all, for ns data symbols
    and a TC of nt symbols."""
    nseg = max(1, -(-ns // SC_D))
    lens = [ns // nseg + (i < ns % nseg) for i in range(nseg)]
    blocks, segs, t = [(0, SC_NA + nt)], [], SC_NA + nt
    for n in lens:
        segs.append((t, t + n))
        t += n
        blocks.append((t, t + nt))
        t += nt
    return blocks, segs, t


def sc_rate(omap, nt=SC_NC, wide=False):
    """Channel bits per second over the data and the training blocks between them."""
    d = SC_SIZES[wide][5]
    return OFDM_BITS[omap] * SC_BANDS[wide][0] * d / (d + nt)


def sc_seconds(ns, nt=SC_NC):
    return SC_GAP + (sc_layout(ns, nt)[2] * SC_SPS + len(SC_G) - 1) / SC_FS


def sc_symbols(bits, omap, nt=SC_NC):
    """Channel bits (padded, scrambled) -> every symbol of the part: training blocks and data."""
    pts, lab = ofdm_points(omap)
    m, inv = OFDM_BITS[omap], np.argsort(lab)
    ns = len(bits) // m
    s = pts[inv[np.asarray(bits).reshape(ns, m) @ (1 << np.arange(m - 1, -1, -1))]]
    ta, tc = sc_training(nt)
    blocks, segs, n = sc_layout(ns, nt)
    x = np.zeros(n, dtype=complex)
    x[:SC_NA], x[SC_NA:SC_NA + nt] = ta, tc
    q = 0
    for (a, b), (c, d) in zip(segs, blocks[1:]):
        x[a:b], x[c:d] = s[q:q + b - a], tc
        q += b - a
    return x


def sc_passband(x):
    u = np.zeros(len(x) * SC_SPS, dtype=complex)
    u[::SC_SPS] = x
    b = np.convolve(u, SC_G)
    k = np.arange(len(b))
    return np.sqrt(2) * np.real(b * np.exp(2j * np.pi * SC_FC / SC_FS * k)) * SC_RMS * math.sqrt(SC_SPS)


def sc_part(bits, omap, nt=SC_NC):
    """The single-carrier part at 48 kHz: the silence, then the symbols (v7.4, the wide band: every other sample of
    SC_FS; the pulses leave nothing above 20.5 kHz to fold back)."""
    p = sc_passband(sc_symbols(bits, omap, nt))[::SC_FS // FS_TX]
    return np.clip(np.r_[np.zeros(round(SC_GAP * FS_TX)), p], -0.9, 0.9)


def sc_channel_bits(code, data, omap):
    """Coded, padded and scrambled bits, and nc (the coded bits)."""
    ch = encode_fec(code, data)
    nc, m = len(ch), OFDM_BITS[omap]
    ns = -(-nc // m)
    b = np.zeros(ns * m, dtype=int)
    b[:nc] = ch
    return b ^ prbs(ns * m), nc


# ---- receiver: linear algebra with FFTs
def nfft(n):
    return 1 << int(math.ceil(math.log2(max(n, 2))))


class SCConv:
    """n symbols -> the rows that see them, y = sum_k w_k (c_k * s): anchor channels c_k (L taps each) blended row
    by row with the weights w_k (one channel: w = 1). fwd and its adjoint, by FFTs."""
    def __init__(self, chs, w, n):
        self.L, self.n = len(chs[0]), n
        self.nr = n + self.L - 1
        self.N = nfft(self.nr)
        self.chs, self.C, self.w = chs, [np.fft.fft(c, self.N) for c in chs], w

    def fwd(self, s):
        S = np.fft.fft(s, self.N)
        return sum(w * np.fft.ifft(S * C)[:self.nr] for w, C in zip(self.w, self.C))

    def adj(self, r):
        return np.fft.ifft(sum(np.fft.fft(w * r, self.N) * np.conj(C) for w, C in zip(self.w, self.C)))[:self.n]

    def spectrum(self, n=4096):
        """|C(f)|^2 averaged over the anchors (the circulant approximation)."""
        return np.mean([np.abs(np.fft.fft(c, max(n, nfft(self.L)))) ** 2 for c in self.chs], axis=0)


class SCMask:
    """An SCConv seen through a row mask w (0 or 1): v5.1, rows lost with the samples the microphone dropped."""
    def __init__(self, op, w):
        self.op, self.w, self.n, self.L, self.nr = op, w, op.n, op.L, op.nr

    def fwd(self, s):
        return self.w * self.op.fwd(s)

    def adj(self, r):
        return self.op.adj(self.w * r)

    def spectrum(self, n=4096):
        return self.op.spectrum(n)


def cg(op, b, x0=None, it=200, tol=1e-6):
    """Conjugate gradients for op(x) = b, op Hermitian positive."""
    x = np.zeros_like(b) if x0 is None else x0.copy()
    r = b - op(x)
    p, rs, b2 = r.copy(), np.vdot(r, r).real, np.vdot(b, b).real
    for _ in range(it):
        if rs <= tol ** 2 * b2:
            break
        Ap = op(p)
        a = rs / np.vdot(p, Ap).real
        x += a * p
        r -= a * Ap
        rn = np.vdot(r, r).real
        p, rs = r + rn / rs * p, rn
    return x


def lbfgs(fg, x, iters=150, mem=10, tol=1e-5):
    """Minimizes f (fg returns f and its gradient) by L-BFGS with a backtracking line search."""
    f, g = fg(x)
    S, Y = [], []
    for _ in range(iters):
        q, al = g.copy(), []
        for s, y in reversed(list(zip(S, Y))):
            a = (s @ q) / (y @ s)
            al.append(a)
            q -= a * y
        if S:
            q *= (S[-1] @ Y[-1]) / (Y[-1] @ Y[-1])
        for (s, y), a in zip(zip(S, Y), reversed(al)):
            q += s * (a - (y @ q) / (y @ s))
        d = -q
        gd = g @ d
        if gd >= 0:
            d, gd = -g, -(g @ g)
        t = 1.0
        while True:
            xn = x + t * d
            fn, gn = fg(xn)
            if fn <= f + 1e-4 * t * gd or t < 1e-10:
                break
            t *= 0.5
        s, yv = xn - x, gn - g
        if s @ yv > 1e-12:
            S, Y = (S + [s])[-mem:], (Y + [yv])[-mem:]
        x, f, g = xn, fn, gn
        if np.linalg.norm(g) < tol * math.sqrt(len(x)):
            break
    return x


def sc_fit(x, y, mask, L, c0, lam, it=100):
    """Taps c (L of them): least squares over the rows where mask = 1 of y ~ x * c, plus lam |c - c0|^2."""
    T = SCConv([x], [1.0], L)                 # x * c = c * x: the taps are the unknown "symbols"
    r0 = mask * (y - T.fwd(c0))
    return c0 + cg(lambda v: T.adj(mask * T.fwd(v)) + lam * v, T.adj(r0), it=it, tol=1e-8)


# ---- receiver: detectors and the turbo decoder's soft output
def qam_llr(z, N, omap):
    """Exact bit LLRs (> 0 means 1) of z = s + noise of variance N, every point equally likely."""
    pts, lab = ofdm_points(omap)
    m = OFDM_BITS[omap]
    d = -np.abs(z[:, None] - pts[None, :]) ** 2 / N
    out = np.empty((len(z), m))
    for k in range(m):
        b = (lab >> (m - 1 - k)) & 1
        out[:, k] = np.logaddexp.reduce(d[:, b == 1], axis=1) - np.logaddexp.reduce(d[:, b == 0], axis=1)
    return out.ravel()


def sc_priors(lp, omap):
    """The decoder's beliefs (bit LLRs, > 0 means 1, m per symbol: I bits then Q bits) -> each axis's level probabilities."""
    amp, g, h = sc_axis(omap)
    Lb, out = lp.reshape(-1, 2 * h), []
    for ax in range(2):
        e = np.zeros((len(Lb), len(amp)))
        for k in range(h):
            bit, l = (g >> (h - 1 - k)) & 1, Lb[:, ax * h + k][:, None]
            e += np.where(bit[None, :] == 1, -np.logaddexp(0, -l), -np.logaddexp(0, l))
        p = np.exp(e - e.max(axis=1, keepdims=True))
        out.append(p / p.sum(axis=1, keepdims=True))
    return out, amp


def sc_lmmse(y, op, N0):
    """Linear MMSE estimate (biased) of the symbols and its bias beta (circulant approximation)."""
    x = cg(lambda v: op.adj(op.fwd(v)) + N0 * v, op.adj(y))
    return x, 1 - np.mean(1 / (1 + op.spectrum() / N0))


def sc_sic(y, op, N0, omap, lprior):
    """Linear MMSE with soft interference cancellation (Tuchler, Koetter and Singer; one prior variance for all
    symbols): the extrinsic LLRs of the bits."""
    (pI, pQ), amp = sc_priors(lprior, omap)
    xb = pI @ amp + 1j * (pQ @ amp)
    vb = max(float(np.mean(pI @ amp ** 2 - (pI @ amp) ** 2 + pQ @ amp ** 2 - (pQ @ amp) ** 2)), 1e-4)
    q = cg(lambda u: vb * op.fwd(op.adj(u)) + N0 * u, y - op.fwd(xb), it=300)
    lam = op.spectrum()
    kap = np.mean(lam / (vb * lam + N0))
    return qam_llr(xb + op.adj(q) / kap, max(1 / kap - vb, 1e-6), omap)


def sc_smap(y, op, N0, omap, lprior, x0, iters=150):
    """Soft SMAP: minimizes |y - c * s|^2 / N0 - sum ln p(s), p a mixture of Gaussians of variance tau^2 around each
    axis's levels, weighted by the decoder's beliefs; then each coordinate's local Gaussian (the data's precision
    at the minimum, from the circulant approximation with the prior's mean curvature) gives the extrinsic LLRs.
    Returns them and the minimum (the next pass starts there)."""
    (pI, pQ), amp = sc_priors(lprior, omap)
    _, g, h = sc_axis(omap)
    n = op.n
    lw = np.log(np.maximum(np.vstack([pI, pQ]), 1e-30))
    tau2, sz2 = SC_WIDTH[omap] * (amp[1] - amp[0]) ** 2, N0 / 2

    def prior(u):
        e = lw - (u[:, None] - amp[None, :]) ** 2 / (2 * tau2)
        mx = e.max(axis=1, keepdims=True)
        w = np.exp(e - mx)
        Z = w.sum(axis=1)
        w /= Z[:, None]
        mu = w @ amp
        return -(np.log(Z) + mx[:, 0]), mu, w @ amp ** 2 - mu ** 2

    def fg(v):
        r = op.fwd(v[:n] + 1j * v[n:]) - y
        gr = op.adj(r) / sz2
        pv, mu, _ = prior(v)
        return np.vdot(r, r).real / (2 * sz2) + pv.sum(), np.r_[gr.real, gr.imag] + (v - mu) / tau2

    x = lbfgs(fg, np.r_[x0.real, x0.imag], iters)
    _, mu, var = prior(x)
    D = 1 / tau2 - var / tau2 ** 2                        # the prior's own curvature at the minimum
    lam = op.spectrum() / sz2
    Db = D.mean()
    delta = max(0.0, 1e-3 - (lam.min() + Db))             # v7: Db, not D.min() (one symbol between two far levels raised
    qv = np.mean(1 / (lam + Db + delta))                  # every symbol's precision, ~2.4x from the second pass on)
    lamx = np.full(2 * n, max(1 / qv - Db - delta, 1e-3))   # the data's precision on each coordinate, the prior taken out
    mloc = x + (x - mu) / tau2 / lamx                     # where the data alone would put it
    d2 = -(mloc[:, None] - amp[None, :]) ** 2 * lamx[:, None] / 2
    Lb = np.empty((2 * n, h))
    for k in range(h):
        b = (g >> (h - 1 - k)) & 1
        Lb[:, k] = np.logaddexp.reduce(d2[:, b == 1], axis=1) - np.logaddexp.reduce(d2[:, b == 0], axis=1)
    return np.hstack([Lb[:n], Lb[n:]]).ravel(), x[:n] + 1j * x[n:]


def sc_soft(lpost, omap):
    """v7.2: the decoder's beliefs (bit LLRs, > 0 means 1, m per symbol) -> each symbol's mean and variance."""
    (pI, pQ), amp = sc_priors(lpost, omap)
    mI, mQ = pI @ amp, pQ @ amp
    return mI + 1j * mQ, pI @ amp ** 2 - mI ** 2 + pQ @ amp ** 2 - mQ ** 2


def sc_relearn(yy, x, w, L):
    """v7.2: one room of L taps for the whole part, by least squares over the rows where w = 1 of yy ~ x * c (x: every
    symbol, the data as the decoder's soft symbols, so that what is still unsure counts as noise), with a ridge of 1e-3
    times the rows' mean power. Returns the room and what it leaves on each row."""
    R = len(x) + L - 1
    P = float(np.sum(w * np.abs(yy[:R]) ** 2) / max(float(w.sum()), 1.0))
    c = sc_fit(x, yy[:R], w, L, np.zeros(L, dtype=complex), 1e-3 * P, it=150)
    return c, w * (yy[:R] - SCConv([x], [1.0], L).fwd(c))


_PRED = [[(s, u) for s in range(8) for u in (0, 1) if TB_NEXT[s, u] == q] for q in range(8)]
_PS, _PU = np.array([[a[0][0], a[1][0]] for a in _PRED]), np.array([[a[0][1], a[1][1]] for a in _PRED])


def _bcjr_app(ls, la, lp, lq=None):
    """As _bcjr, but returns the a-posteriori LLRs (ln P0/P1) of the input and of both parity bits at every step."""
    T = len(ls)
    g = 0.5 * ((ls + la)[:, None, None] * np.array([1, -1]) + lp[:, None, None] * (1 - 2 * TB_PAR))
    if lq is not None:
        g += 0.5 * lq[:, None, None] * (1 - 2 * TB_PAR2)
    Al = np.full((T + 1, 8), -1e9)
    Al[0, 0] = 0
    for t in range(T):
        a = Al[t][_PS] + g[t][_PS, _PU]
        Al[t + 1] = np.logaddexp(a[:, 0], a[:, 1])
        Al[t + 1] -= Al[t + 1].max()
    B, ou, op, oq = np.full(8, -1e9), np.empty(T), np.empty(T), np.empty(T)
    B[0] = 0
    p0, q0 = TB_PAR == 0, TB_PAR2 == 0
    for t in range(T - 1, -1, -1):
        gb = g[t] + B[TB_NEXT]
        m = Al[t][:, None] + gb
        ou[t] = np.logaddexp.reduce(m[:, 0]) - np.logaddexp.reduce(m[:, 1])
        op[t] = np.logaddexp.reduce(m[p0]) - np.logaddexp.reduce(m[~p0])
        oq[t] = np.logaddexp.reduce(m[q0]) - np.logaddexp.reduce(m[~q0])
        B = np.logaddexp(gb[:, 0], gb[:, 1])
        B -= B.max()
    return ou, op, oq


def turbo_siso(llr, n0, r=2, iters=8):
    """turbo_decode that also gives the a-posteriori LLR of every coded bit (> 0 means 1, in the order encode_fec
    sends them): what turbo equalization feeds back to the detector."""
    k, pi = n0, turbo_perm(n0)
    L = -np.clip(deinterleave(np.asarray(llr, dtype=float)), -30, 30)
    part = lambda j: L[j * k:(j + 1) * k]
    ls, tl, z3, t2 = L[:k], L[r * k:], np.zeros(3), 9 if r == 4 else 6
    odd = np.arange(k) % 2 == 1
    lp1, lp2 = (np.where(odd, 0, part(1)), np.where(odd, part(1), 0)) if r == 2 else (part(1), part(2))
    lq1 = np.r_[np.where(odd, 0, part(3)), tl[6:9]] if r == 4 else None
    lq2 = np.r_[np.where(odd, part(3), 0), tl[15:18]] if r == 4 else None
    ls2, le2, x = ls[pi], np.zeros(k), None
    dec1 = lambda le2: _bcjr_app(np.r_[ls, tl[0:3]], np.r_[le2, z3], np.r_[lp1, tl[3:6]], lq1)
    dec2 = lambda le1: _bcjr_app(np.r_[ls2, tl[t2:t2 + 3]], np.r_[le1[pi], z3], np.r_[lp2, tl[t2 + 3:t2 + 6]], lq2)
    for _ in range(iters):
        le1 = dec1(le2)[0][:k] - ls - le2
        a2 = dec2(le1)[0][:k]
        le2, xn = np.empty(k), np.empty(k, dtype=int)
        le2[pi], xn[pi] = a2 - ls2 - le1[pi], a2 < 0
        same = x is not None and np.array_equal(xn, x)
        x = xn
        if same:
            break
    u1, p1, q1 = dec1(le2)
    u2, p2, q2 = dec2(le1)
    s = ls + le1 + le2
    if r == 2:
        out = [s, np.where(odd, p2[:k], p1[:k]), u1[k:], p1[k:], u2[k:], p2[k:]]
    elif r == 3:
        out = [s, p1[:k], p2[:k], u1[k:], p1[k:], u2[k:], p2[k:]]
    else:
        out = [s, p1[:k], p2[:k], np.where(odd, q2[:k], q1[:k]), u1[k:], p1[k:], q1[k:], u2[k:], p2[k:], q2[k:]]
    return x, interleave(-np.concatenate(out))


def fec_siso(code, llr, n0):
    """The decided message and the a-posteriori LLR of every coded bit (> 0 means 1, in the order encode_fec sends
    them), for a turbo, (v7) an LDPC or (v7.3) a polar code: what turbo equalization feeds back to the detector; and the
    decoder's own a-posteriori LLRs, whose signs are its decisions (for LDPC and polar codes, app adds back what the
    decoder's input clipping took off)."""
    if code in POLAR_R:
        x = deinterleave(np.asarray(llr, dtype=float))
        bits, _, post = polar_decode(x, n0, POLAR_R[code], soft=True)
        return bits, interleave(post + x - np.clip(x, -POLAR_CLIP, POLAR_CLIP)), interleave(post)
    if code in LDPC_R:
        bits, app, _, _, post = ldpc_decode(deinterleave(np.asarray(llr, dtype=float)), n0, LDPC_R[code])
        return bits, interleave(app), interleave(post)
    bits, app = turbo_siso(llr, n0, TURBO_R[code])
    return bits, app, app


# ---- receiver: timing, the room, then turbo equalization
def sc_read(y, t, a=8):
    """Band-limited values of y at the (fractional) sample times t (Lanczos kernel)."""
    k0, out = np.floor(t).astype(int), np.zeros(len(t))
    for j in range(-a + 1, a + 1):
        k = k0 + j
        u = t - k
        ok = (k >= 0) & (k < len(y))
        out[ok] += (np.sinc(u) * np.sinc(u / a))[ok] * y[k[ok]]
    return out


def _env_xcorr(z, p):
    """|analytic cross-correlation| of z with p at lags 0 .. len(z) - len(p)."""
    N = nfft(len(z) + len(p))
    X = np.fft.fft(z, N) * np.conj(np.fft.fft(p, N))
    X[N // 2:] = 0
    X[1:N // 2] *= 2
    return np.abs(np.fft.ifft(X))[:len(z) - len(p) + 1]


def _peak(v, i):
    if 0 < i < len(v) - 1:
        a, b, c = v[i - 1], v[i], v[i + 1]
        den = a - 2 * b + c
        return i + (0.5 * (a - c) / den if den < 0 else 0.0)
    return float(i)


def sc_click(y, lo, hi, order=SC_GAP_AR, W=SC_GAP_W):
    """v5.1: the recording sample in [lo, hi) where y jumps. Samples lost or repeated by the microphone join two
    stretches of sound that do not follow each other; the sound itself fills only 3.2 to 16 kHz of the 24 kHz band,
    so a linear predictor fitted on the stretch follows it closely everywhere but there (the click detector of audio
    restoration). Predicted from the past, the W samples after the join go wrong; predicted from the future, the W
    before it: their squared errors summed. Returns the first sample after the join, the runner-up place (more than
    2 SC_GAP_SPAN away) and its score over the best's (near 1: unsure)."""
    lo, hi = max(int(lo), order), min(int(math.ceil(hi)), len(y) - order)
    seg = y[lo - order:hi + order]
    N = nfft(2 * len(seg))
    R = np.fft.irfft(np.abs(np.fft.rfft(seg, N)) ** 2, N)[:order + 1]
    T = R[np.abs(np.arange(order)[:, None] - np.arange(order)[None, :])] + 1e-9 * R[0] * np.eye(order)
    h = np.r_[1.0, -np.linalg.solve(T, R[1:order + 1])]
    m = hi - lo
    ef = np.convolve(seg, h)[order:order + m] ** 2               # sample lo + i from the order samples before it
    eb = np.convolve(seg[::-1], h)[order:order + m][::-1] ** 2   # from the order samples after it
    sc = np.convolve(ef, np.ones(W))[W - 1:W - 1 + m]            # lo + i .. lo + i + W - 1
    sc[1:] += np.convolve(eb, np.ones(W))[:m - 1]                # lo + i - W .. lo + i - 1
    i = int(np.argmax(sc))
    far = np.abs(np.arange(m) - i) > 2 * SC_GAP_SPAN
    j = int(np.argmax(np.where(far, sc, -1.0))) if far.any() else i
    return lo + i, lo + j, float(sc[j] / sc[i]) if j != i and sc[i] > 0 else 0.0


def sc_clock0(y, fs, T0, nt, L):
    """v5.1: the clock (delay gained per transmitted sample) from the first training block alone, for a jump in the
    only interval. The room fitted on each half of the block (rows that see TA and TC only), on all 5 sampling
    phases, predicts the other half; the delay between that prediction and the recording (cross-spectrum phase, as
    the time map) is the clock's drift between the halves, the room cancels. Its error on Sheng's recordings: 0.5 to
    5 ppm, up to 20 ppm when the phone or the speaker moved (sc_receive then tries the other size)."""
    rho, P5 = fs / SC_FS, SC_SPS
    x0 = np.concatenate(sc_training(nt))
    nA = len(x0)
    nrows = nA + L - 1
    n = -P5 * SC_PRE + np.arange(P5 * (nrows + 2 * SC_SPAN + 2))
    v = math.sqrt(2) * sc_read(y, T0 + rho * n) * np.exp(-2j * np.pi * SC_FC / SC_FS * n)
    mf = np.convolve(v, SC_G)
    mid = (L + nA) // 2
    hal = [(L, mid), (mid, nA)]
    Y, M = np.zeros((nrows, P5), complex), [np.zeros((nrows, P5), complex) for _ in hal]
    for p in range(P5):
        Y[:, p] = mf[2 * SC_DLY + p::P5][:nrows]
        for h, (a, b) in enumerate(hal):
            w = np.zeros(nrows)
            w[a:b] = 1
            M[h][:, p] = SCConv([x0], [1.0], L).fwd(sc_fit(x0, Y[:, p], w, L, np.zeros(L, dtype=complex), 0.0))
    nf = nfft(P5 * (mid - L))
    f = np.fft.fftfreq(nf, 1 / SC_FS)
    band = np.abs(f) < round(SC_RS * (1 + SC_BETA) / 2) - 400

    def delay(ys, ms):                                     # ys late on ms by this many samples
        X = (np.fft.fft(ys, nf) * np.conj(np.fft.fft(ms, nf)))[band]
        ph = lambda d: abs(np.sum(X * np.exp(2j * np.pi * f[band] * d / SC_FS)))
        u = SC_FS // FS_TX                                 # v7.4: as much time on the wide band's grid
        grid = np.linspace(-1.5, 1.5, 301) * u
        best = grid[int(np.argmax([ph(d) for d in grid]))]
        grid = np.linspace(best - 0.01 * u, best + 0.01 * u, 41)
        return grid[int(np.argmax([ph(d) for d in grid]))]
    (a0, b0), (a1, b1) = hal
    d10 = delay(Y[a1:b1].ravel(), M[0][a1:b1].ravel())
    d01 = delay(Y[a0:b0].ravel(), M[1][a0:b0].ravel())
    return (d10 - d01) / 2 / (P5 * (mid - L))


def sc_split(y, rho, t, p, u):
    """v5.1: a jump at sample u of the waveform p sent, which arrives near sample t of the recording y. The envelope of
    the correlation of the recording with p's part before u, against that of its part after (cross-correlated, so
    the room cancels): the jump in samples of the recording; and the first strong arrival of the part before."""
    G = SC_GAP_MAX
    lo, n = t - rho * (round(SC_T0_BACK * SC_FS) + G), round(2 * SC_T0_BACK * SC_FS) + 2 * G
    pb, pq = p[:u - SC_GAP_SPAN], p[u + SC_GAP_SPAN:]
    cb = _env_xcorr(sc_read(y, lo + rho * np.arange(n + len(pb))), pb)
    cq = _env_xcorr(sc_read(y, lo + rho * (u + SC_GAP_SPAN + np.arange(n + len(pq)))), pq)
    xc = np.correlate(np.r_[np.zeros(G), cq, np.zeros(G)], cb, "valid")   # lags -G .. G
    return (_peak(xc, int(np.argmax(xc))) - G) * rho, lo + rho * _peak(cb, sc_first_arrival(cb))


SC_T0_FRAC, SC_T0_BACK = 0.5, 0.015   # timing: the earliest arrival with half the strongest's correlation, up to 15 ms before it


def sc_first_arrival(env):
    """Index of the first strong arrival in the envelope correlation env: the earliest local peak at least SC_T0_FRAC
    of the strongest, within SC_T0_BACK seconds before it. Where a later echo is as strong as the direct sound, the
    strongest peak would put the early paths before the room's SC_PRE taps, where nothing is fitted."""
    i = int(np.argmax(env))
    for j in range(max(1, i - round(SC_T0_BACK * SC_FS)), i):
        if env[j] >= SC_T0_FRAC * env[i] and env[j] >= env[j - 1] and env[j] > env[j + 1]:
            return j
    return i


def sc_mmse_first(xA, y1, nA, L):
    """MMSE room on the first block, on the rows that see the training only (as least squares; the rows before it
    would add the control part's late echoes). Pilot: least squares on L taps; its residual gives N0 (echoes beyond
    L taps count as noise). Prior power of tap l: the pilot's |c_l|^2 averaged over SC_SMOOTH taps, minus the pilot's
    estimation floor f = N0 / (nA - L) (at least f / 10). Then a ridge with the weight N0 / p_l on tap l. Returns the
    taps, N0 of the fit, and N0 for detection: plus the taps' remaining error, sum p_l f / (p_l + f) (white training)."""
    mask = np.zeros(nA + L - 1)
    mask[L:nA] = 1
    c1 = sc_fit(xA, y1, mask, L, np.zeros(L, dtype=complex), 0.0, it=200)
    r = mask * (y1 - SCConv([xA], [1.0], L).fwd(c1))
    N0b = np.vdot(r, r).real / (nA - 2 * L)
    f = N0b / (nA - L)
    p = np.maximum(np.convolve(np.abs(c1) ** 2, np.ones(SC_SMOOTH) / SC_SMOOTH, mode="same") - f, f / 10)
    c0 = sc_fit(xA, y1, mask, L, np.zeros(L, dtype=complex), N0b / p, it=200)
    return c0, N0b, N0b + float(np.sum(p * f / (p + f)))


def sc_receive(y, fs, t0, omap, code, n0, nc, det="smap", verbose=False, nt=SC_NC, effort="fast", est="ls", swap=None):
    """The single-carrier part of a frame whose control part ends near sample t0 of the recording y.
    1. Timing: the first strong arrival of the first training block (its correlation with the waveform sent): the
       earliest peak with at least half the strongest's, up to 15 ms before it (a later echo can be the strongest).
    2. Time scale: the clean stretch of each block's TC (rows that see no data) is the same waveform each time; the
       delay of each against the first (the peak of their band-limited cross-correlation, to 1/1000 of a sample)
       gives the time map, linear between blocks: clock offset, a moving phone. The part is read on that map.
       v5.1: the delay is also searched SC_GAP_MAX samples farther; a block that agrees clearly better there marks
       samples lost or repeated by the microphone since the block before. Their number is the jump less the clock's
       share (the clock from the nearest interval without a jump, else from the first block alone, sc_clock0),
       rounded to whole samples of the recording. Where they went missing: where the prediction errors forward and
       backward peak (sc_click), from the middle of the block before to the middle of the block after; if the frame
       then fails, it is tried again with the other whole number or the runner-up place, whichever was less
       certain. The map then jumps there, and the rows the jump
       smears are masked: the room is fitted without them, and the detector neither uses them nor trusts the
       symbols they carry (their LLRs are 0, the code fills them in).
       A jump inside the first block, whose TC is the reference: the arrival of the block's part before the
       prediction error's peak against that of the part after (step 1, the envelopes of their correlations
       cross-correlated); the timing is then that of the part before. The same in the last block, which no later
       block checks.
    3. Down-conversion, matched filter, symbol-rate samples for each of the 5 sampling phases; the phase whose room
       (fitted on the first block) would carry the most, mean of log2(1 + |C(f)|^2 / N0), is kept: the 5 fits leave
       the same error to 0.1 dB, but the room folded at the symbol rate differs from phase to phase by 1 to 1.5 dB
       for detection.
    4. The room: 256, 512 or 1024 taps (SC_TAPS, as the TC of nt symbols allows) by least squares on the first
       block (rows that see no data), the noise from what they leave; with the long TC, 1024 taps if they leave
       less, their own error counted; each later block by ridge regression from the one before (rows that see TC
       only); over the data between two blocks, the two blended linearly row by row. With est 'mmse', the first
       block's room is instead the most taps the TC allows, by MMSE (sc_mmse_first).
    5. Each data segment's rows, the training's echoes removed; detector (soft SMAP or linear MMSE) on each segment,
       turbo decoder over the whole message, its beliefs back to the detector; stops once the decisions form a
       codeword or after the passes of the effort chosen (SC_EFFORT).
    6. v7.2, the second try: if the passes end without a codeword, one room of SC_RELEARN_TAPS taps for the whole part,
       by least squares over every row the recording holds, the data symbols as the decoder's soft symbols (the mean
       under its last beliefs); the noise of each segment from what it leaves there, less the soft symbols' variance,
       the taps' degrees of freedom counted. Each segment then sees every other symbol removed (the other segments'
       as soft symbols), and at most SC_RELEARN_PASSES more passes follow, the room learned again after each."""
    import time as _time
    t_start = _time.time()
    m = OFDM_BITS[omap]
    ns = -(-nc // m)
    blocks, segs, ntot = sc_layout(ns, nt)
    rho, L, P5 = fs / SC_FS, SC_TAPS[nt][0], SC_SPS
    ta, tc = sc_training(nt)
    x = sc_symbols(np.zeros(ns * m, dtype=int), omap, nt)
    steps, npass = SC_EFFORT[effort]
    known = np.ones(ntot, dtype=bool)
    for a, b in segs:
        x[a:b], known[a:b] = 0, False
    # 1. the first block's first strong arrival: frame time 0 (the part's first sample after the silence) arrives at T0
    pa = sc_passband(np.r_[ta, tc])
    g0 = t0 + SC_GAP * fs
    lo = g0 - 0.03 * fs
    nz = len(pa) + round(0.09 * SC_FS)
    env = _env_xcorr(sc_read(y, lo + np.arange(nz) * rho), pa)
    T0 = lo + _peak(env, sc_first_arrival(env)) * rho
    # v5.1: a jump inside the first block (its TC is every later block's reference): the arrival of the block's part
    # before the prediction error's peak against that of the part after (the envelopes of their correlations with
    # what was sent, cross-correlated; the room cancels)
    gaps0, unsure = [], []                                 # unsure: (how unsure, (gap, "size" or "place")) for the second try
    ub = SC_SPS * (SC_NA + nt) + 2 * SC_DLY - SC_GAP_EDGE
    if ub > SC_SPS * L:
        tg = sc_click(y, T0 + rho * SC_SPS * L, T0 + rho * ub)[0]
        raw, Tb = sc_split(y, rho, T0, pa, int((tg - T0) / rho))
        if abs(raw) > SC_GAP_MIN:
            T0 = Tb                                            # the arrival of the part before it
            sz = round(raw)
            if swap == (0, "size"):                            # the second try: the other whole number
                sz += 1 if raw > sz else -1
            unsure.append((2 * abs(raw - round(raw)), (0, "size")))
            gaps0.append(((tg - T0) / rho, sz / rho))
    J0 = lambda n: sum(sz * (np.asarray(n, dtype=float) >= ng) for ng, sz in gaps0)
    # 2. the time map T(n) = T0 + rho (n + D(n)), D through the blocks' delays
    starts = [SC_NA] + [b for b, _ in blocks[1:]]           # where each TC starts (symbols)
    w0 = SC_DLY * 2 + P5 * (L - SC_PRE)                    # clean stretch of a TC, in samples after its start
    w1 = P5 * (nt - SC_PRE)
    win = lambda k: P5 * starts[k] + np.arange(w0, w1)
    f = np.fft.rfftfreq(nfft(w1 - w0), 1 / SC_FS)
    band = (f > round(SC_FC - SC_RS * (1 + SC_BETA) / 2) + 200) & (f < round(SC_FC + SC_RS * (1 + SC_BETA) / 2) - 200)
    W0 = sc_read(y, T0 + rho * (win(0) + J0(win(0))))
    F0 = np.conj(np.fft.rfft(W0, nfft(w1 - w0))[band])
    def fine(k, at):                                       # sub-sample delay of block k near at, and the coherence
        Wk = sc_read(y, T0 + rho * (win(k) + at))
        X = np.fft.rfft(Wk, nfft(w1 - w0))[band] * F0
        ph = lambda d: abs(np.sum(X * np.exp(2j * np.pi * f[band] * d / SC_FS)))
        u = SC_FS // FS_TX                                 # v7.4: as much time on the wide band's grid
        grid = np.linspace(-1, 1, 201) * u
        best = grid[int(np.argmax([ph(d) for d in grid]))]
        grid = np.linspace(best - 0.01 * u, best + 0.01 * u, 41)
        v = [ph(d) for d in grid]
        i = int(np.argmax(v))
        return grid[i], v[i] / max(float(np.sum(np.abs(X))), 1e-30)

    mk, dk, cand = [float(win(0).mean())], [float(J0(win(0).mean()))], []   # cand (v5.1): (block after a jump, its sample in the recording, runner-up, unsure)
    dc = lambda i: dk[i] - J0(mk[i])                       # the clock alone, but for the jumps found after it
    for k in range(1, len(starts)):
        c = float(win(k).mean())
        free = [i for i in range(1, len(dk)) if i not in [g[0] for g in cand]]   # intervals without a jump: the clock
        slope = (dc(free[-1]) - dc(free[-1] - 1)) / (mk[free[-1]] - mk[free[-1] - 1]) if free else 0.0
        pred = dk[-1] + slope * (c - mk[-1]) + J0(c) - J0(mk[-1])
        S = round(SC_MAX_PPM * 1e-6 * (c - mk[-1]) * (1 if len(dk) > 1 else 2)) + 8   # the lags searched
        zl = sc_read(y, T0 + rho * (P5 * starts[k] + np.arange(w0 - S, w1 + S) + pred))
        e = _env_xcorr(zl, W0)
        lag = _peak(e, int(np.argmax(e))) - S
        best, coh = fine(k, pred + lag)
        Sw = S + SC_GAP_MAX                                # v5.1: do the blocks agree best farther away? a gap?
        e = _env_xcorr(sc_read(y, T0 + rho * (P5 * starts[k] + np.arange(w0 - Sw, w1 + Sw) + pred)), W0)
        j = int(np.argmax(e))
        if abs(j - Sw) > S:
            lag2 = _peak(e, j) - Sw
            best2, coh2 = fine(k, pred + lag2)
            if coh2 > SC_GAP_COH * coh:
                mid, ex = (w0 + w1) / 2, (w1 - w0) / 5         # a jump early in a block moves the block, a late one the next
                lo = T0 + rho * (P5 * starts[k - 1] + mid - ex + dk[-1])
                hi = T0 + rho * (P5 * starts[k] + mid + ex + pred + lag2 + best2)
                cand.append((k,) + sc_click(y, lo, hi))
                lag, best = lag2, best2
        mk.append(c)
        dk.append(pred + lag + best)
    mk, dk = np.array(mk), np.array(dk)
    # v5.1: a jump's size is the delay it adds over its interval, less the clock's share, in whole samples of the
    # recording; the clock from the nearest interval without a jump, else from the first block alone (sc_clock0)
    gaps, ks, clk0 = list(gaps0), [g[0] for g in cand], None
    for k, tg, alt, un in cand:
        near = sorted((i for i in range(1, len(mk)) if i not in ks), key=lambda i: abs(i - k))
        if near:
            slope = (dc(near[0]) - dc(near[0] - 1)) / (mk[near[0]] - mk[near[0] - 1])
        else:
            clk0 = sc_clock0(y, fs, T0, nt, L) if clk0 is None else clk0
            slope = clk0
        raw = (dc(k) - dc(k - 1) - slope * (mk[k] - mk[k - 1])) * rho
        sz = round(raw)
        if swap == (len(gaps), "size"):                    # the second try: the other whole number
            sz += 1 if raw > sz else -1
        if swap == (len(gaps), "place"):                   # or the runner-up place
            tg = alt
        unsure += [(2 * abs(raw - round(raw)), (len(gaps), "size")), (un, (len(gaps), "place"))]
        u = (tg - T0) / rho                                # the jump's transmitted sample, on the map before it
        ng = u - (dk[k - 1] + slope * (u - mk[k - 1]))
        ng = u - (dk[k - 1] + slope * (ng - mk[k - 1]))
        gaps.append((float(ng), sz / rho))
    # v5.1: a jump inside the last block's clean stretch, which no later block sees
    tb = T0 + rho * (P5 * blocks[-1][0] + dk[-1])          # where the last block arrives
    lo, hi = tb + rho * w0, tb + rho * (P5 * nt + 2 * SC_DLY - SC_GAP_EDGE)
    if len(mk) > 1 and hi > lo and len(mk) - 1 not in ks:
        tg = sc_click(y, lo, hi)[0]
        raw, _ = sc_split(y, rho, tb, sc_passband(tc), int((tg - tb) / rho))
        if abs(raw) > SC_GAP_MIN:
            sz = round(raw)
            if swap == (len(gaps), "size"):
                sz += 1 if raw > sz else -1
            unsure.append((2 * abs(raw - round(raw)), (len(gaps), "size")))
            gaps.append(((tg - T0) / rho - dk[-1], sz / rho))
    # the clock alone: less the jumps each block's delay holds (those of the blocks that saw them, wherever they fell)
    dcl = dk - np.array([J0(mm) + sum(g[1] for g, b in zip(gaps[len(gaps0):], ks) if b <= i) for i, mm in enumerate(mk)])

    def D(n):
        n = np.asarray(n, dtype=float)
        if len(mk) == 1:
            out = np.zeros_like(n)
        else:
            out = np.interp(n, mk, dcl)
            s0, s1 = (dcl[1] - dcl[0]) / (mk[1] - mk[0]), (dcl[-1] - dcl[-2]) / (mk[-1] - mk[-2])
            out = np.where(n < mk[0], dcl[0] + s0 * (n - mk[0]), np.where(n > mk[-1], dcl[-1] + s1 * (n - mk[-1]), out))
        for ng, sz in gaps:
            out = out + sz * (n >= ng)
        return out
    ppm = float((dcl[-1] - dcl[0]) / (mk[-1] - mk[0]) * 1e6) if len(mk) > 1 else 0.0
    lin = np.polyval(np.polyfit(mk, dcl, 1), mk) if len(mk) > 2 else dcl
    bend = float(np.max(np.abs(dcl - lin)))                 # departure from a steady clock (samples)
    # 3. symbol-rate samples, row i <-> symbol i - SC_PRE at its first strong arrival
    nrows = ntot + max(max(SC_TAPS[nt]), SC_RELEARN_TAPS) - 1   # v7.2: as far as the second try's room reaches
    n = -P5 * SC_PRE + np.arange(P5 * (nrows + 2 * SC_SPAN + 2))
    tpos = T0 + rho * (n + D(n.astype(float)))
    z = sc_read(y, tpos)
    v = math.sqrt(2) * z * np.exp(-2j * np.pi * SC_FC / SC_FS * n)
    mf = np.convolve(v, SC_G)
    phases = [mf[2 * SC_DLY + p::P5][:nrows] for p in range(P5)]
    good = np.ones(nrows)                                  # v5.1: 0 on the rows a gap smears (row r: sample 5 (r - 64) + 60)
    nr_ = P5 * (np.arange(nrows) - SC_PRE) + SC_DLY
    for ng, sz in gaps:
        good[(nr_ > ng - SC_GAP_SPAN - P5) & (nr_ < ng + max(0.0, -sz) + SC_GAP_SPAN + P5)] = 0
    nA = SC_NA + nt
    mask0 = np.zeros(nA + L - 1)
    mask0[L:nA] = 1
    mask0 *= good[:nA + L - 1]
    best = None
    for p, yy in enumerate(phases):
        c0 = sc_fit(x[:nA], yy[:nA + L - 1], mask0, L, np.zeros(L, dtype=complex), 0.0)
        r = mask0 * (yy[:nA + L - 1] - SCConv([x[:nA]], [1.0], L).fwd(c0))
        e = np.vdot(r, r).real
        rate = float(np.mean(np.log2(1 + np.abs(np.fft.fft(c0, 4096)) ** 2 / (e / (nA - 2 * L)))))
        if best is None or rate > best[0]:
            best = (rate, p, c0, e)
    _, phase, c0, e = best
    yy = phases[phase]
    nfit = nA - L
    N0b = e / (nfit - L)
    N0 = N0b * (1 + L / nfit)                              # plus the error of the taps' estimate
    for L2 in SC_TAPS[nt][1:]:                             # a longer room, if it leaves less (its own error counted)
        m2 = np.zeros(nA + L2 - 1)
        m2[L2:nA] = 1
        m2 *= good[:nA + L2 - 1]
        c2 = sc_fit(x[:nA], yy[:nA + L2 - 1], m2, L2, np.zeros(L2, dtype=complex), 0.0, it=200)
        r2 = m2 * (yy[:nA + L2 - 1] - SCConv([x[:nA]], [1.0], L2).fwd(c2))
        e2, nf2 = np.vdot(r2, r2).real, nA - L2
        if e2 / (nf2 - L2) * (1 + L2 / nf2) < N0:
            L, c0, e, nfit = L2, c2, e2, nf2
            N0b = e / (nfit - L)
            N0 = N0b * (1 + L / nfit)
    if est == "mmse":
        L = max(SC_TAPS[nt])
        c0, N0b, N0 = sc_mmse_first(x[:nA], yy[:nA + L - 1], nA, L)
    # 4. the room at each block
    chs, tk, change = [c0], [(L + nA) / 2], []
    for k, (b0, b1) in enumerate(blocks[1:], 1):
        last = k == len(blocks) - 1
        mk_ = np.zeros(nt + L - 1)
        mk_[L - 1:nt + (L - 1 if last else 0)] = 1         # the last block: the rows after it see silence
        mk_ *= good[b0:b0 + nt + L - 1]
        cp = chs[-1]
        ck = sc_fit(tc, yy[b0:b0 + nt + L - 1], mk_, L, cp, N0b / (SC_REL * np.vdot(cp, cp).real / L), it=100 if L <= 512 else 200)
        change.append(float(10 * np.log10(np.vdot(ck - cp, ck - cp).real / np.vdot(cp, cp).real)))
        chs.append(ck)
        tk.append(b0 + (L - 1 + nt + (L - 1 if last else 0)) / 2)
    # the noise of each training block: its room fitted on its own (least squares on its rows that see training only),
    # the first block's as above; each data segment's detector takes the mean of its two blocks', and the SNR shown is
    # the median over the blocks, so a disturbed block (a hand moving) sets only its neighbours
    n0k, snrs = [N0], [float(10 * np.log10(np.vdot(c0, c0).real / N0))]
    for k, (b0, b1) in enumerate(blocks[1:], 1):
        last = k == len(blocks) - 1
        mk_ = np.zeros(nt + L - 1)
        mk_[L - 1:nt + (L - 1 if last else 0)] = 1
        mk_ *= good[b0:b0 + nt + L - 1]
        nk = float(mk_.sum())
        if nk < L + 100:                                   # too few rows left (a gap): the block before's
            n0k.append(n0k[-1])
            continue
        yk = yy[b0:b0 + nt + L - 1]
        ck = sc_fit(tc, yk, mk_, L, np.zeros(L, dtype=complex), 0.0, it=100 if L <= 512 else 200)
        r = mk_ * (yk - SCConv([tc], [1.0], L).fwd(ck))
        n0k.append(float(np.vdot(r, r).real / (nk - L) * (1 + L / nfit)))   # plus the taps' error, as the first block's
        snrs.append(float(10 * np.log10(np.vdot(ck, ck).real / n0k[-1])))
    N0s = [(n0k[j] + n0k[j + 1]) / 2 for j in range(len(segs))]
    # 5. each segment: rows S .. E + L - 2, the training's part removed; its two anchors blended
    obs, ops = [], []
    for j, (S, E) in enumerate(segs):
        rows = np.arange(S, E + L - 1)
        wb = np.clip((rows - tk[j]) / (tk[j + 1] - tk[j]), 0, 1)
        w = [1 - wb, wb]
        xk = x[S - L + 1:E + L - 1]
        yk = sum(wi * SCConv([xk], [1.0], L).fwd(c)[L - 1:L - 1 + len(rows)] for wi, c in zip(w, chs[j:j + 2]))
        gw = good[S:E + L - 1]
        obs.append(gw * (yy[S:E + L - 1] - yk))
        ops.append(SCConv(chs[j:j + 2], w, E - S) if gw.all() else SCMask(SCConv(chs[j:j + 2], w, E - S), gw))
    erased = np.repeat(good[SC_PRE:SC_PRE + ntot] == 0, m)   # symbols whose first arrival is on a masked row
    key = prbs(ns * m)
    flip = np.where(key == 1, -1.0, 1.0)
    lprior, xs, first, passes, hist = np.zeros(ns * m), [None] * len(segs), None, 0, []
    cw, relearn = None, None
    for it in range(npass + SC_RELEARN_PASSES):
        ext, q = [], 0
        for j, (S, E) in enumerate(segs):
            lp = lprior[q * m:(q + E - S) * m]
            if det == "lmmse":
                ext.append(sc_sic(obs[j], ops[j], N0s[j], omap, lp))
            else:
                if xs[j] is None:
                    xs[j] = sc_lmmse(obs[j], ops[j], N0s[j])[0]
                le, xs[j] = sc_smap(obs[j], ops[j], N0s[j], omap, lp, xs[j], steps[0] if it == 0 else steps[1])
                ext.append(le)
            q += E - S
        le_ = np.concatenate(ext)
        le_[erased[np.concatenate([np.arange(S * m, E * m) for S, E in segs])]] = 0.0
        lin_ = (flip * le_)[:nc]
        if first is None:
            first = lin_ > 0
        bits, app, dec = fec_siso(code, lin_, n0)
        passes = it + 1
        cw = encode_fec(code, bits)
        hist.append(int(np.sum((dec > 0) != (cw > 0))))          # coded bits that disagree with the decided message
        if hist[-1] == 0:
            break
        e = np.clip(app - lin_, -20, 20)
        lprior = flip * np.r_[e, np.zeros(ns * m - nc)]
        if npass <= it + 1 < npass + SC_RELEARN_PASSES:  # v7.2: the second try (6.)
            L2 = SC_RELEARN_TAPS
            xm, xv = sc_soft(flip * np.r_[np.clip(app, -30, 30), np.full(ns * m - nc, -30.0)], omap)
            xsoft, vs, q = x.copy(), np.zeros(ntot), 0
            for S, E in segs:
                xsoft[S:E], vs[S:E] = xm[q:q + E - S], xv[q:q + E - S]
                q += E - S
            ins = tpos < len(y) - 8                        # the rows the recording holds
            w2 = good[:ntot + L2 - 1] * ins[np.minimum(2 * SC_DLY + phase + P5 * np.arange(ntot + L2 - 1), len(ins) - 1)]
            c2, res2 = sc_relearn(yy, xsoft, w2, L2)
            P2, nw = float(np.vdot(c2, c2).real), float(w2.sum())
            obs, ops, N0s = [], [], []
            for S, E in segs:
                xo = xsoft.copy()
                xo[S:E] = 0
                gw = w2[S:E + L2 - 1]
                obs.append(gw * (yy[S:E + L2 - 1] - SCConv([xo], [1.0], L2).fwd(c2)[S:E + L2 - 1]))
                ops.append(SCConv([c2], [1.0], E - S) if gw.all() else SCMask(SCConv([c2], [1.0], E - S), gw))
                e2 = float(np.sum(np.abs(res2[S:E + L2 - 1]) ** 2)) / max(float(gw.sum()), 1.0) * nw / max(nw - L2, nw / 4)
                N0s.append(max(e2 - P2 * float(np.mean(vs[S:E])), 0.3 * e2))
            xs = [None] * len(segs)
            relearn = dict(taps=L2, snr=float(10 * np.log10(P2 / np.mean(N0s))))
    if relearn is not None:
        relearn["passes"] = passes - npass
    out = dict(bits=bits, passes=passes, codeword=bool(np.array_equal((dec > 0).astype(int), cw)),
               chan_errors=int(np.sum(first != (cw > 0))), snr=float(np.median(snrs)), snrs=snrs,
               ppm=ppm, bend=bend, phase=phase, change=change, taps=L, train=nt, effort=effort, est=est, nseg=len(segs), T0=T0,
               gaps=[(ng / SC_FS, sz) for ng, sz in gaps], erased=int(sum(np.sum(good[SC_PRE + a:SC_PRE + b] == 0) for a, b in segs)), swap=swap, sec=_time.time() - t_start, det=det, hist=hist, relearn=relearn)
    if gaps and not out["codeword"] and swap is None:     # v5.1: the other size, or the runner-up place, of the least sure gap
        o2 = sc_receive(y, fs, t0, omap, code, n0, nc, det, False, nt, effort, est, swap=max(unsure)[1])
        if o2["codeword"]:
            out = dict(o2, sec=_time.time() - t_start)
    if verbose:
        print(out)
    return out


def sc_desc(omap, nt=SC_NC, wide=False):
    rs, _, fc, _ = SC_BANDS[wide]
    train, taps, d = (SC_SIZES[wide][k] for k in (2, 4, 5))
    name = [k for k, v in train.items() if v == nt][0]
    return (f"single carrier at {fc / 1000:g} kHz, {rs:.0f} symbols/s (roll-off {SC_BETA:.2f}, {fc - rs * (1 + SC_BETA) / 2:.0f} to "
            f"{fc + rs * (1 + SC_BETA) / 2:.0f} Hz), {OFDM_NAMES[omap]}, {name} training ({nt} symbols, rooms up to "
            f"{max(taps[nt]) - SC_PRE} taps after the first strong arrival, {(max(taps[nt]) - SC_PRE) / rs * 1000:.0f} ms) every "
            f"{d / rs:.1f} s, {sc_rate(omap, nt, wide):.0f} bit/s before the code")


# ---------------------------------------------------------------- receiver
def tone_energies(y, fs, rate, tones="low", M=2):
    """|sum of y(n) exp(-j 2 pi F n / fs)| over a window of one symbol starting at each sample, for each tone F."""
    Lw = round(fs * math.log2(M) / rate)
    n = np.arange(len(y))
    out = []
    for F in tone_set(tones, M):
        c = np.concatenate([[0], np.cumsum(y * np.exp(-2j * np.pi * F * n / fs))])
        out.append(np.abs(c[Lw:] - c[:-Lw]))      # window [m, m+Lw-1], m = 0 .. len(y)-Lw
    return out


def pilot_correlation(d, Lr, P):
    off = np.round(np.arange(len(P)) * Lr).astype(int)
    M = max(0, len(d) - off[-1])                  # a recording shorter than the preamble has no start for it
    c = np.zeros(M)
    for j in range(len(P)):
        c += (2 * P[j] - 1) * d[off[j]:off[j] + M]
    return c / len(P)


# ---------------------------------------------------------------- v6: the short control part (OFDM and single carrier)
# The 400 bit/s control part lasts 0.65 s, more than an OFDM or single-carrier message of a few words: it was made for
# FSK, whose frames must be found far away and in reverberant rooms. On Sheng's recordings it survives noise 7 to 10 dB
# above its own level (in-band, 3.2 to 16 kHz), 25 dB more than the OFDM and single-carrier messages there need. A
# wideband control part is 3.4 times shorter:
#   P P P | H, 0.192 s: three periods of P (C6_N samples at 48 kHz: OFDM's reference symbol for N = 2048 with its phases
#   negated, so that its sweep runs the other way and OFDM's own reference symbols do not pass for it), then the header
#   symbol H, its guard of C6_CP samples, then C6_N.
# H carries the 92 coded header bits by BPSK on P's subcarriers: subcarrier j (in increasing frequency) is P's value
# times +1 for bit 0, -1 for bit 1, bit (37 j) mod 92: each bit has 5 or 6 copies, 2.2 kHz apart, and consecutive bits
# fall far apart. The first period lets the speaker and the microphone settle and fills the room; the last two, alike
# once the room is full, give the receiver the room, against which it reads H (no channel estimate to learn first,
# and the room's phases and the timing error cancel). Played in the rooms of Sheng's recordings (the room fitted on
# each frame's 400 bit/s control part), it was found and read in all 30 frames tried. It survived about as much noise
# as the 400 bit/s one in Impasse Carnot (5 dB less to 3 dB more), 3 to 8 dB less in Rue Joliot Curie (3 dB of it
# because it plays 3 dB quieter: rms 0.45 like OFDM, against FSK's 0.64), still 17 dB more than 16QAM needs there.
C6_N, C6_CP, C6_REP, C6_ILV = 2048, 1024, 3, 37
C6_LEN = C6_REP * C6_N + C6_CP + C6_N            # 9216 samples, 0.192 s
C6_THR, C6_BAL = 0.02, 0.5                       # detection: rho above C6_THR and bal above C6_BAL (see c6_detect)
C6_K = ofdm_carriers(C6_N)
C6_PH = -ofdm_reference(len(C6_K))


def c6_symbol(X):
    """Values on the subcarriers C6_K -> one symbol of C6_N samples, at OFDM's level (rms 0.45 for unit values)."""
    F = np.zeros(C6_N // 2 + 1, dtype=complex)
    F[C6_K] = X
    return np.fft.irfft(F, C6_N) * (C6_N / 2 * OFDM_RMS * math.sqrt(2 / len(C6_K)))


def c6_index(n=NHC):
    """The coded header bit that each subcarrier of H carries."""
    return (C6_ILV * np.arange(len(C6_K))) % n


def c6_wave(hc):
    """The 92 coded header bits -> the short control part at 48 kHz. P's peaks are 0.84; H's are clipped at 0.9 like OFDM's."""
    S = np.exp(1j * C6_PH)
    p, h = c6_symbol(S), c6_symbol((1 - 2 * np.asarray(hc)[c6_index(len(hc))]) * S)
    return np.clip(np.concatenate([np.tile(p, C6_REP), h[-C6_CP:], h]), -0.9, 0.9)


def c6_detect(y):
    """y at 48 kHz. rho[n], for P starting at sample n: the energy that the matched filter of one period collects at n,
    n + N and n + 2N, over that of P and over 3 times the in-band energy of the loudest of the three periods. It is the
    share of the energy that a single path explains: 1 for P alone at any level, about 1/1100 for noise, and in a room
    the share of its strongest path (0.05 in Impasse Carnot, 0.17 to 0.35 in Rue Joliot Curie); 0.013 at most over
    Sheng's recordings of v5 frames. A burst followed by silence gets no credit for the silence. bal[n] = c c' /
    mean(c^2, c'^2), c and c' the matched filter at n + N and n + 2N: 1 for P (the same room for both periods), random
    for noise."""
    if len(y) < C6_LEN:
        return np.zeros(0), np.zeros(0)
    F = np.fft.rfft(y)
    f = np.fft.rfftfreq(len(y), 1 / FS_TX)
    F[(f < OFDM_F0) | (f >= OFDM_F1)] = 0
    yb = np.fft.irfft(F, len(y))
    p = c6_symbol(np.exp(1j * C6_PH))
    L = nfft(len(yb) + C6_N)
    c = np.fft.irfft(np.fft.rfft(yb, L) * np.conj(np.fft.rfft(p, L)), L)[:len(yb) - C6_N + 1]   # c[n] = sum_m yb[n + m] p[m]
    q = c * c
    M = len(q) - (C6_REP - 1) * C6_N
    D = sum(q[j * C6_N:j * C6_N + M] for j in range(C6_REP))
    e = np.concatenate([[0.0], np.cumsum(yb * yb)])
    Emax = np.maximum.reduce([e[(j + 1) * C6_N:(j + 1) * C6_N + M] - e[j * C6_N:j * C6_N + M] for j in range(C6_REP)])
    rho = D / (np.sum(p * p) * C6_REP * Emax + 1e-30)
    c1, c2 = c[C6_N:C6_N + M], c[2 * C6_N:2 * C6_N + M]
    return rho, c1 * c2 / (0.5 * (c1 * c1 + c2 * c2) + 1e-30)


def c6_read(y, n):
    """y at 48 kHz, P detected at sample n. The room from P's last two periods (cyclic, N lags), its first strong
    arrival, then H read against them. Returns dict(hdr, hb, t = P's first arrival, drift = samples per period,
    late = the room's energy that reaches H from P), or None."""
    S = np.exp(1j * C6_PH)
    if n < 0 or n + C6_REP * C6_N > len(y):
        return None
    Y2, Y3 = (np.fft.rfft(y[n + j * C6_N:n + (j + 1) * C6_N])[C6_K] for j in (1, 2))
    G = np.zeros(C6_N, dtype=complex)
    G[C6_K] = 0.5 * (Y2 + Y3) / S * np.hanning(len(C6_K) + 2)[1:-1]
    g = np.abs(np.fft.ifft(G)) ** 2                          # the room's envelope, cyclic over N lags
    # the first strong arrival, as the single carrier's (sc_first_arrival): the earliest local peak with half the
    # strongest's energy, up to SC_T0_BACK before it. The energy outside 64 lags before it to the guard's end reaches H
    # from P.
    i = int(np.argmax(g))
    s = i
    for j in range(i - round(SC_T0_BACK * FS_TX), i):
        if g[j % C6_N] >= SC_T0_FRAC * g[i] and g[j % C6_N] >= g[(j - 1) % C6_N] and g[j % C6_N] > g[(j + 1) % C6_N]:
            s = j % C6_N
            break
    late = float(1 - np.sum(g[(s + np.arange(-64, C6_CP + 1)) % C6_N]) / np.sum(g))
    t = n + (s if s < C6_N // 2 else s - C6_N)
    if t < 0 or t + C6_LEN > len(y):
        return None
    Y2, Y3, YH = (np.fft.rfft(y[t + o:t + o + C6_N])[C6_K] for o in (C6_N, 2 * C6_N, 3 * C6_N + C6_CP))
    # the clock: Y3 = Y2 exp(-j 2 pi k d / N), d samples per period
    dd = np.linspace(-1.5, 1.5, 301)
    d = float(dd[int(np.argmax(np.abs(np.exp(2j * np.pi * np.outer(dd, C6_K) / C6_N) @ (Y3 * np.conj(Y2)))))])
    rot = np.exp(-2j * np.pi * C6_K * d / C6_N)
    ker = np.ones(33) / 33
    v = np.convolve(np.abs(Y3 - Y2 * rot) ** 2 / 2, ker, "same") + 1e-30     # the noise of one period's FFT
    # the room from both periods, its impulse response kept from 64 lags before the first arrival to the end of the
    # guard (the noise of this reference is v / 2 (CP + 129) / N), carried to H's window by the clock
    F = np.zeros(C6_N // 2 + 1, dtype=complex)
    F[C6_K] = 0.5 * (Y3 + Y2 * rot) / S
    h = np.fft.irfft(F, C6_N)
    h[C6_CP + 65:C6_N - 64] = 0
    ref = np.fft.rfft(h)[C6_K] * S * np.exp(-2j * np.pi * C6_K * d * (C6_N + C6_CP) / C6_N ** 2)
    isi = late * np.convolve(np.abs(ref) ** 2, ker, "same")
    llr = -4 * np.real(YH * np.conj(ref)) / (v * (1 + (C6_CP + 129) / (2 * C6_N)) + isi)   # > 0: bit 1
    hb = [int(b) for b in decode_fec("conv", HBITS, np.bincount(c6_index(), weights=llr, minlength=NHC))[0]]
    return dict(hdr=parse_header(hb), hb=hb, t=t, drift=d, late=late)


def parse_header(hb):
    """hb: the 40 decoded header bits. None unless the CRC and the fields are valid."""
    crc, rx = crc8(hb[:32]), [int(x) for x in hb[32:40]]
    x = int("".join(str(c ^ r) for c, r in zip(crc, rx)), 2)   # the mask: 0xFF for codes 4 to 7, xor that of a coherent OFDM mapping
    ext, mi = next(((e, i) for i, mk in enumerate(OFDM_MASKS) for e in (0, 1) if x == mk ^ (0xFF * e)), (None, None))
    if ext is None:
        return None
    val = lambda s, n: int("".join(map(str, hb[s:s + n])), 2)
    four = val(28, 2) == 3                       # tones 3: 4 tones per symbol, the tone set in the last two bits
    M = 4 if four else 2
    typ, a, b, c, rate = val(0, 2), val(2, 11), val(13, 11), val(24, 2) + 4 * ext, (RATES4 if four else RATES)[val(26, 2)]
    tones = val(30, 2) if four else val(28, 2)
    ofdm, sc = None, False
    if four and tones == 3:                      # tone set 3: 6400 bit/s (rate field 3), or OFDM (v3, rate field 0 or 1)
        if val(26, 2) == 3:
            tones, rate = TONE_SETS.index("hop"), RATE_WIDE
        elif val(26, 2) < len(OFDM_MODES):
            ofdm, tones = OFDM_MODES[val(26, 2)], TONE_SETS.index("hop")
        elif val(26, 2) == SC_FIELD:             # v5: single carrier
            sc, tones = True, TONE_SETS.index("hop")
        else:
            return None
    if mi and not (ofdm or sc):                  # a mapping mask on anything but OFDM or the single carrier
        return None
    wide = sc and typ == 2                       # v7.4: the wide band, read as type 3
    codes = TYPE3_CODES if typ == 3 or wide else CODES    # v7: type 3, an LDPC code (v7.3: or a polar code)
    code = codes[c] if c < len(codes) else None
    if typ < 2 and (c == 7 or sc and c == 2):    # v7.6: LDPC 5/6, type 0 or 1 with the code 7 (the single carrier's wide band: 2)
        code, wide = LDPC56, sc and c == 2
    omap = OFDM_MAPS[mi]
    if sc:                                       # v7.5: 8PSK's mask = 1024QAM
        omap = next((k for k, v in SC_SLOT.items() if v == omap), omap)
    if code is None or (sc and (omap not in SC_MAPS or code not in SC_CODES)):
        return None
    if tones >= len(TONE_SETS) or (not four and val(30, 2)):
        return None
    nt = SC_NC
    if sc:                                       # v5: the training length in the top bits of a and b
        sc_set_band(wide)                        # v7.4: in symbols of the frame's band
        T, a, b = 2 * (a >> 10) + (b >> 10), a & 1023, b & 1023
        if T == 3:
            return None
        nt = [k for k, v in SC_TCODE.items() if v == T][0]
    if typ == 3 or typ == 2 and wide:            # v7: a text if b = 0, else an image
        typ = 0 if b == 0 else 1
    h = None
    if typ == 0 and a >= 1 and b == 0:
        h = dict(type="text", n0=8 * a, len=a)
    if typ == 1 and 1 <= a <= 512 and 1 <= b <= 512 and a * b <= 70000:   # v4: up to 4x images (the demo: 216 x 320, 69120 pixels; before: 20000 pixels)
        h = dict(type="image", n0=a * b, h=a, w=b)
    if typ == 1 and a >= 513 and 1 <= b <= 512 and 512 * (a - 513) + b <= JPEG_MAX:   # v7.1: a JPEG file of n bytes
        n = 512 * (a - 513) + b
        h = dict(type="jpeg", n0=8 * (n + 4), len=n)
    if typ == 2 and not four and a in SOUND_TAPS and 1 <= b <= 8:
        h = dict(type="sound", n0=0, order=a, periods=b)
    if h and sc:                                 # ns: data symbols
        nc = coded_length(code, h["n0"])
        h.update(code=code, nc=nc, rate=sc_rate(omap, nt, wide), tones="sc", M=None, sc=True, map=omap, ns=-(-nc // OFDM_BITS[omap]), train=nt,
                 wide=wide)
    elif h and ofdm:                             # ns: data symbols after the reference symbol
        nc, omap = coded_length(code, h["n0"]), OFDM_MAPS[mi]
        h.update(code=code, nc=nc, rate=ofdm_rate(ofdm, omap), tones="ofdm", M=None, ofdm=ofdm["name"], map=omap,
                 ns=ofdm_layout(ofdm, omap, nc)[0])
    elif h:
        nc, K = coded_length(code, h["n0"]), int(math.log2(M))
        h.update(code=code, nc=nc, rate=rate, tones=TONE_SETS[tones], M=M, baud=rate // K, ns=-(-nc // K))   # ns: message symbols
    return h


def message(hdr, bits):
    """The message of a frame: its text, its image (rows x columns, True = white) or its JPEG file (v7.1: as it
    arrived, and whether its CRC-32 holds)."""
    if hdr["type"] == "text":
        return dict(text=bytes(np.packbits(bits).tolist()).decode("latin-1"))
    if hdr["type"] == "jpeg":
        f, ok = jpeg_unpack(bits, hdr["len"])
        return dict(file=f, crc_ok=ok)
    return dict(img=bits.reshape(hdr["w"], hdr["h"]).T.astype(bool))


def frame_what(res):
    if res["type"] == "text":
        return repr(res["text"])
    if res["type"] == "jpeg":
        return f"JPEG file of {res['len']} bytes, {'intact (its CRC-32 holds)' if res['crc_ok'] else 'corrupted (its CRC-32 fails)'}"
    return f"{res['w']}x{res['h']} image"


def decode(y, fs, verbose=True, soft="cal", det="smap", effort="fast", est="ls"):
    """Find and decode every frame in the recording y; no settings needed. Returns a list of dicts.
    soft: 'cal' (LLR calibrated on the postamble) or 'ratio' (ln(E1/E0)). det: the single carrier's detector,
    'smap' or 'lmmse'; effort: its steps and passes, 'fast' or 'full' (SC_EFFORT); est: its estimate of the room,
    'ls' or 'mmse' (SC_EST)."""
    y = np.asarray(y, dtype=float)
    y = y - y.mean()
    Lr = fs / CTRL["rate"]
    Lw = round(Lr)
    z0, z1 = tone_energies(y, fs, CTRL["rate"], CTRL["tones"])
    d = (z1 - z0) / (z1 + z0 + 1e-12)
    # one correlator per preamble length, scaled by its threshold 0.3 sqrt(127 / length): the noise on
    # the correlation falls as 1/sqrt(length), so all three have the same false-alarm rate
    corr = {n: pilot_correlation(d, Lr, P) for n, P in PILOTS.items()}
    hcorr = hop_correlation(y, fs, Lr)                       # the same three for a hopping control part
    keys = [(n, False) for n in PILOTS] + [(n, True) for n in PILOTS]
    cs = [(hcorr if hop else corr)[n] for n, hop in keys]
    M = max(len(c) for c in cs)                              # a longer preamble cannot start near the end
    S = np.stack([np.pad(c, (0, M - len(c))) / (0.3 * math.sqrt(127 / n)) for c, (n, _) in zip(cs, keys)])
    score, which = S.max(axis=0), S.argmax(axis=0)
    # v6: the short control part, on the transmitter's grid (c6_detect)
    r48 = FS_TX / fs
    y48 = y if fs == FS_TX else resample(y, r48)
    rho6, bal6 = c6_detect(y48)
    ok6 = (rho6 > C6_THR) & (bal6 > C6_BAL)
    energies = {}                                             # tone energies at each (rate, tones) of the message part
    frames, start, start6 = [], 0, 0                          # start6: where the search for the short control part goes on
    while True:
        above = np.nonzero(score[start:] > 1)[0]
        a6 = np.nonzero(ok6[start6:])[0]
        if len(above) == 0 and len(a6) == 0:
            break
        if len(a6) and (len(above) == 0 or (start6 + a6[0]) / r48 < start + above[0]):   # v6: the short control part first
            c = start6 + a6[0]
            n6 = c + int(np.argmax(rho6[c:c + 2 * C6_N]))       # rho first passes C6_THR up to a period early (2 periods of 3)
            r = c6_read(y48, n6)
            hdr = None if r is None else r["hdr"]
            if hdr is None or not (hdr.get("sc") or hdr.get("ofdm")):
                start6 = n6 + C6_N // 4
                continue
            p1, n_p, hop = int(round(r["t"] / r48)), "short", None   # p1: P's first arrival, in samples of the recording
            t0c, ctl = (r["t"] + C6_LEN) / r48, "the short control part (v6)"   # t0c: where the message part starts
        else:
            cand = start + above[0]
            p1 = cand + int(np.argmax(score[cand:cand + 2 * Lw + 1]))
            n_p, hop = keys[which[p1]]
            P, NH, c1 = PILOTS[n_p], HREPS[n_p] * NHC, (hcorr if hop else corr)[n_p]
            hm = np.round(p1 + np.arange(n_p + NH) * Lr).astype(int)         # preamble and header
            if hm[-1] >= len(z0):                                 # runs past the end: try the next candidate (a shorter preamble)
                start = p1 + Lw
                continue
            if hop:                                               # each tone's level normalized over preamble and header
                (e0, e1), _ = hop_normalize(CTRL["rate"], hop_energies(y, fs, CTRL["rate"], 0, hm, Lw))
                L = np.log((e1[n_p:] + 1e-9) / (e0[n_p:] + 1e-9))
            else:
                L = np.log((z1[hm[n_p:]] + 1e-9) / (z0[hm[n_p:]] + 1e-9))
            L = L.reshape(HREPS[n_p], NHC).sum(axis=0)           # add up the copies
            hdr = parse_header([int(v) for v in decode_fec("conv", HBITS, L)[0]])
            if hdr is None:
                start = p1 + Lw
                continue
            t0c, ctl = p1 + (n_p + NH) * Lr, f"preamble {n_p} bits"
        if hdr["type"] == "sound":                            # room sounding: the MLS, then the postamble in the control mode
            N = 2 ** hdr["order"] - 1
            Delta = (n_p + NH) * Lr + (hdr["periods"] + 1) * N * fs / FS_TX
            tol = round(0.005 * Delta) + Lw
            lo, hi = round(p1 + Delta) - tol, round(p1 + Delta) + tol
            cp = (hcorr if hdr["tones"] == "hop" else corr)[n_p]   # the postamble hops from its first bit, or 4 / 8 kHz
            if hi >= len(cp):
                break
            p2 = lo + int(np.argmax(cp[lo:hi + 1]))
            found = cp[p2] >= 0.5 * c1[p1] and abs((p2 - p1) / Delta - 1) < 1e-3   # clocks are within 1000 ppm
            rho = (p2 - p1) / Delta if found else 1.0
            noise = y[max(0, round(p1 - 1.1 * fs)):max(0, round(p1 - 0.15 * fs))]   # before the warm-up
            room = analyse_sounding(y, fs, p1 + (n_p + NH) * Lr * rho, fs / FS_TX * rho, hdr["order"], hdr["periods"],
                                    noise if len(noise) >= 4096 else None)
            res = dict(hdr, pre=n_p, hop=hop, t=p1 / fs, ppm=(rho * (1 + room["drift"] / N) - 1) * 1e6 if found or hdr["periods"] >= 2 else None,
                       room=room, corr=score, p1=p1, p2=p2)
            frames.append(res)
            if verbose:
                print_room(res)
            start = round(p1 + Delta * rho + n_p * Lr)
            start6 = max(start6, round(start * r48))
            continue
        if hdr.get("sc"):                                     # v5, single carrier: after the control part and a silence
            sc_set_band(hdr["wide"])                          # v7.4: the usual band or the wide one
            t0, dur = t0c, sc_seconds(hdr["ns"], hdr["train"]) * fs
            if t0 + dur * (1 + 1e-3) + 0.05 * fs > len(y):
                break
            o = sc_receive(y, fs, t0, hdr["map"], hdr["code"], hdr["n0"], hdr["nc"], det, nt=hdr["train"], effort=effort, est=est)
            bits = o["bits"]
            res = dict(hdr, pre=n_p, hop=hop, bits=bits, chan_errors=o["chan_errors"], rs_fails=None, t=p1 / fs,
                       ppm=o["ppm"], sc_rx=o, corr=score, p1=p1)
            res.update(message(hdr, bits))
            frames.append(res)
            if verbose:
                what = frame_what(res)
                print(f"frame at {res['t']:.2f} s: {what}, clock offset {o['ppm']:+.0f} ppm (from the training blocks), "
                      f"code {hdr['code']} corrected {o['chan_errors']} of {hdr['nc']} channel bits")
                print(f"  from the header: {sc_desc(hdr['map'], hdr['train'], hdr['wide'])}, {hdr['ns']} data symbols in {o['nseg']} segment{'s' if o['nseg'] > 1 else ''}, {ctl}")
                ch = ", ".join(f"{c:.0f}" for c in o["change"])
                print(f"  SNR {o['snr']:.1f} dB, the median of the training blocks ({' / '.join(f'{v:.1f}' for v in o['snrs'])} dB: each block's room, {o['taps']} taps, against what it leaves"
                      f"{', MMSE' if est == 'mmse' else ''}), "
                      f"room change at each later block {ch} dB, sampling phase {o['phase']}")
                print(f"  {'SMAP' if det == 'smap' else 'linear MMSE'} ({effort}) and the {'LDPC' if hdr['code'] in LDPC_R else 'polar' if hdr['code'] in POLAR_R else 'turbo'} decoder: {o['passes']} pass{'es' if o['passes'] > 1 else ''}, "
                      f"{'a codeword' if o['codeword'] else 'not a codeword (some bits are wrong)'}, {o['sec']:.1f} s"
                      + (f"; the clock wanders by {o['bend']:.1f} samples (samples lost or repeated?)" if o["bend"] > 1 else ""))
                if o.get("relearn"):
                    r2 = o["relearn"]
                    print(f"  second try (v7.2): the room learned again from the frame, {r2['taps']} taps (SNR {r2['snr']:.1f} dB), "
                          f"{r2['passes']} more pass{'es' if r2['passes'] > 1 else ''}: {'a codeword' if o['codeword'] else 'still not a codeword'}")
                for tg, sz in o["gaps"]:
                    print(f"  {abs(sz):.0f} samples {'lost' if sz < 0 else 'repeated'} by the microphone {tg:.2f} s into the part, "
                          f"{o['erased']} symbols erased (v5.1: read again after the gap)")
            start = round(t0 + dur)
            start6 = max(start6, round(start * r48))
            continue
        if hdr.get("ofdm"):                                   # OFDM: the reference symbol follows the control part
            mode = ofdm_mode(hdr["ofdm"])
            t0, Ls = t0c, (mode["N"] + mode["CP"]) * fs / FS_TX
            if t0 + (hdr["ns"] + 2) * Ls * (1 + 1e-3) > len(y):
                break
            o = ofdm_receive(y, fs, t0, mode, hdr["ns"], hdr["nc"], hdr["map"])
            llr, key = o["llr"], prbs(hdr["nc"])
            bits, rs_fails = decode_fec(hdr["code"], hdr["n0"], np.where(key == 1, -llr, llr), np.inf)
            sent = encode_fec(hdr["code"], bits) ^ key
            chan_errors = int(np.sum((llr > 0).astype(int) != sent))
            res = dict(hdr, pre=n_p, hop=hop, bits=bits, llr=llr, chan_errors=chan_errors, rs_fails=rs_fails, t=p1 / fs,
                       ppm=o["eps"], ofdm_rx=o, corr=score, p1=p1)
            res.update(message(hdr, bits))
            frames.append(res)
            if verbose:
                what = frame_what(res)
                coded = "" if hdr["code"] == "none" else f", code {hdr['code']} corrected {chan_errors} of {hdr['nc']} channel bits"
                print(f"frame at {res['t']:.2f} s: {what}, clock offset {o['eps']:+.0f} ppm (from the data){coded}")
                print(f"  from the header: {ofdm_desc(mode, hdr['map'])}, {hdr['ns']} data symbol{'s' if hdr['ns'] > 1 else ''}, {ctl}")
                print(f"  SNR per subcarrier {np.min(o['snr']):.0f} to {np.max(o['snr']):.0f} dB (median {np.median(o['snr']):.0f}, "
                      f"groups of {OFDM_GROUP}), {100 * (1 - o['inside']):.0f}% of the received energy after the guard, "
                      f"windows {o['tau']:+d} samples from the preamble's timing")
                for i, d in o["slips"]:
                    print(f"  {abs(d):.0f} samples {'lost' if d < 0 else 'repeated'} by the microphone before symbol {i} of {hdr['ns']}"
                          f"{', symbol ' + ', '.join(map(str, o['erased'])) + ' erased' if o['erased'] else ''} (v4.1: read again after the gap)")
            start = round(t0 + (hdr["ns"] + 2) * Ls)
            start6 = max(start6, round(start * r48))
            continue
        n0, nc, rate, tones, M, ns = hdr["n0"], hdr["nc"], hdr["rate"], hdr["tones"], hdr["M"], hdr["ns"]
        S = P if M == 2 else to_symbols(P, M)                 # the postamble's symbols
        NS = len(S)
        Ld = fs / hdr["baud"]                                 # one symbol
        Delta = (n_p + NH) * Lr + ns * Ld                     # nominal time from preamble to postamble
        tol = round(MAX_PPM * 1e-6 * Delta + Ld / 4)
        off = np.round(np.arange(NS) * Ld).astype(int)
        lo, hi = round(p1 + Delta) - tol, round(p1 + Delta) + tol
        if hi + off[-1] + round(Ld) > len(y):
            break
        if tones == "hop":                                    # tone hopping: energies of each symbol's pair (group)
            c2 = hop_postamble(y, fs, rate, ns, S, lo, hi, Ld, M)
        else:
            if (rate, tones, M) not in energies:
                Ef = np.array(tone_energies(y, fs, rate, tones, M))
                energies[rate, tones, M] = Ef, (Ef[1] - Ef[0]) / (Ef[1] + Ef[0] + 1e-12) if M == 2 else None
            Ef, dd = energies[rate, tones, M]
            if M == 2:
                c2 = dd[np.arange(lo, hi + 1)[:, None] + off] @ (2 * P - 1) / n_p
            else:
                c2 = sum(known_contrast(Ef[:, lo + o:hi + 1 + o], v) for o, v in zip(off, S)) / NS
        p2 = lo + int(np.argmax(c2))
        # At 6400 bit/s every symbol uses the same 4 tones, so the early echoes (a few ms, many symbols) cap the
        # postamble's contrast near 0.4 even in a quiet room: a fifth of the preamble's score is enough there (its
        # side lobes, a symbol or more away, stay near 0.02).
        found = c2.max() >= (0.2 if rate == RATE_WIDE else 0.5) * c1[p1]
        rho = (p2 - p1) / Delta if found else 1.0
        mk = np.round(p1 + ((n_p + NH) * Lr + np.arange(ns) * Ld) * rho).astype(int)
        post = np.round((p2 if found else p1 + Delta) + np.arange(NS) * Ld * rho).astype(int)
        shift = np.zeros(ns, dtype=int)                       # fine timing on the message, with hopping (see fine_timing)
        if tones == "hop":
            read = lambda s: hop_energies(y, fs, rate, 0, np.clip(mk + s, 0, len(y) - round(Ld)), round(Ld), M)
            shift = fine_timing(read, ns, min(int(Ld // 2), round(5e-4 * fs)))
        mk = mk + shift
        spread = None
        if tones == "hop":                                    # each tone levelled over the message and the postamble
            A, spread = hop_normalize(rate, hop_energies(y, fs, rate, 0, np.concatenate([mk, post]), round(Ld), M))
            mE, pE = A[:, :ns], A[:, ns:]
        else:
            mE, pE = Ef[:, mk], Ef[:, post]
        # the postamble is at the message's tones and rate; without it, the message itself calibrates the soft output
        cal = calibrate(pE, S) if found else calibrate_blind(mE)
        g = 1.0 if soft == "ratio" else level_track(cal, pE if found else mE, mE)
        llr = soft_llr(soft, cal, mE, g).ravel()[:nc]
        key = prbs(nc)
        bits, rs_fails = decode_fec(hdr["code"], n0, np.where(key == 1, -llr, llr), 4 if soft == "ratio" else np.inf)
        sent = encode_fec(hdr["code"], bits) ^ key            # re-encoded decision
        chan_errors = int(np.sum((llr > 0).astype(int) != sent))
        res = dict(hdr, pre=n_p, hop=hop, bits=bits, llr=llr, chan_errors=chan_errors, rs_fails=rs_fails, cal=cal, t=p1 / fs,
                   ppm=(rho - 1) * 1e6 if found else None, spread=spread,
                   pilot_errors=post_errors(pE, S) if found else None,
                   corr=score, p1=p1, p2=p2, timing=(int(shift[0]), int(shift[-1])) if shift.any() else None)
        res.update(message(hdr, bits))
        frames.append(res)
        if verbose:
            ppm = "postamble not found (soft output calibrated on the message)" if res["ppm"] is None else f"clock offset {res['ppm']:+.0f} ppm"
            what = frame_what(res)
            coded = "" if hdr["code"] == "none" else f", code {hdr['code']} corrected {chan_errors} of {nc} channel bits"
            if rs_fails:
                coded += f" ({rs_fails} Reed-Solomon blocks failed)"
            print(f"frame at {res['t']:.2f} s: {what}, {ppm}, postamble errors {res['pilot_errors']}{coded}")
            four = f"4 tones ({hdr['baud']} symbols/s), " if M == 4 else ""
            print(f"  from the header: {rate} bit/s, {four}{tone_desc(tones, rate, M)}, preamble {n_p} bits")
            snr = " / ".join(f"{t['snr_db']:.0f}" for t in cal)
            if tones == "hop":
                print(f"  SNR per symbol {snr} dB ({'bit-0 / bit-1 tones' if M == 2 else 'tones 1 to 4 of each group'}, levelled), "
                      f"tone levels spread over {spread:.0f} dB, soft output {soft}")
            else:
                print(f"  SNR per symbol {snr} dB ({'low / high tone' if M == 2 else 'tones from low to high'}), "
                      f"{'high' if M == 2 else 'highest'} tone {20 * np.log10(cal[-1]['mu'] / cal[0]['mu']):+.1f} dB relative to "
                      f"{'low' if M == 2 else 'lowest'} tone, soft output {soft}")
            if "alpha" in cal[-1]:
                print(f"  calibrated LLR scaled by {cal[-1]['alpha']:.2f} (fitted on the postamble)")
            if res["timing"]:
                print(f"  symbol timing fine-tuned on the message: {res['timing'][0]:+d} -> {res['timing'][1]:+d} samples (start -> end)")
        start = round(p1 + Delta + n_p * Ld)
        start6 = max(start6, round(start * r48))
    return frames


# ---------------------------------------------------------------- room sounding
def _lanczos_at(y, t, a=6):
    """Band-limited interpolation of y at the (fractional) sample positions t."""
    i0 = np.floor(t).astype(int)
    v = np.zeros(len(t))
    for k in range(-a + 1, a + 1):
        i = i0 + k
        x = t - i
        w = np.sinc(x) * np.sinc(x / a) * (np.abs(x) < a)
        ok = (i >= 0) & (i < len(y))
        v[ok] += w[ok] * y[i[ok]]
    return v


def _oct_smooth(p, df, grid, frac):
    """Average the spectrum p (bins of df Hz) over 1/frac octave around each frequency of the grid."""
    out = []
    for f in grid:
        lo = max(1, int(np.floor(f * 2 ** (-0.5 / frac) / df)))
        hi = min(len(p) - 1, int(np.ceil(f * 2 ** (0.5 / frac) / df)))
        out.append(p[lo:hi + 1].mean())
    return np.array(out)


def _welch(x, nfft=2048):
    """Power spectral density (arbitrary scale), Hann window, 50 % overlap."""
    w = 0.5 - 0.5 * np.cos(2 * np.pi * np.arange(nfft) / nfft)
    segs = [x[o:o + nfft] * w for o in range(0, len(x) - nfft + 1, nfft // 2)]
    return np.mean(np.abs(np.fft.rfft(segs, axis=1)) ** 2, axis=0) if segs else None


def analyse_sounding(y, fs, first, step, order, periods, noise=None):
    """Impulse response of the chain speaker + room + microphone from the received MLS, and the room measures.
    first: position in y of the first chip; step: samples of y per chip (fs / 48000 corrected for the clock offset);
    noise: recording without signal (before the frame), for the SNR per frequency."""
    s = 2.0 * mls(order, SOUND_TAPS[order]) - 1
    N = len(s)
    S = np.conj(np.fft.fft(s))
    # the MLS autocorrelation is N at lag 0 and -1 elsewhere, so h = (circular xcorr + sum r) / (N + 1);
    # the band above 20 kHz only carries interpolation errors: circular low-pass (Blackman-windowed sinc)
    K = 48
    u = np.arange(-K, K + 1)
    fc = min(20000, 0.45 * fs) / FS_TX
    lp = 2 * fc * np.sinc(2 * fc * u) * (0.42 + 0.5 * np.cos(np.pi * u / K) + 0.08 * np.cos(2 * np.pi * u / K))
    lpN = np.zeros(N)
    lpN[u % N] = lp / lp.sum()
    LP = np.fft.fft(lpN)
    period = lambda j, st, o=0.0: _lanczos_at(y, first + (j * N + np.arange(N) + o) * st)
    ir = lambda r, filt=False: (np.real(np.fft.ifft(np.fft.fft(r) * S * (LP if filt else 1))) + r.sum()) / (N + 1)

    def peak(h):
        m = int(np.argmax(np.abs(h)))
        a, b, c = np.abs(h[(m - 1) % N]), np.abs(h[m]), np.abs(h[(m + 1) % N])
        den = a - 2 * b + c
        return m + (0.5 * (a - c) / den if den else 0)
    # the direct-path peak of each period: its drift from one period to the next measures the residual clock offset;
    # a jump (samples lost or repeated by the audio system) is absorbed by realigning each period on its own peak
    drift, glitch, use, off = 0.0, False, list(range(1, periods + 1)), np.zeros(periods + 1)
    if periods >= 2:
        pk = np.array([peak(ir(period(j, step))) for j in range(1, periods + 1)])
        pk = pk[0] + (pk - pk[0] + N / 2) % N - N / 2
        jj = np.arange(len(pk)) - (len(pk) - 1) / 2
        drift = float(np.sum(jj * (pk - pk.mean())) / np.sum(jj ** 2))      # chips per period
        if abs(drift) > 50 or np.any(np.abs(pk - pk.mean() - jj * drift) > 1):
            glitch = True
            dd = [v for v in np.diff(pk) if abs(v) <= 50]
            drift = float(np.median(dd)) if len(dd) else 0.0
            off[1:] = pk - pk[0] - np.arange(periods) * drift
    st = step * (1 + drift / N)
    if glitch:                     # the period cut by the jump holds two partial responses: drop the periods far from the median response
        hs = np.array([ir(period(j, st, off[j])) for j in use])
        dev = np.sum((hs - np.median(hs, axis=0)) ** 2, axis=1)
        if len(use) >= 3:
            use = [j for j, v in zip(use, dev) if v <= 3 * np.sort(dev)[len(dev) // 2]]
        else:
            pe = [np.sum(np.roll(h, 48 - round(peak(h)))[:97] ** 2) for h in hs]
            use = [use[int(np.argmax(pe))]]
        # the clock drift again, from consecutive periods that are both kept and have no jump between them
        dd = [pk[j] - pk[j - 1] for j in use if j + 1 in use and abs(pk[j] - pk[j - 1] - drift) < 5]
        if dd:
            drift = float(np.mean(dd))
            off[1:] = pk - pk[0] - np.arange(periods) * drift
            st = step * (1 + drift / N)
    r = np.mean([period(j, st, off[j]) for j in use], axis=0)
    h0 = ir(r, True)
    pre = round(0.005 * FS_TX)
    h = np.roll(h0, pre - round(peak(h0)))                     # direct sound 5 ms after the start
    e = h ** 2
    n9 = round(0.9 * N)
    floor = e[n9:].mean()                                      # noise floor: the end of the period, just before the direct sound
    W = round(0.002 * FS_TX)
    cs = np.concatenate([[0], np.cumsum(e)])
    i = np.arange(N)
    env = (cs[np.minimum(N, i + W + 1)] - cs[np.maximum(0, i - W)]) / (np.minimum(N, i + W + 1) - np.maximum(0, i - W))
    below = np.nonzero(env[pre + W:n9] < 2 * floor)[0]
    iend = pre + W + int(below[0]) if len(below) else n9
    inr = 10 * np.log10(env[pre] / floor)                      # decay range above the noise: T20 needs ~35 dB, T10 ~25 dB
    # Schroeder backward integration (noise floor removed), and the reverberation time from its slope; the decay
    # hidden under the noise after iend is added back from the fitted slope (as in Lundeby's method), then refitted
    ex = np.maximum(0, e[:iend] - floor)

    def schroeder(comp):
        edc = np.concatenate([np.cumsum(ex[::-1])[::-1], [0]]) + comp
        return edc, 10 * np.log10(edc / edc[0] + 1e-30)

    def fit(edc_db, hi_db, lo_db):
        idx = np.nonzero((np.arange(len(edc_db)) > pre) & (edc_db <= hi_db) & (edc_db >= lo_db))[0]
        if len(idx) < 10:
            return None
        sl, at = np.polyfit(idx, edc_db[idx], 1)               # dB per sample
        return dict(rt60=-60 / sl / FS_TX, slope=sl, at=at) if sl < 0 else None

    def reverb(edc_db):
        f20 = fit(edc_db, -5, -25) if inr >= 35 else None
        f10 = fit(edc_db, -5, -15) if f20 is None and inr >= 25 else None
        return dict(f20, method="T20") if f20 else dict(f10, method="T10") if f10 else None
    comp = 0.0
    edc, edc_db = schroeder(comp)
    rt = reverb(edc_db)
    for _ in range(2):
        if rt is None:
            break
        comp = edc[0] * 10 ** ((rt["at"] + rt["slope"] * iend) / 10)
        edc, edc_db = schroeder(comp)
        rt = reverb(edc_db)
    tot = lambda a, b: cs[min(N, max(0, b))] - cs[min(N, max(0, a))]
    d25, d50 = round(0.0025 * FS_TX), round(0.05 * FS_TX)
    direct = tot(pre - d25, pre + d25)
    late = max(1e-30, tot(pre + d25, iend) - floor * (iend - pre - d25) + comp)
    drr = 10 * np.log10(direct / late)
    c50 = 10 * np.log10(max(1e-30, tot(pre - d25, pre + d50) - floor * (d50 + d25))
                        / max(1e-30, tot(pre + d50, iend) - floor * (iend - pre - d50) + comp))
    v = np.maximum(0, e[pre - d25:iend] - floor)
    t = (np.arange(pre - d25, iend) - pre) / FS_TX
    m0, m1, m2 = v.sum(), (v * t).sum(), (v * t * t).sum()
    if rt and comp:                                            # moments of the exponential decay after iend
        ts, te = -10 / (rt["slope"] * np.log(10)) / FS_TX, (iend - pre) / FS_TX
        m0, m1, m2 = m0 + comp, m1 + comp * (te + ts), m2 + comp * (te * te + 2 * te * ts + 2 * ts * ts)
    tau_rms = float(np.sqrt(max(0, m2 / m0 - (m1 / m0) ** 2)))
    # frequency response of the whole chain, 1/6 octave
    nf = 1 << int(np.ceil(np.log2(iend + 1)))
    pw = np.abs(np.fft.rfft(h[:iend + 1], nf)) ** 2
    fmax = min(20000, 0.45 * fs)
    grid = 100 * 2 ** (np.arange(int(np.log2(fmax / 100) * 24 + 1e-9) + 1) / 24)
    H = _oct_smooth(pw, FS_TX / nf, grid, 6)
    H_db = 10 * np.log10(H / H.max() + 1e-30)
    snr = None                                                 # received sounding against the noise heard before the frame
    if noise is not None and len(noise) >= 4096:
        a = round(first + N * step)
        ps, pn = _welch(y[a:a + round(periods * N * step)]), _welch(np.asarray(noise, dtype=float))
        if ps is not None and pn is not None:
            ps, pn = _oct_smooth(ps, fs / 2048, grid, 3), _oct_smooth(pn, fs / 2048, grid, 3)
            snr = 10 * np.log10(np.maximum(1e-12, ps - pn) / pn)
    at = lambda arr, f: None if arr is None or f > fmax else float(arr[np.argmin(np.abs(np.log(grid / f)))])
    tones = [dict(f=f, gain=at(H_db, f), snr=at(snr, f)) for f in (4000, 8000, 16000)]
    return dict(fs_ir=FS_TX, h=h, pre=pre, iend=iend, inr=inr, env=env, edc_db=edc_db, rt=rt, drr=drr, c50=c50,
                tau_rms=tau_rms, grid=grid, H_db=H_db, snr=snr, tones=tones, drift=drift, glitch=glitch, used=len(use), N=N)


# Es/N0 that a tone of amplitude 0.9 at 400 bit/s gets, relative to the sounding's SNR per Hz (MLS spread over 0 .. 24 kHz)
ESN0_400 = 10 * np.log10(0.405 / SOUND_AMP ** 2 * (FS_TX / 2) / 400)


def print_room(res):
    s = res["room"]
    ppm = "postamble not found" if res["ppm"] is None else f"clock offset {res['ppm']:+.0f} ppm"
    print(f"room sounding at {res['t']:.2f} s: MLS of {s['N']} samples ({s['N'] / FS_TX:.2f} s) x {res['periods']} periods, "
          f"preamble {res['pre']}{', hopping' if res['tones'] == 'hop' else ''}, {ppm}")
    if s["glitch"]:
        print(f"  the audio skipped or repeated samples: periods realigned, {s['used']} of {res['periods']} used")
    rt = "not measurable (the decay is too short above the noise)" if s["rt"] is None else f"{s['rt']['rt60']:.2f} s ({s['rt']['method']})"
    print(f"  reverberation time RT60 {rt}, decay range {s['inr']:.0f} dB above the noise")
    print(f"  direct-to-reverberant ratio {s['drr']:+.1f} dB, clarity C50 {s['c50']:+.1f} dB, rms delay spread {1000 * s['tau_rms']:.1f} ms")
    g = {t["f"]: t for t in s["tones"]}
    line = f"  gain of 8 vs 4 kHz {g[8000]['gain'] - g[4000]['gain']:+.1f} dB"
    if g[16000]["gain"] is not None:
        line += f", 16 vs 8 kHz {g[16000]['gain'] - g[8000]['gain']:+.1f} dB"
    print(line + " (speaker, room, microphone)")
    if s["snr"] is not None:
        tt = [t for t in s["tones"] if t["snr"] is not None]
        print(f"  Es/N0 of a tone at 400 bit/s (noise only): " + " / ".join(f"{t['snr'] + ESN0_400:.0f}" for t in tt)
              + " dB at " + " / ".join(f"{t['f'] // 1000}" for t in tt) + " kHz")


def plot_room(res):
    import matplotlib.pyplot as plt
    s = res["room"]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8, 7))
    n = min(s["N"], max(round(0.1 * FS_TX), round(1.25 * s["iend"])))
    t = (np.arange(n) - s["pre"]) / FS_TX * 1000
    a1.plot(t, 10 * np.log10(s["env"][:n] / s["env"].max() + 1e-12), lw=0.8, label="energy of h(t), 4 ms average")
    a1.plot(t[:len(s["edc_db"])], s["edc_db"][:n], lw=2, label="energy decay curve (Schroeder)")
    if s["rt"]:
        i = np.array([s["pre"], min(n, (-65 - s["rt"]["at"]) / s["rt"]["slope"])])
        a1.plot((i - s["pre"]) / FS_TX * 1000, s["rt"]["at"] + s["rt"]["slope"] * i, "k--", label=f"RT60 = {s['rt']['rt60']:.2f} s")
    a1.set(ylim=(-70, 3), xlabel="time after the direct sound (ms)", ylabel="dB")
    a1.legend(); a1.grid(alpha=0.3)
    a2.semilogx(s["grid"], s["H_db"], lw=2, label="frequency response (1/6 octave)")
    if s["snr"] is not None:
        a2.semilogx(s["grid"], s["snr"], lw=2, label="SNR of the sounding")
    for f in (4000, 8000, 16000):
        a2.axvline(f, color="k", ls=":", lw=1)
    a2.set(xlabel="frequency (Hz)", ylabel="dB", ylim=(-40, None))
    a2.legend(); a2.grid(alpha=0.3, which="both")
    fig.tight_layout()
    plt.show()


# ---------------------------------------------------------------- room simulator
def resample(x, ratio, A=24):
    """Windowed-sinc resampling: output sample i is x at time i / ratio (ratio = output rate / input rate)."""
    x = np.asarray(x, dtype=float)
    fc = 0.95 * min(1.0, ratio)
    n = int(len(x) * ratio)
    y = np.zeros(n)
    half = int(np.ceil(A / fc))
    for c in range(0, n, 65536):
        t = np.arange(c, min(n, c + 65536)) / ratio
        k0 = np.floor(t).astype(int)
        for k in range(-half + 1, half + 1):
            i = k0 + k
            u = (t - i) * fc
            w = fc * np.sinc(u) * (0.5 + 0.5 * np.cos(np.pi * u / A)) * (np.abs(u) < A)
            ok = (i >= 0) & (i < len(x))
            y[c:c + len(t)][ok] += w[ok] * x[i[ok]]
    return y


def synthetic_room(rt60, drr, fs=FS_TX, seed=1):
    """Direct sound + diffuse tail: Gaussian noise decaying by 60 dB in rt60 seconds, scaled for the direct-to-
    reverberant ratio drr (dB, direct = within 2.5 ms of the direct sound)."""
    rng = np.random.default_rng(seed)
    d0, n = round(0.005 * fs), round((0.005 + 1.5 * rt60) * fs)
    h = np.zeros(n)
    h[d0] = 1.0
    t = (np.arange(n) - d0) / fs
    tail = rng.standard_normal(n) * np.exp(-6.9078 * t / rt60) * np.clip((t - 0.001) / 0.005, 0, 1)
    h += tail * np.sqrt(10 ** (-drr / 10) / np.sum(tail ** 2))
    return h


def simulate(x, fs_x, h, fs_h, snr_db, ppm=0.0, fs_out=44100, seed=1, lead=1.0):
    """What a microphone at fs_out hears: x through the impulse response h, its clock offset by ppm, white noise
    at snr_db below the received signal, and lead seconds of noise alone first."""
    rng = np.random.default_rng(seed + 1)
    if fs_h != fs_x:
        h = resample(h, fs_x / fs_h) * fs_h / fs_x
    L = 1 << int(np.ceil(np.log2(len(x) + len(h))))
    y = np.fft.irfft(np.fft.rfft(x, L) * np.fft.rfft(h, L), L)[:len(x) + len(h) - 1]
    y = np.concatenate([np.zeros(round(lead * fs_x)), y])
    y = resample(y, fs_out / fs_x * (1 + ppm * 1e-6))
    on = np.abs(y) > 1e-3 * np.abs(y).max()
    p = np.mean(y[on] ** 2)
    return y + rng.standard_normal(len(y)) * np.sqrt(p / 10 ** (snr_db / 10))


# ---------------------------------------------------------------- command line
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("tx", "rx", "sim"):
        p = sub.add_parser(name)
        if name == "tx":
            p.add_argument("--text")
            p.add_argument("--rate", type=int, choices=sorted(set(RATES + RATES4 + [RATE_WIDE])), default=400,
                           help="200 with 2 tones only, 3200 and 6400 with 4 tones only (6400 with --tones hop: one group of 4 tones over the whole band)")
            p.add_argument("--mod", choices=("2", "4", "ofdm", "sc"), default="2",
                           help="FSK with 2 tones (1 bit per symbol) or 4 (2 bits), OFDM, or a single carrier (v5)")
            p.add_argument("--guard", choices=[m["name"] for m in OFDM_MODES], default="long",
                           help="OFDM: long (N 2048, guard 21 ms) or short (N 1024, guard 5 ms)")
            p.add_argument("--map", choices=OFDM_MAPS + ["1024qam"],
                           help="OFDM: DQPSK (no pilots, the default), or QPSK, 8PSK, 16QAM, 64QAM, 256QAM with pilots (2, 3, 4, 6, 8 bits per subcarrier); "
                                "single carrier: QPSK, 16QAM, 64QAM, 256QAM (the default) or 1024QAM (v7.5, phones close together in a quiet room)")
            p.add_argument("--tones", choices=TONE_SETS, default="low", help="hop: tone hopping, for reverberant rooms")
            p.add_argument("--train", choices=list(SC_TRAIN), default="medium",
                           help="single carrier: training of 640, 1280 or 2560 symbols (wide band: 896, 1792, 3584), echoes up to 20, 47 or 100 ms (as OFDM's guard)")
            p.add_argument("--wide", action="store_true",
                           help="v7.4, single carrier: the wide band, 13 714 symbols/s from 2.3 to 20.5 kHz (an LDPC or a polar code; receivers before v7.4 ignore it)")
            p.add_argument("--pre", type=int, choices=sorted(PILOTS), default=127, help="preamble length, longer reaches farther")
            p.add_argument("--control", choices=("short", "long"), default="short",
                           help="OFDM and single carrier: the short control part (v6, 0.19 s; v5 receivers do not read it) "
                                "or the 400 bit/s one (0.65 s with --pre 127; it reaches farther)")
            p.add_argument("--image")
            p.add_argument("--step", type=int, default=1, help="image: every step-th pixel")
            p.add_argument("--up", type=int, choices=(1, 2, 4), default=1, help="image: up times the width and height (smooth edges)")
            p.add_argument("--jpeg", help="v7.1: a JPEG file, sent as it is (see --quality to make one with restart markers)")
            p.add_argument("--quality", type=int, help="v7.1: send the --image as a JPEG of this quality (1 to 100), in colour, "
                           "a restart marker on every row of blocks, saved next to the WAV")
            p.add_argument("--width", type=int, help="with --quality: scale the picture down to this width")
            p.add_argument("--grey", action="store_true", help="with --quality: a grey JPEG")
            p.add_argument("--code", choices=[c for c in CODES if c not in (None, "rs")] + LDPC_CODES + [LDPC56] + POLAR_CODES, default="none",
                           help="turbo, turbo3, turbo4: turbo code of rate 1/2, 1/3, 1/4; ldpc4, ldpc3, ldpc, ldpc23, ldpc34, ldpc56: "
                                "LDPC code (v7) of rate 1/4, 1/3, 1/2, 2/3, 3/4 (v7.6: 5/6); polar3, polar, polar23: polar code (v7.3) of rate 1/3, 1/2, 2/3, "
                                "with an 11-bit CRC per block and list decoding; Reed-Solomon (v0) is still decoded but no longer sent")
            p.add_argument("--sound", action="store_true", help="room sounding frame instead of a message")
            p.add_argument("--order", type=int, choices=sorted(SOUND_TAPS), default=15, help="sounding period 2^order - 1 samples at 48 kHz")
            p.add_argument("--periods", type=int, choices=range(1, 9), default=4, help="sounding periods averaged")
            p.add_argument("-o", "--out", default="tx.wav")
            p.add_argument("--play", action="store_true")
        elif name == "sim":
            p.add_argument("wav", help="what the speaker plays (e.g. from tx)")
            p.add_argument("--ir", help="impulse response WAV (e.g. saved by rx from a sounding); default: a model room")
            p.add_argument("--rt60", type=float, default=0.6, help="model room: reverberation time in s")
            p.add_argument("--drr", type=float, default=0.0, help="model room: direct-to-reverberant ratio in dB")
            p.add_argument("--snr", type=float, default=30.0, help="signal-to-noise ratio at the microphone in dB")
            p.add_argument("--ppm", type=float, default=0.0, help="clock offset of the microphone")
            p.add_argument("--fs", type=int, default=44100, help="sample rate of the microphone")
            p.add_argument("--seed", type=int, default=1)
            p.add_argument("-o", "--out", default="heard.wav")
        else:
            p.add_argument("wav", nargs="?")
            p.add_argument("--listen", type=float, help="record this many seconds from the microphone")
            p.add_argument("--fs", type=int, default=48000)
            p.add_argument("--soft", choices=("cal", "ratio"), default="cal")
            p.add_argument("--det", choices=("smap", "lmmse"), default="smap", help="single carrier: soft SMAP or linear MMSE detector")
            p.add_argument("--effort", choices=list(SC_EFFORT), default="fast",
                           help="single carrier: fast (SMAP 30 steps, at most 4 passes) or full (150 steps, 8 passes)")
            p.add_argument("--est", choices=SC_EST, default="ls",
                           help="single carrier: the room estimate, least squares or MMSE (prior power per tap from a pilot)")
            p.add_argument("--ir", default="room_ir.wav", help="where to save the impulse response of a room sounding")
            p.add_argument("--plot", action="store_true", help="plot the room sounding (needs matplotlib)")
            p.add_argument("--jpg", default="received", help="v7.1: JPEG files received are saved as JPG_1.jpg, JPG_2.jpg...")
    a = ap.parse_args()
    if a.cmd == "sim":
        x, fs_x = read_wav(a.wav)
        if a.ir:
            h, fs_h = read_wav(a.ir)
        else:
            h, fs_h = synthetic_room(a.rt60, a.drr, fs_x, a.seed), fs_x
        y = simulate(x, fs_x, h, fs_h, a.snr, a.ppm, a.fs, a.seed)
        write_wav(a.out, y / max(1.0, np.abs(y).max()), a.fs)
        room = a.ir or f"model room RT60 {a.rt60} s, DRR {a.drr:+.1f} dB"
        print(f"{a.wav} through {room}, SNR {a.snr:.0f} dB, {a.ppm:+.0f} ppm, {a.fs} Hz -> {a.out}")
        return
    if a.cmd == "tx":
        jpeg = None
        if a.jpeg:                                             # v7.1: a JPEG file as it is
            jpeg = open(a.jpeg, "rb").read()
            if b"\xff\xdd" not in jpeg[:jpeg.find(b"\xff\xda") + 1]:
                print("note: this JPEG has no restart markers, so a wrong bit garbles the rest of the picture (see --quality)")
        elif a.image and a.quality:                            # v7.1: the picture made into a JPEG
            jpeg = jpeg_encode(load_rgb(a.image, a.width), a.quality, a.grey)
            with open(os.path.splitext(a.out)[0] + ".jpg", "wb") as f:
                f.write(jpeg)
        img = load_image(a.image, a.step, a.up) if a.image and jpeg is None else None
        sound = (a.order, a.periods) if a.sound else None
        mod = a.mod if a.mod in ("ofdm", "sc") else int(a.mod)
        omap = a.map or ("256qam" if mod == "sc" else "dqpsk")
        code = ("ldpc" if a.wide else "turbo") if mod == "sc" and a.code == "none" else a.code   # single carrier: turbo, LDPC or polar codes (wide: LDPC or polar)
        segs, data = build_frame(text=a.text if img is None and jpeg is None else None, img=img, code=code, rate=a.rate, tones=a.tones, pre=a.pre, sound=sound,
                                 mod=mod, guard=a.guard, omap=omap, train=a.train, control=a.control, jpeg=jpeg, wide=a.wide and mod == "sc")
        x = modulate(segs)
        write_wav(a.out, x, FS_TX)
        if sound:
            print(f"room sounding: MLS of {2 ** a.order - 1} samples x {a.periods + 1}, preamble {a.pre}, {len(x) / FS_TX:.1f} s -> {a.out}")
        else:
            how = ofdm_desc(ofdm_mode(a.guard), omap) if mod == "ofdm" else sc_desc(omap, SC_SIZES[a.wide][2][a.train], a.wide) if mod == "sc" else f"{mod} tones"
            if mod in ("ofdm", "sc"):
                how += ", the short control part (v6)" if a.control == "short" else f", preamble {a.pre}"
            what = f"a JPEG file of {len(jpeg)} bytes ({os.path.splitext(a.out)[0]}.jpg) and its CRC-32: " if jpeg is not None and not a.jpeg else ""
            print(f"{what}{len(data)} message bits, code {code}, {how}, {len(x) / FS_TX:.1f} s -> {a.out}")
        if a.play:
            import sounddevice as sd
            sd.play(x, FS_TX, blocking=True)
    else:
        if a.listen:
            import sounddevice as sd
            print(f"listening for {a.listen} s...")
            y, fs = sd.rec(int(a.listen * a.fs), samplerate=a.fs, channels=1, blocking=True)[:, 0], a.fs
            write_wav("recording.wav", y, fs)
        else:
            y, fs = read_wav(a.wav)
        frames = decode(y, fs, soft=a.soft, det=a.det, effort=a.effort, est=a.est)
        if not frames:
            print("no frame found")
        for f in frames:
            if f["type"] == "sound":
                s = f["room"]
                n = min(s["N"], s["iend"] + round(0.05 * FS_TX))
                write_wav(a.ir, 0.9 * s["h"][:n] / np.abs(s["h"]).max(), FS_TX)
                print(f"impulse response ({n / FS_TX:.2f} s at {FS_TX} Hz) -> {a.ir}")
                if a.plot:
                    plot_room(f)
            elif f["type"] == "text":
                print("decoded:", f["text"])
            elif f["type"] == "jpeg":                          # v7.1: the file as it arrived; corrupted, also read again stripe by stripe
                k = sum(g["type"] == "jpeg" for g in frames[:frames.index(f) + 1])
                with open(f"{a.jpg}_{k}.jpg", "wb") as out:
                    out.write(f["file"])
                more = ""
                if not f["crc_ok"]:
                    with open(f"{a.jpg}_{k}_resync.jpg", "wb") as out:
                        out.write(jpeg_resync(f["file"]))
                    more = f"; {a.jpg}_{k}_resync.jpg: the same without the markers that wrong bits made (see jpeg_resync)"
                print(f"JPEG file {'intact' if f['crc_ok'] else 'CORRUPTED'} -> {a.jpg}_{k}.jpg{more}")
            else:
                for row in f["img"][::2]:
                    print("".join(" " if p else "#" for p in row))


if __name__ == "__main__":
    main()
