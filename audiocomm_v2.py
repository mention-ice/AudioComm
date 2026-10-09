"""audiocomm_v2.py - acoustic modem (FSK with 2 or 4 tones) and room sounding, Python version for the labs.

Same frames as acoustic_modem_v2.html, so a WAV produced or recorded by either can be
decoded by the other. The receiver needs no settings: they travel in the header.

  control part, at 400 bit/s, 4 / 8 kHz or tone hopping (two modes known to every receiver in advance):
    0.1 s warm-up | MLS preamble (127, 255 or 511 bits) | header (sent 1, 2 or 4 times)
  message part, in the tones and rate announced by the header:
    coded, scrambled data | MLS postamble (same sequence as the preamble)

  header : type 2 | a 11 | b 11 | code 2 | rate 2 | tones 2 | 0 0  (MSB first), then CRC-8 (x^8 + x^2 + x + 1),
           all through the convolutional code;  type 0 text: a = number of characters, b = 0;
           type 1 image: a = rows, b = columns;  type 2 room sounding: a = MLS order (14, 15, 16), b = periods (1..8);
           rate: 200 400 800 1600;  tones: low (4/8 kHz), high (8/16 kHz), hop (tone hopping, 3.2 to 16 kHz)
           v2: tones 3 = 4 tones per symbol (4-FSK), the tone set (0 low, 1 high, 2 hop) in the last two bits, and the
           rate field read in 400 800 1600 3200; receivers older than v2 find no tone set 3 and ignore the frame;
           tone set 3 with rate field 3 = 6400 bit/s, one group of 4 tones over the whole band (3.2 / 6.4 / 9.6 /
           12.8 kHz, tone hopping with a single group); receivers older than 6400 bit/s ignore it
  code   : 0 none, 1 Hamming (7,4) + interleaver, 2 convolutional (171,133) K = 7 + interleaver,
           3 Reed-Solomon (255,223) over bytes (received only: v1 no longer sends it), 4 turbo (13, 15 octal,
           rate 1/2, as in 4G) + interleaver, 5 and 6 the same turbo code at rates 1/3 and 1/4 (1/4 adds cdma2000's second parity);
           codes 4, 5, 6 are sent as 0, 1, 2 with the CRC inverted (as 4G masks a CRC), so receivers that predate
           them ignore them
  data   : text = 8 bits per character (Latin-1, MSB first);
           image = pixels in column-major order (MATLAB img(:)), 1 = white
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
  The MATLAB scripts audioComm_tx.m / audioComm_rx.m use the earlier frame and do not read these frames.

Usage
  python audiocomm_v2.py tx --text "Hello" -o hello.wav     # write a WAV (add --play to play it)
  python audiocomm_v2.py tx --image bartS.png --step 2 -o img.wav
  python audiocomm_v2.py tx --text "Hello" --code conv -o coded.wav   # codes: none hamming conv turbo turbo3 turbo4
  python audiocomm_v2.py rx recording.wav                    # decode a recording
  python audiocomm_v2.py rx recording.wav --soft ratio       # soft output ln(E1/E0) instead of the calibrated LLR
  python audiocomm_v2.py tx --text "Hello" --code conv --rate 200 --pre 511 -o far.wav   # long range
  python audiocomm_v2.py tx --text "Hello" --tones high --rate 800 -o fast.wav
  python audiocomm_v2.py tx --text "Hello" --tones hop --code conv -o hop.wav   # tone hopping, for reverberant rooms
  python audiocomm_v2.py tx --text "Hello" --tones hop --code turbo -o turbo.wav   # turbo code (turbo3, turbo4: rates 1/3, 1/4)
  python audiocomm_v2.py tx --text "Hello" --mod 4 --tones hop --rate 3200 --code turbo -o fast4.wav   # 4 tones, 2 bits per symbol
  python audiocomm_v2.py tx --text "Hello" --mod 4 --tones hop --rate 6400 --code turbo -o fast6.wav   # the whole band, no hopping
  python audiocomm_v2.py rx --listen 20                      # record 20 s from the mic, then decode
  python audiocomm_v2.py tx --sound -o sound.wav             # room sounding: MLS of 0.68 s, 4 periods averaged
  python audiocomm_v2.py tx --sound --order 16 --periods 8 --pre 511 -o sound.wav   # large or reverberant room
  python audiocomm_v2.py rx recorded_sound.wav --ir room_ir.wav --plot   # room measures; saves the impulse response
  python audiocomm_v2.py sim hello.wav --ir room_ir.wav --snr 30 -o heard.wav   # play a frame through a measured room
  python audiocomm_v2.py sim hello.wav --rt60 0.6 --drr 0 --snr 30 --ppm 50 -o heard.wav   # or through a model room
  python audiocomm_v2.py sim hop.wav --rt60 0.6 --drr -6 -o far.wav && python audiocomm_v2.py rx far.wav   # beyond the critical distance
--play / --listen need the `sounddevice` package, --plot needs matplotlib.
"""
import argparse
import math
import wave
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
CODES = ["none", "hamming", "conv", "rs", "turbo", "turbo3", "turbo4"]
TURBO_R = {"turbo": 2, "turbo3": 3, "turbo4": 4}   # turbo codes: rate 1/r
RS_K, RS_P = 223, 32


def coded_length(code, n0):
    if code == "hamming":
        return 7 * math.ceil(n0 / 4)
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


# Hamming (7,4), systematic: [d0 d1 d2 d3 p0 p1 p2]
HAM_P = np.array([[1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1]])
HAM_COLS = np.vstack([HAM_P, np.eye(3, dtype=int)])        # the 7 columns of H, one per row


def ham_encode(d):
    u = np.concatenate([np.asarray(d, dtype=int), np.zeros(-len(d) % 4, dtype=int)]).reshape(-1, 4)
    return np.hstack([u, u @ HAM_P % 2]).ravel()


def ham_decode(c, n0):
    w = np.asarray(c, dtype=int).reshape(-1, 7).copy()
    s = w @ HAM_COLS % 2                                      # syndrome = column of H at the error
    for j in range(7):
        w[(s == HAM_COLS[j]).all(axis=1), j] ^= 1
    return w[:, :4].ravel()[:n0]


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


def encode_fec(code, d):
    if code == "hamming":
        return interleave(ham_encode(d))
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
    if code == "hamming":
        return ham_decode((deinterleave(llr) > 0).astype(int), n0), None
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


def header_bits(typ, a, b, code, rate, tones, M=2):
    c, t, four = CODES.index(code), TONE_SETS.index(tones), M == 4   # codes 4 to 7: field c mod 4 and the CRC inverted
    r = 3 if rate == RATE_WIDE else (RATES4 if four else RATES).index(rate)
    if rate == RATE_WIDE:
        t = 3                                    # 6400 bit/s: tone set 3, rate field 3
    v = (msb(typ, 2) + msb(a, 11) + msb(b, 11) + msb(c & 3, 2) + msb(r, 2)
         + msb(3 if four else t, 2) + msb(t if four else 0, 2))   # 4 tones: tones 3, the tone set in the last two bits
    return encode_fec("conv", np.array(v + [x ^ (c >> 2) for x in crc8(v)]))


def build_frame(text=None, img=None, code="none", rate=400, tones="low", pre=127, sound=None, mod=2):
    """Return (segments, data): a list of (bits, tones, rate[, index of the first hop[, M]]), then the message bits;
    with 4 tones (mod=4) the message segment holds tones 0..3 instead of bits, and M = 4.
    img: 2-D boolean array (True = white). code: none, hamming, conv, turbo, turbo3 or turbo4.
    sound = (order, periods): room sounding frame; its MLS segment is (samples at 48 kHz, None, None)."""
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
    else:
        typ, (a, b) = 1, img.shape
        data = img.astype(int).flatten(order="F")
    M = 4 if mod == 4 else 2
    if rate not in (RATES4 + [RATE_WIDE] if M == 4 else RATES):
        raise ValueError(f"rate {rate} bit/s not available with {M} tones: {RATES4 + [RATE_WIDE] if M == 4 else RATES}")
    if rate == RATE_WIDE and tones != "hop":
        raise ValueError(f"{RATE_WIDE} bit/s uses one group of 4 tones over the whole band: choose tones hop")
    warm = [(i + 1) % 2 for i in range(round(0.1 * CTRL["rate"]))]
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


def load_image(path, step=1):
    from PIL import Image
    I = np.asarray(Image.open(path).convert("RGB"), dtype=float).mean(axis=2)[::step, ::step]
    return I > I.max() / 2


def write_wav(path, x, fs):
    with wave.open(path, "wb") as w:
        w.setnchannels(1); w.setsampwidth(2); w.setframerate(fs)
        w.writeframes((np.clip(x, -1, 1) * 32767).astype("<i2").tobytes())


def read_wav(path):
    """16- or 32-bit PCM WAV -> (first channel as float, fs)."""
    with wave.open(path, "rb") as w:
        fs, ch, sw = w.getframerate(), w.getnchannels(), w.getsampwidth()
        raw = w.readframes(w.getnframes())
    y = np.frombuffer(raw, dtype={2: "<i2", 4: "<i4"}[sw]).reshape(-1, ch)[:, 0]
    return y / float(2 ** (8 * sw - 1)), fs


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


def parse_header(hb):
    """hb: the 40 decoded header bits. None unless the CRC and the fields are valid."""
    crc, rx = crc8(hb[:32]), [int(x) for x in hb[32:40]]
    if crc == rx:
        ext = 0
    elif all(c != r for c, r in zip(crc, rx)):   # inverted CRC: codes 4 to 7
        ext = 1
    else:
        return None
    val = lambda s, n: int("".join(map(str, hb[s:s + n])), 2)
    four = val(28, 2) == 3                       # tones 3: 4 tones per symbol, the tone set in the last two bits
    M = 4 if four else 2
    typ, a, b, c, rate = val(0, 2), val(2, 11), val(13, 11), val(24, 2) + 4 * ext, (RATES4 if four else RATES)[val(26, 2)]
    tones = val(30, 2) if four else val(28, 2)
    if four and tones == 3:                      # tone set 3: 6400 bit/s, one hopping group over the whole band
        if val(26, 2) != 3:
            return None
        tones, rate = TONE_SETS.index("hop"), RATE_WIDE
    if tones >= len(TONE_SETS) or c >= len(CODES) or (not four and val(30, 2)):
        return None
    code = CODES[c]
    h = None
    if typ == 0 and a >= 1 and b == 0:
        h = dict(type="text", n0=8 * a, len=a)
    if typ == 1 and 1 <= a <= 256 and 1 <= b <= 256 and a * b <= 20000:
        h = dict(type="image", n0=a * b, h=a, w=b)
    if typ == 2 and not four and a in SOUND_TAPS and 1 <= b <= 8:
        h = dict(type="sound", n0=0, order=a, periods=b)
    if h:
        nc, K = coded_length(code, h["n0"]), int(math.log2(M))
        h.update(code=code, nc=nc, rate=rate, tones=TONE_SETS[tones], M=M, baud=rate // K, ns=-(-nc // K))   # ns: message symbols
    return h


def decode(y, fs, verbose=True, soft="cal"):
    """Find and decode every frame in the recording y; no settings needed. Returns a list of dicts.
    soft: 'cal' (LLR calibrated on the postamble) or 'ratio' (ln(E1/E0))."""
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
    energies = {}                                             # tone energies at each (rate, tones) of the message part
    frames, start = [], 0
    while True:
        above = np.nonzero(score[start:] > 1)[0]
        if len(above) == 0:
            break
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
        if hdr["type"] == "text":
            res["text"] = bytes(np.packbits(bits).tolist()).decode("latin-1")
        else:
            res["img"] = bits.reshape(hdr["w"], hdr["h"]).T.astype(bool)
        frames.append(res)
        if verbose:
            ppm = "postamble not found (soft output calibrated on the message)" if res["ppm"] is None else f"clock offset {res['ppm']:+.0f} ppm"
            what = repr(res["text"]) if hdr["type"] == "text" else f"{hdr['w']}x{hdr['h']} image"
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
            p.add_argument("--mod", type=int, choices=(2, 4), default=2, help="tones per symbol: 2 (1 bit) or 4 (2 bits)")
            p.add_argument("--tones", choices=TONE_SETS, default="low", help="hop: tone hopping, for reverberant rooms")
            p.add_argument("--pre", type=int, choices=sorted(PILOTS), default=127, help="preamble length, longer reaches farther")
            p.add_argument("--image")
            p.add_argument("--step", type=int, default=1)
            p.add_argument("--code", choices=[c for c in CODES if c != "rs"], default="none",
                           help="turbo, turbo3, turbo4: turbo code of rate 1/2, 1/3, 1/4; Reed-Solomon (v0) is still decoded but no longer sent")
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
            p.add_argument("--ir", default="room_ir.wav", help="where to save the impulse response of a room sounding")
            p.add_argument("--plot", action="store_true", help="plot the room sounding (needs matplotlib)")
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
        img = load_image(a.image, a.step) if a.image else None
        sound = (a.order, a.periods) if a.sound else None
        segs, data = build_frame(text=a.text if img is None else None, img=img, code=a.code, rate=a.rate, tones=a.tones, pre=a.pre, sound=sound, mod=a.mod)
        x = modulate(segs)
        write_wav(a.out, x, FS_TX)
        if sound:
            print(f"room sounding: MLS of {2 ** a.order - 1} samples x {a.periods + 1}, preamble {a.pre}, {len(x) / FS_TX:.1f} s -> {a.out}")
        else:
            print(f"{len(data)} message bits, code {a.code}, {a.mod} tones, {len(x) / FS_TX:.1f} s -> {a.out}")
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
        frames = decode(y, fs, soft=a.soft)
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
            else:
                for row in f["img"][::2]:
                    print("".join(" " if p else "#" for p in row))


if __name__ == "__main__":
    main()
