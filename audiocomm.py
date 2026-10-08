"""audiocomm.py - acoustic modem (binary FSK), Python version for the labs.

Same frame as acoustic_modem.html, so a WAV produced or recorded by either can be
decoded by the other. The receiver needs no settings: they travel in the header.

  control part, always 4 / 8 kHz at 400 bit/s (known to every receiver in advance):
    0.1 s warm-up | MLS preamble (127, 255 or 511 bits) | header (sent 1, 2 or 4 times)
  message part, in the tones and rate announced by the header:
    coded, scrambled data | MLS postamble (same sequence as the preamble)

  header : type 2 | a 11 | b 11 | code 2 | rate 2 | tones 2 | 0 0  (MSB first), then CRC-8 (x^8 + x^2 + x + 1),
           all through the convolutional code;  type 0 text: a = number of characters, b = 0;
           type 1 image: a = rows, b = columns;  rate: 200 400 800 1600;  tones: low (4/8 kHz), high (8/16 kHz)
  code   : 0 none, 1 Hamming (7,4) + interleaver, 2 convolutional (171,133) K = 7 + interleaver,
           3 Reed-Solomon (255,223) over bytes
  data   : text = 8 bits per character (Latin-1, MSB first);
           image = pixels in column-major order (MATLAB img(:)), 1 = white
  scrambling: XOR with the LFSR x^15 + x^14 + 1 sequence, all-ones seed
  bit 0 -> sin(2 pi F0 t), bit 1 -> sin(2 pi F1 t), one symbol = 1/rate seconds
  The MATLAB scripts audioComm_tx.m / audioComm_rx.m use the earlier frame and do not read these frames.

Usage
  python audiocomm.py tx --text "Hello" -o hello.wav     # write a WAV (add --play to play it)
  python audiocomm.py tx --image bartS.png --step 2 -o img.wav
  python audiocomm.py tx --text "Hello" --code conv -o coded.wav   # codes: none hamming conv rs
  python audiocomm.py rx recording.wav                    # decode a recording
  python audiocomm.py rx recording.wav --soft ratio       # soft output ln(E1/E0) instead of the calibrated LLR
  python audiocomm.py tx --text "Hello" --code conv --rate 200 --pre 511 -o far.wav   # long range
  python audiocomm.py tx --text "Hello" --tones high --rate 800 -o fast.wav
  python audiocomm.py rx --listen 20                      # record 20 s from the mic, then decode
--play / --listen need the `sounddevice` package.
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
RATES, TONE_SETS = [200, 400, 800, 1600], ["low", "high"]
TONES = {"low": (4000, 8000), "high": (8000, 16000)}
FS_TX = 48000


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


# ---------------------------------------------------------------- channel codes
# Same codes as acoustic_modem.html; the code is announced in the header, so the receiver finds it there.
CODES = ["none", "hamming", "conv", "rs"]
RS_K, RS_P = 223, 32


def coded_length(code, n0):
    if code == "hamming":
        return 7 * math.ceil(n0 / 4)
    if code == "conv":
        return 2 * (n0 + 6)
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


def calibrate(pe0, pe1):
    """Level A (mu) and noise variance per real dimension (s2) of each tone, from its energies on the postamble:
    E[E^2] = A^2 + 2 s2 when the tone is sent (Rician), 2 s2 when it is not (Rayleigh)."""
    def tone(e, on):
        P = PILOTS[len(e)]
        p_on, p_off = np.mean(e[P == on] ** 2), np.mean(e[P != on] ** 2)
        a2, s2 = max(p_on - p_off, 1e-6 * p_on), max(p_off / 2, 1e-12 * p_on)
        return dict(mu=np.sqrt(a2), s2=s2, snr_db=10 * np.log10(a2 / (2 * s2)))
    return tone(np.asarray(pe0), 0), tone(np.asarray(pe1), 1)


def level_track(cal, pe0, pe1, e0, e1, W=24):
    """Signal level (whichever tone is on) averaged over +-W symbols, relative to the preamble."""
    t0, t1 = cal
    q = lambda a, b: np.maximum(a ** 2 / t0["mu"] ** 2, b ** 2 / t1["mu"] ** 2)
    qref = np.mean(q(np.asarray(pe0), np.asarray(pe1)))
    cs = np.concatenate([[0], np.cumsum(q(np.asarray(e0), np.asarray(e1)))])
    j = np.arange(len(e0))
    a, b = np.maximum(0, j - W), np.minimum(len(e0), j + W + 1)
    return np.sqrt(np.clip((cs[b] - cs[a]) / (b - a) / qref, 0.01, 4))


def soft_llr(soft, cal, e0, e1, g=1.0):
    """'ratio': ln(E1/E0).  'cal': LLR of non-coherent FSK with the tone levels and noise measured on the preamble.
    In a room most of the disturbance is the echo of the signal, so it follows the signal level g: the energies
    are divided by g before using the preamble's levels and noise."""
    e0, e1 = np.asarray(e0, dtype=float), np.asarray(e1, dtype=float)
    if soft == "ratio":
        return np.log((e1 + 1e-9) / (e0 + 1e-9))
    t0, t1 = cal
    e0, e1 = e0 / g, e1 / g
    return (ln_i0(t1["mu"] * e1 / t1["s2"]) - t1["mu"] ** 2 / (2 * t1["s2"])
            - ln_i0(t0["mu"] * e0 / t0["s2"]) + t0["mu"] ** 2 / (2 * t0["s2"]))


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


def header_bits(typ, a, b, code, rate, tones):
    v = (msb(typ, 2) + msb(a, 11) + msb(b, 11) + msb(CODES.index(code), 2) + msb(RATES.index(rate), 2)
         + msb(TONE_SETS.index(tones), 2) + [0, 0])
    return encode_fec("conv", np.array(v + crc8(v)))


def build_frame(text=None, img=None, code="none", rate=400, tones="low", pre=127):
    """Return (segments, data): a list of (bits, tones, rate), then the message bits.
    img: 2-D boolean array (True = white). code: none, hamming, conv or rs."""
    if text is not None:
        codes = [ord(c) if ord(c) < 256 else 63 for c in text[:1023]]
        typ, a, b = 0, len(codes), 0
        data = np.array([x for c in codes for x in msb(c, 8)], dtype=int)
    else:
        typ, (a, b) = 1, img.shape
        data = img.astype(int).flatten(order="F")
    warm = [(i + 1) % 2 for i in range(round(0.1 * CTRL["rate"]))]
    hdr = header_bits(typ, a, b, code, rate, tones)
    ch = encode_fec(code, data)
    ctrl = np.concatenate([warm, PILOTS[pre], np.tile(hdr, HREPS[pre])]).astype(int)
    body = np.concatenate([ch ^ prbs(len(ch)), PILOTS[pre]]).astype(int)
    return [(ctrl, CTRL["tones"], CTRL["rate"]), (body, tones, rate)], data


def modulate(segs, fs=FS_TX):
    out = [np.zeros(round(0.2 * fs))]
    for bits, tones, rate in segs:
        F0, F1 = TONES[tones]
        t = np.arange(round(fs / rate)) / fs
        w = 0.9 * np.stack([np.sin(2 * np.pi * F0 * t), np.sin(2 * np.pi * F1 * t)])
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


# ---------------------------------------------------------------- receiver
def tone_energies(y, fs, rate, tones="low"):
    """|sum of y(n) exp(-j 2 pi F n / fs)| over a window of one symbol starting at each sample."""
    F0, F1 = TONES[tones]
    Lw = round(fs / rate)
    n = np.arange(len(y))
    out = []
    for F in (F0, F1):
        c = np.concatenate([[0], np.cumsum(y * np.exp(-2j * np.pi * F * n / fs))])
        out.append(np.abs(c[Lw:] - c[:-Lw]))      # window [m, m+Lw-1], m = 0 .. len(y)-Lw
    return out


def pilot_correlation(d, Lr, P):
    off = np.round(np.arange(len(P)) * Lr).astype(int)
    M = len(d) - off[-1]
    c = np.zeros(M)
    for j in range(len(P)):
        c += (2 * P[j] - 1) * d[off[j]:off[j] + M]
    return c / len(P)


def parse_header(hb):
    """hb: the 40 decoded header bits. None unless the CRC and the fields are valid."""
    if crc8(hb[:32]) != list(hb[32:40]):
        return None
    val = lambda s, n: int("".join(map(str, hb[s:s + n])), 2)
    typ, a, b, code, rate, tones = val(0, 2), val(2, 11), val(13, 11), CODES[val(24, 2)], RATES[val(26, 2)], val(28, 2)
    if tones >= len(TONE_SETS):
        return None
    h = None
    if typ == 0 and a >= 1 and b == 0:
        h = dict(type="text", n0=8 * a, len=a)
    if typ == 1 and 1 <= a <= 256 and 1 <= b <= 256 and a * b <= 20000:
        h = dict(type="image", n0=a * b, h=a, w=b)
    if h:
        h.update(code=code, nc=coded_length(code, h["n0"]), rate=rate, tones=TONE_SETS[tones])
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
    M = max(len(c) for c in corr.values())                   # a longer preamble cannot start near the end
    S = np.stack([np.pad(corr[n], (0, M - len(corr[n]))) / (0.3 * math.sqrt(127 / n)) for n in PILOTS])
    score, which = S.max(axis=0), np.array(list(PILOTS))[S.argmax(axis=0)]
    energies = {}                                             # tone energies at each (rate, tones) of the message part
    frames, start = [], 0
    while True:
        above = np.nonzero(score[start:] > 1)[0]
        if len(above) == 0:
            break
        cand = start + above[0]
        p1 = cand + int(np.argmax(score[cand:cand + 2 * Lw + 1]))
        n_p = int(which[p1])
        P, NH = PILOTS[n_p], HREPS[n_p] * NHC
        hm = np.round(p1 + (n_p + np.arange(NH)) * Lr).astype(int)
        if hm[-1] >= len(z0):
            break
        L = np.log((z1[hm] + 1e-9) / (z0[hm] + 1e-9)).reshape(HREPS[n_p], NHC).sum(axis=0)   # add up the copies
        hdr = parse_header([int(v) for v in decode_fec("conv", HBITS, L)[0]])
        if hdr is None:
            start = p1 + Lw
            continue
        n0, nc, rate, tones = hdr["n0"], hdr["nc"], hdr["rate"], hdr["tones"]
        if (rate, tones) not in energies:
            e = tone_energies(y, fs, rate, tones)
            energies[rate, tones] = e[0], e[1], (e[1] - e[0]) / (e[1] + e[0] + 1e-12)
        w0, w1, dd = energies[rate, tones]
        Ld = fs / rate
        Delta = (n_p + NH) * Lr + nc * Ld                     # nominal time from preamble to postamble
        tol = round(0.005 * Delta) + round(Ld)
        off = np.round(np.arange(n_p) * Ld).astype(int)
        lo, hi = round(p1 + Delta) - tol, round(p1 + Delta) + tol
        if hi + off[-1] >= len(dd):
            break
        c2 = dd[np.arange(lo, hi + 1)[:, None] + off] @ (2 * P - 1) / n_p
        p2 = lo + int(np.argmax(c2))
        found = c2.max() >= 0.5 * corr[n_p][p1]
        rho = (p2 - p1) / Delta if found else 1.0
        mk = np.round(p1 + ((n_p + NH) * Lr + np.arange(nc) * Ld) * rho).astype(int)
        post = np.round((p2 if found else p1 + Delta) + np.arange(n_p) * Ld * rho).astype(int)
        cal = calibrate(w0[post], w1[post])                    # the postamble is at the message's tones and rate
        g = 1.0 if soft == "ratio" else level_track(cal, w0[post], w1[post], w0[mk], w1[mk])
        llr = soft_llr(soft, cal, w0[mk], w1[mk], g)
        key = prbs(nc)
        bits, rs_fails = decode_fec(hdr["code"], n0, np.where(key == 1, -llr, llr), 4 if soft == "ratio" else np.inf)
        sent = encode_fec(hdr["code"], bits) ^ key            # re-encoded decision
        chan_errors = int(np.sum((llr > 0).astype(int) != sent))
        res = dict(hdr, pre=n_p, bits=bits, llr=llr, chan_errors=chan_errors, rs_fails=rs_fails, cal=cal, t=p1 / fs,
                   ppm=(rho - 1) * 1e6 if found else None,
                   pilot_errors=int(np.sum((w1[post] > w0[post]) != P)) if found else None,
                   corr=score, p1=p1, p2=p2)
        if hdr["type"] == "text":
            res["text"] = bytes(np.packbits(bits).tolist()).decode("latin-1")
        else:
            res["img"] = bits.reshape(hdr["w"], hdr["h"]).T.astype(bool)
        frames.append(res)
        if verbose:
            ppm = "postamble not found" if res["ppm"] is None else f"clock offset {res['ppm']:+.0f} ppm"
            what = repr(res["text"]) if hdr["type"] == "text" else f"{hdr['w']}x{hdr['h']} image"
            coded = "" if hdr["code"] == "none" else f", code {hdr['code']} corrected {chan_errors} of {nc} channel bits"
            if rs_fails:
                coded += f" ({rs_fails} Reed-Solomon blocks failed)"
            print(f"frame at {res['t']:.2f} s: {what}, {ppm}, postamble errors {res['pilot_errors']}{coded}")
            print(f"  from the header: {rate} bit/s, tones {tones} {TONES[tones]}, preamble {n_p} bits")
            print(f"  SNR per symbol {cal[0]['snr_db']:.0f} dB (low tone) / {cal[1]['snr_db']:.0f} dB (high tone), "
                  f"high tone {20 * np.log10(cal[1]['mu'] / cal[0]['mu']):+.1f} dB relative to low tone, soft output {soft}")
        start = round(p1 + Delta + n_p * Ld)
    return frames


# ---------------------------------------------------------------- command line
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("tx", "rx"):
        p = sub.add_parser(name)
        if name == "tx":
            p.add_argument("--text")
            p.add_argument("--rate", type=int, choices=RATES, default=400)
            p.add_argument("--tones", choices=TONES, default="low")
            p.add_argument("--pre", type=int, choices=sorted(PILOTS), default=127, help="preamble length, longer reaches farther")
            p.add_argument("--image")
            p.add_argument("--step", type=int, default=1)
            p.add_argument("--code", choices=CODES, default="none")
            p.add_argument("-o", "--out", default="tx.wav")
            p.add_argument("--play", action="store_true")
        else:
            p.add_argument("wav", nargs="?")
            p.add_argument("--listen", type=float, help="record this many seconds from the microphone")
            p.add_argument("--fs", type=int, default=48000)
            p.add_argument("--soft", choices=("cal", "ratio"), default="cal")
    a = ap.parse_args()
    if a.cmd == "tx":
        img = load_image(a.image, a.step) if a.image else None
        segs, data = build_frame(text=a.text if img is None else None, img=img, code=a.code, rate=a.rate, tones=a.tones, pre=a.pre)
        x = modulate(segs)
        write_wav(a.out, x, FS_TX)
        print(f"{len(data)} message bits, code {a.code}, {len(x) / FS_TX:.1f} s -> {a.out}")
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
            if f["type"] == "text":
                print("decoded:", f["text"])
            else:
                for row in f["img"][::2]:
                    print("".join(" " if p else "#" for p in row))


if __name__ == "__main__":
    main()
