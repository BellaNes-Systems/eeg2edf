"""Synthetic recordings in every supported format, built from the layouts the
readers document. Nothing here imports the package, so the same files can be
fed to any version of the converters."""
import datetime as dt
import struct
import uuid
from pathlib import Path

import numpy as np

# --- Nihon Kohden (EEG-1100 legacy layout) ---------------------------------

NK_START = dt.datetime(2024, 2, 6, 19, 55, 47)
NK_SFREQ = 100
NK_E21 = [0, 1, 2, 31, 42]  # Fp1, Fp2, C3 (renamed), EKG1, DC01
NK_BLOCK_SECONDS = (2.0, 3.0)


def _bcd(n):
    return (n // 10) << 4 | n % 10


def _nk_samples(block, n_samples, n_ch):
    """Offset-binary channel words plus a raw event word per frame."""
    rng = np.random.default_rng(block)
    ch = rng.integers(32768 - 3000, 32768 + 3000, size=(n_samples, n_ch), dtype=np.int64)
    ev = np.zeros((n_samples, 1), dtype=np.int64)
    ev[10, 0] = 1
    ev[20, 0] = 0xFFFF  # reads back as int16 -1
    return np.concatenate([ch, ev], 1).astype("<u2")


def write_nk(directory, stem="SYNTH"):
    """A two-block EEG-1100 .EEG with .21E, .PNT, .LOG and a .PTN folder.
    Returns {'eeg': path, 'blocks': [uint16 frames per block]}."""
    directory = Path(directory)
    n_ch = len(NK_E21)
    ctl = 0x100
    first = 0x200
    blocks, frames = [], []
    at = first
    start = NK_START
    for b, secs in enumerate(NK_BLOCK_SECONDS):
        dur10 = int(secs * 10)
        n = dur10 * NK_SFREQ // 10
        head = bytearray(0x27)
        head[0] = 1
        head[1:17] = b"TIME195547000000"
        head[0x14:0x1A] = bytes(_bcd(v) for v in (start.year % 100, start.month, start.day,
                                                  start.hour, start.minute, start.second))
        struct.pack_into("<H", head, 0x1A, NK_SFREQ)
        struct.pack_into("<I", head, 0x1C, dur10)
        head[0x26] = n_ch
        table = bytearray(n_ch * 10)
        for i, e in enumerate(NK_E21):
            table[i * 10] = e
        data = _nk_samples(b, n, n_ch)
        frames.append(data)
        blocks.append((at, bytes(head) + bytes(table) + data.tobytes()))
        at += len(blocks[-1][1])
        start += dt.timedelta(seconds=secs + 60)

    buf = bytearray(first)
    buf[:16] = b"EEG-1100A V01.00"
    buf[0x91] = 1
    struct.pack_into("<I", buf, 0x92, ctl)
    buf[ctl + 0x11] = len(blocks)
    for j, (addr, _) in enumerate(blocks):
        struct.pack_into("<I", buf, ctl + 0x12 + j * 20, addr)
    for _, raw in blocks:
        buf += raw

    eeg = directory / f"{stem}.EEG"
    eeg.write_bytes(bytes(buf))
    (directory / f"{stem}.21E").write_text(
        "[ELECTRODE]\n0000=Fp1\n0001=Fp2\n0002=C3\n0031=EKG1\n0003=\n"
        "[REFERENCE]\n0041=$0V\n"
        "[SYSTEM_SETUP]\nSystemReference=C3C4\nDeviceName=EEG-1200A\n",
        encoding="latin-1",
    )

    pnt = bytearray(0x680)
    pnt[0x64A:0x64A + 4] = b"Male"
    pnt[0x660:0x660 + 10] = b"2010/05/01"
    pnt[0x672:0x672 + 2] = b"13"
    (directory / f"{stem}.PNT").write_bytes(bytes(pnt))

    entries = [("REC START EMU1 EEG", NK_START),
               ("Eyes closed", NK_START + dt.timedelta(seconds=1)),
               ("Seizure", NK_START + dt.timedelta(seconds=63.0))]
    log = bytearray(0x300)
    log[0x91] = 1
    struct.pack_into("<I", log, 0x92, 0x100)
    log[0x100 + 0x12] = len(entries)
    for j, (text, when) in enumerate(entries):
        rec = text.encode().ljust(20, b"\0")
        rec += f"000000({when:%y%m%d%H%M%S})".encode()
        rec = rec.ljust(45, b"\0")
        lo = 0x100 + 0x14 + j * 45
        log[lo:lo + 45] = rec
    (directory / f"{stem}.LOG").write_bytes(bytes(log))

    ptn = directory / f"{stem}.PTN"
    ptn.mkdir(exist_ok=True)
    pat = bytearray(0x410 + 64 * 0x50)
    pat[0x80:0x84] = b"EMU1"
    for i, (a, b) in enumerate([(0, 1), (1, 2), (41, 31)]):  # Fp1-Fp2, Fp2-C3, 0V-EKG1
        o = 0x410 + i * 0x50
        pat[o], pat[o + 1], pat[o + 7] = a, b, 1
        struct.pack_into("<H", pat, o + 12, i)
    (ptn / "Pattern_001.PTN").write_bytes(bytes(pat))
    return {"eeg": eeg, "blocks": frames}


# --- Nicolet / Nervus .e ----------------------------------------------------

NIC_MAGIC = bytes.fromhex("0fc37b1ccf2d6d4b8aea1f64ced2b917")
NIC_TS_GUID = "{A271CCCB-515D-4590-B6A1-DC170C8D6EE2}"
NIC_DERIV_GUID = "{8A19AA48-BEA0-40D5-B89F-667FC578D635}"
NIC_ANNOTATION = "{A5A95612-A7F8-11CF-831A-0800091B5BDA}"
NIC_SEIZURE = "{A5A95608-A7F8-11CF-831A-0800091B5BDA}"
NIC_START = dt.datetime(2022, 9, 21, 3, 50, 11)
NIC_SEGMENTS = ((NIC_START, 2.0), (NIC_START + dt.timedelta(seconds=30), 1.0))
# label, active, reference, sfreq, resolution
NIC_CHANNELS = (("Fp1", "Fp1", "REF", 256.0, 0.25),
                ("Fp2", "Fp2", "REF", 256.0, 0.25),
                ("EKG", "EKG", "REF", 256.0, 0.5),
                ("Rate", "Rate", "", 1.0, 1.0))


def _w(text, size):
    return text.encode("utf-16-le").ljust(size, b"\0")[:size]


def _ole(when):
    return (when - dt.datetime(1899, 12, 30)).total_seconds() / 86400.0


def _event(when, duration, guid, channel, label):
    tail = (label + "\0").encode("utf-16-le")
    size = max(272, 0x108 + len(tail))
    rec = bytearray(size)
    struct.pack_into("<Q", rec, 0x10, size)
    struct.pack_into("<dd", rec, 0x20, _ole(when), 0.0)
    struct.pack_into("<d", rec, 0x30, duration)
    rec[0x88:0x98] = uuid.UUID(guid).bytes_le
    rec[0xA8:0xA8 + 64] = _w(channel, 64)
    rec[0x108:0x108 + len(tail)] = tail
    return bytes(rec)


def write_nicolet(directory, stem="SAMPLE"):
    """Returns {'path': .e path, 'streams': {channel: stored int16 stream}}."""
    ts = bytearray(760 + len(NIC_CHANNELS) * 552)
    struct.pack_into("<I", ts, 752, len(NIC_CHANNELS))
    for i, (label, active, ref, rate, res) in enumerate(NIC_CHANNELS):
        at = 760 + i * 552
        ts[at:at + 64] = _w(label, 64)
        ts[at + 0x80:at + 0xC0] = _w(active, 64)
        ts[at + 0xC0:at + 0x100] = _w(ref, 64)
        struct.pack_into("<dddd", ts, at + 0x100, 0.5 if ref else 0.0,
                         70.0 if ref else 0.0, rate, res)
        struct.pack_into("<HH", ts, at + 0x120, 0, 60 if ref else 0)

    traces = (("Fp1-Fp2", "Fp1", "Fp2"), ("EKG-Bipolar", "EKG", ""))
    deriv = bytearray(752 + len(traces) * 520)
    deriv[0x28:0x28 + 64] = _w("Bipolar", 64)
    struct.pack_into("<I", deriv, 744, len(traces))
    for i, (label, active, ref) in enumerate(traces):
        at = 752 + i * 520
        deriv[at:at + 64] = _w(label, 64)
        deriv[at + 0x80:at + 0xC0] = _w(active, 64)
        deriv[at + 0xC0:at + 0x100] = _w(ref, 64)

    seg = bytearray()
    for start, duration in NIC_SEGMENTS:
        rec = bytearray(152)
        struct.pack_into("<d", rec, 0, _ole(start))
        struct.pack_into("<d", rec, 0x10, duration)
        seg += rec

    events = (_event(NIC_START + dt.timedelta(seconds=1), 0.5, NIC_ANNOTATION, "", "onset")
              + _event(NIC_START + dt.timedelta(seconds=30.5), 0.0, NIC_SEIZURE, "Fp1", "big one")
              + _event(NIC_START + dt.timedelta(seconds=10), 0.0, NIC_ANNOTATION, "", "in gap"))

    total = sum(d for _, d in NIC_SEGMENTS)
    rng = np.random.default_rng(7)
    streams = {}
    for i, (_l, _a, _r, rate, _res) in enumerate(NIC_CHANNELS):
        n = int(round(total * rate))
        v = rng.integers(-2000, 2000, size=n).astype("<i2")
        v[0] = -32768  # has no positive counterpart once inverted
        streams[i] = v

    # section key -> list of payload chunks; channel streams are split in two
    # and interleaved so a reader has to follow the index.
    names = {1: NIC_TS_GUID, 2: "SegmentStream", 3: "Events", 4: NIC_DERIV_GUID}
    payload = {1: [bytes(ts)], 2: [bytes(seg)], 3: [events], 4: [bytes(deriv)]}
    for i in streams:
        key = 10 + i
        names[key] = str(i)
        raw = streams[i].tobytes()
        cut = (len(raw) // 2) & ~1
        payload[key] = [raw[:cut], raw[cut:]]
    order = [(k, 0) for k in payload] + [(k, 1) for k in payload if len(payload[k]) > 1]

    tag_at = 0xB0
    index_at = tag_at + len(names) * 84
    data_at = index_at + 8 + len(order) * 24
    head = bytearray(index_at)
    head[:16] = NIC_MAGIC
    struct.pack_into("<I", head, 0x18, index_at)
    struct.pack_into("<I", head, 0xAC, len(names))
    for j, (key, name) in enumerate(names.items()):
        at = tag_at + j * 84
        head[at:at + 80] = _w(name, 80)
        struct.pack_into("<I", head, at + 80, key)

    index = bytearray(struct.pack("<Q", len(order)))
    body = bytearray()
    for key, part in order:
        chunk = payload[key][part]
        index += struct.pack("<QQII", key, data_at + len(body), len(chunk), len(chunk))
        body += chunk
    path = Path(directory) / f"{stem}.e"
    path.write_bytes(bytes(head) + bytes(index) + bytes(body))
    return {"path": path, "streams": streams}


# --- Micromed VWR -----------------------------------------------------------

VWR_HEADER_SIZE = 640
VWR_SEGMENT_FIRST = 176
VWR_SEGMENT_SIZE = 16
VWR_ORDER_SIZE = 512
VWR_LABCOD_SIZE = 128
VWR_NOTE_COUNT, VWR_NOTE_SIZE = 1000, 44
VWR_MONTAGE_SIZE, VWR_MONTAGE_NAME, VWR_MONTAGE_INPUTS = 4096, 264, 328
VWR_HISTORY_TIMES = 512
VWR_TRIGGER_SIZE, VWR_TRONCA_SIZE, VWR_FLAGS_SIZE = 6, 8, 8
VWR_EVENT_NAME_SIZE, VWR_EVENT_COUNT = 64, 100

# LABCOD 0 is the recording ground, as in a real file: a montage input of 0
# means "this channel's own ground".
LABCOD = (
    (b"G2", b"", 0, 0, 256),
    (b"A", b"REF", 100, 0, 256),
    (b"B", b"REF", 50, 0, 256),
)
MONTAGES = (
    ("Bipolar", ((2, 1), (0, 2))),   # A-B, then B against its own ground
    ("Ground", ((0, 1),)),
)


def _montage_block(name, traces):
    block = bytearray(VWR_MONTAGE_SIZE)
    struct.pack_into("<4H", block, 0, len(traces), 0, 10, 1)
    block[VWR_MONTAGE_NAME:VWR_MONTAGE_NAME + len(name)] = name.encode()
    for index, (reference, active) in enumerate(traces):
        struct.pack_into("<2H", block, VWR_MONTAGE_INPUTS + index * 4, reference, active)
    return bytes(block)


def write_vwr(path, width=2):
    order, rate = 2, 4
    sizes = [
        ("ORDER", VWR_ORDER_SIZE),
        ("MONTAGE", len(MONTAGES) * VWR_MONTAGE_SIZE),
        ("LABCOD", len(LABCOD) * VWR_LABCOD_SIZE),
        ("NOTE", VWR_NOTE_COUNT * VWR_NOTE_SIZE),
        ("HISTORY", VWR_HISTORY_TIMES + VWR_MONTAGE_SIZE),
        ("TRIGGER", 8 * VWR_TRIGGER_SIZE),
        ("TRONCA", 2 * VWR_TRONCA_SIZE),
        ("FLAGS", VWR_FLAGS_SIZE),
        ("EVENT A", VWR_EVENT_NAME_SIZE + VWR_EVENT_COUNT * 8),
    ]
    offsets, at = {}, VWR_HEADER_SIZE
    for name, size in sizes:
        offsets[name] = at
        at += size
    data_offset = at

    header = bytearray(VWR_HEADER_SIZE)
    header[:16] = b"Micromed VWR\x1a".ljust(16, b"\0")
    header[128:134] = bytes((2, 1, 124, 3, 4, 5))
    struct.pack_into("<I", header, 138, data_offset)
    struct.pack_into("<H", header, 142, order)
    struct.pack_into("<H", header, 146, rate)
    struct.pack_into("<H", header, 148, width)
    header[175] = 4
    for index, (name, size) in enumerate(sizes):
        slot = VWR_SEGMENT_FIRST + index * VWR_SEGMENT_SIZE
        header[slot:slot + 8] = name.encode().ljust(8, b"\0")
        struct.pack_into("<II", header, slot + 8, offsets[name], size)

    order_table = struct.pack("<256H", 2, 1, *([0] * 254))
    montage = b"".join(_montage_block(name, traces) for name, traces in MONTAGES)
    labcod = bytearray(len(LABCOD) * VWR_LABCOD_SIZE)
    for index, (label, ground, lground, pmin, pmax) in enumerate(LABCOD):
        slot = index * VWR_LABCOD_SIZE
        struct.pack_into("<BB", labcod, slot, 1, 0)
        labcod[slot + 2:slot + 8] = label.ljust(6, b"\0")
        labcod[slot + 8:slot + 14] = ground.ljust(6, b"\0")
        struct.pack_into("<iiiii", labcod, slot + 14, 0, 255, lground, pmin, pmax)
        struct.pack_into("<h", labcod, slot + 34, 0)
    notes = bytearray(VWR_NOTE_COUNT * VWR_NOTE_SIZE)
    notes[4:9] = b"start"  # frame 0, which the reader used to drop
    struct.pack_into("<I", notes, VWR_NOTE_SIZE, 2)
    notes[VWR_NOTE_SIZE + 4:VWR_NOTE_SIZE + 9] = b"onset"
    # HISTORY names the montage the recording was made with.
    history = b"\xff" * VWR_HISTORY_TIMES + _montage_block(*MONTAGES[1])
    trigger = struct.pack("<IH", 1, 7) + b"\xff\xff\xff\xff\x00\x00" * 7
    tronca = struct.pack("<4I", 10, 0, 20, 2)
    flags = struct.pack("<2I", 0, 2)
    event_a = b"Seizure".ljust(VWR_EVENT_NAME_SIZE, b"\0")
    event_a += struct.pack(f"<{VWR_EVENT_COUNT}I", 1, *([0] * (VWR_EVENT_COUNT - 1)))
    event_a += struct.pack(f"<{VWR_EVENT_COUNT}I", 3, *([0] * (VWR_EVENT_COUNT - 1)))

    values = np.array([[51, 101], [52, 102], [53, 103], [54, 104]], dtype=f"<u{width}")
    blocks = {"ORDER": order_table, "MONTAGE": montage, "LABCOD": bytes(labcod),
              "NOTE": bytes(notes), "HISTORY": history, "TRIGGER": trigger,
              "TRONCA": tronca, "FLAGS": flags, "EVENT A": event_a}
    with open(path, "wb") as fh:
        fh.write(header)
        for name, size in sizes:
            assert len(blocks[name]) == size, name
            fh.write(blocks[name])
        fh.write(values.tobytes())
    return values
