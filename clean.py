import io

SEGMENTS = {
    "text": ((10, 17), (18, 25)),
    "data": ((26, 33), (34, 41)),
    "analysis": ((42, 49), (50, 57)),
    "other": ((58, None), None)
}

filepath = "./fcs1.lmd"

fileout = "./fcs1_cleaned.lmd"

with open(filepath, "rb") as f:
    data = f.read()


def get_segments(data):
    SEGMENTS = {
        "text": ((10, 17), (18, 25)),
        "data": ((26, 33), (34, 41)),
        "analysis": ((42, 49), (50, 57)),
        "other": ((58, None), None)
    }

    def parse_pos(t, next_segment):
        if t is None:
            return t
        a, b = t
        if b is None:
            b = next_segment

        raw_pos = data[a:b+1]
        try:
            pos = int(raw_pos.decode())
        except ValueError:
            pos = None
        return pos

    positions = {}
    next_segment = None
    all_pos = []
    for name, (start, end) in SEGMENTS.items():
        positions[name] = (
            parse_pos(start, next_segment),
            parse_pos(end, next_segment),
        )
        all_pos += [p for p in positions[name] if p is not None]
        next_segment = min(all_pos)
    return positions


def get_segment_data(data, positions=None):
    if positions is None:
        positions = get_segments(data)

    segment_data = {}
    for name, (start, end) in positions.items():
        if start is None or end is None:
            segment_data[name] = None
            continue
        segment_data[name] = data[start:end + 1]
    return segment_data


def get_text_keys(textdata):
    delim = textdata[0:1]
    parts = textdata[1:].split(delim)
    return dict(zip(parts[::2], parts[1::2]))


def pack_meta(meta, encoding="utf-8", separator=b"\\"):
    encoded = {}
    for key, value in meta.items():
        if not isinstance(value, bytes):
            value.encode(encoding)
        encoded[key] = value

    parts = [t for p in encoded.items() for t in p]
    return separator + separator.join(parts) + separator


def generate_header(segments, start_offset=256, tag="FCS2.0", encoding="utf-8", padding=128):
    header = bytearray(b" "*start_offset)
    header[:len(tag)] = tag.encode(encoding)

    def insert_value(data, val, position):
        start, end = position
        val_enc = str(val).encode(encoding)
        val_len = len(val_enc)

        if len(val_enc) > end - start:
            raise ValueError("Size is too large.")
        data[end-val_len:end + 1] = val_enc

    last_pos = start_offset
    positions = {}
    for name, (start, end) in SEGMENTS.items():
        segment_length = len(segments[name]) if segments[name] is not None else None
        if segment_length is None:
            continue

        if segment_length > 0:
            segment_start = last_pos
            segment_end = last_pos + segment_length - 1
        else:
            segment_start = 0
            segment_end = 0

        if start is not None and start[1] is not None:
            insert_value(header, segment_start, start)

        if end is not None:
            insert_value(header, segment_end, end)

        last_pos = segment_end + padding
    return header, last_pos


def get_key(meta, keyname, typefun=lambda d: d):
    return typefun(meta[keyname.encode()].decode())


def calc_crc(_):
    # no crc implementation yet
    return b"00000000"


DIRTY_KEYS = [
    b"$INSTADDRESS",
    b"@LOCATION",
    b"$RUNNUMBER",
    b"@FILEGUID",
    b"$DATE",
    b"@Y2KDATE",
    b"@SETTINGSFILE",
    b"@SETTINGSFILEDATETIME",
    b"@SAMPLEID1",
    b"@SAMPLEID2",
    b"@SAMPLEID3",
    b"@SAMPLEID4",
    b"@CYTOMETERID",
    b"@BUILDNUMBER",
    b"$FIL",
    b"@Acquisition Protocol Offset"
]


class FCS:
    def __init__(self, data):
        self.segments = get_segment_data(data)
        meta = get_text_keys(self.segments["text"])
        next_data_index = get_key(meta, "$NEXTDATA", int)

        if next_data_index:
            self.next = FCS(data[next_data_index:])
        else:
            self.next = None

    @property
    def meta(self):
        return get_text_keys(self.segments["text"])

    def set_text_key(self, key, value):
        rkey = key.encode()
        rvalue = str(value).encode()
        m = self.meta
        m[rkey] = rvalue
        self.set_text(m)

    def set_text(self, meta):
        self.segments["text"] = pack_meta(meta)

    def compile(self):
        if self.next is None:
            self.set_text_key("$NEXTDATA", 0)
            header, last_pos = generate_header(self.segments)
        else:
            _, last_pos = generate_header(self.segments)

            # set next segment to some offset after last segment
            self.set_text_key("$NEXTDATA", last_pos + 32)

            header, _ = generate_header(self.segments)

        positions = get_segments(header)

        out_data = bytearray(last_pos)
        out_data[:len(header)] = header
        print(header)

        for name, (start, end) in positions.items():
            if start is None or end is None:
                continue
            out_data[start:end + 1] = self.segments[name]

        out_data += calc_crc(out_data)
        if self.next:
            out_data += bytes(24)
            out_data += self.next.compile()

        return out_data


fcsdata = FCS(data)
fcsdata.next = None
test = fcsdata.compile()

tseg = get_segment_data(bytes(test))
tsegh = get_segments(bytes(test))
t2 = FCS(bytes(test))

with open(fileout, "wb") as f:
    f.write(test)
