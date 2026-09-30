"""Frequency-only redshift hypotheses, not posterior probabilities."""
from dataclasses import dataclass, asdict
import csv
import math
from pathlib import Path

C_KMS = 299792.458

@dataclass(frozen=True)
class Transition:
    name: str
    rest_ghz: float

@dataclass(frozen=True)
class Pair:
    pair_id: str
    frequency_a_ghz: float
    frequency_b_ghz: float
    frame_a: str
    frame_b: str

@dataclass(frozen=True)
class Hypothesis:
    pair_id: str
    line_a: str
    line_b: str
    rest_a_ghz: float
    rest_b_ghz: float
    z: float
    z_a: float
    z_b: float
    velocity_mismatch_kms: float
    observed_a_ghz: float
    observed_b_ghz: float
    predicted_a_ghz: float
    predicted_b_ghz: float
    frame_status: str
    status: str = 'frequency_consistent_hypothesis'


def positive(value, field):
    try:
        x = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f'{field}: expected a number') from exc
    if not math.isfinite(x) or x <= 0:
        raise ValueError(f'{field}: expected a finite positive number')
    return x


def read_lines(path, unit='MHz'):
    """Legacy whitespace format: rest_frequency line_name [# comment]."""
    scales = {'Hz': 1e-9, 'kHz': 1e-6, 'MHz': 1e-3, 'GHz': 1.0}
    if unit not in scales:
        raise ValueError('line frequency unit must be Hz, kHz, MHz or GHz')
    result = []
    names = set()
    for i, raw in enumerate(Path(path).read_text(encoding='utf-8').splitlines(), 1):
        text = raw.split('#', 1)[0].strip()
        if not text:
            continue
        fields = text.split()
        if len(fields) != 2:
            raise ValueError(f'{path}:{i}: expected frequency and one whitespace-free name')
        frequency, name = fields
        if name in names:
            raise ValueError(f'{path}:{i}: duplicate transition name {name}')
        names.add(name)
        result.append(Transition(name, positive(frequency, f'line {i}') * scales[unit]))
    if not result:
        raise ValueError('empty line catalogue')
    return sorted(result, key=lambda r: (r.rest_ghz, r.name))


def normalize_frame(value):
    text = str(value or '').strip().upper()
    return {'BARY': 'BARYCENT', 'BARYCENTRIC': 'BARYCENT', 'TOPO': 'TOPOCENT',
            'TOPOCENTRIC': 'TOPOCENT', '': 'UNKNOWN'}.get(text, text)


def read_pairs(path):
    required = {'pair_id', 'frequency_a_ghz', 'frequency_b_ghz'}
    result = []
    ids = set()
    with open(path, newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        if not required.issubset(reader.fieldnames or []):
            raise ValueError('pair CSV requires: ' + ', '.join(sorted(required)))
        for row in reader:
            ident = row['pair_id'].strip()
            if not ident or ident in ids:
                raise ValueError('pair_id must be nonempty and unique')
            ids.add(ident)
            result.append(Pair(ident, positive(row['frequency_a_ghz'], ident + ' frequency A'),
                               positive(row['frequency_b_ghz'], ident + ' frequency B'),
                               normalize_frame(row.get('frame_a')), normalize_frame(row.get('frame_b'))))
    if not result:
        raise ValueError('empty pair catalogue')
    return result


def frame_status(pair, allow_mixed=False):
    a, b = normalize_frame(pair.frame_a), normalize_frame(pair.frame_b)
    if a == 'UNKNOWN' or b == 'UNKNOWN':
        return 'unverified_unknown_frame'
    if a != b:
        if not allow_mixed:
            raise ValueError(f'{pair.pair_id}: mixed frames {a}/{b}; transform frequencies first or explicitly use --allow-mixed-frames')
        return 'unverified_mixed_frames'
    return 'same_declared_frame_' + a


def match_pair(pair, lines, tolerance_kms=100., z_min=0., z_max=10., allow_mixed=False):
    """Enumerate distinct-transition assignments; retain every compatible one.

    z_i = rest_i / observed_i - 1; z = (z_a+z_b)/2.
    Acceptance: c*abs(z_a-z_b)/(1+z) <= tolerance_kms.
    This is a pair-wise velocity mismatch, NOT a per-line +/- tolerance.
    It is not a line-width, centroid uncertainty, or statistical significance.
    """
    positive(pair.frequency_a_ghz, 'frequency A')
    positive(pair.frequency_b_ghz, 'frequency B')
    tolerance = positive(tolerance_kms, 'tolerance_kms')
    if not all(math.isfinite(x) for x in [z_min, z_max]) or not 0 <= z_min <= z_max:
        raise ValueError('require finite 0 <= z_min <= z_max')
    status = frame_status(pair, allow_mixed)
    out = []
    for a in lines:
        for b in lines:
            if a.name == b.name or math.isclose(a.rest_ghz, b.rest_ghz, rel_tol=1e-10):
                continue  # Same transition / alias is not a second redshift constraint.
            za = a.rest_ghz / pair.frequency_a_ghz - 1
            zb = b.rest_ghz / pair.frequency_b_ghz - 1
            z = (za + zb) / 2
            if not z_min <= z <= z_max:
                continue
            dv = C_KMS * abs(za-zb) / (1+z)
            if dv <= tolerance + 1e-8:
                out.append(Hypothesis(pair.pair_id, a.name, b.name, a.rest_ghz, b.rest_ghz,
                                      z, za, zb, dv, pair.frequency_a_ghz, pair.frequency_b_ghz,
                                      a.rest_ghz/(1+z), b.rest_ghz/(1+z), status))
    return sorted(out, key=lambda h: (h.velocity_mismatch_kms, h.z, h.line_a, h.line_b))


def anchor_hypotheses(anchor_ghz, lines, z_min=0., z_max=10.):
    positive(anchor_ghz, 'anchor_ghz')
    return [{'anchor_line': r.name, 'z': r.rest_ghz / anchor_ghz - 1}
            for r in lines if z_min <= r.rest_ghz / anchor_ghz - 1 <= z_max]


def predicted_lines(z, lines, coverage):
    """Coverage is supplied intervals, never an inferred sensitivity limit."""
    return [dict(line=r.name, observed_ghz=r.rest_ghz/(1+z),
                 in_coverage=any(lo <= r.rest_ghz/(1+z) <= hi for lo, hi in coverage)) for r in lines]
