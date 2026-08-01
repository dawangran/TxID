"""Small pure implementation of the normative exact identity algorithms."""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Callable

from .models import IdentityBundle, StructuralIdentity, TranscriptModel

ALGORITHM_VERSION = "1"
PUBLIC_DIGEST_LENGTH = 24
_HEX_256 = re.compile(r"^[0-9a-f]{64}$")
DigestFunction = Callable[[bytes], str]


def canonical_json(value: dict[str, object]) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _make_identity(
    family: str,
    canonical_object: dict[str, object],
    digest_function: DigestFunction,
    public_digest_length: int,
) -> StructuralIdentity:
    serialized = canonical_json(canonical_object)
    full_digest = digest_function(serialized.encode("utf-8")).lower()
    if not _HEX_256.fullmatch(full_digest):
        raise ValueError("digest function must return exactly 64 lowercase hexadecimal characters")
    if not 1 <= public_digest_length <= 64:
        raise ValueError("public digest length must be between 1 and 64 hexadecimal characters")
    public_digest = full_digest[:public_digest_length]
    return StructuralIdentity(
        family=family,
        public_id=f"txid:{family}.{public_digest}",
        public_digest=public_digest,
        full_digest=full_digest,
        canonical_json=serialized,
    )


def identify(
    transcript: TranscriptModel,
    assembly_fingerprint: str,
    *,
    digest_function: DigestFunction = _sha256,
    public_digest_length: int = PUBLIC_DIGEST_LENGTH,
) -> IdentityBundle:
    common: dict[str, object] = {
        "assembly": assembly_fingerprint,
        "contig": transcript.contig,
        "strand": transcript.strand,
    }
    if transcript.exon_count == 1:
        single_object = {
            **common,
            "algorithm": "SE1",
            "start": transcript.start,
            "end": transcript.end,
        }
        return IdentityBundle(
            form=_make_identity(
                "SE1", single_object, digest_function, public_digest_length
            )
        )

    introns = [[start, end] for start, end in transcript.introns]
    splice_object = {**common, "algorithm": "SC1", "introns": introns}
    form_object = {
        **common,
        "algorithm": "TF1",
        "introns": introns,
        "tss": transcript.tss,
        "tes": transcript.tes,
    }
    return IdentityBundle(
        splice_chain=_make_identity(
            "SC1", splice_object, digest_function, public_digest_length
        ),
        form=_make_identity("TF1", form_object, digest_function, public_digest_length),
    )
