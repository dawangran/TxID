"""Deterministic complete-linkage fuzzy grouping of exact transcript forms."""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Iterable

from .errors import TxIDError
from .models import StructuralIdentity


@dataclass(frozen=True, slots=True)
class FuzzyCluster:
    members: tuple[str, ...]
    bridge_status: str = "unambiguous"


def _coordinates(identity: StructuralIdentity) -> dict[str, object]:
    value = json.loads(identity.canonical_json)
    if identity.family not in {"TF1", "SE1"}:
        raise TxIDError(f"fuzzy grouping requires TF1 or SE1 forms, found {identity.family}")
    return value


def compatible(
    left: StructuralIdentity,
    right: StructuralIdentity,
    *,
    splice_tolerance: int,
    end_tolerance: int,
) -> bool:
    if splice_tolerance < 0 or end_tolerance < 0:
        raise ValueError("fuzzy tolerances must be non-negative")
    if left.family != right.family:
        return False
    a = _coordinates(left)
    b = _coordinates(right)
    if any(a[key] != b[key] for key in ("assembly", "contig", "strand")):
        return False
    if left.family == "SE1":
        return (
            abs(int(a["start"]) - int(b["start"])) <= end_tolerance
            and abs(int(a["end"]) - int(b["end"])) <= end_tolerance
        )
    introns_a = a["introns"]
    introns_b = b["introns"]
    if not isinstance(introns_a, list) or not isinstance(introns_b, list):
        return False
    if len(introns_a) != len(introns_b):
        return False
    junctions_match = all(
        abs(int(pair_a[0]) - int(pair_b[0])) <= splice_tolerance
        and abs(int(pair_a[1]) - int(pair_b[1])) <= splice_tolerance
        for pair_a, pair_b in zip(introns_a, introns_b)
    )
    return (
        junctions_match
        and abs(int(a["tss"]) - int(b["tss"])) <= end_tolerance
        and abs(int(a["tes"]) - int(b["tes"])) <= end_tolerance
    )


def cluster_forms(
    identities: Iterable[StructuralIdentity],
    *,
    splice_tolerance: int,
    end_tolerance: int,
) -> list[FuzzyCluster]:
    if splice_tolerance < 0 or end_tolerance < 0:
        raise ValueError("fuzzy tolerances must be non-negative")
    ordered = sorted(identities, key=lambda item: item.public_id)
    by_id = {item.public_id: item for item in ordered}
    clusters: list[list[str]] = []
    statuses: list[str] = []
    for identity in ordered:
        matches = [
            index
            for index, members in enumerate(clusters)
            if all(
                compatible(
                    identity,
                    by_id[member],
                    splice_tolerance=splice_tolerance,
                    end_tolerance=end_tolerance,
                )
                for member in members
            )
        ]
        if len(matches) == 1:
            clusters[matches[0]].append(identity.public_id)
        elif len(matches) > 1:
            clusters.append([identity.public_id])
            statuses.append("ambiguous_bridge")
        else:
            clusters.append([identity.public_id])
            statuses.append("unambiguous")
    results = [
        FuzzyCluster(tuple(sorted(members)), status)
        for members, status in zip(clusters, statuses)
    ]
    return sorted(results, key=lambda cluster: cluster.members)

