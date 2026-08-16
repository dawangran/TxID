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


@dataclass(frozen=True, slots=True)
class _FuzzyVector:
    group: tuple[str, str, str, str, int]
    coordinates: tuple[int, ...]
    tolerances: tuple[int, ...]


@dataclass(slots=True)
class _ClusterState:
    members: list[str]
    minimum: list[int]
    maximum: list[int]

    @classmethod
    def create(cls, public_id: str, vector: _FuzzyVector) -> "_ClusterState":
        coordinates = list(vector.coordinates)
        return cls([public_id], coordinates.copy(), coordinates)

    def accepts(self, vector: _FuzzyVector) -> bool:
        # Pairwise complete linkage is equivalent to requiring the new value to
        # fall within tolerance of both extrema on every coordinate.
        return all(
            maximum - tolerance <= value <= minimum + tolerance
            for value, minimum, maximum, tolerance in zip(
                vector.coordinates,
                self.minimum,
                self.maximum,
                vector.tolerances,
            )
        )

    def add(self, public_id: str, vector: _FuzzyVector) -> None:
        self.members.append(public_id)
        for index, value in enumerate(vector.coordinates):
            self.minimum[index] = min(self.minimum[index], value)
            self.maximum[index] = max(self.maximum[index], value)


def _coordinates(identity: StructuralIdentity) -> dict[str, object]:
    value = json.loads(identity.canonical_json)
    if identity.family not in {"TF1", "SE1"}:
        raise TxIDError(f"fuzzy grouping requires TF1 or SE1 forms, found {identity.family}")
    return value


def _vector(
    identity: StructuralIdentity,
    *,
    splice_tolerance: int,
    end_tolerance: int,
) -> _FuzzyVector:
    value = _coordinates(identity)
    common = (
        identity.family,
        str(value["assembly"]),
        str(value["contig"]),
        str(value["strand"]),
    )
    if identity.family == "SE1":
        return _FuzzyVector(
            (*common, 0),
            (int(value["start"]), int(value["end"])),
            (end_tolerance, end_tolerance),
        )
    introns = value["introns"]
    if not isinstance(introns, list):
        raise TxIDError("TF1 canonical object has a non-list intron chain")
    flattened_introns = tuple(
        int(coordinate) for pair in introns for coordinate in pair
    )
    return _FuzzyVector(
        (*common, len(introns)),
        (int(value["tss"]), int(value["tes"]), *flattened_introns),
        (
            end_tolerance,
            end_tolerance,
            *(splice_tolerance for _ in flattened_introns),
        ),
    )


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
    a = _vector(
        left,
        splice_tolerance=splice_tolerance,
        end_tolerance=end_tolerance,
    )
    b = _vector(
        right,
        splice_tolerance=splice_tolerance,
        end_tolerance=end_tolerance,
    )
    return a.group == b.group and all(
        abs(left_coordinate - right_coordinate) <= tolerance
        for left_coordinate, right_coordinate, tolerance in zip(
            a.coordinates, b.coordinates, a.tolerances
        )
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
    vectors = [
        _vector(
            identity,
            splice_tolerance=splice_tolerance,
            end_tolerance=end_tolerance,
        )
        for identity in ordered
    ]
    clusters: list[_ClusterState] = []
    statuses: list[str] = []
    buckets: dict[
        tuple[str, str, str, str, int], dict[int, list[int]]
    ] = {}
    # TSS (or SE start) is always governed by end_tolerance. A compatible
    # cluster's first member must therefore be in the same or an adjacent bucket.
    bucket_width = end_tolerance + 1
    for identity, vector in zip(ordered, vectors):
        bucket = vector.coordinates[0] // bucket_width
        group_buckets = buckets.setdefault(vector.group, {})
        candidates = sorted(
            {
                index
                for candidate_bucket in (bucket - 1, bucket, bucket + 1)
                for index in group_buckets.get(candidate_bucket, ())
            }
        )
        matches = [
            index
            for index in candidates
            if clusters[index].accepts(vector)
        ]
        if len(matches) == 1:
            clusters[matches[0]].add(identity.public_id, vector)
        elif len(matches) > 1:
            clusters.append(_ClusterState.create(identity.public_id, vector))
            statuses.append("ambiguous_bridge")
            group_buckets.setdefault(bucket, []).append(len(clusters) - 1)
        else:
            clusters.append(_ClusterState.create(identity.public_id, vector))
            statuses.append("unambiguous")
            group_buckets.setdefault(bucket, []).append(len(clusters) - 1)
    results = [
        FuzzyCluster(tuple(sorted(cluster.members)), status)
        for cluster, status in zip(clusters, statuses)
    ]
    return sorted(results, key=lambda cluster: cluster.members)
