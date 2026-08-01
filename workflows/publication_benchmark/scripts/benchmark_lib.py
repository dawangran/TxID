"""Shared helpers for the publication benchmark wrappers."""

from __future__ import annotations

import hashlib
import json
import os
import random
import resource
import subprocess
import time
from collections import Counter
from pathlib import Path
from typing import Iterable, Mapping, Sequence


def _percentile(values: Sequence[float], probability: float) -> float:
    """Return a deterministic linearly interpolated sample percentile."""

    if not values:
        raise ValueError("cannot calculate a percentile from no values")
    if not 0.0 <= probability <= 1.0:
        raise ValueError("percentile probability must be between zero and one")
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] + fraction * (ordered[upper] - ordered[lower])


def clustered_ratio_bootstrap(
    strata: Sequence[tuple[int, int]],
    *,
    replicates: int,
    seed: int,
    confidence: float = 0.95,
) -> dict[str, object]:
    """Bootstrap a ratio of summed counts by resampling independent strata."""

    if replicates <= 0:
        raise ValueError("bootstrap replicates must be positive")
    if not 0.0 < confidence < 1.0:
        raise ValueError("bootstrap confidence must be between zero and one")
    checked = []
    for numerator, denominator in strata:
        if numerator < 0 or denominator < 0 or numerator > denominator:
            raise ValueError(
                "bootstrap strata require 0 <= numerator <= denominator"
            )
        checked.append((numerator, denominator))
    if not checked:
        return {
            "lower": 0.0,
            "upper": 0.0,
            "effective_replicates": 0,
        }
    if all(numerator == 0 for numerator, _ in checked):
        return {
            "lower": 0.0,
            "upper": 0.0,
            "effective_replicates": (
                replicates if any(denominator for _, denominator in checked) else 0
            ),
        }

    rng = random.Random(seed)
    estimates = []
    grouped = list(Counter(checked).items())
    binomial = getattr(rng, "binomialvariate", None)
    if binomial is None:
        for _ in range(replicates):
            sampled = [checked[rng.randrange(len(checked))] for _ in checked]
            denominator = sum(value[1] for value in sampled)
            if denominator:
                estimates.append(sum(value[0] for value in sampled) / denominator)
    else:
        # Resampling n strata with replacement is a multinomial draw. Grouping
        # identical strata and drawing the category counts conditionally keeps
        # the ordinary cluster bootstrap exact while avoiding O(n * replicates)
        # work on large public-data partitions.
        population = len(checked)
        for _ in range(replicates):
            draws_left = population
            population_left = population
            sampled_numerator = 0
            sampled_denominator = 0
            for index, ((numerator, denominator), category_size) in enumerate(
                grouped
            ):
                if index == len(grouped) - 1:
                    category_draws = draws_left
                else:
                    category_draws = binomial(
                        draws_left, category_size / population_left
                    )
                sampled_numerator += category_draws * numerator
                sampled_denominator += category_draws * denominator
                draws_left -= category_draws
                population_left -= category_size
            if sampled_denominator:
                estimates.append(sampled_numerator / sampled_denominator)
    if not estimates:
        return {
            "lower": 0.0,
            "upper": 0.0,
            "effective_replicates": 0,
        }
    tail = (1.0 - confidence) / 2.0
    return {
        "lower": _percentile(estimates, tail),
        "upper": _percentile(estimates, 1.0 - tail),
        "effective_replicates": len(estimates),
    }


def binary_precision_recall_bootstrap(
    truth_hits: Sequence[bool],
    prediction_hits: Sequence[bool],
    *,
    replicates: int,
    seed: int,
    confidence: float = 0.95,
) -> dict[str, float | int]:
    """Bootstrap transcript-level precision, recall, and F1 intervals."""

    if replicates <= 0:
        raise ValueError("bootstrap replicates must be positive")
    if not 0.0 < confidence < 1.0:
        raise ValueError("bootstrap confidence must be between zero and one")
    truth = [int(value) for value in truth_hits]
    predictions = [int(value) for value in prediction_hits]
    rng = random.Random(seed)
    estimates = {"precision": [], "recall": [], "f1": []}
    for _ in range(replicates):
        recall = (
            sum(truth[rng.randrange(len(truth))] for _ in truth) / len(truth)
            if truth
            else 0.0
        )
        precision = (
            sum(
                predictions[rng.randrange(len(predictions))]
                for _ in predictions
            )
            / len(predictions)
            if predictions
            else 0.0
        )
        f1 = (
            2.0 * precision * recall / (precision + recall)
            if precision + recall
            else 0.0
        )
        estimates["precision"].append(precision)
        estimates["recall"].append(recall)
        estimates["f1"].append(f1)
    tail = (1.0 - confidence) / 2.0
    result: dict[str, float | int] = {
        "replicates": replicates,
        "seed": seed,
    }
    for metric, values in estimates.items():
        result[f"{metric}_ci95_lower"] = _percentile(values, tail)
        result[f"{metric}_ci95_upper"] = _percentile(values, 1.0 - tail)
    return result


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_records(paths: Iterable[str | Path]) -> list[dict[str, object]]:
    records = []
    for raw_path in paths:
        path = Path(raw_path)
        records.append(
            {
                "path": str(path.resolve()),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    return records


def run_command(
    command: Sequence[str | Path],
    *,
    log_prefix: str | Path,
    env: Mapping[str, str] | None = None,
    cwd: str | Path | None = None,
) -> dict[str, object]:
    """Run one argv without a shell and retain complete stdout/stderr."""

    argv = [str(value) for value in command]
    prefix = Path(log_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    started_at = time.time()
    started = time.perf_counter()
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    result = subprocess.run(
        argv,
        cwd=cwd,
        env=dict(os.environ, **dict(env or {})),
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=False,
    )
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    record = {
        "argv": argv,
        "cwd": str(Path(cwd).resolve()) if cwd else str(Path.cwd().resolve()),
        "started_unix": started_at,
        "ended_unix": time.time(),
        "wall_seconds": time.perf_counter() - started,
        "child_user_seconds": after.ru_utime - before.ru_utime,
        "child_system_seconds": after.ru_stime - before.ru_stime,
        "max_rss_kib": after.ru_maxrss,
        "returncode": result.returncode,
        "stdout": str(prefix.with_suffix(".stdout.log").resolve()),
        "stderr": str(prefix.with_suffix(".stderr.log").resolve()),
    }
    prefix.with_suffix(".stdout.log").write_text(result.stdout, encoding="utf-8")
    prefix.with_suffix(".stderr.log").write_text(result.stderr, encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError(
            f"command failed with exit {result.returncode}: {argv!r}; "
            f"see {record['stderr']}"
        )
    return record


def write_json(path: str | Path, value: object) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f".{target.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, target)
