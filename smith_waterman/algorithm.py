# -*- coding: utf-8 -*-

from dataclasses import dataclass, field
from typing import Tuple

import numpy as np


GAP = "-"


@dataclass(frozen=True)
class ScoringScheme:
    """Scoring parameters for Smith-Waterman local alignment."""

    match: int = 3
    mismatch: int = -3
    gap: int = -2

    def __post_init__(self) -> None:
        if self.match <= 0:
            raise ValueError("match must be positive")
        if self.mismatch >= 0:
            raise ValueError("mismatch must be negative")
        if self.gap >= 0:
            raise ValueError("gap must be negative")


@dataclass(frozen=True)
class AlignmentResult:
    """The best local alignment and its coordinates."""

    query_alignment: str
    reference_alignment: str
    query_start: int
    query_end: int
    reference_start: int
    reference_end: int
    score: int
    scoring_matrix: np.ndarray = field(repr=False, compare=False)

    @property
    def query_span(self) -> Tuple[int, int]:
        return self.query_start, self.query_end

    @property
    def reference_span(self) -> Tuple[int, int]:
        return self.reference_start, self.reference_end

    @property
    def matches(self) -> int:
        return sum(
            query_base == reference_base != GAP
            for query_base, reference_base in zip(self.query_alignment, self.reference_alignment)
        )

    @property
    def mismatches(self) -> int:
        return sum(
            query_base != reference_base and query_base != GAP and reference_base != GAP
            for query_base, reference_base in zip(self.query_alignment, self.reference_alignment)
        )

    @property
    def gaps(self) -> int:
        return sum(
            query_base == GAP or reference_base == GAP
            for query_base, reference_base in zip(self.query_alignment, self.reference_alignment)
        )


def smith_waterman(
    query: str,
    reference: str,
    scoring: ScoringScheme = ScoringScheme(),
) -> AlignmentResult:
    """Return the highest-scoring Smith-Waterman local alignment."""

    matrix = _build_scoring_matrix(query, reference, scoring)
    score, end_i, end_j = _find_best_cell(matrix)
    query_alignment, reference_alignment, start_i, start_j = _traceback(
        query,
        reference,
        matrix,
        scoring,
        end_i,
        end_j,
    )
    matrix.setflags(write=False)

    return AlignmentResult(
        query_alignment=query_alignment,
        reference_alignment=reference_alignment,
        query_start=start_i,
        query_end=end_i,
        reference_start=start_j,
        reference_end=end_j,
        score=score,
        scoring_matrix=matrix,
    )


def _build_scoring_matrix(query: str, reference: str, scoring: ScoringScheme) -> np.ndarray:
    matrix = np.zeros((len(query) + 1, len(reference) + 1), dtype=np.int64)

    for i, query_base in enumerate(query, start=1):
        for j, reference_base in enumerate(reference, start=1):
            substitution = scoring.match if query_base == reference_base else scoring.mismatch
            diagonal = matrix[i - 1, j - 1] + substitution
            up = matrix[i - 1, j] + scoring.gap
            left = matrix[i, j - 1] + scoring.gap
            matrix[i, j] = max(0, diagonal, up, left)

    return matrix


def _find_best_cell(matrix: np.ndarray) -> Tuple[int, int, int]:
    row, column = np.unravel_index(int(matrix.argmax()), matrix.shape)
    return int(matrix[row, column]), int(row), int(column)


def _traceback(
    query: str,
    reference: str,
    matrix: np.ndarray,
    scoring: ScoringScheme,
    end_i: int,
    end_j: int,
) -> Tuple[str, str, int, int]:
    query_alignment = []
    reference_alignment = []
    i = end_i
    j = end_j

    while i > 0 and j > 0 and matrix[i, j] > 0:
        query_base = query[i - 1]
        reference_base = reference[j - 1]
        substitution = scoring.match if query_base == reference_base else scoring.mismatch

        if matrix[i, j] == matrix[i - 1, j - 1] + substitution:
            query_alignment.append(query_base)
            reference_alignment.append(reference_base)
            i -= 1
            j -= 1
        elif matrix[i, j] == matrix[i - 1, j] + scoring.gap:
            query_alignment.append(query_base)
            reference_alignment.append(GAP)
            i -= 1
        elif matrix[i, j] == matrix[i, j - 1] + scoring.gap:
            query_alignment.append(GAP)
            reference_alignment.append(reference_base)
            j -= 1
        else:
            break

    return (
        "".join(reversed(query_alignment)),
        "".join(reversed(reference_alignment)),
        i,
        j,
    )
