# Smith-Waterman

A small Python implementation of the Smith-Waterman local sequence alignment
algorithm.

## Usage

```python
from algorithm import ScoringScheme, smith_waterman

result = smith_waterman(
    "AACCGGTT",
    "TTCCGGAA",
    ScoringScheme(match=3, mismatch=-3, gap=-2),
)

print(result.query_alignment)      # CCGG
print(result.reference_alignment)  # CCGG
print(result.query_span)           # (2, 6)
print(result.reference_span)       # (2, 6)
print(result.score)                # 12
```

`smith_waterman` returns an immutable `AlignmentResult` with:

- aligned query and reference strings
- 0-based half-open spans for both input sequences
- alignment score
- match, mismatch, and gap counts
- the NumPy scoring matrix

## Development

Install dependencies:

```bash
python3 -m pip install -r requirements.txt
```

Run tests:

```bash
python3 -m unittest -v
```
