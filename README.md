# Smith-Waterman

A small Python implementation of the Smith-Waterman local sequence alignment
algorithm, plus browser-based teaching demos.

## Project Layout

```text
smith-waterman/
├── smith_waterman/          # Python package
├── tests/                   # Unit tests
├── demos/
│   ├── classic.html         # Single-file teaching demo
│   └── claude-design/       # Multi-style Claude Design demo
├── prompts/                 # Prompt sources used to generate demos
├── requirements.txt
└── README.md
```

## Python Usage

```python
from smith_waterman import ScoringScheme, smith_waterman

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

`smith_waterman` returns an immutable `AlignmentResult` with aligned strings,
0-based half-open spans, score, match/mismatch/gap counts, and the NumPy
scoring matrix.

## Demos

Open the classic standalone demo:

```text
demos/classic.html
```

Open the Claude Design random-style entry:

```text
demos/claude-design/index.html
```

The Claude Design entry randomly routes to one of the paper, terminal, or
editorial visual styles. Those pages share `demos/claude-design/sw-core.js`.

## Development

Install dependencies:

```bash
python3 -m pip install -r requirements.txt
```

Run tests:

```bash
python3 -m unittest -v
```
