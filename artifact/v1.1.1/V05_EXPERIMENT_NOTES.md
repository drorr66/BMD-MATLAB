# BMD-MATLAB v0.5 — focused cold-multiplication crossover

## Question

v0.4 showed that for the repeated-block family, canonical BMD equality and exact repeated multiplication strongly beat a conventional sparse term list.  The remaining scientific question is narrower:

> Does a **first, uncached multiplication** of two compact BMD operands become faster than sparse multiplication once structural sharing is large enough?

v0.5 tests only that question.

## Fixed polynomial family

The operands are fixed to

- `inner_n = 255`, i.e. `1 + x + ... + x^255`, and
- block spacing `2^10 = 1024`.

The second operand is

`1 + x^1024 + x^(2*1024) + ... + x^((blocks-1)*1024)`.

Their product has `(255+1)*blocks = 256*blocks` explicit terms, while the BMD remains very small because low exponent bits and block-index bits form shared decision structure.

Only `blocks` changes:

`96, 128, 160, 192, 256, 384, 512, 1024`.

## Timing protocol

- Operands are built before timing.
- Each BMD cold trial starts from a fresh compact manager, so the multiplication cache is empty.
- Sparse operands are also pre-built; only `sparse_terms_multiply` is timed.
- Nine independent trials are retained for every block count.
- Measurement order alternates BMD-first / sparse-first to reduce order bias.
- An untimed warm-up is performed before recording trials to reduce JIT/allocation startup effects.

The primary ratio is:

`ratio_median = sparse_median_time / BMD_median_time`.

Thus a ratio above 1 means BMD is faster.

## Robustness criterion

A median crossing alone can be noisy. v0.5 also reports a conservative ratio:

`robust_ratio_low = sparse_Q25 / BMD_Q75`.

If this is above 1, BMD is faster even when comparing a relatively fast sparse quartile against a relatively slow BMD quartile. That is labeled a robust BMD win.

The companion upper bound is:

`robust_ratio_high = sparse_Q75 / BMD_Q25`.

If it is below 1, sparse is a robust winner. Otherwise the timing distributions overlap.

## Outputs

- `results/cold_crossover_results_v05.csv` — one summary row per block count.
- `results/cold_crossover_trials_v05.csv` — all raw timing trials.
- `results/cold_crossover_summary_v05.csv` — bracketed median/robust crossover estimates and sustained-win points.
- `results/run_metadata_v05.txt` — MATLAB version and experiment parameters.

## Interpretation rule

The strongest desired result is not merely one point with `ratio_median > 1`. It is:

1. a bracketed crossing near the v0.4 prediction,
2. ratios that remain above 1 at larger block counts, and ideally
3. `robust_ratio_low > 1` for a sustained range.

That would support the claim that structural sharing can make canonical BMD multiplication faster than an explicit sparse term-list multiplication even on the **first** operation, rather than only after memoization.
