# FracFixD — User-facing Changelog

This file records the **public-facing** release history of
FracFixD.  Run `fracfixd --changelog` (or `--changelog-full`) for
the detailed in-binary changelog with per-session deliverables.

---

## v2.0.0 "Quokka" — 2026-05-19

First public release.

FracFixD ships as a binary-only subfolder inside the
[`Arnaroo/FracFixR`](https://github.com/Arnaroo/FracFixR)
umbrella repository, alongside the open-source FracFixR R
package (CRAN-ready v1.1.0).

### Highlights

- **Linux x86_64 release artefacts** (three microarchitecture-
  tuned GUI+CLI binaries plus a static CLI variant):
  - `fracfixd-linux-znver2-x86_64`     (AMD Zen 2 / 3 / 4)
  - `fracfixd-linux-broadwell-x86_64`  (Intel Broadwell or newer)
  - `fracfixd-linux-generic-x86_64`    (x86-64-v3 baseline)
  - `fracfixd-cli-linux-x86_64-static` (CLI only, no runtime deps
    beyond libc + libblas)
- **Single binary, GUI + CLI** — no command-line flags → GTK3
  GUI; `--cli` or any subcommand → headless console mode.
- **Maximally optimised release builds** — `-O3 --flto=full
  -boundscheck=off` with `--link-defaultlib-shared=false` to
  statically link the D runtime.  Per-microarch `-mcpu=znver2 /
  broadwell / x86-64-v3` tuning.
- **macOS arm64 and Windows x86_64 builds planned** for
  follow-up releases; source-build walk-throughs are shipped
  today in `docs/BUILD_MACOS.md` and `docs/BUILD_WINDOWS.md`.

### Method (carried forward from v1.5.x)

The statistical method has been stable since v1.5.0 and is
**not changed** for the public v2.0.0 cut.  Highlights:

- **Compositional fixup** via per-replicate non-negative least
  squares (`fracfix` subcommand): plain, ridge-penalised, and
  κ(X)-triggered-auto variants.
- **Differential-proportion testing** (`diffprop` subcommand)
  with seven test backends — GLM-LRT (`glm`), binomial Wald
  (`logit`), beta-binomial Wald (`wald`), Rao score (`score`),
  HC0 / HC3 sandwich, and quasi-binomial.
- **Multi-condition global tests** (`--multi-cond`) — LRT or
  joint Wald χ² across K ≥ 2 conditions with optional per-pair
  contrasts.
- **Three dispersion modes** — `global`, `trend`, `per-transcript`
  — with optional Bayesian shrinkage prior on log-φ
  (`--phi-prior trended`, df controlled by `--phi-prior-df`).
- **Three FDR procedures** — Benjamini-Hochberg, Benjamini-
  Yekutieli, Storey single-λ q-values — with optional
  independent filtering on a `mean` or `count` covariate.
- **Three posterior-shrinkage estimators** — empirical-Bayes
  Gaussian (`normal`), apeglm (Cauchy-prior posterior mode),
  ashr (adaptive-shrinkage scale mixture).
- **Permutation p-values** (`--permutation N`) with cluster-
  friendly deterministic per-transcript seeding.
- **Resumable pipeline** (`--cache-fits` / `--from-cache`) —
  save FracFix fits once, reuse from any number of `diffprop`
  invocations.
- **Native FFXD1BIN proportions format** — ~30× faster than TSV
  to parse on 100k-transcript fixtures.
- **Native SVG volcano plots** (`fracfixd plot`) with optional
  EnhancedVolcano-style R reproduction script export.

### Equivalence

Default-flag output is verified against FracFixR on every
release via the in-tree equivalence harness:

- `--test glm` / `--test logit`: Spearman ρ(log₂FC) ≥ 0.99,
  ρ(−log₁₀ p) ≥ 0.95.
- `--test wald`: same Spearman targets; `|Δlog₁₀ p|` envelope
  is 0.20 (the beta-binomial MLE has a wider convergence basin
  than the GLM LRT route, so small drift in z-scores is
  expected).
- `--bb-step-rule deterministic` is **EXPERIMENTAL**: on the
  N=1000 dev fixture it reaches log₂FC Spearman 0.98 vs
  classic; the residual gap is a small number of φ-boundary
  transcripts.  Default classic mode is unchanged.

### Notes

- Version jump v1.5.4 → v2.0.0 marks the public-release
  milestone.  No statistical-behaviour changes between v1.5.4
  "Basin-Tighten" and v2.0.0 "Quokka" — purely a packaging /
  version-bump cut.  Default-flag classic-mode output is
  bit-identical to v1.5.4.
- Binaries are released under **CC-BY-NC-ND-4.0**; D source is
  closed and confidential.  For commercial licensing or
  source-code access see [`README.md` → Source code](README.md#source-code).

---

## Earlier development history

The full per-session development history (35+ "sessions" from
the initial scaffold in May 2026 through to v1.5.4 "Basin-
Tighten") is embedded in the binary itself and accessible via:

```bash
fracfixd --changelog-full
```

This includes the v1.0.0 (Release-1.0), v1.1.x (Cache-Resume,
Equivalence-Harness), v1.2.x (Shrinkage, Sandwich, Score),
v1.3.x (FDR, Dispersion, Permutation), v1.4.x (Multi-condition,
Contrasts), and v1.5.x (Bayes-CI, Robust-NNLS, Resume,
Basin-SIMD, Basin-Tighten) pre-public-release milestones.
