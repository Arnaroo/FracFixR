# FracFixD User-facing Changelog

This file records the release history of FracFixD.  Run
`fracfixd --changelog` for the same notes from the binary
itself.

---

## v2.0.0 "Quokka" 2026-05-19

First public release.

FracFixD ships as a binary-only subfolder inside the
[`Arnaroo/FracFixR`](https://github.com/Arnaroo/FracFixR)
umbrella repository, alongside the open-source FracFixR R
package.

### Highlights

- **Linux x86_64 release artefacts**, three microarchitecture-
  tuned GUI+CLI binaries plus a static CLI variant:
  - `fracfixd-linux-znver2-x86_64`     (AMD Zen 2 / 3 / 4)
  - `fracfixd-linux-broadwell-x86_64`  (Intel Broadwell or newer)
  - `fracfixd-linux-generic-x86_64`    (x86-64-v3 baseline)
  - `fracfixd-cli-linux-x86_64-static` (CLI only, no runtime deps
    beyond libc + libblas)
- **macOS arm64 release artefacts**, a drag-to-Applications
  `.dmg` plus a relocatable `.tar.gz` for CLI / pipeline use.
  GTK 3, OpenBLAS and gfortran runtimes are bundled into the
  `.app`; the launcher is a native Mach-O so Tahoe Gatekeeper
  accepts it.
- **Windows x86_64 release artefact**, a portable ZIP
  (`fracfixd-v2.0.0-windows-x86_64.zip`) containing
  `fracfixd.exe` plus the full GTK 3 + OpenBLAS + gfortran
  runtime closure (around 70 DLLs).  Drop the unzipped folder
  on any Windows 10/11 host and double-click `fracfixd.exe`,
  no installer or admin rights required.  A signed Inno Setup
  `.exe` installer can be produced from the source tree.
- **Single binary, GUI + CLI**: no command-line flags launches
  the GTK3 GUI; `--cli` or any subcommand drops into headless
  console mode.
- **Maximally optimised release builds** with `-O3 --flto=full
  -boundscheck=off` and a statically-linked D runtime.  Per-
  microarch `-mcpu=znver2 / broadwell / x86-64-v3` tuning.
- **Comprehensive statistical surface**, carried forward
  unchanged into the public release:
  - NNLS compositional fixup (plain / ridge / auto).
  - Seven differential-proportion test backends: GLM-LRT,
    logit-Wald, beta-binomial Wald, Rao score, HC0 / HC3
    sandwich, quasi-binomial.
  - Multi-condition global LRT or joint Wald across K ≥ 2
    conditions, optional per-pair contrasts.
  - Three dispersion modes (global, trend, per-transcript) with
    optional Bayesian shrinkage prior on log-φ.
  - Three FDR procedures (BH, BY, Storey) with optional
    independent filtering.
  - Three posterior log₂FC shrinkage estimators (normal,
    apeglm, ashr).
  - Per-transcript permutation p-values with deterministic
    seeding.
- **Resumable pipeline** via `--cache-fits` / `--from-cache`,
  save FracFix fits once and reuse them from any number of
  `diffprop` invocations.
- **Native FFXD1BIN proportions format**, ~30× faster than TSV
  to parse on 100k-transcript fixtures.
- **Native SVG volcano plots** with optional EnhancedVolcano-
  style R reproduction script export.

Source-build instructions for all three platforms ship in
[`docs/BUILD_LINUX.md`](docs/BUILD_LINUX.md),
[`docs/BUILD_MACOS.md`](docs/BUILD_MACOS.md) and
[`docs/BUILD_WINDOWS.md`](docs/BUILD_WINDOWS.md).

### Equivalence

Default-flag output is verified against FracFixR on every
release via an in-tree equivalence harness:

- `--test glm` / `--test logit`: Spearman ρ(log₂FC) ≥ 0.99,
  ρ(−log₁₀ p) ≥ 0.95 on the standard synthetic fixtures.
- `--test wald`: same Spearman targets; the `|Δlog₁₀ p|`
  envelope is 0.20 because the beta-binomial MLE has a wider
  convergence basin than the GLM LRT route.

### Notes

- Binaries are released under **CC-BY-NC-ND-4.0**.  D source
  code is closed and confidential; for commercial licensing or
  source-tree access see [`README.md` → Source code](README.md#source-code).
- The accompanying FracFixR R package (sibling [`../CRAN/`](../CRAN/))
  is the canonical scientific reference and is released
  open-source under CC BY 4.0.
