# FracFixD User-facing Changelog

This file records the release history of FracFixD.  Run
`fracfixd --changelog` for the same notes from the binary
itself.

---

## v2.0.6 "Quokka-5" 2026-05-31

Documentation and in-app help release.  No numeric-kernel changes
versus 2.0.5; analysis output is unchanged on all valid inputs.

### Highlights

- **New in-app Help tab**: the full user guide (tool introduction,
  GUI tour with screenshots, every CLI subcommand and option, and
  worked CLI examples) now renders inside the application, baked
  into the binary, so GUI-only users get the complete reference
  without a terminal or a browser.
- **About / License credits updated**: COMPASS (a division of
  Biocodecs) added to the About credits; the binary licence text
  now attributes copyright to Biocodecs and Arnaroo Ribologicals.
- **No surface changes** to commands, file formats, or output;
  numeric results are byte-identical to v2.0.5.
- All three platforms rebuilt for v2.0.6: Linux x86_64 microarch
  variants + static CLI, macOS arm64 .dmg + relocatable tarball,
  Windows x86_64 portable ZIP.

---

## v2.0.5 "Quokka-4" 2026-05-29

Audit-closeout and equivalence-revalidation release.  Output is
byte-identical to v2.0.4 on all valid inputs.

### Highlights

- **FracFixR equivalence reconfirmed end-to-end** against FracFixR
  1.1.0 at manuscript-grade thresholds: `glm` / `logit` Spearman
  ρ(log₂FC) = 1.0000, and corrected-proportion / normalized-count
  min Spearman = 1.0000.  Confirms the v2.0.4 pooled-denominator
  fix collapses the FracFixD-vs-FracFixR proportion scatter onto
  the y = x diagonal.
- **Three boundary sanity-checks added** in the differential-
  proportion path (`extractConditionMatrix`,
  `extractConditionMatrixMulti`, `runPermutationSweep`): each
  validates matrix / annotation index bounds eagerly.  No-ops on
  consistent inputs, so all numeric outputs stay byte-identical to
  v2.0.4.
- **No surface changes**: CLI subcommands and options, GTK GUI
  layout, FFXD1BIN / FFXD1ICA cache formats, and output filenames
  are unchanged.
- All three platforms rebuilt for v2.0.5: Linux x86_64 microarch
  variants + static CLI, macOS arm64 `.dmg` + relocatable tarball,
  Windows x86_64 portable ZIP + installer.

---

## v2.0.4 "Quokka-3" 2026-05-29

Numeric-correctness fix in the per-transcript correction path,
identified by Alice Cleynen.  No CLI, GUI, IO, or output-format
surface changes.

### Highlights

- **Pooled denominator in per-transcript proportions**: the
  per-transcript proportion now uses the condition-pooled Total
  (`TotalSum`, summed over all replicate Total columns) as the
  denominator ceiling, matching FracFixR 1.1.0's
  `ProcessReplicate()` semantics, instead of the single-replicate
  Total used previously.  Eliminates a ~3× systematic
  over-estimate of source-column proportions and the resulting
  triangular FracFixD-vs-FracFixR scatter.
- **Total-column proportions unchanged** (filled by a separate
  stage); NNLS, GLM-IRLS, beta-binomial Wald, FDR, shrinkage, CI
  and multi-condition paths were already pooled-correct and are
  bytewise unchanged.

---

## v2.0.3 "Quokka-Static" 2026-05-27

Static-linkage fix release for the Linux GUI binaries.  No changes
to the numeric kernels, the GUI surface, or the CLI subcommand
interface — equivalence-harness outputs are byte-identical to v2.0.2.

### Highlights

- **Fully static Linux GUI binaries**: `fracfixd-linux-{znver2,
  broadwell,generic}-x86_64` are now truly statically linked
  against OpenBLAS (with bundled netlib LAPACK).  `ldd` drops to
  `libc / libm / libgcc_s` only, matching the long-standing
  CLI-static variant.  Resolves the v2.0.2 `libopenblas.so.0:
  cannot open shared object file` failure on hosts without system
  BLAS / LAPACK installed.
- **Full LTO + per-microarch tuning preserved**: `--flto=full`,
  `-O3`, `-boundscheck=off`, `--mcpu={znver2,broadwell,x86-64-v3}`
  carried across all four Linux release artefacts.  CLI-static
  rebuilt against the same `DYNAMIC_ARCH=1` OpenBLAS archive for
  parity.
- **Numerical regression verified**: same equivalence-harness
  fixture run with v2.0.2 dynamic and v2.0.3 static CLI binaries
  produced byte-identical TSV outputs (MD5-equal).  Static
  linkage is numerically transparent.
- **`dub.json` build types added**: new `linux-static` config
  plus four combined build types
  (`release-{znver2,broadwell,generic}-static` and
  `release-cli-static-flto`) capture the working static-link
  line for reproducible re-cuts.

macOS arm64 and Windows x86_64 artefacts from v2.0.2 are carried
forward unchanged; this fix is Linux-only.

---

## v2.0.2 "Quokka-2" 2026-05-26

GUI usability + multi-condition visualisation release.  No
changes to the numeric kernels — equivalence-harness numbers and
CLI output bytes carry over from v2.0.1 unchanged.

### Highlights

- **Multi-condition plots**: four new visualisations for the
  K ≥ 3 / global-test pipeline, surfaced via a Pairwise /
  Multi-Cond sub-notebook inside the [Plots] tab:
  - **Global test (test stat)** — log₁₀(χ²/Wald) vs −log₁₀(padj),
    coloured by significance threshold, top-N labelled.
  - **p-value histogram** — 50-bin distribution with the cutoff
    drawn in.  Standard diagnostic for any multi-test pipeline:
    uniform = no signal, spike near 0 = real effects, spike near
    1 = test misspecification.
  - **Per-contrast volcano grid** — one pairwise volcano per
    contrast pair (e.g. `Mix1_vs_Mix2`, `Mix1_vs_Mix3`) tiled in
    a grid so you can see which contrasts drive the global hit.
  - **Top-N condition-means heatmap** — row z-score (or raw mean)
    proportion per condition for the top-N most-significant
    transcripts; blue/red colour scale, configurable row count.
- **Plot zoom**: Ctrl+wheel zooms whichever plot pane the cursor
  is over; Ctrl++ / Ctrl+- / Ctrl+0 zooms every visible pane.
  Re-rasterised from source SVG at each zoom step for vector-
  sharp output.
- **GUI load-choke fix**: a hand-rolled TSV parser plus an
  on-disk parsed-input cache (FFXD1ICA format, sidecar
  `<input>.ffxdcache` or XDG cache-dir fallback) cuts the
  205 000-row reference load from many seconds to ~95 ms cold /
  ~80 ms warm.  Determinate byte-progress bar and a working
  Cancel button now appear on long jobs (rooted GLib Idle /
  Timeout sources fix the GUI-freeze-on-load regression v2.0.1
  shipped with).
- **GUI window-resize polish**: minimum window size dropped to
  360×240 (fits split-screen, tablet, VNC); each tab grows
  horizontal and vertical scrollbars on demand so no widget
  becomes unreachable when the user shrinks the window.
- **Volcano legend XML-escape fix**: the literal `padj<cut`
  legend text was being parsed by librsvg as the start of a
  `<cut>` tag, which silently blanked the GUI volcano preview.
  Legend strings now route through the existing `escapeXml()`
  helper (`padj&lt;cut`, `|log2FC|&gt;cut`).

### Compatibility

- Default-flag CLI invocations and GUI runs that leave the new
  selectors at their defaults produce byte-identical output to
  v2.0.1.
- The on-disk FFXD1ICA cache is opt-in via the GUI's Data tab;
  the CLI is unaffected.  Cache files are validated against the
  source TSV's mtime + size before reuse, so editing the input
  silently invalidates the cache.
- All three platforms (Linux x86_64 microarch variants, macOS
  arm64 .dmg + relocatable tarball, Windows x86_64 portable
  ZIP) rebuilt for v2.0.2.  macOS built on the macincloud
  Apple Silicon host (macOS 26.2 Tahoe, Homebrew GTK 3.24.52,
  OpenBLAS 0.3.33, LDC 1.42.0).

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
