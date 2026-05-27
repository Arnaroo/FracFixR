# FracFixD v2.0.3 "Quokka-Static" — 2026-05-27

Static-linkage fix release for the Linux GUI binaries.

**No numeric-kernel changes** from v2.0.2.  Equivalence-harness
outputs are byte-identical to v2.0.2 on the synthetic and real-data
fixtures; the manuscript figures, the CLI default-flag output bytes
and the FFXD1BIN container layout all carry over unchanged.

## Background

v2.0.2 shipped three Linux GUI binaries
(`fracfixd-linux-{generic,broadwell,znver2}-x86_64`) that were
mis-described in the release notes and the manuscript as
"statically linked".  In fact they were dynamically linked against
the host's BLAS / LAPACK / gfortran / gomp stack, and would fail on
any Linux host without those system libraries installed.  Affected
users saw, on launch:

```
./fracfixd-linux-<arch>-x86_64: error while loading shared libraries:
libopenblas.so.0: cannot open shared object file
```

Only the CLI-static variant (`fracfixd-cli-linux-x86_64-static`) was
genuinely statically linked.

## What v2.0.3 fixes

All three Linux GUI binaries — and the CLI-static binary — are now
**fully statically linked** against OpenBLAS (with bundled netlib
LAPACK) and the gfortran / gomp / pthread runtimes.  `ldd` on every
Linux release artefact now shows only:

```
linux-vdso.so.1, libm.so.6, libgcc_s.so.1, libc.so.6,
/lib64/ld-linux-x86-64.so.2
```

— the standard glibc system trio that cannot realistically be
statically linked on Linux.  No `libopenblas`, `liblapack`,
`libblas`, `libgfortran`, `libgomp` or `libquadmath` references
remain in any `DT_NEEDED` entry.

GTK3 is loaded at runtime by the `fracfixd gui` subcommand via
`dlopen` and does not appear as a link-time dependency, matching
the v2.0.0 / v2.0.1 / v2.0.2 behaviour.

## Verification

- **Static linkage**: `ldd` per-binary confirmation; see "Linux
  binary downloads" below for `SHA256SUMS`.
- **Numerical regression**: the equivalence harness
  (`tests/equivalence/full/`) was run against FracFixR 1.1.0
  (CRAN tarball) with both the v2.0.2 dynamic and v2.0.3 static
  CLI binaries on the same synthetic 1000-transcript fixture.
  Both binaries produced **byte-identical** `d_diffprop_glm.tsv`
  and `d_diffprop_logit.tsv` outputs (matching MD5).  Static
  linkage is numerically transparent — it changes how OpenBLAS
  symbols are resolved (in-archive vs at runtime), not what they
  compute.
- **Per-microarch tuning preserved**: each GUI variant compiled
  with `--mcpu=` set to `znver2`, `broadwell`, or `x86-64-v3`
  respectively, `--flto=full`, `-O3`, `-boundscheck=off`,
  matching the v2.0.2 microarch surface.

## Linux binary downloads

| File | Microarch | Size (stripped) |
|---|---|---|
| `fracfixd-linux-znver2-x86_64` | AMD Zen 2 (AVX2 + BMI2 tuned) | 7.9 MB |
| `fracfixd-linux-broadwell-x86_64` | Intel Broadwell (AVX2) | 8.0 MB |
| `fracfixd-linux-generic-x86_64` | x86-64-v3 (generic AVX2 baseline) | 8.0 MB |
| `fracfixd-cli-linux-x86_64-static` | x86-64-v2 (SSE4.2, CLI-only) | 2.3 MB |

See `SHA256SUMS` attached to this release.

## macOS arm64 and Windows x86_64

The macOS arm64 `.dmg` / tarball and the Windows x86_64 portable
ZIP from v2.0.2 are not affected by this fix; they are carried
forward unchanged.  Verify with `otool -L` on macOS and the
bundled `OpenBLAS.dll` on Windows.

## Build provenance

- Host: Manjaro Linux, AMD Ryzen 5 4600G, 30 GB RAM
- Toolchain: LDC 1.42.0 (LLVM 21.1.8, DMD frontend 2.112.1)
- OpenBLAS: v0.3.27 built with `NO_SHARED=1 USE_THREAD=0
  NUM_THREADS=1 DYNAMIC_ARCH=1 NO_AFFINITY=1 NO_WARMUP=1
  NO_LAPACKE=1 BUILD_LAPACK_DEPRECATED=0` via
  `scripts/build_static_openblas.sh`
- Build types: `release-{znver2,broadwell,generic}-static` (GUI)
  and `release-cli-static-flto` (CLI-only); see `dub.json`.

## Carry-overs from v2.0.2

Everything from v2.0.2 stays in place:

- Multi-condition diagnostic plot suite (global-test volcano,
  p-value histogram, per-contrast volcano grid, top-N heatmap)
- FFXD1ICA fast-load cache for the 205k-row reference TSV
- Determinate progress + working Cancel on long-running jobs
- Ctrl-wheel / Ctrl-+/-/0 plot zoom
- 360×240 window minimum + per-tab scrollbars on shrink
- The full v2.0.0 statistical surface: seven differential-proportion
  test backends, three dispersion modes, three FDR procedures,
  three log₂ fold-change shrinkage estimators, permutation
  inference, Bayesian credible intervals, multi-condition global
  LRT / Wald.

## Citation

If you use FracFixD in published work please cite:

- the FracFixR method paper (Cleynen, Hauling & Shirokikh, 2026 —
  see the repository `CITATION.cff`), and
- this release via its Zenodo concept DOI
  [10.5281/zenodo.20234583](https://doi.org/10.5281/zenodo.20234583)
  or v2.0.3 version DOI (to be minted with this release).
