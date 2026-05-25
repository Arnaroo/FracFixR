# FracFixD v2.0.2 "Quokka-2" — Build Provenance

This document records exactly how the binaries shipped in
[`../bin/`](../bin/) were produced.  It is intended to give
downstream users full transparency on toolchain, optimisation
flags, runtime dependencies, and verification recipes.

The D source code that produced these binaries is **closed and
confidential**; see [`../README.md` → Source code](../README.md#source-code)
for licensing enquiries.  This document describes the
**build environment**, not the source.

---

## Compiler & build tools

| Tool | Version | Notes |
|---|---|---|
| LDC2 | 1.42.0 (LLVM 21.1.8, frontend DMD 2.112.1) | Self-hosted; built with LDC 1.41.0 |
| DUB  | 1.41.0 (Jan 2026) | D package + build manager |
| Host CPU | AMD Zen 2 (znver2) | Development workstation; binaries are tested cross-microarch |
| Host OS | Manjaro Linux 7.0.3 (rolling) | glibc target compatibility ≥ 2.34 |

---

## Optimisation flags (per microarch)

All three Linux x86_64 release builds share the following flag
set, varying only on `-mcpu`:

```
--link-defaultlib-shared=false   (static D runtime: phobos2-ldc, druntime-ldc)
-O3                              (full LLVM optimisation)
--mcpu=<target>                  (microarchitecture-specific)
--flto=full                      (link-time optimisation across modules)
-boundscheck=off                 (no array-bounds runtime checks in release)
```

| Artefact | `-mcpu=` value | Target |
|---|---|---|
| `fracfixd-linux-znver2-x86_64`    | `znver2`      | AMD Zen 2 / 3 / 4 (Ryzen 3000/4000/5000, EPYC Rome/Milan/Genoa) |
| `fracfixd-linux-broadwell-x86_64` | `broadwell`   | Intel Broadwell or newer (i5/i7 5th-gen+, Xeon E5 v4+) |
| `fracfixd-linux-generic-x86_64`   | `x86-64-v3`   | x86-64-v3 baseline (AVX2 + BMI2 + F16C; most 2013+ CPUs) |
| `fracfixd-cli-linux-x86_64-static` | (LDC default) | CLI-only build with the static D runtime; no GUI; no LAPACK |

Resource files (logos, icons) are embedded into the binary at
compile time via the LDC `-J=resources` flag.

The corresponding `dub build` invocations:

```bash
# znver2
LIBRARY_PATH=build/lib-shim dub build \
    --config=linux --build=release-znver2 --compiler=ldc2 --force

# broadwell
LIBRARY_PATH=build/lib-shim dub build \
    --config=linux --build=release-broadwell --compiler=ldc2 --force

# generic x86-64-v3 baseline
LIBRARY_PATH=build/lib-shim dub build \
    --config=linux --build=release-generic --compiler=ldc2 --force

# Static CLI
LIBRARY_PATH=build/lib-shim dub build \
    --config=cli-only-static --build=release-static --compiler=ldc2 --force
```

After each build the artefact is stripped:

```bash
strip --strip-unneeded fracfixd-linux-<microarch>-x86_64
```

---

## Runtime dependencies

### GUI builds (`fracfixd-linux-{znver2,broadwell,generic}-x86_64`)

```
linux-vdso.so.1                  (kernel-provided virtual DSO)
liblapack.so.3                   (system LAPACK)
libblas.so.3                     (system BLAS)
libm.so.6                        (libc math)
libgcc_s.so.1                    (GCC runtime)
libc.so.6                        (GNU libc)
libgfortran.so.5                 (Fortran runtime; linked transitively via LAPACK)
ld-linux-x86-64.so.2             (dynamic linker)
```

GTK 3 (`libgtk-3-0`) is **NOT** linked at compile time — it is
loaded dynamically via `dlopen` at GUI-launch time.  This means:

- CLI mode (`fracfixd --cli ...` or any subcommand) runs on
  headless hosts with no GTK3 installed at all.
- The same binary becomes a working CLI when its GUI mode
  cannot start (e.g. over SSH with no X-forwarding).

### Static CLI (`fracfixd-cli-linux-x86_64-static`)

```
linux-vdso.so.1
libblas.so.3
libm.so.6
libgcc_s.so.1
libc.so.6
libgfortran.so.5
ld-linux-x86-64.so.2
```

LAPACK is statically resolved out of the CLI's call set (the
core compositional fixup uses BLAS only; LAPACK is consumed
by the GUI's interactive recompute path).  The D runtime is
statically linked.

---

## File-format summary (release artefacts)

| File | Type | Size (stripped) |
|---|---|---|
| `fracfixd-linux-znver2-x86_64`     | ELF 64-bit LSB pie | ~7.9 MB |
| `fracfixd-linux-broadwell-x86_64`  | ELF 64-bit LSB pie | ~7.9 MB |
| `fracfixd-linux-generic-x86_64`    | ELF 64-bit LSB pie | ~7.9 MB |
| `fracfixd-cli-linux-x86_64-static` | ELF 64-bit LSB pie | ~1.9 MB |

All four are stripped (no debug symbols) and use position-
independent code (PIE).

---

## Verification

Each release ships a `SHA256SUMS` sidecar in
[`../bin/SHA256SUMS`](../bin/SHA256SUMS).  Verify a downloaded
binary:

```bash
cd FracFixD/bin
sha256sum -c SHA256SUMS
```

Or check a single file:

```bash
sha256sum fracfixd-linux-znver2-x86_64
# Compare against the entry in SHA256SUMS
```

After install, the binary's own self-report should agree:

```bash
fracfixd --version
# +----------------------------------------------------------------+
# |   Version:   2.0.2                                             |
# |   Codename:  Quokka-2                                           |
# |   Status:    Stable                                            |
# +----------------------------------------------------------------+
```

---

## CPU microarchitecture selection guide

If you don't know which binary to run:

1. **`uname -m`** says `x86_64` → you can run at least one of
   the three Linux variants.
2. **AMD Ryzen 3000 series or later** (3000/4000/5000/7000) /
   **AMD EPYC Rome / Milan / Genoa** → use the `znver2` build.
3. **Intel 5th-gen Core (Broadwell, ~2015) or later** → use the
   `broadwell` build.
4. **Anything else**, **unsure**, or **mixed-CPU cluster** →
   use the `generic` (x86-64-v3) build.  It runs everywhere
   the other two do; you give up a few % of microarch-specific
   gain.
5. **Headless / containerised / no GTK install** → use the
   `static` CLI build.

If you run a binary tuned for a microarch your CPU does not
support, Linux's dynamic linker will refuse to load it
(`SIGILL` on the first AVX-512 instruction etc.); pick a
broader build.

---

## Reproducibility

The build is **not bit-reproducible across machines** (LDC's
LTO incorporates timestamps and host paths).  It **is**
deterministic on a single machine across re-runs given a
fixed source tree, a fixed `dub.selections.json`, and the
same LDC version.  Downstream users who need bit-identical
re-builds should pin all three.

The `--bb-step-rule deterministic` flag on the `diffprop`
subcommand is the relevant *runtime* SIMD-reproducibility
knob: it routes the beta-binomial fitter through a basin-
desensitised L-BFGS-B variant that produces the same MLE
across CPU microarchitectures and SIMD widths.  Default
classic mode is deterministic on a single machine given a
fixed binary.
