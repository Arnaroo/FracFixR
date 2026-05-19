# Building FracFixD on Linux

> **NOTE.** The FracFixD source code is closed and confidential
> (see [`../README.md` → Source code](../README.md#source-code)).
> This document describes the **build environment and recipe**
> as a provenance reference and for licensed source-tree
> recipients.  End users should download the pre-built binaries
> from [`../bin/`](../bin/) instead.

---

## 1. Prerequisites

### Toolchain

- **LDC** 1.42.0 (LLVM 21.1.8 backend, DMD 2.112.1 frontend) or
  newer.  Earlier 1.4x versions usually work; 1.4x is the
  reference target.
- **DUB** 1.41.0 or newer (ships with LDC on most distros).
- **GCC** 11+ for the linker driver and binutils.

### Distro packages

```bash
# Arch / Manjaro
sudo pacman -S ldc dub gcc make pkg-config gtk3 openblas lapack gfortran

# Ubuntu 22.04+
sudo apt-get install ldc dub build-essential pkg-config \
                     libgtk-3-dev libgtkd-dev \
                     libopenblas-dev liblapack-dev libgfortran-13-dev

# Fedora 38+
sudo dnf install ldc dub gcc make pkg-config \
                 gtk3-devel openblas-devel lapack-devel libgfortran
```

### Source-tree access

The D source is held in a private development tree.  Licensed
source-tree recipients have direct access; see
[`../README.md` → Source code](../README.md#source-code) for
licensing enquiries.  Once you have the source tree the
following steps reproduce the public-release binaries.

---

## 2. Build the public-release variants

From the source-tree root:

```bash
# Build the lib-shim (one-time, used by all configurations)
make -C build/lib-shim                         # or similar; see the source-tree README

# AMD Zen 2 / 3 / 4 (Ryzen 3000+, EPYC Rome / Milan / Genoa)
LIBRARY_PATH=build/lib-shim dub build \
    --config=linux --build=release-znver2 --compiler=ldc2 --force
mv fracfixd dist/fracfixd-linux-znver2-x86_64

# Intel Broadwell or newer (i5/i7 5th-gen+, Xeon E5 v4+)
LIBRARY_PATH=build/lib-shim dub build \
    --config=linux --build=release-broadwell --compiler=ldc2 --force
mv fracfixd dist/fracfixd-linux-broadwell-x86_64

# Generic x86-64-v3 baseline (AVX2 + BMI2; most 2013+ CPUs)
LIBRARY_PATH=build/lib-shim dub build \
    --config=linux --build=release-generic --compiler=ldc2 --force
mv fracfixd dist/fracfixd-linux-generic-x86_64

# Static CLI (no GUI, minimal runtime deps)
LIBRARY_PATH=build/lib-shim dub build \
    --config=cli-only-static --build=release-static --compiler=ldc2 --force
mv fracfixd-cli-static dist/fracfixd-cli-linux-x86_64-static
```

Strip and checksum:

```bash
cd dist
strip --strip-unneeded fracfixd-linux-{znver2,broadwell,generic}-x86_64 \
                       fracfixd-cli-linux-x86_64-static
sha256sum * > SHA256SUMS
```

---

## 3. Smoke tests

```bash
./dist/fracfixd-linux-generic-x86_64 --version
# Should report:  Version 2.0.0, Codename Quokka, Status Stable

./dist/fracfixd-linux-generic-x86_64 --cli help diffprop
# Should print diffprop help WITHOUT initialising GTK

./dist/fracfixd-cli-linux-x86_64-static --version
ldd ./dist/fracfixd-cli-linux-x86_64-static
# Should list ≤ 8 dynamic libs (libc, libm, libblas, libgfortran, etc.)
```

Unit-test suite:

```bash
LIBRARY_PATH=build/lib-shim dub test --config=cli-only --compiler=ldc2
# 19/19 modules pass
```

Performance gate (in-tree fixture):

```bash
python3 tests/benchmark/perf_gate.py \
    --binary $(realpath ./dist/fracfixd-linux-generic-x86_64) \
    --workdir tests/benchmark/work_perf_gate_local \
    --gen-script tests/benchmark/gen_sizes.py \
    --time-script tests/benchmark/time_fracfixd.sh \
    --baseline tests/benchmark/baseline.json
# All cells PASS
```

---

## 4. Packaging

Drop the four binaries plus the `SHA256SUMS` into
[`FracFixD/bin/`](../bin/).  No installer is needed on Linux —
the binaries run directly.

For desktop integration (optional), install a `.desktop` file:

```ini
[Desktop Entry]
Type=Application
Name=FracFixD
GenericName=Compositional RNA Fractionation Analysis
Comment=Compositional fixup + differential proportion testing for fractional RNA-seq
Exec=/usr/local/bin/fracfixd
Icon=fracfixd
Terminal=false
Categories=Science;Biology;
```

---

## 5. Troubleshooting

- **"`ldc2: command not found`"** — install LDC from your
  distro or from the official tarball at
  https://github.com/ldc-developers/ldc/releases.
- **"`Illegal instruction (core dumped)`"** when running a
  microarch-tuned binary — the binary was built for a newer CPU
  than yours.  Use the `generic` build.
- **"`error while loading shared libraries: libgfortran.so.5`"** —
  install `libgfortran5` (Debian / Ubuntu) or `libgfortran`
  (Fedora) or use the **static CLI** build which avoids the
  dependency.
- **GUI fails to launch but `--cli` works** — GTK 3 is not
  installed.  `sudo apt-get install libgtk-3-0` (Debian /
  Ubuntu), `sudo dnf install gtk3` (Fedora), or use CLI mode.
