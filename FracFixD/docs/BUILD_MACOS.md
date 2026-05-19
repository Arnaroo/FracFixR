# Building FracFixD on macOS (arm64) — verified recipe

> **NOTE.** The FracFixD source code is closed and confidential
> (see [`../README.md` → Source code](../README.md#source-code)).
> This document is the **verified macOS build recipe** that
> produced the `FracFixD-2.0.0-macos-arm64.dmg` shipped in
> [`../bin/`](../bin/).  It was validated end-to-end on a
> macincloud Apple Silicon build host running macOS 26 (Tahoe).
>
> If you have access to the pre-built `.dmg` you can skip this
> document entirely.

The build is two steps: `dub build` produces a Mach-O binary,
then [`../installer/package-macos.sh`](../installer/package-macos.sh)
bundles the GTK 3 / OpenBLAS / gfortran dylibs into a relocatable
tree, wraps it in an `.app` bundle (with a native Mach-O launcher
from [`../installer/applauncher.c`](../installer/applauncher.c)
for Tahoe-Gatekeeper compatibility), and packages the whole
thing as a `.dmg`.

---

## Step 1 — Install prerequisites

### Hardware

- Apple Silicon (M1 / M2 / M3 / M4) Mac running macOS 13
  (Ventura) or newer.  macincloud "Apple Silicon Build"
  tier is sufficient.
- 8 GB RAM, 10 GB free disk.

### Apple toolchain

```bash
xcode-select --install     # Command Line Tools (clang, ld, codesign, install_name_tool)
```

### Homebrew + GTK 3 + OpenBLAS + LAPACK + dylibbundler

```bash
# 1. Homebrew (skip if installed)
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# 2. Add brew to PATH for the session
eval "$(/opt/homebrew/bin/brew shellenv)"

# 3. Build dependencies
brew install dub gtk+3 openblas lapack dylibbundler pkg-config
```

`lapack` is keg-only on Homebrew (macOS provides LAPACK via the
Accelerate framework).  We still install the keg-only formula so
that `dub`'s `-lopenblas -llapack` link flags resolve against
brew's BLAS/LAPACK pair rather than against Accelerate — keeps
the dylib closure walkable by `dylibbundler`.

### LDC 1.42 for arm64-darwin

Brew's `ldc` formula lags occasionally; use the official tarball
for predictability:

```bash
cd ~
curl -fsSLO https://github.com/ldc-developers/ldc/releases/download/v1.42.0/ldc2-1.42.0-osx-arm64.tar.xz
tar -xJf ldc2-1.42.0-osx-arm64.tar.xz
mv ldc2-1.42.0-osx-arm64 ldc2
export PATH=$HOME/ldc2/bin:$PATH
ldc2 --version       # should report 1.42.0
```

(Add `export PATH=$HOME/ldc2/bin:$PATH` to `~/.zshrc` for
future shells.)

---

## Step 2 — Get the FracFixD source

The source tree is closed and provided under a separate licence.
Extract it to your home directory:

```bash
cd ~
# Extract the source-tree archive here.
# Result:  ~/FracFixD/dub.json
# Result:  ~/FracFixD/source/
# Result:  ~/FracFixD/installer/
```

---

## Step 3 — Build the binary

```bash
eval "$(/opt/homebrew/bin/brew shellenv)"
export PATH=$HOME/ldc2/bin:$PATH
export LDFLAGS="-L/opt/homebrew/opt/openblas/lib -L/opt/homebrew/opt/lapack/lib"

cd ~/FracFixD
dub build --config=macos --build=release-static --compiler=ldc2 --force
```

The `macos` configuration in `dub.json` uses:

```json
"libs": ["openblas", "lapack"],
"excludedSourceFiles": ["source/benchmark.d"],
"lflags": [
    "-L/opt/homebrew/lib",
    "-framework", "CoreFoundation"
]
```

The `release-static` buildType statically links the D runtime
(`--link-defaultlib-shared=false`).  GTK 3 / OpenBLAS / gfortran
stay dynamic and get bundled by the packaging step below.

Verify the build:

```bash
file ./fracfixd
# Should report: Mach-O 64-bit executable arm64

DYLD_LIBRARY_PATH=/opt/homebrew/lib ./fracfixd --cli --version
# Should print the v2.0.0 "Quokka" banner.
```

---

## Step 4 — Bundle dylibs + build the `.app` and `.dmg`

```bash
cd ~/FracFixD
chmod +x installer/package-macos.sh
bash installer/package-macos.sh
```

The script does five things, in order:

1. **Relocatable tree.**  Copies the binary to
   `dist/fracfixd-macos/bin/fracfixd` and the top-level GTK 3 /
   OpenBLAS / gfortran dylibs (including `libgcc_s.1.1.dylib`,
   `libgfortran.5.dylib`, `libquadmath.0.dylib`) to
   `dist/fracfixd-macos/lib/`.
2. **`dylibbundler` walk.**  Walks every transitive dependency
   of the binary AND each top-level dylib, copies all referenced
   dylibs to `lib/`, and rewrites their install names to
   `@executable_path/../lib/`.  After this pass there are ~50
   dylibs in `lib/`.
3. **LC_RPATH dedupe.**  `dylibbundler` appends `LC_RPATH
   @executable_path/../lib/` every time it touches a binary,
   producing N duplicates after N passes.  Duplicate LC_RPATHs
   cause `dyld` to refuse to load the dylib with an error of
   the form `tried: <path> (duplicate LC_RPATH '<rpath>')`.
   The script counts duplicates with `otool -l`, strips them
   with `install_name_tool -delete_rpath` until only one
   remains, re-adds the canonical entry, and re-signs the
   dylib.  This is the only FracFixD-specific deviation from
   the TagGen recipe — it is needed because OpenBLAS pulls in
   many transitive dylibs that dylibbundler processes multiple
   times.
4. **`.app` bundle with native Mach-O launcher.**  macOS
   Sequoia / Tahoe Gatekeeper rejects shell-script launchers in
   `Contents/MacOS/` (error `-10669` even after ad-hoc
   codesign).  The script compiles
   [`../installer/applauncher.c`](../installer/applauncher.c)
   to `Contents/MacOS/FracFixD` — a 51 KB native Mach-O that
   sets `DYLD_LIBRARY_PATH` to the bundle's `lib/` and execs
   the real binary.
5. **`.dmg` packaging.**  Stages the `.app` plus a `/Applications`
   symlink in a temp dir and builds a UDZO-compressed `.dmg`
   via `hdiutil`.  Ad-hoc signs the `.app` so `codesign --verify`
   reports a usable signature.

Expected output on success:

```
Bundle ready:    /Users/admin/FracFixD/dist/fracfixd-macos  (55M)
.app ready:      /Users/admin/FracFixD/dist/FracFixD.app     (56M)
.dmg ready:      /Users/admin/FracFixD/dist/FracFixD-2.0.0-macos-arm64.dmg     (21M)
```

---

## Step 5 — Verify the build

```bash
# Via the .app's native launcher (what end-users will run)
open ~/FracFixD/dist/FracFixD.app
# (GUI launches; close it)

# Via the .app's internal CLI binary
DYLD_LIBRARY_PATH=~/FracFixD/dist/FracFixD.app/Contents/Resources/lib \
    ~/FracFixD/dist/FracFixD.app/Contents/Resources/bin/fracfixd --cli --version

# Via the relocatable bundle's launcher script (CLI use)
~/FracFixD/dist/fracfixd-macos/bin/fracfixd-launcher.sh --cli --version

# Verify the .dmg
hdiutil verify ~/FracFixD/dist/FracFixD-2.0.0-macos-arm64.dmg

# Hash for SHA256SUMS
shasum -a 256 ~/FracFixD/dist/FracFixD-2.0.0-macos-arm64.dmg
```

The `--cli --version` output should display the standard banner:

```
+----------------------------------------------------------------+
|   Version:   2.0.0                                             |
|   Codename:  Quokka                                            |
|   Status:    Stable                                            |
+----------------------------------------------------------------+
```

---

## Step 6 — Pull the artefacts back to the release host

If you're building remotely (e.g. on macincloud), pull both the
`.dmg` and the relocatable `fracfixd-macos/` tree back to the
release host:

```bash
# On the macincloud host (where the build ran)
cd ~/FracFixD
tar -czf fracfixd-2.0.0-macos-arm64.tar.gz dist/fracfixd-macos

# On the release host
scp builder@DXSnnn.macincloud.com:FracFixD/dist/FracFixD-2.0.0-macos-arm64.dmg ./
scp builder@DXSnnn.macincloud.com:FracFixD/fracfixd-2.0.0-macos-arm64.tar.gz ./
```

Drop both into `FracFixD/bin/` of the release repo and update
`SHA256SUMS`:

```bash
cd FracFixD/bin
sha256sum * > SHA256SUMS
```

---

## End-user experience

The user mounts the `.dmg` and drags `FracFixD.app` into the
`/Applications` symlink.  On first launch macOS Gatekeeper
warns that the app is from an unidentified developer (ad-hoc
signed, not Developer-ID notarised).  The standard one-time
workaround is right-click → Open, which whitelists the
signature.  Or via Terminal:

```bash
xattr -dr com.apple.quarantine /Applications/FracFixD.app
```

CLI users can either invoke the `.app`'s internal binary
directly:

```bash
/Applications/FracFixD.app/Contents/Resources/bin/fracfixd-launcher.sh --cli --version
```

Or extract the relocatable tarball:

```bash
tar -xzf fracfixd-2.0.0-macos-arm64.tar.gz
./fracfixd-macos/bin/fracfixd-launcher.sh --cli --version
```

---

## Troubleshooting

### `dyld[NNN]: Library not loaded: ... (duplicate LC_RPATH ...)`

The LC_RPATH dedupe step did not run.  Re-run
`installer/package-macos.sh` — the script's RPATH-dedupe loop
will fix the affected dylib.

### `Library not loaded: libopenblas.0.dylib`

OpenBLAS isn't installed, or the build linker didn't see
`-L/opt/homebrew/opt/openblas/lib`.  Set the export:

```bash
export LDFLAGS="-L/opt/homebrew/opt/openblas/lib -L/opt/homebrew/opt/lapack/lib"
```

and re-run `dub build --force`.

### `_LSOpenURLsWithCompletionHandler() failed with error -10669`

macOS Tahoe (15+) rejects shell-script launchers in
`Contents/MacOS/`.  The script handles this by compiling
`applauncher.c` to a native Mach-O — but if you're running an
older version of the script, the launcher might still be a
shell wrapper.  Re-pull the latest `installer/package-macos.sh`
and re-run.

### `codesign --verify` reports "no usable signature"

The `.app` was created but ad-hoc signing failed.  Run
`codesign --force --deep --sign - dist/FracFixD.app` manually.

### Mac is M-series but `file` reports `Mach-O ... x86_64`

You built on a Rosetta x86_64 shell.  Open a native arm64
terminal (`arch` should report `arm64`) and rebuild.

---

## Reference: macincloud build-host snapshot

The reference build host used for the v2.0.0 release-engineering
cycle:

| Component | Version |
|---|---|
| macOS | 26.2 (Tahoe) |
| Xcode CLT | (clang 21.0.0, ld-prime) |
| Homebrew | 5.x |
| GTK 3 | 3.24.52 |
| OpenBLAS | 0.3.33 |
| LAPACK | (Apple Accelerate via the brew keg-only formula) |
| dylibbundler | 1.0.5 |
| DUB | 1.41.0 |
| LDC | 1.42.0 (DMD 2.112.1, LLVM 21.1.8) |

The recipe is known to work end-to-end on this exact stack.
Newer macOS releases and brew versions usually work as well;
report regressions to Nick.
