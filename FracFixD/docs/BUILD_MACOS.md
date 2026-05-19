# Building FracFixD on macOS (arm64)

> **NOTE.** The FracFixD source code is closed and confidential
> (see [`../README.md` → Source code](../README.md#source-code)).
> This document describes the **macOS build recipe** as a
> provenance reference and for licensed source-tree recipients.
> A pre-built macOS arm64 binary is **planned** for a follow-up
> release of FracFixD; until it lands, this walk-through is the
> path for adventurous source-tree users.

The recipe below mirrors the TagGen v1.2.5 macOS release
process and has been validated against that companion project
on macincloud arm64 hosts.

---

## 1. Prerequisites

### Hardware

- Apple Silicon (M1 / M2 / M3 / M4) Mac running macOS 13
  (Ventura) or newer.
- Recommended: a macincloud "Apple Silicon Build" tier or a
  physical Mac with at least 8 GB RAM.

### Toolchain

```bash
# 1. Homebrew (skip if installed)
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# 2. Build tools, GTK 3, dylib bundler, pkg-config
brew install dub gtk+3 dylibbundler pkg-config

# 3. LDC 1.42 (or matching) for arm64-darwin
# Use the official tarball for predictability (Homebrew's ldc lags occasionally).
curl -fsSLO https://github.com/ldc-developers/ldc/releases/download/v1.42.0/ldc2-1.42.0-osx-arm64.tar.xz
tar -xJf ldc2-1.42.0-osx-arm64.tar.xz
sudo mv ldc2-1.42.0-osx-arm64 /opt/ldc

# Add to PATH for the session:
export PATH=/opt/ldc/bin:$PATH
ldc2 --version       # confirm 1.42.0
```

### Source-tree access

See [`../README.md` → Source code](../README.md#source-code).
The build commands below assume you have the source tree
checked out and that `dub.json` ships the `macos`
configuration (it does in the public-release v2.0.0 tree).

---

## 2. Build the GUI binary

Build configuration `macos` is the GTK3 + GUI configuration
for Apple Silicon.  The `release-static` buildType is the
right one — it adds the `-Wl,-force_load` flags for the D
runtime archives which Apple's `ld-prime` linker needs.

```bash
LIBRARY_PATH=build/lib-shim dub build \
    --config=macos --build=release-static --compiler=/opt/ldc/bin/ldc2 --force
# Produces `./fracfixd` (Mach-O 64-bit executable, arm64).
```

The `macos` configuration in `dub.json` already includes:

```json
"lflags": [
    "-L/opt/homebrew/lib",
    "-framework", "CoreFoundation"
]
```

For full reproducibility, also add the static D-runtime
force-load (already present in the public-release `dub.json`
`macos` variant; mention here so source-tree forks know to
keep them):

```json
"-Wl,-force_load,/opt/ldc/lib/libdruntime-ldc.a",
"-Wl,-force_load,/opt/ldc/lib/libphobos2-ldc.a"
```

These are **essential**: Apple's new linker (`ld-prime`) will
not pull weak external symbols out of static D archives
unless forced.

---

## 3. Bundle the GTK dylibs

Run the in-tree packaging script (lives in
[`../installer/package-macos.sh`](../installer/package-macos.sh)
on a configured release-host source tree; structure mirrors
TagGen's `installer/package-macos.sh`).  Logical steps:

```bash
# Output dir layout:
#   fracfixd-macos-arm64/
#     bin/fracfixd                    (the executable)
#     bin/fracfixd-launcher.sh        (sets DYLD_LIBRARY_PATH)
#     lib/*.dylib                     (transitive GTK3 closure, ~48 dylibs)

mkdir -p fracfixd-macos-arm64/{bin,lib}
cp ./fracfixd fracfixd-macos-arm64/bin/

# Copy the top-level GTK / GLib / GdkPixbuf / Cairo / Pango dylibs
cp /opt/homebrew/lib/libgtk-3.0.dylib       fracfixd-macos-arm64/lib/
cp /opt/homebrew/lib/libgdk-3.0.dylib       fracfixd-macos-arm64/lib/
cp /opt/homebrew/lib/libgdk_pixbuf-2.0.0.dylib fracfixd-macos-arm64/lib/
cp /opt/homebrew/lib/libcairo.2.dylib       fracfixd-macos-arm64/lib/
cp /opt/homebrew/lib/libpango-1.0.0.dylib   fracfixd-macos-arm64/lib/
cp /opt/homebrew/lib/libglib-2.0.0.dylib    fracfixd-macos-arm64/lib/
cp /opt/homebrew/lib/libgobject-2.0.0.dylib fracfixd-macos-arm64/lib/
# ... (full list of 12-15 top-level dylibs)

# Walk transitive deps and rewrite install names to @executable_path/../lib/
dylibbundler -of -b \
    -x fracfixd-macos-arm64/bin/fracfixd \
    -d fracfixd-macos-arm64/lib/ \
    -p @executable_path/../lib/
# ~48 dylibs end up in lib/ after the walk.
```

Wrapper script (`bin/fracfixd-launcher.sh`):

```sh
#!/bin/sh
DIR="$(cd "$(dirname "$0")" && pwd)"
export DYLD_LIBRARY_PATH="$DIR/../lib:${DYLD_LIBRARY_PATH:-}"
exec "$DIR/fracfixd" "$@"
```

CLI users can invoke `./fracfixd-macos-arm64/bin/fracfixd-launcher.sh`
directly; the `.app` bundle and `.dmg` below are for GUI users.

---

## 4. `.app` bundle for GUI users

Tahoe Gatekeeper rejects shell scripts in `CFBundleExecutable`
— `.app` must contain a native Mach-O launcher.  Compile the
launcher from [`../installer/applauncher.c`](../installer/applauncher.c):

```bash
clang -O2 -arch arm64 -o FracFixD.app/Contents/MacOS/FracFixD \
    installer/applauncher.c \
    -DTARGET_BIN=\"../Resources/bin/fracfixd\" \
    -DTARGET_LIB=\"../Resources/lib\"
```

The `Info.plist` should declare:

```xml
<key>CFBundleIdentifier</key>          <string>com.arnaroo.fracfixd</string>
<key>CFBundleName</key>                <string>FracFixD</string>
<key>CFBundleDisplayName</key>         <string>FracFixD</string>
<key>CFBundleVersion</key>             <string>2.0.0</string>
<key>CFBundleShortVersionString</key>  <string>2.0.0</string>
<key>CFBundleExecutable</key>          <string>FracFixD</string>
<key>CFBundleIconFile</key>            <string>fracfixd-icon.icns</string>
<key>LSMinimumSystemVersion</key>      <string>13.0</string>
```

Then drop the relocatable tree into the bundle:

```bash
mkdir -p FracFixD.app/Contents/Resources/bin FracFixD.app/Contents/Resources/lib
cp fracfixd-macos-arm64/bin/fracfixd      FracFixD.app/Contents/Resources/bin/
cp fracfixd-macos-arm64/lib/*.dylib       FracFixD.app/Contents/Resources/lib/
cp resources/fracfixd-icon.icns           FracFixD.app/Contents/Resources/
```

Ad-hoc sign:

```bash
codesign --force --deep --sign - FracFixD.app
```

---

## 5. `.dmg` packaging

```bash
# Create a temp staging dir with the .app + Applications symlink
STAGING=$(mktemp -d)
cp -R FracFixD.app "$STAGING/"
ln -s /Applications "$STAGING/Applications"

# Build the DMG
hdiutil create -volname "FracFixD 2.0.0" \
               -srcfolder "$STAGING" \
               -ov -format UDZO \
               -fs HFS+ \
               FracFixD-2.0.0-macos-arm64.dmg
```

End-user experience: double-click `.dmg` → drag FracFixD to
Applications.  First launch requires right-click → Open to
clear the one-time Gatekeeper confirmation (ad-hoc signed,
not Developer-ID notarised).

---

## 6. Verification

```bash
codesign --verify --verbose FracFixD.app
spctl --assess --verbose FracFixD.app          # Will report "unsigned by trusted authority" — expected for ad-hoc

# Smoke test
./FracFixD.app/Contents/Resources/bin/fracfixd --version
./FracFixD.app/Contents/Resources/bin/fracfixd --cli help diffprop

# Hash the artefacts
shasum -a 256 FracFixD-2.0.0-macos-arm64.dmg
```

---

## 7. Known macOS-specific issues

- **Gatekeeper "FracFixD is damaged" message** — the ad-hoc
  sign is rejected by Gatekeeper on some MDM-managed Macs.
  Right-click → Open is the standard workaround.
- **GTK 3 deprecation warnings on stderr** — harmless; GTK 3
  is in long-term support on macOS via Homebrew but emits
  deprecation noise.  We do not surface this to the user.
- **HiDPI scaling** — the GTK 3 backend on Apple Silicon
  reports the logical (post-HiDPI) screen size; pixel-perfect
  SVG output is unaffected because rendering is vector.

---

## 8. Sign-and-notarise (optional, not used in v2.0.0)

For a fully Gatekeeper-clean distribution, swap the ad-hoc
sign for a Developer ID Application certificate and run
through `notarytool`:

```bash
codesign --force --deep --sign "Developer ID Application: <Your Name>" \
         --options runtime --timestamp \
         --entitlements installer/fracfixd.entitlements \
         FracFixD.app

ditto -c -k --keepParent FracFixD.app FracFixD.app.zip

xcrun notarytool submit FracFixD.app.zip \
    --apple-id you@example.com \
    --team-id ABCDE12345 \
    --password "@keychain:notary-password" \
    --wait

xcrun stapler staple FracFixD.app
```

This is currently outside the v2.0.0 scope; the ad-hoc-signed
DMG is sufficient for academic / research distribution.
