#!/bin/bash
# package-macos.sh - Bundle FracFixD binary with GTK+3 dylibs for relocatable
# distribution on macOS arm64. Produces:
#   dist/fracfixd-macos/{bin,lib}     -- relocatable bundle (CLI use)
#   dist/FracFixD.app                  -- .app bundle (Finder use)
#   dist/FracFixD-<ver>-macos-arm64.dmg -- drag-to-Applications installer
#
# Usage: bash installer/package-macos.sh   (run on macOS build host)
# Adapted from TagGen v1.2.5 installer/package-macos.sh.

set -e

SRC="${SRC:-$HOME/FracFixD}"
INSTALLER="$SRC/installer"
DIST="$SRC/dist"
BUNDLE="$DIST/fracfixd-macos"
BIN="$BUNDLE/bin"
LIB="$BUNDLE/lib"
APP="$DIST/FracFixD.app"
HBLIB="${HBLIB:-/opt/homebrew/lib}"
DYLIBBUNDLER="${DYLIBBUNDLER:-/opt/homebrew/bin/dylibbundler}"

# Pull version from version_.d (single-source-of-truth in the D tree).
VERSION="${VERSION:-$(grep -m1 -E '^enum VERSION_MAJOR' "$SRC/source/version_.d" 2>/dev/null | grep -oE '[0-9]+')}"
VMINOR=$(grep -m1 -E '^enum VERSION_MINOR' "$SRC/source/version_.d" 2>/dev/null | grep -oE '[0-9]+')
VPATCH=$(grep -m1 -E '^enum VERSION_PATCH' "$SRC/source/version_.d" 2>/dev/null | grep -oE '[0-9]+')
VERSION="${VERSION}.${VMINOR}.${VPATCH}"
[ -z "$VERSION" ] || [ "$VERSION" = ".." ] && VERSION="0.0.0"
DMG="$DIST/FracFixD-${VERSION}-macos-arm64.dmg"

# Locate the freshly-built binary.
if   [ -x "$SRC/fracfixd" ];     then SRC_EXE="$SRC/fracfixd"
elif [ -x "$SRC/bin/fracfixd" ]; then SRC_EXE="$SRC/bin/fracfixd"
else echo "ERROR: fracfixd executable not found in $SRC or $SRC/bin"; exit 1
fi
echo "Using executable: $SRC_EXE  (version $VERSION)"

#
# 1. Relocatable bin/lib bundle
#
rm -rf "$BUNDLE"
mkdir -p "$BIN" "$LIB"
cp "$SRC_EXE" "$BIN/fracfixd"

GTK_LIBS=(
  libatk-1.0.dylib libcairo.dylib libgdk-3.0.dylib libgdk_pixbuf-2.0.dylib
  libglib-2.0.dylib libgmodule-2.0.dylib libgobject-2.0.dylib libgio-2.0.dylib
  libgthread-2.0.dylib libgtk-3.0.dylib libpango-1.0.dylib libpangocairo-1.0.dylib
  libharfbuzz.dylib libfribidi.dylib libepoxy.dylib
)
for lib in "${GTK_LIBS[@]}"; do
  src=$(readlink -f "$HBLIB/$lib" 2>/dev/null || true)
  [ -f "$src" ] || { echo "  skip (missing): $lib"; continue; }
  base=$(basename "$src")
  cp -f "$src" "$LIB/$base"
  chmod +w "$LIB/$base"
  [ "$lib" != "$base" ] && ln -sf "$base" "$LIB/$lib"
done

# Also bundle openblas + gfortran runtime so the binary is portable
# (Mac LAPACK ships via Accelerate, but openblas was the build-time
# choice; copy its dylib so the user doesn't need brew installed).
for lib in libopenblas.dylib libgfortran.5.dylib libgcc_s.1.1.dylib \
           libquadmath.0.dylib; do
  src=$(readlink -f "$HBLIB/$lib" 2>/dev/null || \
        readlink -f /opt/homebrew/opt/openblas/lib/$lib 2>/dev/null || \
        readlink -f /opt/homebrew/opt/gcc/lib/gcc/*/$lib 2>/dev/null | head -1 || true)
  [ -f "$src" ] || continue
  base=$(basename "$src")
  cp -f "$src" "$LIB/$base"
  chmod +w "$LIB/$base"
  [ "$lib" != "$base" ] && ln -sf "$base" "$LIB/$lib"
done

XARGS=("-x" "$BIN/fracfixd")
for f in "$LIB"/*.dylib; do
  [ -f "$f" ] && [ ! -L "$f" ] && XARGS+=("-x" "$f")
done
"$DYLIBBUNDLER" -of -b "${XARGS[@]}" \
  -d "$LIB" \
  -p "@executable_path/../lib/" \
  -s "$HBLIB" \
  -s "/opt/homebrew/opt/openblas/lib" \
  -s "/opt/homebrew/opt/gcc/lib"

# Dedupe LC_RPATH entries — dylibbundler appends @executable_path/../lib/
# every time it touches a binary, producing N duplicates after N passes.
# Duplicate LC_RPATHs cause dyld to refuse to load the dylib with an
# error of the form: "tried: <path> (duplicate LC_RPATH '<rpath>')".
echo "Deduplicating LC_RPATH entries..."
for f in "$BIN/fracfixd" "$LIB"/*.dylib; do
  [ -f "$f" ] && [ ! -L "$f" ] || continue
  # Count occurrences of the rpath; if >1, delete all but one.
  count=$(otool -l "$f" | awk '/LC_RPATH/{getline; getline; print}' \
                       | grep -c "@executable_path/../lib/" || true)
  if [ "$count" -gt 1 ]; then
    # Remove all duplicates, then re-add one canonical entry.
    while [ "$count" -gt 0 ]; do
      install_name_tool -delete_rpath "@executable_path/../lib/" "$f" 2>/dev/null || break
      count=$((count - 1))
    done
    install_name_tool -add_rpath "@executable_path/../lib/" "$f"
    # Re-sign after install_name_tool surgery.
    codesign --force --sign - "$f" 2>/dev/null || true
  fi
done

# Shell-script launcher kept for CLI / scripted use (Terminal works fine).
cat > "$BIN/fracfixd-launcher.sh" <<'WRAP'
#!/bin/bash
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LIBDIR="$(cd "$HERE/../lib" && pwd)"
export DYLD_LIBRARY_PATH="$LIBDIR${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}"
exec "$HERE/fracfixd" "$@"
WRAP
chmod +x "$BIN/fracfixd-launcher.sh"

echo "Relocatable bundle:"
du -sh "$BUNDLE"

#
# 2. .app bundle with NATIVE Mach-O launcher
#
rm -rf "$APP"
mkdir -p "$APP/Contents/MacOS" "$APP/Contents/Resources"
cp -R "$BUNDLE/bin" "$APP/Contents/Resources/"
cp -R "$BUNDLE/lib" "$APP/Contents/Resources/"

cc -O2 -arch arm64 -mmacosx-version-min=13.0 \
   -DTARGET_BIN=\"fracfixd\" \
   "$INSTALLER/applauncher.c" -o "$APP/Contents/MacOS/FracFixD"
chmod +x "$APP/Contents/MacOS/FracFixD"

cat > "$APP/Contents/Info.plist" <<PLIST
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>CFBundleExecutable</key>            <string>FracFixD</string>
    <key>CFBundleIdentifier</key>            <string>org.biocodecs.fracfixd</string>
    <key>CFBundleName</key>                  <string>FracFixD</string>
    <key>CFBundleDisplayName</key>           <string>FracFixD</string>
    <key>CFBundleVersion</key>               <string>${VERSION}</string>
    <key>CFBundleShortVersionString</key>    <string>${VERSION}</string>
    <key>CFBundlePackageType</key>           <string>APPL</string>
    <key>CFBundleSignature</key>             <string>FFXD</string>
    <key>CFBundleInfoDictionaryVersion</key> <string>6.0</string>
    <key>LSMinimumSystemVersion</key>        <string>13.0</string>
    <key>NSHighResolutionCapable</key>       <true/>
    <key>LSApplicationCategoryType</key>     <string>public.app-category.education</string>
    <key>NSHumanReadableCopyright</key>      <string>Copyright (c) 2026 Biocodecs / Arnaroo Ribologicals / RMODEL. CC-BY-NC-ND-4.0.</string>
</dict>
</plist>
PLIST

# Ad-hoc sign so codesign --verify is happy.
codesign --force --deep --sign - "$APP" >/dev/null

echo ".app bundle: $APP"
DYLD_LIBRARY_PATH="$APP/Contents/Resources/lib" "$APP/Contents/Resources/bin/fracfixd" --cli --version 2>&1 | head -3 || true

#
# 3. .dmg installer
#
DMG_STAGE="$DIST/dmg_stage"
rm -rf "$DMG_STAGE" "$DMG"
mkdir -p "$DMG_STAGE"
cp -R "$APP" "$DMG_STAGE/"
ln -s /Applications "$DMG_STAGE/Applications"
cat > "$DMG_STAGE/README.txt" <<RM
FracFixD ${VERSION} - macOS arm64

Drag "FracFixD.app" into the "Applications" folder.

First launch: right-click FracFixD.app -> Open
(Gatekeeper warns since the build is ad-hoc-signed; one-time confirmation.)

Or remove the quarantine attribute from Terminal:
  xattr -dr com.apple.quarantine /Applications/FracFixD.app

CLI use (from Terminal, no GUI):
  /Applications/FracFixD.app/Contents/Resources/bin/fracfixd-launcher.sh --cli --version
RM

hdiutil create -volname "FracFixD ${VERSION}" -srcfolder "$DMG_STAGE" -ov \
  -format UDZO -fs HFS+ "$DMG" >/dev/null
hdiutil verify "$DMG" >/dev/null
rm -rf "$DMG_STAGE"

echo
echo "Bundle ready:    $BUNDLE  ($(du -sh "$BUNDLE" | awk '{print $1}'))"
echo ".app ready:      $APP     ($(du -sh "$APP"    | awk '{print $1}'))"
echo ".dmg ready:      $DMG     ($(du -sh "$DMG"    | awk '{print $1}'))"
