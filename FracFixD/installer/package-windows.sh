#!/bin/bash
# Package FracFixD for Windows installer / portable ZIP
# Run from MSYS2 MINGW64 terminal after building fracfixd.exe
#
# Usage:
#   ./installer/package-windows.sh           # from project root
#   ./installer/package-windows.sh --zip     # also create ZIP archive
#
# This script collects fracfixd.exe, GTK DLLs, themes, icon loaders,
# and resources into dist/fracfixd-windows/ — ready for Inno Setup
# or standalone distribution.
#
# Adapted from TagGen v1.2.5 installer/package-windows.sh.
# License: MIT

set -e

# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------
DIST_DIR="dist/fracfixd-windows"
MINGW_PREFIX="/mingw64"
PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"

cd "$PROJECT_ROOT"

# Version: prefer env var, fall back to version_.d, then CHANGELOG.md, then default
if [ -z "$VERSION" ]; then
    if [ -f "source/version_.d" ]; then
        VMAJ=$(grep -m1 -E '^enum VERSION_MAJOR' source/version_.d 2>/dev/null | grep -oE '[0-9]+')
        VMIN=$(grep -m1 -E '^enum VERSION_MINOR' source/version_.d 2>/dev/null | grep -oE '[0-9]+')
        VPAT=$(grep -m1 -E '^enum VERSION_PATCH' source/version_.d 2>/dev/null | grep -oE '[0-9]+')
        [ -n "$VMAJ" ] && [ -n "$VMIN" ] && [ -n "$VPAT" ] && VERSION="${VMAJ}.${VMIN}.${VPAT}"
    fi
fi
if [ -z "$VERSION" ] && [ -f "CHANGELOG.md" ]; then
    VERSION=$(grep -m1 -oE 'v[0-9]+\.[0-9]+\.[0-9]+' CHANGELOG.md | tr -d 'v' | head -n 1)
fi
[ -z "$VERSION" ] && VERSION="0.0.0"
echo "Packaging FracFixD $VERSION for Windows"

# ------------------------------------------------------------------
# Verify executable exists
# ------------------------------------------------------------------
if [ ! -f "fracfixd.exe" ]; then
    echo "ERROR: fracfixd.exe not found in project root."
    echo "Build first:  dub build --compiler=ldc2 --config=windows-gui -b release-static"
    exit 1
fi

echo "=== FracFixD Windows Packaging ==="
echo ""

# ------------------------------------------------------------------
# Clean and create staging directory
# ------------------------------------------------------------------
rm -rf "$DIST_DIR"
mkdir -p "$DIST_DIR"

# ------------------------------------------------------------------
# 1. Copy executable
# ------------------------------------------------------------------
echo "[1/6] Copying fracfixd.exe..."
cp fracfixd.exe "$DIST_DIR/"

# ------------------------------------------------------------------
# 2. Copy GTK DLLs (all transitive dependencies)
# ------------------------------------------------------------------
echo "[2/6] Collecting GTK + scientific DLLs..."

# Roots from which transitive dependencies are walked.
ROOT_DLLS=(
    libgtk-3-0.dll libgdk-3-0.dll libgdk_pixbuf-2.0-0.dll
    libgio-2.0-0.dll libglib-2.0-0.dll libgobject-2.0-0.dll libgmodule-2.0-0.dll
    libcairo-2.dll libcairo-gobject-2.dll
    libpango-1.0-0.dll libpangocairo-1.0-0.dll libpangowin32-1.0-0.dll libpangoft2-1.0-0.dll
    libatk-1.0-0.dll libharfbuzz-0.dll libfribidi-0.dll libepoxy-0.dll
    librsvg-2-2.dll
    # Scientific stack — FracFixD-specific
    libopenblas.dll liblapack.dll libgfortran-5.dll libquadmath-0.dll
)

# Method 1: Walk transitive deps via ldd.
for root_dll in "${ROOT_DLLS[@]}"; do
    if [ -f "$MINGW_PREFIX/bin/$root_dll" ]; then
        cp -n "$MINGW_PREFIX/bin/$root_dll" "$DIST_DIR/" 2>/dev/null || true
        ldd "$MINGW_PREFIX/bin/$root_dll" 2>/dev/null | grep -i mingw64 | awk '{print $3}' | sort -u | while read dll; do
            if [ -f "$dll" ]; then
                cp -n "$dll" "$DIST_DIR/" 2>/dev/null || true
            fi
        done
    fi
done

# Walk transitive deps of fracfixd.exe itself for any DLLs we missed.
ldd fracfixd.exe 2>/dev/null | grep -i mingw64 | awk '{print $3}' | sort -u | while read dll; do
    if [ -f "$dll" ]; then
        cp -n "$dll" "$DIST_DIR/" 2>/dev/null || true
    fi
done

# Method 2: Belt-and-braces — explicit copy of essential support DLLs.
EXTRA_DLLS=(
    libgcc_s_seh-1.dll libwinpthread-1.dll libstdc++-6.dll
    libintl-8.dll libiconv-2.dll libffi-8.dll libpcre2-8-0.dll
    libpng16-16.dll zlib1.dll libbz2-1.dll libexpat-1.dll
    libbrotlidec.dll libbrotlicommon.dll libfontconfig-1.dll libfreetype-6.dll
    libpixman-1-0.dll libgraphite2.dll libthai-0.dll libdatrie-1.dll
    libjpeg-8.dll libtiff-6.dll libdeflate.dll libLerc.dll
    liblzma-5.dll libzstd.dll libjbig-0.dll libwebp-7.dll libsharpyuv-0.dll
    # Scientific stack
    libopenblas.dll liblapack.dll libgfortran-5.dll libquadmath-0.dll
)
for dll in "${EXTRA_DLLS[@]}"; do
    if [ -f "$MINGW_PREFIX/bin/$dll" ] && [ ! -f "$DIST_DIR/$dll" ]; then
        cp -n "$MINGW_PREFIX/bin/$dll" "$DIST_DIR/" 2>/dev/null || true
    fi
done

DLL_COUNT=$(ls -1 "$DIST_DIR"/*.dll 2>/dev/null | wc -l)
echo "       Copied $DLL_COUNT DLLs"
if [ "$DLL_COUNT" -lt 20 ]; then
    echo ""
    echo "  WARNING: Only $DLL_COUNT DLLs collected (expected ~70-80)."
    echo "  Possible causes:"
    echo "    - This script must run inside MSYS2 MINGW64 bash (\$MSYSTEM=$MSYSTEM)"
    echo "    - MSYS2 GTK3 install may be incomplete: pacman -S mingw-w64-x86_64-gtk3"
    echo "    - OpenBLAS missing: pacman -S mingw-w64-x86_64-openblas mingw-w64-x86_64-lapack"
    echo "    - \$MINGW_PREFIX=$MINGW_PREFIX (should be /mingw64)"
    echo ""
fi

# ------------------------------------------------------------------
# 3. Copy GdkPixbuf loaders (needed for image rendering)
# ------------------------------------------------------------------
echo "[3/6] Copying GdkPixbuf loaders..."
PIXBUF_DIR="$MINGW_PREFIX/lib/gdk-pixbuf-2.0"
if [ -d "$PIXBUF_DIR" ]; then
    mkdir -p "$DIST_DIR/lib/gdk-pixbuf-2.0"
    cp -r "$PIXBUF_DIR"/* "$DIST_DIR/lib/gdk-pixbuf-2.0/"
    # Update the loaders cache to use relative paths
    if [ -f "$DIST_DIR/lib/gdk-pixbuf-2.0/2.10.0/loaders.cache" ]; then
        sed -i "s|$MINGW_PREFIX/lib/gdk-pixbuf-2.0/2.10.0/loaders/|lib/gdk-pixbuf-2.0/2.10.0/loaders/|g" \
            "$DIST_DIR/lib/gdk-pixbuf-2.0/2.10.0/loaders.cache"
    fi
    echo "       Done"
else
    echo "       WARNING: GdkPixbuf loaders not found at $PIXBUF_DIR"
fi

# ------------------------------------------------------------------
# 4. Copy GTK themes and icons
# ------------------------------------------------------------------
echo "[4/6] Copying GTK themes and icons..."
mkdir -p "$DIST_DIR/share"

if [ -d "$MINGW_PREFIX/share/themes/Default" ]; then
    mkdir -p "$DIST_DIR/share/themes"
    cp -r "$MINGW_PREFIX/share/themes/Default" "$DIST_DIR/share/themes/"
fi
if [ -d "$MINGW_PREFIX/share/themes/MS-Windows" ]; then
    cp -r "$MINGW_PREFIX/share/themes/MS-Windows" "$DIST_DIR/share/themes/"
fi

if [ -d "$MINGW_PREFIX/share/icons/Adwaita" ]; then
    mkdir -p "$DIST_DIR/share/icons"
    cp -r "$MINGW_PREFIX/share/icons/Adwaita" "$DIST_DIR/share/icons/"
fi
if [ -d "$MINGW_PREFIX/share/icons/hicolor" ]; then
    cp -r "$MINGW_PREFIX/share/icons/hicolor" "$DIST_DIR/share/icons/"
fi

if [ -d "$MINGW_PREFIX/share/glib-2.0/schemas" ]; then
    mkdir -p "$DIST_DIR/share/glib-2.0"
    cp -r "$MINGW_PREFIX/share/glib-2.0/schemas" "$DIST_DIR/share/glib-2.0/"
fi

echo "       Done"

# ------------------------------------------------------------------
# 5. Copy application resources and docs
# ------------------------------------------------------------------
echo "[5/6] Copying resources and documentation..."
[ -d "resources" ] && cp -r resources "$DIST_DIR/"
[ -f "README.md" ] && cp README.md "$DIST_DIR/"
[ -f "LICENSE" ]   && cp LICENSE   "$DIST_DIR/LICENSE.txt"
[ -f "LICENSE.txt" ] && cp LICENSE.txt "$DIST_DIR/"
[ -f "CHANGELOG.md" ] && cp CHANGELOG.md "$DIST_DIR/"
echo "       Done"

# ------------------------------------------------------------------
# 6. Summary
# ------------------------------------------------------------------
echo "[6/6] Calculating sizes..."
EXE_SIZE=$(du -h "$DIST_DIR/fracfixd.exe" | cut -f1)
TOTAL_SIZE=$(du -sh "$DIST_DIR" | cut -f1)
DLL_COUNT=$(ls -1 "$DIST_DIR"/*.dll 2>/dev/null | wc -l)

echo ""
echo "=== Packaging complete ==="
echo "  Staging directory: $DIST_DIR/"
echo "  Executable:        $EXE_SIZE"
echo "  DLLs:              $DLL_COUNT files"
echo "  Total size:        $TOTAL_SIZE"
echo ""

# ------------------------------------------------------------------
# Optional: Create ZIP archive
# ------------------------------------------------------------------
if [[ "$1" == "--zip" ]]; then
    echo "Creating ZIP archive..."
    cd dist
    ZIP_NAME="fracfixd-v${VERSION}-windows-x86_64.zip"
    rm -f "$ZIP_NAME"
    zip -qr "$ZIP_NAME" "fracfixd-windows/"
    ZIP_SIZE=$(du -h "$ZIP_NAME" | cut -f1)
    echo "  ZIP: dist/$ZIP_NAME ($ZIP_SIZE)"
    cd ..
    echo ""
fi

# ------------------------------------------------------------------
# Next steps
# ------------------------------------------------------------------
echo "Next steps:"
echo ""
echo "  Create installer (requires Inno Setup):"
echo "    iscc installer/fracfixd-installer.iss"
echo ""
echo "  Or create ZIP only:"
echo "    $0 --zip"
echo ""
echo "  Test the staging directory directly:"
echo "    $DIST_DIR/fracfixd.exe"
echo ""
