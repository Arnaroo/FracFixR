# FracFixD installer / packaging scripts

This directory holds the per-platform packaging scripts and
installer source files for FracFixD.  They are released under
the **MIT** licence so the same recipes can be reused for
downstream redistributions.

| File | Purpose | Status (v2.0.0) |
|---|---|---|
| `package-linux.sh` | Linux strip + sha256 + staging script | **Active** — used to stage the four Linux binaries shipped in `../bin/` |
| `package-macos.sh` | macOS dylib-bundle + `.app` + `.dmg` builder | **Active** — produced `FracFixD-2.0.0-macos-arm64.dmg` shipped in `../bin/` |
| `applauncher.c`    | Native Mach-O launcher for the `.app` bundle (Tahoe-Gatekeeper-compatible) | **Active** — compiled in by `package-macos.sh` |
| `package-windows.sh` | MSYS2 + DLL-collection + ZIP builder | **Active** — invoked by `build-windows.bat` |
| `build-windows.bat` | All-in-one Windows build + ZIP + (optional) Inno Setup driver | **Active** — see `../docs/BUILD_WINDOWS.md` for step-by-step instructions |
| `fracfixd-installer.iss` | Inno Setup script (creates `.exe` installer) | **Active** — invoked when `build-windows.bat --installer` is passed |

For the **canonical build walk-throughs** see:

- [`../docs/BUILD_LINUX.md`](../docs/BUILD_LINUX.md) — produces
  three microarch-tuned Linux binaries plus a static CLI variant.
- [`../docs/BUILD_MACOS.md`](../docs/BUILD_MACOS.md) —
  macincloud-verified end-to-end recipe; produces a relocatable
  bundle, an `.app`, and a `.dmg`.
- [`../docs/BUILD_WINDOWS.md`](../docs/BUILD_WINDOWS.md) —
  step-by-step recipe for licensed source-tree recipients;
  produces a portable ZIP and an Inno Setup `.exe`.

These scripts were adapted from the TagGen v1.2.5 release
pipeline (`https://github.com/Arnaroo/taggen`) — the macOS
script gained a `LC_RPATH` dedupe step (needed because OpenBLAS
pulls in many transitive dylibs that `dylibbundler` processes
multiple times), the Windows scripts grew explicit OpenBLAS /
LAPACK / gfortran DLL lists, and `applauncher.c` is reused
verbatim.
