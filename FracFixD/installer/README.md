# FracFixD installer / packaging templates

This directory holds the per-platform packaging scripts and
installer source files for FracFixD.  They are released under
the **MIT** licence so the same recipes can be reused for
downstream redistributions.

| File | Purpose | Status (v2.0.0) |
|---|---|---|
| `package-linux.sh` | Linux strip + sha256 + staging script | Stub — see `../docs/BUILD_LINUX.md` for the canonical workflow |
| `package-macos.sh` | macOS dylib-bundle + `.app` + `.dmg` builder | Template — adapt from TagGen v1.2.5 `installer/package-macos.sh` |
| `applauncher.c`    | Native Mach-O launcher for the `.app` bundle (Tahoe-Gatekeeper-compatible) | Template — adapt from TagGen v1.2.5 `installer/applauncher.c` |
| `package-windows.sh` | MSYS2 + DLL-collection + ZIP builder | Template — adapt from TagGen v1.2.5 `installer/package-windows.sh` |
| `build-windows.bat` | All-in-one Windows build + ZIP + (optional) Inno Setup driver | Template — adapt from TagGen v1.2.5 `build-windows.bat` |
| `fracfixd-installer.iss` | Inno Setup script (creates `.exe` installer) | Template — adapt from TagGen v1.2.5 `installer/taggen-installer.iss` |

For **v2.0.0 "Quokka"** the Linux release is the only platform
shipped pre-built.  macOS and Windows scaffolds in this
directory are provided as a starting point for follow-up
releases (and for source-tree licensees who want to produce
their own builds today).

The complete TagGen v1.2.5 reference implementations are at
`https://github.com/Arnaroo/taggen` (specifically the
`installer/` directory and `build-windows.bat`).  Adapting
them to FracFixD is mechanical: change `taggen` → `fracfixd`,
update the version string, switch the resource filenames, and
keep the dylib/DLL collection logic untouched.

For the canonical build walk-throughs see:

- [`../docs/BUILD_LINUX.md`](../docs/BUILD_LINUX.md)
- [`../docs/BUILD_MACOS.md`](../docs/BUILD_MACOS.md)
- [`../docs/BUILD_WINDOWS.md`](../docs/BUILD_WINDOWS.md)
