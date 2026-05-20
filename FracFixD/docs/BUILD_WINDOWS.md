# Building FracFixD on Windows (x86_64) — step-by-step

> **NOTE.** The FracFixD source code is closed and confidential
> (see [`../README.md` → Source code](../README.md#source-code)).
> This document is the **step-by-step Windows build recipe** for
> licensed source-tree recipients (Seva and others).  Following
> the steps below produces both a portable ZIP and an Inno Setup
> `.exe` installer with all GTK 3 + OpenBLAS / LAPACK / gfortran
> DLLs bundled.
>
> If you do not have access to the source tree, download the
> pre-built Linux / macOS binaries from
> [`../bin/`](../bin/) and skip this document.

The build is a single end-to-end script — [`../installer/build-windows.bat`](../installer/build-windows.bat)
— that handles everything: prerequisites check, DUB build, DLL
collection (via the MSYS2-side [`../installer/package-windows.sh`](../installer/package-windows.sh)
helper), GdkPixbuf-loader / GTK-theme staging, portable ZIP
creation, and optional Inno Setup `.exe` installer compilation.

---

## Step 1 — Install one-time prerequisites

You only need to do this once per machine.

### 1.1 MSYS2 + GTK 3

MSYS2's MinGW openblas ships in `.a` format, which the MSVC
linker LDC uses cannot consume.  We use MSYS2 only for the GTK
3 runtime stack; OpenBLAS comes from the official Windows
pre-build in Step 1.5 below (auto-downloaded by the build
script if missing).

1. Download MSYS2 from **https://www.msys2.org/** and run the
   installer.  Accept the default install path (`C:\msys64`).
2. Open the **"MSYS2 MINGW64"** shortcut from the Start Menu
   (NOT the plain "MSYS2 MSYS", they are different).
3. In the MINGW64 terminal, run:

```bash
pacman -Syu                                    # initial sync; close + reopen shell after
pacman -Syu                                    # second pass
pacman -S mingw-w64-x86_64-gtk3 \
          mingw-w64-x86_64-pkg-config \
          mingw-w64-x86_64-gcc \
          mingw-w64-x86_64-binutils \
          unzip wget zip
```

Verify GTK 3 is installed:

```bash
ls /mingw64/bin/libgtk-3-0.dll          # should exist (~10 MB)
```

If missing, re-run `pacman -S mingw-w64-x86_64-gtk3`.

### 1.5 OpenBLAS Windows pre-build (auto-installed)

The build script will auto-download the official OpenBLAS
Windows pre-build (https://github.com/OpenMathLib/OpenBLAS/releases)
to `C:\OpenBLAS\` if it isn't already there.  No manual install
needed.  If you want to pre-populate it yourself (e.g. on a
network-restricted machine):

1. Download
   `OpenBLAS-0.3.30-x64.zip` from the OpenBLAS releases page.
2. Extract to `C:\OpenBLAS\` so that
   `C:\OpenBLAS\lib\libopenblas.lib` and
   `C:\OpenBLAS\bin\libopenblas.dll` both exist.

The OpenBLAS Windows pre-build includes LAPACK in the same DLL,
so no separate LAPACK install is needed.

### 1.2 Visual Studio Build Tools

LDC's `-mtriple=x86_64-windows-msvc` requires the MSVC linker.

1. Download **Visual Studio 2022 Build Tools** from
   https://visualstudio.microsoft.com/visual-cpp-build-tools/.
2. Run the installer and select the **"Desktop development with
   C++"** workload.  The default-checked sub-components are fine.
3. After install, confirm by opening
   **"x64 Native Tools Command Prompt for VS 2022"** from the
   Start Menu and running:

```cmd
cl
```

You should see `Microsoft (R) C/C++ Optimizing Compiler Version ...`.
If you see "command not found" you opened the wrong terminal.

### 1.3 LDC2 D compiler

1. Download `ldc2-1.42.0-windows-multilib.7z` (or newer) from
   https://github.com/ldc-developers/ldc/releases.
2. Extract to `C:\D\ldc2` so that `C:\D\ldc2\bin\ldc2.exe` exists.
   (The build script also searches `C:\D\ldc2-*\bin\`.)
3. Verify by opening the **x64 Native Tools Command Prompt** and
   running:

```cmd
set PATH=%PATH%;C:\D\ldc2\bin
ldc2 --version
```

You should see `LDC - the LLVM D compiler (1.42.0): ...`.

### 1.4 Inno Setup 6 (optional — only for `.exe` installer)

If you only need the portable ZIP, skip this step.

```cmd
winget install JRSoftware.InnoSetup
```

Or download from https://jrsoftware.org/isinfo.php and run the
default installer.  Default install path is
`C:\Program Files (x86)\Inno Setup 6\`.

---

## Step 2 — Get the FracFixD source

The source tree is closed and provided under a separate licence
agreement.  Once you have the tree (typically as a tarball or
git bundle), extract it to a convenient location:

```cmd
cd C:\Repos
:: Extract the source-tree archive here
:: Result:  C:\Repos\FracFixD\dub.json
:: Result:  C:\Repos\FracFixD\source\...
:: Result:  C:\Repos\FracFixD\installer\
```

The directory structure should look like:

```
C:\Repos\FracFixD\
+-- dub.json
+-- dub.selections.json
+-- source\
|   +-- app.d
|   +-- cli.d
|   +-- version_.d
|   +-- gui_main.d
|   +-- gui_worker.d
|   +-- gui_state.d
|   +-- embedded_resources.d
|   +-- fracfix\
+-- resources\
|   +-- logo.svg
|   +-- logo-256.png
|   +-- ...
+-- installer\
|   +-- build-windows.bat
|   +-- package-windows.sh
|   +-- fracfixd-installer.iss
|   +-- package-macos.sh
|   +-- applauncher.c
+-- ...
```

---

## Step 3 — Run the build

1. Open **"x64 Native Tools Command Prompt for VS 2022"** from
   the Start Menu.  This is critical — a regular `cmd.exe` will
   not have the MSVC linker on `PATH`.
2. `cd` to the source-tree root:

```cmd
cd C:\Repos\FracFixD
```

3. (Optional) Add LDC to `PATH` if not already done in your user
   environment:

```cmd
set PATH=%PATH%;C:\D\ldc2\bin
```

4. Run the build script:

```cmd
installer\build-windows.bat
```

The script will:

- **Step 0/6**: Check prerequisites (MSVC, LDC, DUB, MSYS2 + GTK3
  + OpenBLAS, Inno Setup if `--installer` was passed).
- **Step 1/6**: Run `dub build --config=windows-gui --build=release-static
  --compiler=ldc2 --force` to produce `fracfixd.exe`.
- **Step 2/6**: Stage `fracfixd.exe` in `dist\fracfixd-windows\`.
- **Step 3/6**: Delegate to `installer\package-windows.sh` (running
  inside MSYS2 MINGW64) which walks `ldd` over the binary and
  GTK / OpenBLAS roots to collect ~70-80 DLLs.  Copies them all
  into `dist\fracfixd-windows\` alongside the executable.
- **Step 4/6**: Copy GdkPixbuf loaders, GTK Default + MS-Windows
  themes, Adwaita icon theme, GLib schemas — everything GTK 3
  needs at runtime — into `dist\fracfixd-windows\lib\` and
  `dist\fracfixd-windows\share\`.
- **Step 5/6**: Wrap the whole staging tree in a portable ZIP:
  `dist\fracfixd-v2.0.0-windows-x86_64.zip` (~40 MB).
- **Step 6/6**: (skipped without `--installer`).

Expected output on success:

```
============================================================
 BUILD COMPLETE
============================================================

  Output files in dist\:

    fracfixd.exe                                  6500000 bytes
    fracfixd-v2.0.0-windows-x86_64.zip            38000000 bytes

  Staging directory: dist\fracfixd-windows\

  To test the build:
    dist\fracfixd-windows\fracfixd.exe                  [GTK + BLAS DLLs bundled in this dir]
    dist\fracfixd-windows\fracfixd.exe --cli --version
```

---

## Step 4 — Build the `.exe` installer (optional)

If Inno Setup 6 is installed, re-run with `--installer`:

```cmd
installer\build-windows.bat --installer
```

This re-runs the build (or `--skip-build` to reuse the existing
`fracfixd.exe` if you don't want to recompile) and then calls
`iscc.exe` on `installer\fracfixd-installer.iss` to produce:

```
dist\FracFixD-2.0.0-windows-x86_64-setup.exe         (~30 MB)
```

End-user experience: double-click `.exe`, accept the licence,
choose install location (default `C:\Program Files\FracFixD`),
optionally tick "add to PATH" and "create desktop shortcut",
click Install.  The installer writes an uninstaller entry under
"Add or Remove Programs".

---

## Step 5 — Verify the build

From the same Command Prompt (still in source root):

```cmd
:: GUI smoke (launches the GTK3 window)
dist\fracfixd-windows\fracfixd.exe

:: CLI smoke (no GUI — uses --cli short-circuit)
dist\fracfixd-windows\fracfixd.exe --cli --version
dist\fracfixd-windows\fracfixd.exe --cli help diffprop

:: SHA-256 hash for the release-notes SHA256SUMS file
powershell -NoProfile -Command "Get-FileHash dist\fracfixd-v2.0.0-windows-x86_64.zip -Algorithm SHA256"
powershell -NoProfile -Command "Get-FileHash dist\FracFixD-2.0.0-windows-x86_64-setup.exe -Algorithm SHA256"
```

The `--version` output should display the banner:

```
+----------------------------------------------------------------+
|   Version:   2.0.0                                             |
|   Codename:  Quokka                                            |
|   Status:    Stable                                            |
+----------------------------------------------------------------+
```

---

## Step 6 — Ship the artefacts

Send Nick (or whoever is collating the release):

1. `dist\fracfixd-v2.0.0-windows-x86_64.zip`
2. `dist\FracFixD-2.0.0-windows-x86_64-setup.exe` (if you ran
   `--installer`)
3. The two SHA-256 hashes from Step 5.

They will land in `FracFixD/bin/` of the release repo alongside
the existing Linux + macOS artefacts.

---

## Troubleshooting

### `ldc2: command not found`

LDC is not on `PATH`.  Either:

- Set it for the session: `set PATH=%PATH%;C:\D\ldc2\bin`
- Or set it permanently: Control Panel → System → Advanced
  → Environment Variables → User PATH → Edit → New →
  `C:\D\ldc2\bin`.

### `cl: command not found`

You opened a regular `cmd.exe` instead of the **"x64 Native
Tools Command Prompt for VS 2022"**.  Close and reopen from
the correct Start Menu shortcut.

### `LINK : fatal error LNK1181: cannot open input file 'openblas.lib'`

The MSVC linker cannot find the OpenBLAS import library.  The
build script auto-downloads OpenBLAS to `C:\OpenBLAS\` if
missing; this error means either the download did not run or
it failed.  Manual recovery:

1. Download
   `OpenBLAS-0.3.30-x64.zip` from
   https://github.com/OpenMathLib/OpenBLAS/releases.
2. Extract to `C:\OpenBLAS\` so that
   `C:\OpenBLAS\lib\libopenblas.lib` exists.
3. Re-run `installer\build-windows.bat`.

Do NOT use MSYS2's `mingw-w64-x86_64-openblas` package: it
ships in `.a` (MinGW) format, which the MSVC linker cannot
consume.

### `No DLLs were collected`

The MSYS2-side helper (`package-windows.sh`) couldn't find GTK 3
in `/mingw64/bin/`.  Either:

- MSYS2 isn't installed at the default `C:\msys64` — edit
  `installer\build-windows.bat` and change the `MSYS2_ROOT`
  variable.
- GTK 3 is missing — install via
  `pacman -S mingw-w64-x86_64-gtk3` from the MINGW64 terminal.

### `Defender SmartScreen warning on first launch`

The unsigned `.exe` triggers a one-time SmartScreen confirmation.
Click "More info" → "Run anyway".  Code signing requires an
Authenticode certificate and is outside the v2.0.0 scope.

### `vcruntime140.dll not found` on end-user machine

The end user is missing the Visual C++ Redistributable 2015-2022
runtime.  They can install it from
https://aka.ms/vs/17/release/vc_redist.x64.exe (free download).
Most Windows installs have it already; this only affects very
clean systems.

### `Build succeeds but DLL count < 50`

Run `installer\package-windows.sh` manually from an MSYS2 MINGW64
terminal so you can see all the `ldd` warnings:

```bash
cd /c/Repos/FracFixD
./installer/package-windows.sh
```

Likely cause: a transitive dependency in the `ROOT_DLLS` or
`EXTRA_DLLS` lists couldn't be found.  Add the missing DLL to
the `EXTRA_DLLS` array in `package-windows.sh` and rerun.

---

## Reference: build-host snapshot

The reference Windows build host used during the v2.0.0
release-engineering cycle:

| Component | Version |
|---|---|
| Windows | 11 (build 22631) |
| MSYS2 | 2026-01 rolling release |
| GTK 3 (mingw-w64) | 3.24.43+ |
| OpenBLAS (official Windows pre-build at `C:\OpenBLAS\`) | 0.3.30 |
| Visual Studio Build Tools | 2022 v17.10+ |
| LDC | 1.42.0 |
| Inno Setup | 6.3.3+ |

Building on older MSYS2 / GTK versions should work but is
untested.  If you hit a compatibility issue let Nick know
which versions you're running.
