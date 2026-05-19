# Building FracFixD on Windows (x86_64)

> **NOTE.** The FracFixD source code is closed and confidential
> (see [`../README.md` → Source code](../README.md#source-code)).
> This document describes the **Windows build recipe** as a
> provenance reference and for licensed source-tree recipients.
> A pre-built Windows x86_64 binary is **planned** for a
> follow-up release; until it lands, this walk-through is the
> path for source-tree users.

The recipe below mirrors the TagGen v1.2.5 Windows release
process.

---

## 1. Prerequisites

### Hardware

- 64-bit Windows 10 or Windows 11.
- At least 4 GB RAM, 10 GB disk for the build environment.

### Toolchain

Install the following, in this order:

#### 1.1 MSYS2 (for GTK 3 + dependencies)

Download from https://www.msys2.org/ and run the installer.
Once installed, open the **MSYS2 MINGW64** shell and:

```bash
pacman -Syu                                    # initial sync; close + reopen shell after
pacman -Syu                                    # second pass
pacman -S mingw-w64-x86_64-gtk3 \
          mingw-w64-x86_64-pkg-config \
          mingw-w64-x86_64-gcc \
          mingw-w64-x86_64-binutils \
          unzip wget
```

GTK 3 should now be at `C:\msys64\mingw64\bin\` (around 70
DLLs).

#### 1.2 Visual Studio Build Tools (for the MSVC linker)

LDC's `-mtriple=x86_64-windows-msvc` requires the MSVC linker.

1. Download Visual Studio 2022 Build Tools from
   https://visualstudio.microsoft.com/downloads/.
2. Install the **"Desktop development with C++"** workload.
3. Confirm `link.exe` is on PATH from the **x64 Native Tools
   Command Prompt for VS 2022**.

#### 1.3 LDC 1.42 (Windows release)

```cmd
:: Download LDC 1.42.0-windows-x64.7z from
::   https://github.com/ldc-developers/ldc/releases
:: Extract to C:\D\ldc2

set PATH=%PATH%;C:\D\ldc2\bin
ldc2 --version
```

#### 1.4 Inno Setup 6 (optional, for the `.exe` installer)

Download from https://jrsoftware.org/isinfo.php.  Default
install location is `C:\Program Files (x86)\Inno Setup 6\`.

---

## 2. Source-tree access

See [`../README.md` → Source code](../README.md#source-code).
The build commands below assume you have the source tree
checked out and that `dub.json` ships the `windows`
configuration (it does in the public-release v2.0.0 tree).

Open the **x64 Native Tools Command Prompt for VS 2022** and
`cd` to the source-tree root.

---

## 3. Build the executable

```cmd
set PATH=%PATH%;C:\D\ldc2\bin

:: Default Windows config is CLI-only (no GTK at compile time).
:: For the GUI build use the macos-style config adapted for win,
:: which is wired up in `windows-gui` (custom; see source-tree dub.json).

dub build --config=windows --build=release-static --compiler=ldc2 --force

:: Produces fracfixd-cli.exe (~4 MB)
:: For the GUI build:
dub build --config=windows-gui --build=release-static --compiler=ldc2 --force
:: Produces fracfixd.exe (~5 MB; needs GTK 3 DLLs at runtime)
```

For the GUI binary, LDC must target the MSVC ABI:

```json
"dflags-ldc": ["-mtriple=x86_64-windows-msvc"]
```

(This is already present in the `windows` and `windows-gui`
configurations.)

---

## 4. Run `build-windows.bat`

The in-tree all-in-one builder script
[`../installer/build-windows.bat`](../installer/build-windows.bat)
wraps prerequisite checks + DUB build + GTK-DLL collection +
optional Inno Setup invocation:

```cmd
cd FracFixD\installer
build-windows.bat                  :: Build + portable ZIP
build-windows.bat --installer      :: Build + ZIP + Inno Setup .exe
```

The script's logical steps:

1. Locate LDC on PATH (search `C:\D\ldc2\bin\` if not found).
2. Run `dub build --config=windows-gui --build=release-static
   --compiler=ldc2 --force` (or `--config=windows` for the
   CLI-only build).
3. Stage the executable in `dist\fracfixd-windows\`.
4. Walk `ldd` (msys2's `ntldd` is the equivalent) on
   `fracfixd.exe` and copy every DLL transitively required
   into `dist\fracfixd-windows\`.  Around 70 DLLs end up here
   (GTK 3, GLib, GObject, GIO, Pango, Cairo, GdkPixbuf,
   HarfBuzz, FreeType, libintl, libpcre, fontconfig, etc.).
5. Copy GdkPixbuf loaders to
   `dist\fracfixd-windows\lib\gdk-pixbuf-2.0\2.10.0\` and
   patch `loaders.cache` to use relative paths.
6. Copy GTK 3 themes (Adwaita), Adwaita icon theme, GLib
   schemas.
7. Drop `README.md`, `LICENSE` into the staging dir.
8. ZIP the staging dir:
   `fracfixd-v2.0.0-windows-x86_64.zip` (~30-40 MB).
9. If `--installer` was passed, invoke Inno Setup on
   `installer\fracfixd-installer.iss` to produce
   `FracFixD-2.0.0-windows-x86_64-setup.exe` (~25 MB).

---

## 5. Inno Setup installer (optional)

The shipped script
[`../installer/fracfixd-installer.iss`](../installer/fracfixd-installer.iss)
declares:

```iss
[Setup]
AppName=FracFixD
AppVersion=2.0.0
AppPublisher=Arnaroo Ribologicals
AppPublisherURL=https://github.com/Arnaroo/FracFixR
DefaultDirName={autopf}\FracFixD
DefaultGroupName=FracFixD
LicenseFile=..\LICENSE
OutputBaseFilename=FracFixD-2.0.0-windows-x86_64-setup
SetupIconFile=..\resources\fracfixd-icon.ico
Compression=lzma2
SolidCompression=yes

[Files]
Source: "..\..\dist\fracfixd-windows\*"; DestDir: "{app}"; Flags: recursesubdirs

[Icons]
Name: "{group}\FracFixD"; Filename: "{app}\fracfixd.exe"
Name: "{commondesktop}\FracFixD"; Filename: "{app}\fracfixd.exe"; Tasks: desktopicon

[Tasks]
Name: desktopicon; Description: "Create a desktop shortcut"; GroupDescription: "Optional:"
Name: addtopath;  Description: "Add to system PATH";          GroupDescription: "Optional:"

[Code]
{ ... PATH modification helpers ... }
```

Build with Inno Setup's `iscc.exe`:

```cmd
"C:\Program Files (x86)\Inno Setup 6\iscc.exe" installer\fracfixd-installer.iss
```

Output: `FracFixD-2.0.0-windows-x86_64-setup.exe`.

---

## 6. Verification

```cmd
:: Smoke test
dist\fracfixd-windows\fracfixd.exe --version
dist\fracfixd-windows\fracfixd.exe --cli help diffprop

:: Hash the artefacts (PowerShell)
Get-FileHash dist\fracfixd-v2.0.0-windows-x86_64.zip -Algorithm SHA256
Get-FileHash FracFixD-2.0.0-windows-x86_64-setup.exe -Algorithm SHA256
```

End-user experience for the ZIP: extract the folder anywhere,
double-click `fracfixd.exe`.  For the installer: run the
`.exe`, accept the license, optionally add to PATH and create
a desktop shortcut.

---

## 7. Known Windows-specific issues

- **Console window appears alongside the GUI** — by default
  LDC produces console-subsystem executables.  Add
  `-Wl,/SUBSYSTEM:WINDOWS` to suppress the console for the
  GUI variant (a no-console build is planned; v2.0.0 uses the
  console-subsystem default to keep the build script
  simple).
- **"`vcruntime140.dll` not found"** — install the Visual C++
  Redistributable 2015-2022 from Microsoft.  Most users
  already have it.
- **Defender SmartScreen warning on first launch** — the
  unsigned `.exe` triggers a one-time SmartScreen confirmation.
  Code signing is not yet in scope for v2.0.0.
- **GTK 3 theme appears wrong** — confirm Adwaita theme and
  schemas were copied to `lib\share\` in the staging dir.  The
  installer script handles this automatically; manual builds
  may forget the `share\glib-2.0\schemas\` directory.

---

## 8. Build-host snapshot (for full reproducibility)

The reference Windows build host used during the v2.0.0
release-engineering cycle:

| Component | Version |
|---|---|
| Windows | 11 (build 22631) |
| MSYS2 | 2026-01-12 release |
| GTK 3 (mingw-w64) | 3.24.43 |
| LDC | 1.42.0 |
| Visual Studio Build Tools | 2022 v17.10 |
| Inno Setup | 6.3.3 |

A future release MAY ship a CI workflow (`.github/workflows/`)
that automates this on GitHub Actions; for v2.0.0 the build is
manual on a Windows host.
