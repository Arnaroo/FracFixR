@echo off
setlocal enabledelayedexpansion
REM ============================================================================
REM  FracFixD Windows Build Script
REM  All-in-one: build, collect GTK + scientific DLLs, package, and optionally
REM  create an Inno Setup installer.
REM
REM  PREREQUISITES (install once):
REM    1. MSYS2           -- https://www.msys2.org/  (default: C:\msys64)
REM       Then in MSYS2 MINGW64 terminal:
REM         pacman -Syu && pacman -Su
REM         pacman -S mingw-w64-x86_64-gtk3 \
REM                   mingw-w64-x86_64-openblas \
REM                   mingw-w64-x86_64-lapack
REM    2. Visual Studio Build Tools (Desktop C++ workload)
REM       -- https://visualstudio.microsoft.com/visual-cpp-build-tools/
REM    3. LDC2 + DUB (via VisualD or standalone)
REM       -- https://github.com/ldc-developers/ldc/releases
REM       Recommended: ldc2-1.42.0-windows-multilib.7z, extract to C:\D\ldc2
REM    4. Inno Setup (optional, for .exe installer)
REM       -- winget install JRSoftware.InnoSetup
REM          or https://jrsoftware.org/isdl.php
REM
REM  HOW TO RUN:
REM    Open "x64 Native Tools Command Prompt for VS 2022"
REM    cd C:\path\to\fracfixd-source
REM    installer\build-windows.bat
REM
REM  OPTIONS:
REM    build-windows.bat              Build + package (portable ZIP)
REM    build-windows.bat --installer  Build + package + create .exe installer
REM    build-windows.bat --skip-build Package only (reuse existing fracfixd.exe)
REM    build-windows.bat --cli-only   Build only the CLI variant (no GUI)
REM    build-windows.bat --help       Show this help
REM ============================================================================

set "VERSION=2.0.0"
set "MSYS2_ROOT=C:\msys64"
set "DIST_DIR=dist\fracfixd-windows"
set "DO_INSTALLER=0"
set "SKIP_BUILD=0"
set "BUILD_CONFIG=windows-gui"

REM --- Parse arguments ---
:parse_args
if "%~1"=="" goto :done_args
if /I "%~1"=="--installer"  set "DO_INSTALLER=1"
if /I "%~1"=="--skip-build" set "SKIP_BUILD=1"
if /I "%~1"=="--cli-only"   set "BUILD_CONFIG=windows"
if /I "%~1"=="--help"       goto :show_help
shift
goto :parse_args
:done_args

echo.
echo  ============================================================
echo   FracFixD v%VERSION% - Windows Build Script
echo   Config: %BUILD_CONFIG%
echo  ============================================================
echo.

REM ============================================================================
REM  STEP 0: Verify prerequisites
REM ============================================================================
echo [0/6] Checking prerequisites...
echo.

REM Check we're in the x64 Native Tools environment
where cl >nul 2>&1
if errorlevel 1 (
    echo  ERROR: MSVC compiler 'cl' not found.
    echo  You must run this script from the x64 Native Tools Command Prompt.
    echo  Find it in: Start Menu, Visual Studio 2022, x64 Native Tools...
    echo.
    exit /b 1
)
echo   [OK] MSVC environment detected

REM Check LDC
where ldc2 >nul 2>&1
if errorlevel 1 (
    echo   [!!] ldc2 not on PATH. Searching common locations...
    if exist "C:\D\ldc2\bin\ldc2.exe" (
        set "PATH=%PATH%;C:\D\ldc2\bin"
        echo   [OK] Found LDC2 at C:\D\ldc2\bin
    ) else (
        for /d %%D in ("C:\D\ldc2-*") do (
            if exist "%%D\bin\ldc2.exe" (
                set "PATH=%PATH%;%%D\bin"
                echo   [OK] Found LDC2 at %%D\bin
                goto :found_ldc
            )
        )
        echo   ERROR: LDC2 not found. Install from:
        echo     https://github.com/ldc-developers/ldc/releases
        echo   Then add its bin\ directory to PATH, e.g.:
        echo     set PATH=%%PATH%%;C:\D\ldc2-1.42.0-windows-multilib\bin
        exit /b 1
    )
)
:found_ldc
for /f "tokens=*" %%v in ('ldc2 --version 2^>^&1') do (
    echo   [OK] %%v
    goto :ldc_ver_done
)
:ldc_ver_done

REM Check DUB
where dub >nul 2>&1
if errorlevel 1 (
    echo   ERROR: dub not found. It should be in the same directory as ldc2.
    exit /b 1
)
echo   [OK] dub found

REM Check MSYS2 + GTK3
if not exist "%MSYS2_ROOT%\mingw64\bin\libgtk-3-0.dll" (
    echo   ERROR: GTK3 not found at %MSYS2_ROOT%\mingw64\bin\
    echo   Install MSYS2 from https://www.msys2.org/ then run:
    echo     pacman -S mingw-w64-x86_64-gtk3
    exit /b 1
)
echo   [OK] MSYS2 + GTK3 found at %MSYS2_ROOT%

REM Check OpenBLAS + LAPACK (FracFixD specific)
if not exist "%MSYS2_ROOT%\mingw64\bin\libopenblas.dll" (
    echo   WARNING: OpenBLAS not found at %MSYS2_ROOT%\mingw64\bin\
    echo   FracFixD needs OpenBLAS at runtime. Install with:
    echo     pacman -S mingw-w64-x86_64-openblas mingw-w64-x86_64-lapack
    echo   Build will likely fail at link time.
)
if not exist "%MSYS2_ROOT%\mingw64\bin\liblapack.dll" (
    echo   WARNING: LAPACK not found at %MSYS2_ROOT%\mingw64\bin\
)

REM Check dub.json exists
if not exist "dub.json" (
    echo   ERROR: dub.json not found. Run this script from the FracFixD source root.
    exit /b 1
)
echo   [OK] dub.json found

REM Check the package-windows.sh helper exists
if not exist "installer\package-windows.sh" (
    echo   ERROR: installer\package-windows.sh not found.
    echo   This script delegates DLL collection to that helper.
    exit /b 1
)
echo   [OK] installer\package-windows.sh found
echo.

REM ============================================================================
REM  STEP 1: Build
REM ============================================================================
if "%SKIP_BUILD%"=="1" (
    echo [1/6] Skipping build [--skip-build]
    if not exist "fracfixd.exe" (
        echo   ERROR: fracfixd.exe not found. Cannot skip build without existing binary.
        exit /b 1
    )
    goto :step2
)

echo [1/6] Building FracFixD...
echo   Config: %BUILD_CONFIG% ^| Build: release-static ^| Compiler: ldc2
echo.

REM Clean stale artifacts
dub clean 2>nul

REM Build with MSVC-targeting LDC. release-static gives a self-contained .exe
REM (statically-linked phobos/druntime; only GTK + BLAS DLLs needed at runtime).
dub build --compiler=ldc2 --config=%BUILD_CONFIG% -b release-static --force
if errorlevel 1 (
    echo.
    echo  BUILD FAILED. Common fixes:
    echo    - Ensure you're in the x64 Native Tools Command Prompt
    echo    - Try: dub clean ^&^& rd /s /q %%LOCALAPPDATA%%\dub
    echo    - Check that ldc2 --version shows a working compiler
    echo    - Check OpenBLAS + LAPACK installed in MSYS2
    exit /b 1
)

if not exist "fracfixd.exe" (
    if exist "fracfixd-cli.exe" (
        echo   Renaming fracfixd-cli.exe -^> fracfixd.exe for packaging
        move /y fracfixd-cli.exe fracfixd.exe >nul
    ) else (
        echo   ERROR: Build completed but fracfixd.exe not produced.
        exit /b 1
    )
)

for %%F in (fracfixd.exe) do echo   [OK] fracfixd.exe built [%%~zF bytes]
echo.

REM ============================================================================
REM  STEP 2: Clean and create staging directory
REM ============================================================================
:step2
echo [2/6] Preparing staging directory...
if exist "%DIST_DIR%" rd /s /q "%DIST_DIR%"
mkdir "%DIST_DIR%" 2>nul
copy /y fracfixd.exe "%DIST_DIR%\" >nul
echo   [OK] %DIST_DIR%\fracfixd.exe staged
echo.

REM ============================================================================
REM  STEP 3 + 4: Collect GTK + BLAS DLLs, themes, icons, resources via the
REM  installer\package-windows.sh helper (run inside MSYS2 MINGW64).
REM ============================================================================
echo [3/6] Collecting GTK + scientific DLLs + themes + icons + resources...
echo   Delegating to installer\package-windows.sh [in MINGW64 environment]

set "MSYS2_BASH=%MSYS2_ROOT%\usr\bin\bash.exe"
if not exist "%MSYS2_BASH%" (
    echo   ERROR: %MSYS2_BASH% not found.
    exit /b 1
)
set "MSYSTEM=MINGW64"
set "CHERE_INVOKING=1"
"%MSYS2_BASH%" -lc "VERSION=%VERSION% ./installer/package-windows.sh"
if errorlevel 1 (
    echo   ERROR: package-windows.sh failed inside MSYS2.
    echo   You can also run it manually from an MSYS2 MINGW64 terminal:
    echo     cd /c/path/to/fracfixd-source
    echo     VERSION=%VERSION% ./installer/package-windows.sh
    exit /b 1
)

REM Sanity check: count DLLs landed in the staging dir
set "DLL_COUNT=0"
for %%F in ("%DIST_DIR%\*.dll") do set /a DLL_COUNT+=1
if "%DLL_COUNT%"=="0" (
    echo   ERROR: No DLLs were collected. The .exe will not run standalone.
    exit /b 1
)
echo   [OK] %DLL_COUNT% DLLs in %DIST_DIR%\
echo.

echo [4/6] Resources + docs were copied by package-windows.sh.
if exist "%DIST_DIR%\resources" echo   [OK] %DIST_DIR%\resources\
if exist "%DIST_DIR%\share"     echo   [OK] %DIST_DIR%\share\
if exist "%DIST_DIR%\lib"       echo   [OK] %DIST_DIR%\lib\
echo.

REM ============================================================================
REM  STEP 5: Create portable ZIP
REM ============================================================================
echo [5/6] Creating portable ZIP archive [MSYS2 zip, POSIX paths]...
set "ZIP_NAME=fracfixd-v%VERSION%-windows-x86_64.zip"
if exist "dist\%ZIP_NAME%" del /f /q "dist\%ZIP_NAME%"
"%MSYS2_BASH%" -lc "cd dist && rm -f '%ZIP_NAME%' && zip -qr '%ZIP_NAME%' 'fracfixd-windows/'"
if exist "dist\%ZIP_NAME%" (
    for %%F in ("dist\%ZIP_NAME%") do echo   [OK] dist\%ZIP_NAME% [%%~zF bytes]
) else (
    echo   [WARN] MSYS2 zip not available. Falling back to PowerShell Compress-Archive
    echo          [paths may use backslashes - extract on Windows only].
    powershell -NoProfile -Command "Compress-Archive -Path '%DIST_DIR%' -DestinationPath 'dist\%ZIP_NAME%' -Force"
    if exist "dist\%ZIP_NAME%" (
        for %%F in ("dist\%ZIP_NAME%") do echo   [OK] dist\%ZIP_NAME% [%%~zF bytes]
    ) else (
        echo   [SKIP] ZIP creation failed. You can manually zip %DIST_DIR%\
    )
)
echo.

REM ============================================================================
REM  STEP 6: Create Inno Setup installer (optional)
REM ============================================================================
if "%DO_INSTALLER%"=="0" (
    echo [6/6] Skipping installer [use --installer to create]
    goto :summary
)

echo [6/6] Creating Inno Setup installer...

set "ISCC="
if exist "C:\Program Files (x86)\Inno Setup 6\iscc.exe" set "ISCC=C:\Program Files (x86)\Inno Setup 6\iscc.exe"
if not defined ISCC if exist "C:\Program Files\Inno Setup 6\iscc.exe" set "ISCC=C:\Program Files\Inno Setup 6\iscc.exe"
if not defined ISCC (
    where iscc >nul 2>&1
    if not errorlevel 1 (
        for /f "tokens=*" %%p in ('where iscc') do set "ISCC=%%p"
    )
)

if not defined ISCC (
    echo   ERROR: Inno Setup compiler not found.
    echo   Install with:  winget install JRSoftware.InnoSetup
    echo   Or download from:  https://jrsoftware.org/isdl.php
    echo   After installing, re-run:  build-windows.bat --installer
    goto :summary
)

echo   Using: !ISCC!
"!ISCC!" installer\fracfixd-installer.iss
if errorlevel 1 (
    echo   ERROR: Installer compilation failed.
) else (
    if exist "dist\FracFixD-%VERSION%-windows-x86_64-setup.exe" (
        for %%F in ("dist\FracFixD-%VERSION%-windows-x86_64-setup.exe") do (
            echo   [OK] dist\FracFixD-%VERSION%-windows-x86_64-setup.exe [%%~zF bytes]
        )
    )
)
echo.

REM ============================================================================
REM  Summary
REM ============================================================================
:summary
echo.
echo  ============================================================
echo   BUILD COMPLETE
echo  ============================================================
echo.
echo   Output files in dist\:
echo.
if exist "fracfixd.exe" (
    for %%F in (fracfixd.exe) do echo     fracfixd.exe                                  %%~zF bytes
)
if exist "dist\fracfixd-v%VERSION%-windows-x86_64.zip" (
    for %%F in ("dist\fracfixd-v%VERSION%-windows-x86_64.zip") do echo     fracfixd-v%VERSION%-windows-x86_64.zip            %%~zF bytes
)
if exist "dist\FracFixD-%VERSION%-windows-x86_64-setup.exe" (
    for %%F in ("dist\FracFixD-%VERSION%-windows-x86_64-setup.exe") do echo     FracFixD-%VERSION%-windows-x86_64-setup.exe   %%~zF bytes
)
echo.
echo   Staging directory: %DIST_DIR%\
echo.
echo   To test the build:
echo     %DIST_DIR%\fracfixd.exe                 [GTK + BLAS DLLs bundled in this dir]
echo     %DIST_DIR%\fracfixd.exe --cli --version
echo.
goto :eof

REM ============================================================================
:show_help
echo.
echo  FracFixD v%VERSION% Windows Build Script
echo.
echo  USAGE:
echo    build-windows.bat              Build + package [portable ZIP]
echo    build-windows.bat --installer  Build + package + .exe installer
echo    build-windows.bat --skip-build Package only [reuse existing fracfixd.exe]
echo    build-windows.bat --cli-only   Build only the CLI variant [no GUI]
echo    build-windows.bat --help       Show this help
echo.
echo  PREREQUISITES:
echo    1. MSYS2 with GTK3 + OpenBLAS + LAPACK:
echo       - Install MSYS2 from https://www.msys2.org/
echo       - Open MSYS2 MINGW64 terminal, run:
echo           pacman -Syu ^&^& pacman -Su
echo           pacman -S mingw-w64-x86_64-gtk3 ^
echo                     mingw-w64-x86_64-openblas ^
echo                     mingw-w64-x86_64-lapack
echo.
echo    2. Visual Studio Build Tools:
echo       - https://visualstudio.microsoft.com/visual-cpp-build-tools/
echo       - Select "Desktop development with C++" workload
echo.
echo    3. LDC2 [D compiler]:
echo       - https://github.com/ldc-developers/ldc/releases
echo       - Extract to C:\D\ldc2 [or similar]
echo       - Add bin\ to PATH:  set PATH=%%PATH%%;C:\D\ldc2-1.42.0-windows-multilib\bin
echo.
echo    4. Inno Setup [optional, for --installer]:
echo       - winget install JRSoftware.InnoSetup
echo.
echo  RUN FROM: x64 Native Tools Command Prompt for VS 2022
echo.
goto :eof
