/*
 * applauncher.c -- minimal native Mach-O launcher for the .app bundle
 *
 * Why this exists:
 *   gtk-d's Loader.d calls dlopen() with bare filenames (e.g.
 *   "libatk-1.0.0.dylib"). Apple's dyld does not search @rpath for bare
 *   names, so we must set DYLD_LIBRARY_PATH=<bundle>/Resources/lib before
 *   exec'ing the real binary.
 *
 *   We used to do this with a shell-script launcher in Contents/MacOS/.
 *   macOS Sequoia/Tahoe (>=15) tightened Gatekeeper to require the bundle's
 *   main executable to be a signable Mach-O; shell scripts now trigger
 *   _LSOpenURLsWithCompletionHandler() failed with error -10669 even after
 *   ad-hoc codesign. This native launcher fixes that.
 *
 * Compile with:
 *   cc -O2 -arch arm64 -mmacosx-version-min=13.0 \
 *      -DTARGET_BIN=\"<binary-name>\" \
 *      [-DTARGET_GDK_BACKEND='"quartz,*"'] \
 *      installer/applauncher.c -o <bundle>/Contents/MacOS/<AppName>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <libgen.h>
#include <mach-o/dyld.h>

#ifndef TARGET_BIN
#error "compile with -DTARGET_BIN=\"<name>\""
#endif
#ifndef TARGET_GDK_BACKEND
#define TARGET_GDK_BACKEND ""
#endif

int main(int argc, char *argv[]) {
    char selfpath[4096];
    uint32_t sz = sizeof(selfpath);
    if (_NSGetExecutablePath(selfpath, &sz) != 0) {
        fprintf(stderr, "applauncher: _NSGetExecutablePath failed\n");
        return 1;
    }
    /* selfpath = .../Contents/MacOS/<name>; dirname twice -> .../Contents */
    char dup1[4096];
    strncpy(dup1, selfpath, sizeof(dup1) - 1);
    dup1[sizeof(dup1) - 1] = 0;
    char *macos_dir = dirname(dup1);

    char dup2[4096];
    strncpy(dup2, macos_dir, sizeof(dup2) - 1);
    dup2[sizeof(dup2) - 1] = 0;
    char *contents = dirname(dup2);

    char libdir[4096];
    snprintf(libdir, sizeof(libdir), "%s/Resources/lib", contents);
    char binpath[4096];
    snprintf(binpath, sizeof(binpath), "%s/Resources/bin/%s", contents, TARGET_BIN);

    const char *prev = getenv("DYLD_LIBRARY_PATH");
    char dyld[8192];
    if (prev && *prev) {
        snprintf(dyld, sizeof(dyld), "%s:%s", libdir, prev);
    } else {
        snprintf(dyld, sizeof(dyld), "%s", libdir);
    }
    setenv("DYLD_LIBRARY_PATH", dyld, 1);

    if (strlen(TARGET_GDK_BACKEND) > 0 && getenv("GDK_BACKEND") == NULL) {
        setenv("GDK_BACKEND", TARGET_GDK_BACKEND, 1);
    }

    execv(binpath, argv);
    fprintf(stderr, "applauncher: execv(%s) failed\n", binpath);
    return 1;
}
