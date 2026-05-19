; FracFixD Windows Installer Script (Inno Setup 6+)
; Compile with: iscc fracfixd-installer.iss
;
; Before compiling this script, run package-windows.sh from an MSYS2
; MINGW64 terminal to collect the executable, GTK + scientific DLLs,
; and resources into the dist/fracfixd-windows/ staging directory.
;
; Adapted from TagGen v1.2.5 installer/taggen-installer.iss.

#define MyAppName      "FracFixD"
#define MyAppVersion   "2.0.0"
#define MyAppPublisher "Arnaroo Ribologicals"
#define MyAppURL       "https://github.com/Arnaroo/FracFixR"
#define MyAppExeName   "fracfixd.exe"

; Staging directory produced by package-windows.sh
#define StagingDir     "..\dist\fracfixd-windows"

[Setup]
AppId={{E2A4F1C5-7B6D-4E8A-9F0B-1A2C3D4E5F6A}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
AppPublisherURL={#MyAppURL}
AppSupportURL={#MyAppURL}/issues
DefaultDirName={autopf}\{#MyAppName}
DefaultGroupName={#MyAppName}
AllowNoIcons=yes
OutputDir=..\dist
OutputBaseFilename=FracFixD-{#MyAppVersion}-windows-x86_64-setup
Compression=lzma2/ultra64
SolidCompression=yes
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
ChangesEnvironment=yes
WizardStyle=modern
LicenseFile={#StagingDir}\LICENSE.txt
SetupIconFile={#StagingDir}\resources\fracfixd-icon.ico
UninstallDisplayIcon={app}\resources\fracfixd-icon.ico
MinVersion=10.0

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon";  Description: "Create a &desktop shortcut";       GroupDescription: "Additional shortcuts:"
Name: "addtopath";    Description: "Add FracFixD to the system &PATH"; GroupDescription: "Environment:"; Flags: checkedonce

[Files]
; Main executable
Source: "{#StagingDir}\fracfixd.exe"; DestDir: "{app}"; Flags: ignoreversion

; GTK + scientific DLLs (all .dll files in staging directory)
Source: "{#StagingDir}\*.dll"; DestDir: "{app}"; Flags: ignoreversion

; Resources (icons, logos)
Source: "{#StagingDir}\resources\*"; DestDir: "{app}\resources"; Flags: ignoreversion recursesubdirs createallsubdirs

; GTK runtime data (themes, icons, loaders)
Source: "{#StagingDir}\share\*"; DestDir: "{app}\share"; Flags: ignoreversion recursesubdirs createallsubdirs; Check: DirExists(ExpandConstant('{#StagingDir}\share'))
Source: "{#StagingDir}\lib\*";   DestDir: "{app}\lib";   Flags: ignoreversion recursesubdirs createallsubdirs; Check: DirExists(ExpandConstant('{#StagingDir}\lib'))

; Documentation
Source: "{#StagingDir}\README.md";    DestDir: "{app}"; Flags: ignoreversion
Source: "{#StagingDir}\LICENSE.txt";  DestDir: "{app}"; Flags: ignoreversion
Source: "{#StagingDir}\CHANGELOG.md"; DestDir: "{app}"; Flags: ignoreversion; Check: FileExists(ExpandConstant('{#StagingDir}\CHANGELOG.md'))

[Icons]
Name: "{group}\{#MyAppName}";              Filename: "{app}\{#MyAppExeName}"; WorkingDir: "{app}"
Name: "{group}\Uninstall {#MyAppName}";     Filename: "{uninstallexe}"
Name: "{commondesktop}\{#MyAppName}";       Filename: "{app}\{#MyAppExeName}"; WorkingDir: "{app}"; Tasks: desktopicon

[Registry]
; Add to PATH (user-level, non-destructive)
Root: HKCU; Subkey: "Environment"; \
    ValueType: expandsz; ValueName: "Path"; ValueData: "{olddata};{app}"; \
    Tasks: addtopath; Check: NeedsAddPath(ExpandConstant('{app}'))

[Run]
Filename: "{app}\{#MyAppExeName}"; Description: "Launch {#MyAppName}"; Flags: nowait postinstall skipifsilent

[Code]
// Check whether {app} is already on the user PATH to avoid duplicates
function NeedsAddPath(Param: string): Boolean;
var
  OrigPath: string;
begin
  if not RegQueryStringValue(HKEY_CURRENT_USER,
    'Environment', 'Path', OrigPath)
  then begin
    Result := True;
    exit;
  end;
  Result := Pos(';' + Uppercase(Param) + ';',
                ';' + Uppercase(OrigPath) + ';') = 0;
end;

// Remove {app} from PATH on uninstall
procedure CurUninstallStepChanged(CurUninstallStep: TUninstallStep);
var
  Path: string;
  AppDir: string;
  P: Integer;
begin
  if CurUninstallStep = usPostUninstall then
  begin
    if RegQueryStringValue(HKEY_CURRENT_USER,
      'Environment', 'Path', Path) then
    begin
      AppDir := ExpandConstant('{app}');
      P := Pos(';' + Uppercase(AppDir), Uppercase(Path));
      if P <> 0 then
      begin
        Delete(Path, P, Length(';' + AppDir));
        RegWriteStringValue(HKEY_CURRENT_USER,
          'Environment', 'Path', Path);
      end
      else begin
        P := Pos(Uppercase(AppDir) + ';', Uppercase(Path));
        if P <> 0 then
        begin
          Delete(Path, P, Length(AppDir + ';'));
          RegWriteStringValue(HKEY_CURRENT_USER,
            'Environment', 'Path', Path);
        end
        else if Uppercase(Path) = Uppercase(AppDir) then
        begin
          RegDeleteValue(HKEY_CURRENT_USER,
            'Environment', 'Path');
        end;
      end;
    end;
  end;
end;
