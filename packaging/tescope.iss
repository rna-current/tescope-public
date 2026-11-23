; Inno Setup script for the TEscope Windows installer.
;
; Build order (on Windows, or in CI on windows-latest):
;   1. uv sync --extra gui
;   2. uv pip install pyinstaller
;   3. cd packaging && uv run pyinstaller --noconfirm tescope-gui.spec
;   4. uv run python packaging/collect_licenses.py packaging/dist/TEscope/THIRD-PARTY-NOTICES.txt
;   5. iscc packaging/tescope.iss
;
; Produces packaging/Output/TEscope-<version>-setup.exe (unsigned; code
; signing is a release-phase operator action).

#define AppName "TEscope"
#define AppVersion "0.2.0"
#define AppPublisher "TEscope project"
#define AppURL "https://github.com/rna-current/tescope-public"

[Setup]
AppId={{7E5C0E7A-2B7C-4A48-9D5B-3A1F62C90417}}
AppName={#AppName}
AppVersion={#AppVersion}
AppPublisher={#AppPublisher}
AppSupportURL={#AppURL}
DefaultDirName={autopf}\{#AppName}
DefaultGroupName={#AppName}
DisableProgramGroupPage=yes
LicenseFile=..\LICENSE
OutputBaseFilename=TEscope-{#AppVersion}-setup
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
; per-user install by default: no admin rights needed, matches local-first design
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
ArchitecturesInstallIn64BitMode=x64compatible
; approved TEscope symbol icon (assets/branding/generate_assets.py)
SetupIconFile=tescope.ico
UninstallDisplayIcon={app}\TEscope.exe

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; Flags: unchecked

[Files]
Source: "dist\TEscope\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\{#AppName}"; Filename: "{app}\TEscope.exe"
Name: "{autodesktop}\{#AppName}"; Filename: "{app}\TEscope.exe"; Tasks: desktopicon

[Run]
Filename: "{app}\TEscope.exe"; Description: "{cm:LaunchProgram,{#AppName}}"; Flags: nowait postinstall skipifsilent

[UninstallDelete]
; the reference cache lives under %LOCALAPPDATA%\tescope and survives
; uninstall on purpose (users may reinstall); document, don't delete.
