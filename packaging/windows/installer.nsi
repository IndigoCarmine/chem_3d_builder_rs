; NSIS installer script for Chem 3D Builder.
;
; Compiled by the release workflow with:
;   makensis /DVERSION=<x.y.z> packaging\windows\installer.nsi
; run from the repository root, so the relative paths below resolve against it.

Unicode true
!include "MUI2.nsh"

!ifndef VERSION
  !define VERSION "0.0.0"
!endif

!define APPNAME "Chem 3D Builder"
!define EXENAME "chem_3d_builder_rs.exe"
!define COMPANY "IndigoCarmine"
!define UNINSTKEY "Software\Microsoft\Windows\CurrentVersion\Uninstall\${APPNAME}"

Name "${APPNAME} ${VERSION}"
OutFile "chem_3d_builder-${VERSION}-setup.exe"
InstallDir "$PROGRAMFILES64\${APPNAME}"
InstallDirRegKey HKLM "Software\${APPNAME}" "InstallDir"
RequestExecutionLevel admin

!define MUI_ABORTWARNING
; A custom icon can be dropped in later via !define MUI_ICON / MUI_UNICON.

!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
!define MUI_FINISHPAGE_RUN "$INSTDIR\${EXENAME}"
!insertmacro MUI_PAGE_FINISH

!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES

!insertmacro MUI_LANGUAGE "English"

Section "Install"
  SetOutPath "$INSTDIR"
  File "target\release\${EXENAME}"

  ; OpenBabel runtime. The layout is load-bearing, not incidental:
  ;   - openbabel-3.dll must be beside the exe or it will not start at all,
  ;     and OpenBabel finds its .obf plugins in that DLL's directory.
  ;   - data\ beside the exe is what makes the app resolve BABEL_DATADIR here
  ;     instead of at the absolute path baked in on the build machine. Without
  ;     it the app still launches, and then force fields and 3D generation are
  ;     silently dead — a failure that never reproduces on a dev box.
  ; `cargo build` puts the DLLs and plugins in target\release; the release
  ; workflow stages data\ next to them.
  File "target\release\openbabel-3.dll"
  File "target\release\inchi.dll"
  File "target\release\*.obf"
  File /r "target\release\data"

  WriteUninstaller "$INSTDIR\Uninstall.exe"

  CreateDirectory "$SMPROGRAMS\${APPNAME}"
  CreateShortcut "$SMPROGRAMS\${APPNAME}\${APPNAME}.lnk" "$INSTDIR\${EXENAME}"
  CreateShortcut "$SMPROGRAMS\${APPNAME}\Uninstall ${APPNAME}.lnk" "$INSTDIR\Uninstall.exe"

  WriteRegStr HKLM "Software\${APPNAME}" "InstallDir" "$INSTDIR"

  ; Add/Remove Programs entry.
  WriteRegStr HKLM "${UNINSTKEY}" "DisplayName" "${APPNAME}"
  WriteRegStr HKLM "${UNINSTKEY}" "DisplayVersion" "${VERSION}"
  WriteRegStr HKLM "${UNINSTKEY}" "Publisher" "${COMPANY}"
  WriteRegStr HKLM "${UNINSTKEY}" "DisplayIcon" "$INSTDIR\${EXENAME}"
  WriteRegStr HKLM "${UNINSTKEY}" "UninstallString" "$\"$INSTDIR\Uninstall.exe$\""
  WriteRegDWORD HKLM "${UNINSTKEY}" "NoModify" 1
  WriteRegDWORD HKLM "${UNINSTKEY}" "NoRepair" 1
SectionEnd

Section "Uninstall"
  Delete "$INSTDIR\${EXENAME}"
  Delete "$INSTDIR\openbabel-3.dll"
  Delete "$INSTDIR\inchi.dll"
  Delete "$INSTDIR\*.obf"
  RMDir /r "$INSTDIR\data"
  Delete "$INSTDIR\Uninstall.exe"
  RMDir "$INSTDIR"

  Delete "$SMPROGRAMS\${APPNAME}\${APPNAME}.lnk"
  Delete "$SMPROGRAMS\${APPNAME}\Uninstall ${APPNAME}.lnk"
  RMDir "$SMPROGRAMS\${APPNAME}"

  DeleteRegKey HKLM "${UNINSTKEY}"
  DeleteRegKey HKLM "Software\${APPNAME}"
SectionEnd
