@echo off
if not exist "%~dp0\bin\bash" ( goto installBusyBox )
goto mmseqs

:installBusyBox
"%~dp0\bin\busybox.exe" --install -s "%~dp0\bin" > nul 2>&1
if not exist "%~dp0\bin\bash" ( goto installBusyBoxUAC )
goto mmseqs

:installBusyBoxUAC
echo MMseqs2 will now request administrator permissions to install helper tools through Busybox in a subdirectory.
echo WScript.Sleep 2000 > "%temp%\~ElevateBusyBox.vbs"
echo set UAC = CreateObject^("Shell.Application"^) >> "%temp%\~ElevateBusyBox.vbs"
echo UAC.ShellExecute "%~dp0\bin\busybox.exe", "--install -s ""%~dp0\bin""", , "runas", 0 >> "%temp%\~ElevateBusyBox.vbs"
echo WScript.Sleep 2000 >> "%temp%\~ElevateBusyBox.vbs"
cscript "%temp%\~ElevateBusyBox.vbs" > nul 2>&1
del /Q /F "%temp%\~ElevateBusyBox.vbs"
if not exist "%~dp0\bin\bash" ( goto noBusyBox )
goto mmseqs

:mmseqs
"%~dp0\bin\mmseqs.exe" %*
exit /b 0

:noBusyBox
echo Error: Could not install BusyBox helper tools. Please execute the following command manually: >&2
echo "%~dp0\bin\busybox.exe" --install -s "%~dp0\bin" >&2
exit /b 1
