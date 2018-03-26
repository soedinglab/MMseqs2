@echo off
if not exist "%~dp0\bin\bash" ( "%~dp0\bin\busybox.exe" --install -s "%~dp0\bin" )
"%~dp0\bin\mmseqs.exe" %*
