@echo off
if not exist "%~dp0\bin\bash" ( "%~dp0\bin\busybox.exe" --install -s "%~dp0\bin" )
for /f %%i in ('%~dp0\bin\testcpu.exe') do (
if "%%i" == "avx2" ( "%~dp0\bin\mmseqs_avx2.exe" %* )
if "%%i" == "sse4.1" ( "%~dp0\bin\mmseqs_sse41.exe" %* )
if "%%i" == "fail" ( echo Unsupported CPU! MMseqs2 requires SSE4.1 or AVX2 instruction set support. )
)
