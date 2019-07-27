#!/bin/bash -e
REPO="$(readlink -f $1)"
BUILD="$(readlink -f $2)"
BINARY_NAME="${3:-mmseqs}"

if [ ! -d "$REPO" ]; then
    echo "${BINARY_NAME} repository missing"
    exit 1
fi

mkdir -p "$BUILD/build_sse41" && cd "$BUILD/build_sse41"
cmake -DCMAKE_BUILD_TYPE=Release -DHAVE_TESTS=0 -DHAVE_MPI=0 -DHAVE_SSE4_1=1 -DBUILD_SHARED_LIBS=OFF -DCMAKE_EXE_LINKER_FLAGS="-static -static-libgcc -static-libstdc++" -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" "$REPO"
make -j 4
mkdir -p "$BUILD/${BINARY_NAME}/bin"
cp "src/${BINARY_NAME}" "$BUILD/${BINARY_NAME}/bin/${BINARY_NAME}"
for i in $(ldd "src/${BINARY_NAME}" | awk '{ print $3 }' | grep -v cygdrive | grep -v '???'); do
    cp $i "$BUILD/${BINARY_NAME}/bin";
done

cp -f /usr/libexec/busybox-standalone/bin/busybox.exe "$BUILD/${BINARY_NAME}/bin"

cat << EOF > "$BUILD/${BINARY_NAME}/${BINARY_NAME}.bat"
@echo off
if not exist "%~dp0\\bin\\bash" ( goto installBusyBox )
goto binary

:installBusyBox
"%~dp0\\bin\\busybox.exe" --install -s "%~dp0\\bin" > nul 2>&1
if not exist "%~dp0\\bin\\bash" ( goto installBusyBoxUAC )
goto binary

:installBusyBoxUAC
echo Administrator permissions will now be requested to install helper tools through Busybox.
echo WScript.Sleep 2000 > "%temp%\\~ElevateBusyBox.vbs"
echo set UAC = CreateObject^("Shell.Application"^) >> "%temp%\\~ElevateBusyBox.vbs"
echo UAC.ShellExecute "%~dp0\\bin\\busybox.exe", "--install -s ""%~dp0\\bin""", , "runas", 0 >> "%temp%\\~ElevateBusyBox.vbs"
echo WScript.Sleep 2000 >> "%temp%\\~ElevateBusyBox.vbs"
cscript "%temp%\\~ElevateBusyBox.vbs" > nul 2>&1
del /Q /F "%temp%\\~ElevateBusyBox.vbs"
if not exist "%~dp0\\bin\\bash" ( goto noBusyBox )
goto binary

:binary
"%~dp0\\bin\\${BINARY_NAME}.exe" %*
exit /b 0

:noBusyBox
echo Error: Could not install BusyBox helper tools. Please execute the following command manually: >&2
echo "%~dp0\\bin\\busybox.exe" --install -s "%~dp0\\bin" >&2
exit /b 1
EOF

chmod +x "$BUILD/${BINARY_NAME}/${BINARY_NAME}.bat"

