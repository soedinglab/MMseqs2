#!/bin/sh -e
[ "$#" -ne 2 ] && echo "Please provide <shellcheckBinary> <inputPath>" && exit 1;

SHELLCHECK="$1"
if [ ! -x "$SHELLCHECK" ]; then
    exit 0
fi

INPUT="$2"
INPUT_EXT="${INPUT##*.}"

if [ "${INPUT_EXT}" = "sh" ]; then
    ${SHELLCHECK} "$2"
else
    exit 0
fi


