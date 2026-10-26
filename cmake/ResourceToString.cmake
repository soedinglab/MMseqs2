set(BYTE "\\\\x[0-9a-f][0-9a-f]")
set(GROUP "")
foreach (I RANGE 1 12)
    string(APPEND GROUP "${BYTE}")
endforeach ()

file(READ "${IN}" HEX HEX)
string(LENGTH "${HEX}" LEN)
math(EXPR LEN "${LEN} / 2")
string(REGEX REPLACE "(..)" "\\\\x\\1" BODY "${HEX}")
string(REGEX REPLACE "(${GROUP})" "\"\\1\"\n" BODY "${BODY}")
# The last line holds fewer than 12 bytes unless the size divides evenly.
string(REGEX REPLACE "((${BYTE})+)$" "\"\\1\"" BODY "${BODY}")
string(REGEX REPLACE "\n$" "" BODY "${BODY}")

file(WRITE "${OUT}" "\
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored \"-Woverlength-strings\"
#endif
static const unsigned char ${NAME}[] =
${BODY}
\"\";
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif
static const unsigned int ${NAME}_len = ${LEN};
")
