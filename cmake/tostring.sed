s/^unsigned char \([A-Za-z0-9_]*\)\[\] = {[ 	]*$/#if defined(__GNUC__) || defined(__clang__)\
#pragma GCC diagnostic push\
#pragma GCC diagnostic ignored "-Woverlength-strings"\
#endif\
static const unsigned char \1[] =/
s/^};$/"";\
#if defined(__GNUC__) || defined(__clang__)\
#pragma GCC diagnostic pop\
#endif/
s/^unsigned int /static const unsigned int /
/0x/{
	s/,//g
	s/[ 	]//g
	s/0x/\\x/g
	s/^/"/
	s/$/"/
}
