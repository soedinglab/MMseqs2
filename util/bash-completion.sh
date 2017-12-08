#!/bin/bash
_mmseqs() {
	local cur
	COMPREPLY=()
	cur="${COMP_WORDS[COMP_CWORD]}"

	if [[ ${COMP_CWORD} -eq 1 ]] ; then
		COMPREPLY=( $(LC_COLLATE=C compgen -W "$(mmseqs shellcompletion)" -- "${cur}") )
		return 0
	fi

	if [[ ${COMP_CWORD} -gt 1 ]] ; then
		COMPREPLY=( $(LC_COLLATE=C compgen -f -W "$(mmseqs shellcompletion "${COMP_WORDS[1]}")" -- "${cur}") )
		return 0
	fi

}
complete -o plusdirs -F _mmseqs mmseqs
