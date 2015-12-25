#!/bin/bash
_mmseqs() {
	local cur prev opts
	COMPREPLY=()
	cur="${COMP_WORDS[COMP_CWORD]}"
	prev="${COMP_WORDS[COMP_CWORD-1]}"

	if [[ ${COMP_CWORD} -eq 1 ]] ; then
		COMPREPLY=( $(compgen -W "$(mmseqs shellcompletion)" -- ${cur}) )
		return 0
	fi

	if [[ ${COMP_CWORD} -gt 1 ]] ; then
		COMPREPLY=( $(compgen -f -W "$(mmseqs shellcompletion ${COMP_WORDS[1]})" -- ${cur}) )
		return 0
	fi

}
complete -F _mmseqs mmseqs
