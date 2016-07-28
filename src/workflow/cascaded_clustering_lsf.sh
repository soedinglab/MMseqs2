#!/bin/bash
hasCommand () {
    command -v $1 >/dev/null 2>&1 || { echo >&2 "Please make sure that $1 is in $PATH."; exit 1; }
}
[ -z "$MMDIR" ] && echo "Please set the environment variable $MMDIR to your MMSEQS installation directory." && exit 1;
# check number of input variables
[ "$#" -lt 4 ] && echo "Please provide <sequenceDB> <outDB> <tmpDir> <jobname> [clustermode=0] [minseqid=30] [coverage=0.8] [logDir=.]" && exit 1;
# check if necessary files exist
[ ! -f "$1" ] &&  echo "$1 not found!" && exit 1;
[   -f "$2" ] &&  echo "$2 exists already!" && exit 1;
[ ! -d "$3" ] &&  echo "tmp directory $3 not found!" && exit 1;
hasCommand awk

jobname="$4_$RANDOM"

clustermode=${5:-0}
minseqid=${6:-30}
coverage=${7:-0.8}

# default current folder
logfolder=${8:-.}

# LSF parameters
export LSF_PAM_HOSTLIST_USE=unique

queue=mpi
runtime="48:00"
nodenumber="16"
ncore="16"

bsuboptions="-m hh -R np16 -R haswell -R ngsscratch"
mpioptions="-a openmpi"
ptileoptions="-R span[ptile=$ncore]"

bigjob () {
	if [ -n "$NOMPI" ]; then
		eval $3
		return
	fi  
    if [ "$2" != "START" ]; then
        bsub -J "$1_$jobname" -w "done('$2_${jobname}')" \
            -o "$logfolder/${jobname}_$1" -e "$logfolder/${jobname}_$1.err" \
            -q $queue -W $runtime -n $nodenumber $bsuboptions $mpioptions $ptileoptions \
            <<EOF
#!/bin/bash
mpirun.lsf $3
EOF
    else
        bsub -J "$1_$jobname" \
            -o "$logfolder/${jobname}_$1" -e "$logfolder/${jobname}_$1.err" \
            -q $queue -W $runtime -n $nodenumber $bsuboptions $mpioptions $ptileoptions \
            <<EOF
#!/bin/bash
$3
EOF
    fi
}

smalljob () {
    if [ -n "$NOMPI" ]; then
		eval $4
        return
    fi
	if [ "$2" != "START" ]; then
	    bsub -J "$1_$jobname" -w "done('$2_${jobname}')" \
			-o "$logfolder/${jobname}_$1" -e "$logfolder/${jobname}_$1.err" \
			-q $queue -W $runtime -n $3 $bsuboptions -R 'span[hosts=1]' \
			<<EOF
#!/bin/bash
$4
EOF
	else
		bsub -J "$1_$jobname" \
			-o "$logfolder/${jobname}_$1" -e "$logfolder/${jobname}_$1.err" \
			-q $queue -W $runtime -n $3 $bsuboptions -R 'span[hosts=1]' \
			<<EOF
#!/bin/bash
$4
EOF
	fi
}

########## clusteringworkflow parameters  ##########
PREFILTER1_PAR="-s 1 --search-mode 2 --max-seqs 100 --k-score 130 --z-score 400.0 --threads $ncore -v 3"
ALIGNMENT1_PAR="--min-seq-id 0.$minseqid -c $coverage --max-seqs 100 --threads $ncore -v 3"
CLUSTER1_PAR="--cluster-mode $clustermode --max-seqs 1000 --min-seq-id 0.$minseqid"

PREFILTER2_PAR="-s 2 --search-mode 2 --max-seqs 200 --k-score 110 --z-score 200.0 --threads $ncore -v 3"
ALIGNMENT2_PAR="--min-seq-id 0.$minseqid -c $coverage --max-seqs 300 --threads $ncore -v 3 "
CLUSTER2_PAR="--cluster-mode $clustermode --max-seqs 1000 --min-seq-id 0.$minseqid"

PREFILTER3_PAR="-s 4 --search-mode 1 --max-seqs 1000 --k-score 100 --z-score 50.0 --threads $ncore -v 3"
ALIGNMENT3_PAR="--min-seq-id 0.$minseqid -c $coverage --max-seqs 2000 --threads $ncore -v 3 "
CLUSTER3_PAR="--cluster-mode $clustermode --max-seqs 20000 --min-seq-id 0.$minseqid"

################ clustering step 1 ################
bigjob pref_step1 START \
	"mmseqs prefilter $1 $1 $3/pref_step1 $PREFILTER1_PAR"
bigjob aln_step1 pref_step1 \
	"mmseqs alignment $1 $1 $3/pref_step1 $3/aln_step1 $ALIGNMENT1_PAR"
smalljob clu_step1 aln_step1 $ncore \
	"mmseqs cluster $1 $3/aln_step1 $3/clu_step1 $CLUSTER1_PAR"
smalljob order_step1 clu_step1 1 \
	"cut -f1 $3/clu_step1.index > $3/order_step1"  
smalljob input_step2 order_step1 2 \
	"mmseqs order $3/order_step1 $1 $3/input_step2"

################ clustering step 2 ################
bigjob pref_step2 input_step2 \
	"mmseqs prefilter $3/input_step2 $3/input_step2 $3/pref_step2 $PREFILTER2_PAR"
bigjob aln_step2 pref_step2 \
	"mmseqs alignment $3/input_step2 $3/input_step2 $3/pref_step2 $3/aln_step2 $ALIGNMENT2_PAR"
smalljob clu_step2 aln_step2 $ncore \
	"mmseqs cluster $3/input_step2 $3/aln_step2 $3/clu_step2 $CLUSTER2_PAR"
smalljob order_step2 clu_step2 1 \
    "cut -f1 $3/clu_step2.index > $3/order_step2"
smalljob input_step3 order_step2 2 \
    "mmseqs order $3/order_step2 $1 $3/input_step3"

################ clustering step 3 ################
bigjob pref_step3 input_step3 \
	"mmseqs prefilter $3/input_step3 $3/input_step3 $3/pref_step3 $PREFILTER3_PAR" 
bigjob aln_step3 pref_step3 \
	"mmseqs alignment $3/input_step3 $3/input_step3 $3/pref_step3 $3/aln_step3 $ALIGNMENT3_PAR" 
smalljob clu_step3 aln_step3 $ncore \
	"mmseqs cluster $3/input_step3 $3/aln_step3 $3/clu_step3 $CLUSTER3_PAR" 

# merge cluster results
smalljob merge_step3 clu_step3 $ncore \
	"mmseqs mergecluster $1 $3/clu $3/clu_step1 $3/clu_step2 $3/clu_step3"

# post processing
smalljob postprocess merge_step3 1 \
	"mv -f $3/clu $2 && mv -f $3/clu.index $2.index"
