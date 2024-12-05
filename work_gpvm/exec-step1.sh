#!/bin/bash
DIR_MACRO=$(dirname $(readlink -f $BASH_SOURCE))
DIR_DST=/pnfs/e1039/scratch/users/kenichi/RecoData2024/reco-20241102
##
## Functions
##
function PrintHelp {
    cat <<-EOF
	$(basename $0):  Script to execute the 1st analysis step on local or GRID.

	Typical usage:
	  $0 -n 3
	  Analyze the first three runs on local computer.

	  $0 -n 4-100
	  Analyze the 4th-100th runs on local computer.

	  $0 -g -n 101-
	  Analyze the 101th-last runs on grid.

	Options:
	  -l list_run.txt   | Read the run list from 'list_run.txt'.
	  -n 4-100          | Set the number (range) of the runs to be analyzed.
	  -e 9999           | Analyze only 9999 events per run.
	  -o                | Allow to overwrite existing output files.
	  -g                | Use the Grid computing.
	EOF
}

##
## Main
##
N_EVT=0
RUN_LIST=
RUN_B=0
RUN_E=47
DO_OVERWRITE=yes
USE_GRID=no
OPTIND=1
while getopts ":l:n:e:og" OPT ; do
    case $OPT in
        l ) RUN_LIST=$OPTARG ;;
        n )    RUN_E=$OPTARG ;;
        e )    N_EVT=$OPTARG ;;
        o ) DO_OVERWRITE=yes ;;
        g ) USE_GRID=yes     ;;
        * ) PrintHelp ; exit ;;
    esac
done

##
## Make a list of runs to be analyzed
##
echo "Runs to be analyzed:"
FN_LIST=list_run_spill.txt
LIST_RUN=( $(cut -f 1 $FN_LIST | uniq) )
N_RUN=${#LIST_RUN[*]}

if [ "${RUN_E%-*}" != "$RUN_E" ] ; then # Contain '-'
    RUN_B=${RUN_E%-*} # Before '-'
    RUN_E=${RUN_E#*-} # After '-'
fi
test -z $RUN_B || test $RUN_B -lt 1      && RUN_B=1
test -z $RUN_E || test $RUN_E -gt $N_RUN && RUN_E=$N_RUN

shift $((OPTIND - 1))
echo "N_RUN        = $N_RUN"
echo "N_EVT        = $N_EVT"
echo "RUN_B...E    = $RUN_B...$RUN_E"
echo "DO_OVERWRITE = $DO_OVERWRITE"
echo "USE_GRID     = $USE_GRID"

##
## Set up the working directory
##

if [ $USE_GRID == yes ]; then
    DIR_DATA=/pnfs/e1039/scratch/users/$USER/e906-root-ana
    DIR_WORK=$DIR_DATA/$JOB_NAME
    ln -nfs $DIR_DATA data # for convenience
else
    DIR_WORK=$DIR_MACRO/scratch
fi

echo "DIR_WORK = $DIR_WORK"

cd $DIR_MACRO
mkdir -p $DIR_WORK
rm -f    $DIR_WORK/input.tar.gz
tar czf  $DIR_WORK/input.tar.gz  *.C ../setup.sh ../inst ../opts

for (( RUN_I = $RUN_B ; RUN_I < $RUN_E ; RUN_I++ )) ; do
	
	RUN=${LIST_RUN[((RUN_I-1))]}
    	RUN6=$(printf "%06d" $RUN)
    	DIR_WORK_RUN=$DIR_WORK/run_$RUN6
	
	if [ -e $DIR_WORK_RUN ] ; then
                echo -n "  DIR_WORK_RUN already exists."
                if [ $DO_OVERWRITE = yes ] ; then
                	echo "  Clean up."
                        rm -rf $DIR_WORK_RUN
                else
                        echo "  Skip."
                        continue
                fi
        fi

        mkdir -p $DIR_WORK_RUN/out
	cp -p $DIR_MACRO/{exec-step1-sub.sh,MacroCommon.h} $DIR_WORK_RUN

	##
	## Loop over the spills
	##
	FN_LIST_IN=list_input.txt
	for SPILL in $(awk "{if (\$1==$RUN) print \$2;}" $DIR_MACRO/$FN_LIST) ; do
		FNAME=run_${RUN6}_spill_$(printf '%09d' $SPILL)_spin_reco.root
		echo -e "$SPILL\t$FNAME"
	done >$DIR_WORK_RUN/$FN_LIST_IN

    	echo "----------------------------------------------------------------"

    	if [ $USE_GRID == yes ] ; then
		CMD="/e906/app/software/script/jobsub_submit_spinquest.sh"
		CMD+=" --expected-lifetime='short'" # medium=8h, short=3h, long=23h
		CMD+=" -L $DIR_WORK_RUN/log_exec-step1-sub.txt"
		CMD+=" -f $DIR_WORK/input.tar.gz"
		CMD+=" -f $FN_RAW"
		CMD+=" -f $FN_REC"
		CMD+=" -d OUTPUT $DIR_WORK_RUN/out"
		CMD+=" file://$DIR_WORK_RUN/exec-step1-sub.sh $DATASET $RUN $(basename $FN_RAW) $(basename $FN_REC) $N_EVT"
		$CMD |& tee $DIR_WORK_RUN/log_jobsub_submit.txt
		RET_SUB=${PIPESTATUS[0]}
		test $RET_SUB -ne 0 && exit $RET_SUB
    	else
		export  CONDOR_DIR_INPUT=$DIR_WORK_RUN/in
		export CONDOR_DIR_OUTPUT=$DIR_WORK_RUN/out
		mkdir -p $DIR_WORK_RUN/in
		cp -p $DIR_WORK/input.tar.gz $DIR_WORK_RUN/in
		cp -a $DIR_WORK_RUN/$FN_LIST_IN $DIR_WORK_RUN/in
		while read SPILL FNAME ; do
           		ln -s $DIR_DST/run_$RUN6/spill_00$SPILL/out/$FNAME $DIR_WORK_RUN/in/$FNAME
        	done <$DIR_WORK_RUN/$FN_LIST_IN
		#ln -nfs $FN_RAW $DIR_WORK_RUN/in
		#ln -nfs $FN_REC $DIR_WORK_RUN/in
		mkdir -p $DIR_WORK_RUN/exe
		cd       $DIR_WORK_RUN/exe
		$DIR_WORK_RUN/exec-step1-sub.sh $RUN $FN_LIST_IN $N_EVT |& tee $DIR_WORK_RUN/log_exec-step1-sub.txt
    	fi
done
