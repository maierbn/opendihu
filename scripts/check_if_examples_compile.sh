# get path to examples directory
if [ -z "$OPENDIHU_HOME" ]; then
  OPENDIHU_HOME=$(pwd)/..
fi
EXAMPLE_PATH=$OPENDIHU_HOME/examples
WORKDIR=$(pwd)

# measure duration
START_TOTAL=$(date +%s.%N)

# shortcuts for colorful fonts
RED='\033[0;31m'
YELLOW='\033[0;33m'
RESET='\033[0m'
SUCCESS="1"

# prepare log files
rm -f $WORKDIR/compile_examples_failed.txt
touch $WORKDIR/compile_examples_failed.txt
echo "Testing compilation of examples at $(date)" > $WORKDIR/compile_examples_log.txt

check_example(){

  directory=$1
  
  printf "${YELLOW}Check if example $directory compiles successfully... ${RESET}\n"
  cd $directory
    
  # remove cached scons values
  rm -rf .scon*

  # measure duration
  START=$(date +%s.%N)
    
  # compile example
  python3 $OPENDIHU_HOME/dependencies/scons/scons.py BUILD_TYPE=r
  RESULT=$?
    
  # measured duration
  END=$(date +%s.%N)
  DIFF=$(python -c "print($END - $START)")
    
  # depending on success of compilation output message to log files
  if [ $RESULT -eq 0 ]; then
    echo "Success: Example $directory compiles without errors. Duration: $(date -u -d @$DIFF +%T)" | tee -a ${WORKDIR}/compile_examples_log.txt
  else
    echo "Failed:  Example $directory gives errors. Duration: $(date -u -d @$DIFF +%T)" | tee -a ${WORKDIR}/compile_examples_log.txt
    printf "\n${RED}Example $directory does not compile!${RESET}\n\n\n"
    echo "$directory" >> ${WORKDIR}/compile_examples_failed.txt
    SUCCESS="0"
  fi
}

# number of parallel jobs to build examples
N=$(expr `nproc --all` / 2)
i=0
i_end=$(python -c "print($N-1)")    # i_end = N-1

echo "Compile all examples in $EXAMPLE_PATH, using N=$N processes" | tee -a $WORKDIR/compile_examples_log.txt

# loop over all subdirectories of examples
for directory in `find $EXAMPLE_PATH -type d`; do
  if [ -d "$directory" ]; then
    if [ -f "$directory/SConstruct" ]; then

      # only run with N processes in parallel,
      # wait if i == N-1
      if [ "$i" -eq "$i_end" ]; then 
        wait
      fi

      i=$(python -c "print(($i+1)%$N)")   # advance i by one and wrap around
      
      # check if example in current directory compiles
      check_example "$directory" &
    fi
  fi
done

wait

# prepend the string "failed examples at $(date)" to the compile_examples_failed.txt file, if there were failed examples
if [ $SUCCESS = "0" ]; then
  printf '%s\n%s\n' "failed examples at $(date)" "$(cat $WORKDIR/compile_examples_failed.txt)" > $WORKDIR/compile_examples_failed.txt
fi

END_TOTAL=$(date +%s.%N)
DIFF_TOTAL=$(python -c "print($END_TOTAL - $START_TOTAL)")

# output end message
echo "done after $(date -u -d @$DIFF_TOTAL +%T)"
echo ""
printf "${RED}"
cat ${WORKDIR}/compile_examples_failed.txt
printf "${RESET}"
echo ""
echo "See ${WORKDIR}/compile_examples_log.txt and ${WORKDIR}/compile_examples_failed.txt for further information."

if [ $SUCCESS = "0" ]; then
  exit 1
fi
