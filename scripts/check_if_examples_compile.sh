
# get path to examples directory
EXAMPLE_PATH=$OPENDIHU_HOME/examples
if [ -z "$EXAMPLE_PATH" ]; then
  EXAMPLE_PATH=$(pwd)/../examples
fi

WORKDIR=$(pwd)

# shortcuts for colorful fonts
RED='\033[0;31m'
YELLOW='\033[0;33m'
RESET='\033[0m'
SUCCESS="1"

# prepare log files
rm -f $WORKDIR/compile_examples_failed.txt
touch $WORKDIR/compile_examples_failed.txt
echo "Testing compilation of examples at $(date)" > $WORKDIR/compile_examples_log.txt

# loop over all subdirectories of examples
for directory in $EXAMPLE_PATH/**/*; do
  if [ -d "$directory" ]; then
    if [ -f "$directory/SConstruct" ]; then
      printf "${YELLOW}Check if example $directory compiles successfully... ${RESET}\n"
      cd $directory
      
      # measure duration
      START=$(date +%s.%N)
      
      # compile example
      scons BUILD_TYPE=r
      RESULT=$?
            
      # measured duration
      END=$(date +%s.%N)
      DIFF=$(python -c "print $END - $START")
      
      # depending on success of compilation output message to log files
      if [ $RESULT -eq 0 ]; then
        echo "Success: Example $directory compiles without errors. Duration: $(date -u -d @$DIFF +%T)" | tee -a ${WORKDIR}/compile_examples_log.txt
      else
        echo "Failed:  Example $directory gives errors. Duration: $(date -u -d @$DIFF +%T)" | tee -a ${WORKDIR}/compile_examples_log.txt
        printf "\n${RED}Example $directory does not compile!${RESET}\n\n\n"
        echo "$directory" >> ${WORKDIR}/compile_examples_failed.txt
        SUCCESS="0"
      fi
    fi
  fi
done

# prepend the string "failed examples at $(date)" to the compile_examples_failed.txt file, if there were failed examples
if [ $SUCCESS = "0" ]; then
  printf '%s\n%s\n' "failed examples at $(date)" "$(cat $WORKDIR/compile_examples_failed.txt)" > $WORKDIR/compile_examples_failed.txt
fi

# output end message
echo "done"
echo ""
printf "${RED}"
cat ${WORKDIR}/compile_examples_failed.txt
printf "${RESET}"
echo ""
echo "See ${WORKDIR}/compile_examples_log.txt and ${WORKDIR}/compile_examples_failed.txt for further information."

if [ $SUCCESS = "0" ]; then
  exit 1
fi
