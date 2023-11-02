# .bashrc 
# The content of this file is to be added to the docker container dashfile
export OPENDIHU_HOME=/workspace/opendihu

export PATH=$PATH:$OPENDIHU_HOME/scripts
export PATH=$PATH:$OPENDIHU_HOME/scripts/geometry_manipulation
export PATH=$PATH:$OPENDIHU_HOME/scripts/file_manipulation

# define convenience commands for compilation
alias scons='$OPENDIHU_HOME/dependencies/scons/scons.py'
alias s='scons'
alias sd='$OPENDIHU_HOME/scripts/shortcuts/sd.sh'
alias sdd='$OPENDIHU_HOME/scripts/shortcuts/sdd.sh'
alias sddn='cd .. && scons BUILD_TYPE=d no_tests=yes no_examples=yes; cd -'
alias sdn='scons BUILD_TYPE=d no_tests=yes no_examples=yes'
alias srn='scons BUILD_TYPE=r no_tests=yes no_examples=yes'
alias sr='$OPENDIHU_HOME/scripts/shortcuts/sr.sh'
alias srd='$OPENDIHU_HOME/scripts/shortcuts/srd.sh'
alias srr='$OPENDIHU_HOME/scripts/shortcuts/srr.sh'
alias mkor='$OPENDIHU_HOME/scripts/shortcuts/mkor.sh'
alias mkorn='$OPENDIHU_HOME/scripts/shortcuts/mkorn.sh'
alias mkod='$OPENDIHU_HOME/scripts/shortcuts/mkod.sh'
alias mkodn='$OPENDIHU_HOME/scripts/shortcuts/mkodn.sh'
alias mkordn='$OPENDIHU_HOME/scripts/shortcuts/mkordn.sh'



