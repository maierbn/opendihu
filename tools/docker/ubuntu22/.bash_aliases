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

# optional convenience aliases
alias ..="cd .."
alias ...="cd ../.."
alias ....="cd ../../.."
alias l='ll -h'
alias d='cd build_debug'
alias r='cd build_release/'
alias out='cd out'

# command coloring
parse_git_branch() {
 git branch 2> /dev/null | sed -e '/^[^*]/d' -e 's/* \(.*\)/(\1)/'
}
PS1='\[\e[1;32m\]\u@\h\[\e[m\]\[\e[37m\] $(parse_git_branch)\[\e[m\] \[\e[1;30m\]\w\[\e[m\] \[\e[1;32m\]\n\$\[\e[m\]\[\e[1;37m\]'
trap '[[ -t 1 ]] && tput sgr0' DEBUG
