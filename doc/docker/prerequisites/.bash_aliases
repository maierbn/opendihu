# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific aliases and functions
alias ssh="ssh -Y" 
alias big="du -ah . | sort -n -r | head -n 10"
alias bashrc="vi ~/.bashrc"
alias rs=". ~/.bashrc"
alias ,="cd ."
alias ..="cd .."
alias ,,="cd .."
alias ..y="cd .."
alias ...="cd ../.."
alias ....="cd ../../.."
alias .....="cd ../../../.."
alias ......="cd ../../../../.."
alias .......="cd ../../../../../.."
alias ........="cd ../../../../../../.."
alias .........="cd ../../../../../../../.."
alias ..........="cd ../../../../../../../../.."
alias l='ll -h'
alias pd='popd'
alias gr="grep -iInr"
alias cpr="cp -r"
alias cp-r="cp -r"
alias dc='cd'
alias mk='make clean && make'
alias ccc='rm -rf CMakeFiles/ CMakeCache.txt cmake_install.cmake Makefile CTestTestfile.cmake export mpi_verification packaging Tests support'
alias build='. build.sh'
alias biuld='. build.sh'
alias rmnpy='rm *.npy *.vtr'
alias itg='git'
alias igt='git'
alias gti='git'
alias dusch='du -sck * | sort -n'

alias gdbrun='gdb -ex=run --args '
alias k='kill %%'
alias kk='kill %% && kill %%'
alias kkk='kill %% && kill %% && kill %%'
alias kkkk='kill %% && kill %% && kill %% && kill %%'
alias s='scons'
alias sd='scons BUILD_TYPE=d'
alias sdd='cd .. && scons BUILD_TYPE=d; cd -'
alias sddn='cd .. && scons BUILD_TYPE=d no_tests=yes no_examples=yes; cd -'
alias sdn='scons BUILD_TYPE=d no_tests=yes no_examples=yes'
alias sr='scons BUILD_TYPE=r'
alias srr='cd .. && scons BUILD_TYPE=r; cd -'
alias thop='htop'
alias pdflatex='pdflatex --synctex=1'
alias ptyhon='python'
alias gbd='gdb'
alias memcheck='valgrind --tool=memcheck --suppressions=/workspace/valgrind-python.supp'
alias ubild='cd build_debug'
alias hotp='htop'
alias c='cd build_debug/'
alias b='cd build_debug/'
alias r='cd build_release/'


#---- coloring ------------------------
force_color_prompt=yes
color_prompt=yes

# Add git branch if its present to PS1
parse_git_branch() {
 git branch 2> /dev/null | sed -e '/^[^*]/d' -e 's/* \(.*\)/(\1)/'
}
#färbt Terminal passend ein
if [ "$color_prompt" = yes ]; then
# lbue
PS1='\[\e[1;34m\]\u@\h\[\e[m\]\[\e[1;37m\] $(parse_git_branch)\[\e[m\] \[\e[1;30m\]\w\[\e[m\] \[\e[1;34m\]\n\$\[\e[m\]\[\e[1;37m\]'
# red
# PS1='\[\e[1;32m\]\u@\h\[\e[m\]\[\e[1;37m\] $(parse_git_branch)\[\e[m\] \[\e[1;30m\]\w\[\e[m\] \[\e[1;32m\]\n\$\[\e[m\]\[\e[1;31m\]'
# PS1='\[\e[1;32m\]\u@\h\[\e[m\] $(parse_git_branch)\[\e[m\] \[\e[1;32m\]\n\$\[\e[m\]\[\e[1;37m\]'
else
 PS1='${debian_chroot:+($debian_chroot)}\u@\h:\w$(parse_git_branch)\$ '
fi
unset color_prompt force_color_prompt

#Stellt sicher, dass Ausgabe "normal" formatiert ist
#trap 'printf "\e[0m" "$_"' DEBUG
trap '[[ -t 1 ]] && tput sgr0' DEBUG

#-------- functions ------------
[ -z "$PS1" ] && return,	# if this is in a script, don't use the following aliases
function cd() {
    if [ -d "$@" ]; then
        #echo -n "Stack: "
        pushd "$@" > /dev/null
        #builtin cd "$@"
        pwd
        echo " "
        ls
    else
        builtin cd "$@"
    fi
}
function popd() {
    builtin popd "$@" > /dev/null
    pwd
    echo " "
    ls
}

