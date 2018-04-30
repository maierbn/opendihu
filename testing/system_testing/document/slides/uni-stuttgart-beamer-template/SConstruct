from __future__ import print_function
import os
import subprocess
import sys

# print to stderr
def print_stderr(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# generate handout.tex from talk.tex by adding "handout" to the
# Beamer package options
def make_handout_command(target, source, env):
  with open(source[0].abspath, "r") as f: tex = f.read()
  tex = tex.replace("]{beamer}", ",handout]{beamer}")
  with open(target[0].abspath, "w") as f: f.write(tex)

# action used to compress 4 slides on one sheet for the handout
def pdfnup_action(target, source, env):
  for t in target:
    path = t.get_abspath()
    subprocess.call(["pdfnup", "--nup", "2x2", "--outfile", path, path])



# use PDFLaTeX
env = Environment(tools=["pdftex"], ENV=os.environ)
conf = Configure(env)

# check for pdfnup
if conf.CheckProg("pdfnup"):
  is_nup_installed = True
else:
  print_stderr("Couldn't find pdfnup on your system.",
               "The handout will not be compressed with",
               "multiple slides per sheet.")
  is_nup_installed = False

env = conf.Finish()

# add command line flags
env.Append(PDFLATEXFLAGS="-halt-on-error")
env.Append(PDFLATEXFLAGS="-synctex=1")
env.Append(BIBERFLAGS="-q")

# main presentation (talk.tex -> talk.pdf)
talk = env.PDF("tex/talk.tex")
env.Depends(talk, "tex/beamerthemeStuttgart.sty")

# handout generation (talk.tex -> handout.tex -> handout.pdf)
env.Command("tex/handout.tex", "tex/talk.tex", make_handout_command)
handout = env.PDF("tex/handout.tex")
if is_nup_installed:
  env.AddPostAction(handout, pdfnup_action)
env.Depends(handout, talk)

# proper cleaning with scons -c
env.Clean(talk, "config.log")
for filename in sorted(os.listdir("tex")):
  if filename.endswith(".tex"):
    env.Clean(talk, ["tex/{}.{}".format(filename[:-4], extension)
                     for extension in ["aux", "nav", "out", "run.xml", "snm"]])
