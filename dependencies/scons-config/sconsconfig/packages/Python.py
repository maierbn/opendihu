import subprocess
from Package import Package

##
##  Handles Python interpreter C library
##
class Python(Package):

    def __init__(self, **kwargs):
        super(Python, self).__init__(**kwargs)

    def check(self, ctx):
        env = ctx.env
        ctx.Message('Checking for Python ... ')
        
        try:
          cflags = subprocess.check_output("python-config --cflags", shell=True)
          ldflags = subprocess.check_output("python-config --ldflags", shell=True)
          
          # remove "-Wstrict-prototypes" because it is only valid for c and not c++
          if "-Wstrict-prototypes" in cflags:
            startpos = cflags.index("-Wstrict-prototypes")
            length = len("-Wstrict-prototypes")
            cflags = cflags[0:startpos] + cflags[startpos+length:]
            
        except:
          ctx.Result(False)
          return False
        
        env.MergeFlags(cflags)
        env.MergeFlags(ldflags)

        ctx.Result(True)
        return True
