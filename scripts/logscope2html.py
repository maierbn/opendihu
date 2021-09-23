"""
Transfrom OpenDiHu logs with log scopes into a html document
    cat log | python3 logscope2html.py > log.html
"""

import sys
import re
import html

import sys

re_start = re.compile(r'^Scope Start (.*)$')
re_stop  = re.compile(r'^Scope End (.*)$')

re_highlight_file = re.compile('^(.*/)([^/]*):([0-9]+)(.*)$')
re_highlight_function = re.compile('^(.*::)([a-zA-Z_]*)(\(.*)$')

re_highlight_red    = re.compile('\033\[31m(?P<text>.*?)\033\[0m')
re_highlight_green  = re.compile('\033\[32m(?P<text>.*?)\033\[0m')
re_highlight_orange = re.compile('\033\[33m(?P<text>.*?)\033\[0m')
re_highlight_grey   = re.compile('\033\[90m(?P<text>.*?)\033\[0m')
re_highlight_purple_start = re.compile('^\033\[35m(?P<text>.*?)')
re_highlight_purple_end   = re.compile('(?P<text>.*?)\033\[0m</span>$') # we insert a span beforehand

re_build_information = re.compile('.*This is opendihu ([0-9\.]*), built [a-zA-Z0-9 ]*, C\+\+ [0-9]*, GCC [0-9\.]*, current time: (?P<currenttime>[0-9/: ]*), hostname: (.*), n ranks: [0-9]*$')

file_in = sys.stdin
stack = []

print("""
<!DOCTYPE html>
<html>
  <head>
    <title>"""+"OpenDiHu Log"+"""</title>
    <meta charset="UTF-8"/>
    <style>
        details {
            border-left: solid medium rgba(0,0,0,0.1);
        }
        details > * {
            margin-left: 1cm;
        }
        details > summary {
            margin-left: 0;
            margin-top: 0.1cm;
            margin-bottom: 0.1cm;
            background-color: rgba(0,0,0,0.1);
        }
        div {
            tab-size: 4;
            white-space: pre; /* preserve whitespace */
            font-family: Consolas,Monaco,Lucida Console,Liberation Mono,DejaVu Sans Mono,Bitstream Vera Sans Mono,Courier New, monospace;
        }
        /* empty details elements are converted to a div to avoid clicking on it */
        div.emptydetails {
            margin-top: 0.1cm;
            margin-bottom: 0.1cm;
            background-color: rgba(0,100,0,0.05);
        }
    </style>
    <script>
        function openAllDetails(open) {
            var details = document.getElementsByTagName("details");
            Array.from(details).forEach((d) => {
                d.open=open
            })
        }
    </script>
  </head>
  <body>
""")

purple_count = False

for line in file_in:
    if m := re_start.match(line):
        name = m.group(1)
        stack.append(name)
        html_name = html.escape(name)

        if mm := re_highlight_file.match(html_name):
            html_name = f'{mm.group(1)}<span style="font-weight: bold;">{mm.group(2)}</span>:{mm.group(3)}{mm.group(4)}'
        if mm := re_highlight_function.match(html_name):
            html_name = f'{mm.group(1)}<span style="font-weight: bold;">{mm.group(2)}</span>{mm.group(3)}'

        # open details by default
        print(f"""<details open><summary>{html_name}</summary>""")
    elif m := re_stop.match(line):
        name = m.group(1)
        assert name==stack[-1], "Start/End do not match: "+name
        print("</details>")
        stack.pop()
    elif m := re_build_information.match(line):
        if line[-1] == '\n':
            line = line[:-1]
        html_line = html.escape(line)
        print(f'<div style="position: sticky; top: 0px; background-color:#e0ffff;"><span style="font-weight: bold; text-decoration: underline;">{html_line}</span> <button onclick="openAllDetails(true)">Open all</button><button onclick="openAllDetails(false)">Close all</button> <span id="time"></span></div>')
        # store the log creation time
        current_time = m['currenttime']
        print(f'<script>const time_log = new Date("{current_time    }");</script>')
    else:
        if line[-1] == '\n':
            line = line[:-1]
        html_line = html.escape(line)
        # replace ansi color codes
        html_line = re_highlight_red   .sub('<span style="color:red;"   >\g<text></span>', html_line)
        html_line = re_highlight_green .sub('<span style="color:green;" >\g<text></span>', html_line)
        html_line = re_highlight_orange.sub('<span style="color:orange;">\g<text></span>', html_line)
        html_line = re_highlight_grey  .sub('<span style="color:grey;"  >\g<text></span>', html_line)

        # purple (fatal errors) might span over mutliple lines
        if re_highlight_purple_start.match(html_line):
            html_line = re_highlight_purple_start.sub('\g<text>', html_line)
            purple_count = True
        if purple_count:
            html_line = f'<span style="color:purple;">{html_line}</span>'
        if re_highlight_purple_end.match(html_line):
            html_line = re_highlight_purple_end.sub('\g<text>', html_line)
            purple_count = False

        # highlight solver names
        html_line = re.sub('::setSolverDescription\(&quot;(?P<name>.*?)&quot;\)', '::setSolverDescription(&quot;<span style="background-color:rgba(0,255,0,0.2);">\g<name></span>&quot;)', html_line)
        html_line = re.sub('CouplingOrGodunov\(&quot;(?P<name>.*?)&quot;\)', 'CouplingOrGodunov(&quot;<span style="background-color:rgba(0,255,0,0.2);">\g<name></span>&quot;)', html_line)

        print(f'<div>{html_line}</div>')


if len(stack) != 0:
    print("WARN: stack not empty:\n  "+'\n  '.join(stack), file=sys.stderr)
    for name in stack:
        print(f"</details><!-- Missing End: {name} -->")

print("""
    <script>
        function updateTime() {
            const time = document.getElementById("time")
            const ms = new Date() - time_log
            const days = Math.floor(ms / 1000 / 60 / 60 / 24)
            const hours = Math.floor(ms / 1000 / 60 / 60 - (days * 24))
            const min = Math.floor(ms / 1000 / 60 - (hours * 60))
            if (days > 0) {
                time.innerHTML = `(${days}d ${hours}h ${min}min ago)`
            } else if (hours > 0) {
                time.innerHTML = `(${hours}h ${min}min ago)`
            } else {
                time.innerHTML = `(${min}min ago)`
            }
        }
        updateTime  ()
        setInterval(updateTime, 1 * 60*1000); // every minute
    </script>
  </body>
</html>
""")
