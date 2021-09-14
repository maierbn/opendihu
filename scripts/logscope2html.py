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

file_in = sys.stdin
stack = []

print("""
<!DOCTYPE html>
<html>
  <head>
    <title>"""+"OpenDiHu Log"+"""</title>
    <meta charset="UTF-8"/>
    <style>
        details > * {
            margin-left: 1cm;
        }
        details > summary {
            margin-left: 0;
            margin-top: 0.1cm;
            margin-bottom: 0.1cm;
            background-color: rgba(0,0,0,0.1);
        }
        details > div  {
            white-space: pre; /* preserve whitespace */
            font-family: Consolas,Monaco,Lucida Console,Liberation Mono,DejaVu Sans Mono,Bitstream Vera Sans Mono,Courier New, monospace;
        }
        div {
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
  </head>
  <body>
""")

for line in file_in:
    if m := re_start.match(line):
        name = m.group(1)
        stack.append(name)
        html_name = html.escape(name)

        if mm := re_highlight_file.match(html_name):
            html_name = f'{mm.group(1)}<span style="font-weight: bold;">{mm.group(2)}</span>:{mm.group(3)}{mm.group(4)}'
        if mm := re_highlight_function.match(html_name):
            html_name = f'{mm.group(1)}<span style="font-weight: bold;">{mm.group(2)}</span>{mm.group(3)}'

        print(f"""<details><summary>{html_name}</summary>""")
    elif m := re_stop.match(line):
        name = m.group(1)
        assert name==stack[-1], "Start/End do not match: "+name
        print("</details>")
        stack.pop()
    else:
        if line[-1] == '\n':
            line = line[:-1]
        html_line = html.escape(line)
        # replace ansi color codes
        html_line = re_highlight_red   .sub('<span style="color:red;"   >\g<text></span>', html_line)
        html_line = re_highlight_green .sub('<span style="color:green;" >\g<text></span>', html_line)
        html_line = re_highlight_orange.sub('<span style="color:orange;">\g<text></span>', html_line)
        html_line = re_highlight_grey  .sub('<span style="color:grey;"  >\g<text></span>', html_line)
        # highlight solver names
        html_line = re.sub('::setSolverDescription\(&quot;(?P<name>.*?)&quot;\)', '::setSolverDescription(&quot;<span style="background-color:rgba(0,255,0,0.2);">\g<name></span>&quot;)', html_line)
        html_line = re.sub('CouplingOrGodunov\(&quot;(?P<name>.*?)&quot;\)', 'CouplingOrGodunov(&quot;<span style="background-color:rgba(0,255,0,0.2);">\g<name></span>&quot;)', html_line)

        print(f'<div>{html_line}</div>')


if len(stack) != 0:
    print("WARN: stack not empty:\n  "+'\n  '.join(stack), file=sys.stderr)
    for name in stack:
        print(f"</details><!-- Missing End: {name} -->")

print("""
  </body>
</html>
""")
