"""
Transfrom OpenDiHu logs with log scopes into a html document
    cat log | python3 logscope2html.py > log.html

Useful command to reduce the log size:
    ssh pcsgs03 'cat .../log | grep -v -e "[a-zA-Z_]*[0-9]*[a-zA-Z_]*[0-9]* dof [0-9]* value" -e "create MappingBetweenMeshes" -e "before mapping:" -e "local values on ranks: " -Ee "DEBUG: (jacobian|normal index space|cofactor|geometryValues): "' | python3 ~/opendihu/scripts/logscope2html.py > log.html
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
re_highlight_bright_white_start = re.compile('^\033\[97m(?P<text>.*?)')
re_highlight_bright_white_end   = re.compile('(?P<text>.*?)\033\[0m</span>$') # we insert a span beforehand

re_build_information = re.compile('.*This is opendihu ([0-9\.]*), built [a-zA-Z0-9 ]*, C\+\+ [0-9]*, GCC [0-9\.]*, current time: (?P<currenttime>[0-9/: ]*), hostname: (.*), n ranks: [0-9]*$')

re_disconnected_slot_information = re.compile("^(DEBUG:\s*slotInformation_\[[12](To|to|-&gt;)[12]\]\[[01]\]\[\s*[0-9]+\]\s*=\s*false,)(.*)$")
re_disconnected_slot       = re.compile("^(DEBUG:)?(\s*[0-9]+(\.|\s*-&gt;)\s*)(-1)(|\s+.*)$")
re_disconnected_slot_local = re.compile("^(DEBUG:)?(\s*[0-9]+(\.|\s*-&gt;)\s*)(-2)(|\s+.*)$")
re_disconnected_slot2       = re.compile("^\s*Term[12]\.slot [0-9]+ \(.*\) -&gt; Term[12]\.slot -1 \(.*\)\s*$")
re_disconnected_slot_local2 = re.compile("^\s*Term[12]\.slot [0-9]+ \(.*\) -&gt; Term[12]\.slot -2 .*$")

re_highlight_solver_names = [
    re.compile('(?<=::setSolverDescription\(&quot;)(?P<name>.*?)(?=&quot;\))'),
    re.compile('(?<=CouplingOrGodunov\(&quot;)(?P<name>.*?)(?=&quot;\))'),
]
re_highlight_mesh_names = [
    re.compile('(?<=Mesh configuration for &quot;)(?P<name>.*?)(?=&quot;)'),
    re.compile('(?<= and name &quot;)(?P<name>.*?)(?=&quot;)'),
    re.compile('(?<=Mesh with meshName &quot;)(?P<name>.*?)(?=&quot; requested)'),
    re.compile('(?<=on mesh &quot;)(?P<name>.*?)(?=&quot;)'),
    re.compile('(?<=initialize, use functionSpace_: )(?P<name>.*?)$'),
    re.compile('(?<=functionSpace: )(?P<name>.*?)(?=, nElementsLocal: )'),
]

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

        // show explanation for slot connections on mouse hover
        .hover_uc, .hover_ucl {
            position: relative;
        }
        .hover_uc:after, .hover_ucl:after {
            visibility: hidden;
            opacity: 0;
            border-radius: 5px;
            padding: 2px 2px;
            transition: opacity 500ms ease-in-out;
        }
        .hover_uc:after {
            content: "slot not connected";
            background-color: #444;
            color: #fff;
        }
        .hover_ucl:after {
            content: "slot connected but not present in other term";
            background-color: #335;
            color: #eef;
        }
        .hover_uc:hover:after, .hover_ucl:hover:after {
            opacity: 1;
            visibility: visible;
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
bright_white_count = False

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
        # keep the build information on top of the page
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

        # verbose might span over mutliple lines
        if re_highlight_bright_white_start.match(html_line):
            html_line = re_highlight_bright_white_start.sub('\g<text>', html_line)
            bright_white_count = True
        if bright_white_count:
            html_line = f'<span style="color:#474747;">{html_line}</span>'
        if re_highlight_bright_white_end.match(html_line):
            html_line = re_highlight_bright_white_end.sub('\g<text>', html_line)
            bright_white_count = False

        # highlight names
        subst = '<span style="background-color:rgba(0,255,0,0.2);">\g<name></span>'
        for r in re_highlight_solver_names:
            html_line = r.sub(subst, html_line)
        subst = '<span style="background-color:rgba(0,0,255,0.2);">\g<name></span>'
        for r in re_highlight_mesh_names:
            html_line = r.sub(subst, html_line)

        # deemphesize unconnected slots
        html_line = re_disconnected_slot_information.sub('\g<1><span style="color:grey;">\g<3></span>', html_line)
        html_line = re_disconnected_slot      .sub('\g<1><span class="hover_uc"  style="color:grey;"   >\g<2>\g<4>\g<5></span>', html_line)
        html_line = re_disconnected_slot_local.sub('\g<1><span class="hover_ucl" style="color:#ADD8E6;">\g<2>\g<4>\g<5></span>', html_line)
        html_line = re_disconnected_slot2      .sub('<span class="hover_uc"  style="color:grey;"   >\g<0></span>', html_line)
        html_line = re_disconnected_slot_local2.sub('<span class="hover_ucl" style="color:#ADD8E6;">\g<0></span>', html_line)

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
            const min = Math.floor(ms / 1000 / 60 - (days * 24 * 60) - (hours * 60))
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
