import subprocess, shlex

def process_command_line(cmd):
    if not isinstance(cmd, list):
        cmd = shlex.split(cmd)
    return cmd

def check_output(cmd):
    cmd = process_command_line(cmd)
    return subprocess.check_output(cmd)

def check_call(cmd):
    cmd = process_command_line(cmd)
    return subprocess.check_call(cmd)
