__author__ = 'Connor Morgan-Lang'

import os
import sys
import subprocess


def launch_write_command(cmd_list, collect_all=True):
    """
    Wrapper function for opening subprocesses through subprocess.Popen()
    :param cmd_list: A list of strings forming a complete command call
    :param collect_all: A flag determining whether stdout and stderr are returned
    via stdout or just stderr is returned leaving stdout to be written to the screen
    :return: A string with stdout and/or stderr text and the returncode of the executable
    """
    stdout = ""
    if collect_all:
        proc = subprocess.Popen(' '.join(cmd_list),
                                shell=True,
                                preexec_fn=os.setsid,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT)
        stdout = proc.communicate()[0].decode("utf-8")
    else:
        proc = subprocess.Popen(' '.join(cmd_list),
                                shell=True,
                                preexec_fn=os.setsid)
        proc.wait()

    # Ensure the command completed successfully
    if not validate_success(proc.returncode, cmd_list):
        sys.stderr.write("\nOutput:\n")
        sys.stderr.write(stdout)
        sys.stderr.flush()
        sys.exit(19)

    return stdout, proc.returncode


def validate_success(return_code, cmd_list):
    if return_code != 0:
        sys.stderr.write("ERROR: " + cmd_list[0] +
                         " did not complete successfully! Command used:\n" + ' '.join(cmd_list) + "\n")
        return 0
    else:
        return 1
