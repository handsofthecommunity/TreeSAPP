__author__ = 'Connor Morgan-Lang'

import os
import sys
import logging
import subprocess
from multiprocessing import Process, JoinableQueue


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
    if proc.returncode != 0:
        logging.error(cmd_list[0] + " did not complete successfully! Command used:\n" +
                      ' '.join(cmd_list) + "\nOutput:\n" + stdout)
        sys.exit(19)

    return stdout, proc.returncode


def setup_progress_bar(num_items):
    if num_items > 50:
        progress_bar_width = 50
        step_proportion = float(num_items) / progress_bar_width
    else:
        progress_bar_width = num_items
        step_proportion = 1

    sys.stdout.write("[%s]" % (" " * progress_bar_width))
    sys.stdout.write("\b" * (progress_bar_width+1))
    sys.stdout.flush()

    return step_proportion


class CommandLineWorker(Process):
    def __init__(self, task_queue, commander):
        Process.__init__(self)
        self.task_queue = task_queue
        self.master = commander

    def run(self):
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                # Poison pill means shutdown
                self.task_queue.task_done()
                break
            logging.debug("STAGE: " + self.master + "\n" +
                          "\tCOMMAND:\n" + " ".join(next_task) + "\n")
            launch_write_command(next_task)
            self.task_queue.task_done()
        return


class CommandLineFarmer:
    """
    A worker that will launch command-line jobs using multiple processes in its queue
    """

    def __init__(self, command, num_threads):
        """
        Instantiate a CommandLineFarmer object to oversee multiprocessing of command-line jobs
        :param command:
        :param num_threads:
        """
        self.max_size = 32767  # The actual size limit of a JoinableQueue
        self.task_queue = JoinableQueue(self.max_size)
        self.num_threads = int(num_threads)

        process_queues = [CommandLineWorker(self.task_queue, command) for i in range(int(self.num_threads))]
        for process in process_queues:
            process.start()

    def add_tasks_to_queue(self, task_list):
        """
        Function for adding commands from task_list to task_queue while ensuring space in the JoinableQueue
        :param task_list: List of commands
        :return: Nothing
        """
        num_tasks = len(task_list)

        task = task_list.pop()
        while task:
            if not self.task_queue.full():
                self.task_queue.put(task)
                if num_tasks > 1:
                    task = task_list.pop()
                    num_tasks -= 1
                else:
                    task = None

        i = self.num_threads
        while i:
            if not self.task_queue.full():
                self.task_queue.put(None)
                i -= 1

        return
