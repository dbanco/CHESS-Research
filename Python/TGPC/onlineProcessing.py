# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 12:39:51 2024

@author: dpqb1
"""

import time
import subprocess
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler

# Directory to monitor
directory_to_watch = '/path/to/directory'

# Command to execute for each new file
command = 'qsub your_script.sh'  # Replace 'your_script.sh' with the script to execute

class NewFileHandler(FileSystemEventHandler):
    def on_created(self, event):
        if not event.is_directory:
            # Execute command when a new file is created
            subprocess.run(command, shell=True)

if __name__ == "__main__":
    event_handler = NewFileHandler()
    observer = Observer()
    observer.schedule(event_handler, directory_to_watch, recursive=False)
    observer.start()

    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        observer.stop()

    observer.join()
