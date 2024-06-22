#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 09:48:18 2024

@author: dbanco
"""
import subprocess
import time
 
def collect_job_ids():
    """
    Collects job IDs from the output of the qstat command.

    Returns:
        list: A list of job IDs currently in the queue.
    """
    try:
        # Run the qstat command
        result = subprocess.run(['qstat'], capture_output=True, text=True, check=True)
        
        # Split the output into lines
        lines = result.stdout.strip().split('\n')
        
        # Extract job IDs from the lines (assuming the job ID is in the first column)
        job_ids = []
        for line in lines[2:]:  # Skip the header lines
            parts = line.split()
            if len(parts) > 0:
                job_id = parts[0]
                job_ids.append(job_id)
        
        return job_ids
    
    except subprocess.CalledProcessError as e:
        print(f"Error running qstat: {e}")
        return []
    
def check_job_status(job_id):
    result = subprocess.run(['qstat'], capture_output=True, text=True)
    return job_id in result.stdout

def wait_for_jobs(job_ids):
    while True:
        all_completed = True
        for job_id in job_ids:
            if check_job_status(job_id):
                all_completed = False
                break
        if all_completed:
            break
        time.sleep(5)  # Wait for 1 minute before checking again