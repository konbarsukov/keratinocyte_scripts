#!/usr/bin/python

import subprocess
import sys

# S3 bucket url 
s3_macs_enriched_bucket = 's3://aif-crc-input/data/macsEnriched'
# macsEnriched directory 
macs_enriched_folder = '../macsEnriched'

def main():
    subprocess.call(['mkdir', macs_enriched_folder])
    
    bash_command = "aws s3 sync {} {}".format(s3_macs_enriched_bucket, macs_enriched_folder)
    subprocess.check_call(bash_command.split())

if __name__=="__main__":
    main()
