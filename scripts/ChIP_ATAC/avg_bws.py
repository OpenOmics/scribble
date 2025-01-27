#!/usr/bin/env python3
import json
import argparse
import subprocess

parser = argparse.ArgumentParser(description='Script to grab bigwigs that are associated with a single group and average them to create a single bigwig. Make sure to module load deeptools before use.')
parser.add_argument('-f', required=True, help="Full path to the chrom-seek working directory")
parser.add_argument('-g', required=True, help="Name of the group as listed in the peakcall file from chrom-seek")
parser.add_argument('-r', required=True, help="Part of the individual bigwigs that tell the file format. Options are: '.sorted.RPGC.bw', '.Q5DD.RPGC.bw', or '.Q5DD.RPGC.inputnorm.bw'")

args = parser.parse_args()
folder = args.f
group = args.g
rootname = args.r

config = json.load(open(folder + "config.json"))
groupdata = config['project']['groups']

inBWs = [ folder + "bigwig/" + sample + rootname for sample in groupdata[group] ]
outBW = folder + "bigwig/" + group + rootname

bashCommand = "bigwigAverage -b " + (" ").join(inBWs) + " -o " + outBW
output = subprocess.check_output(['bash','-c', bashCommand])
