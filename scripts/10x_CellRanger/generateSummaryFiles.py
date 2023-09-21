#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jul 31 2018

@author: Vicky

usage: python generateSummaryFiles.py metric_summary

"""


import glob, xlsxwriter, csv, sys, ntpath, os
from shutil import copyfile

metricsPath = 'finalreport/'
summaryPath = 'finalreport/summaries/'
def main(arg1='metric_summary'):
    createMetricsSummary(arg1)
    copyWebSummary()

def createMetricsSummary(arg1):
    try:
        os.makedirs(metricsPath)
    except OSError:
        if not os.path.isdir(metricsPath):
            raise Exception('Error: provided metricsPath, {}, is not a path!'.format(metricsPath))
    files = glob.glob('./*/outs/metrics_summary.csv')

    workbook = xlsxwriter.Workbook(metricsPath + arg1+'.xlsx')
    worksheet = workbook.add_worksheet("metrics_summary")
    worksheet.set_column(0, 2, 10.1)
    worksheet.set_column(3, 29, 12.2)

    formatNum = workbook.add_format({'num_format': '#,###'})
    formatFlt = workbook.add_format({'num_format': '#,###.#'})
    formatPer = workbook.add_format({'num_format': '0.00%'})
    formatHead = workbook.add_format({'bold': True, 'italic': True, 'text_wrap': True, 'align': 'center'})

    summaryDict = dict()
    finalheaders = list()
    for filename in files:
        with open(filename, 'r') as csvfile:
            f = csv.reader(csvfile, delimiter=',', quotechar='"')
            header = next(f)
            line = next(f)
            for index in range(len(header)):
                category = header[index]
                catdict = summaryDict.get(category, dict())
                catdict[filename] = line[index]
                summaryDict[category] = catdict
            if len(set(finalheaders).difference(header)) == 0:
                finalheaders = header
            else:
                finalheaders = header + list(set(finalheaders).difference(header))

    row = 1
    samples = list()
    for filename in files:
        worksheet.write(row, 0, filename.split('/')[1])
        samples.append(filename.split('/')[1])
        col = 1
        for category in finalheaders:
            if filename in summaryDict[category]:
                i = summaryDict[category][filename].strip('"')
                if '%' in i:
                    worksheet.write(row, col, float(i.strip('%'))/100, formatPer)
                elif '.' in i:
                    if float(i.replace(',','')).is_integer():
                        worksheet.write(row, col, float(i.replace(',','')), formatNum)
                    elif i.count('.') > 1:
                        worksheet.write(row, col, i.replace(',', ''))
                    else:
                        worksheet.write(row, col, float(i.replace(',','')), formatFlt)
                else:
                    try:
                        worksheet.write(row, col, int(i.replace(',','')), formatNum)
                    except:
                        worksheet.write(row, col, i)
            col += 1
        row += 1

    col = 1
    row = 0
    worksheet.write(0, 0, "Sample", formatHead)
    for i in finalheaders:
        worksheet.write(row, col, i, formatHead)
        col += 1

    workbook.close()

def copyWebSummary():
    try:
        os.makedirs(summaryPath)
    except OSError:
        if not os.path.isdir(summaryPath):
            raise Exception('Error: provided summaryPath, {}, is not a path!'.format(summaryPath))
    files = glob.glob('./*/outs/web_summary.html')
    for filename in files:
    	copyfile(filename, '%s/%s_web_summary.html' % (summaryPath, filename.split('/')[1]))

if __name__ == "__main__":
    if len(sys.argv) == 1:
        main()
    else:
        main(sys.argv[1])
