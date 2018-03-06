#!/usr/bin/env python
# encoding: utf-8

# File        : change.py
# Author      : Zhenbin Wu
# Contact     : zhenbin.wu@gmail.com
# Date        : 2018 Mar 01
#
# Description :

import argparse
import re
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--globalTag', dest='GT', help='input globalTag')
    parser.add_argument('--sqlite', dest='SQ', help='input sqlite for Ecal laser corretion')
    args = parser.parse_args()
    inputfile = open("l1Ntuple_%s.py" % args.GT, 'r')
    outfile = open("l1Ntuple_%s_%s.py" % (args.GT, args.SQ), 'w')
    for line  in inputfile.readlines():
        outfile.write(line)
        if re.match("^process.GlobalTag\s=\sGlobalTag.*%s.*" % args.GT, line) is not None:
            outfile.write('process.GlobalTag.toGet = cms.VPSet(\n')
            outfile.write('   cms.PSet(record = cms.string("EcalTPGLinearizationConstRcd"),\n')
            outfile.write('            tag = cms.string("EcalTPGLinearizationConst_IOV_%s_beginning_at_1"),\n' % args.SQ)
            outfile.write('            connect =cms.string("sqlite_file:%s/EcalTPG_%s_moved_to_1.db"),\n' % (os.getcwd(), args.SQ))
            outfile.write('            ),\n')
            outfile.write('   cms.PSet(record = cms.string("EcalTPGPedestalsRcd"),\n')
            outfile.write('            tag = cms.string("EcalTPGPedestals_%s_beginning_at_1"),\n' % args.SQ)
            outfile.write('            connect =cms.string("sqlite_file:%s/EcalTPG_%s_moved_to_1.db"),\n' % (os.getcwd(), args.SQ))
            outfile.write('            ),\n')
            outfile.write(')\n')
            # outfile.wirte('process.MessageLogger.cerr.FwkReport.reportEvery = 1000\n')

    outfile.write("process.TFileService.fileName= cms.string('L1Ntuple_%s_%s.root')\n" % (args.GT, args.SQ) )
