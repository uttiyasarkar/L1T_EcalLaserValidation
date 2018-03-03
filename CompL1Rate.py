#!/usr/bin/env python
# encoding: utf-8

# File        : CompL1Rate.py
# Author      : Zhenbin Wu
# Contact     : zhenbin.wu@gmail.com
# Date        : 2018 Mar 02
#
# Description : 


import pandas
import argparse
import re
import os

def ParseCSV(GT, SQ):
    df_menu = pandas.read_csv("results/L1Menu_%s_%s_emu.csv" % (GT, SQ))
    df_seed = pandas.read_csv("results/L1Seed_%s_%s_emu.csv" % (GT, SQ))
    df_comp = df_seed[df_seed.L1Bit < 1000][["L1SeedName","rate0"]]
    df_comp.loc[len(df_comp)] =["L1Menu", df_menu.loc[df_menu.L1Bit == 9999, "rate0"].values[0]]
    colname = [SQ if x == "rate0" else x for x in list(df_comp.columns)  ]
    df_comp.columns=colname
    return df_comp
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--globalTag', dest='GT', help='input globalTag')
    parser.add_argument('--sqlite1', dest='SQ1', help='input sqlite1 for Ecal laser corretion')
    parser.add_argument('--sqlite2', dest='SQ2', help='input sqlite2 for Ecal laser corretion')
    args = parser.parse_args()

    df_comp1 = ParseCSV(args.GT, args.SQ1)
    df_comp2 = ParseCSV(args.GT, args.SQ2)
    comp = df_comp1.merge(df_comp2)
    comp.loc[:, "diff"]= comp[args.SQ1] / comp[args.SQ2]
    print(comp)

