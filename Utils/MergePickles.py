# should run as is if you've kept MakePickleFile and PickleScripts/*.sh are all unchanged
# python3 MergePickles.py

import os
import uproot3
import pandas as pd
import numpy as np
import time
import sys

def combine_files(start,end,increment,folder = "./Data/Merged_charge_122420/", side = 'A'):
    output = []
    start_time = time.time()
    for i in range(start, end, increment):
        if i % 999999 == 0:
            print("event", i, "time", time.time() - start_time)
            start_time = time.time()
        print("reading file: ",f"{side}_{i}.pickle")
        output.append(pd.read_pickle(folder + f"{side}_{i}.pickle"))
    data = pd.concat(output).astype(float)
    print("data: ", data)
    return data


if __name__ == '__main__':

    print("Program starts")
    start = 9999
    end = 1000000
    increment = 10000
    folder =  "/projects/engrit/jzcapa/Users/Mike/condorJobs/ToyFermi_qqFibers_LHC_noPedNoise/pickle/"
    for side in ["A","B"]:
        merged_file = combine_files(int(start),int(end),int(increment),folder,side)
        print(f"merged_file_{side}: ", merged_file)
        merged_file.to_pickle("./ToyFermi_qqFibers_LHC_noPedNoise"+side+".pickle")
