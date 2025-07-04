# Copyright (c) 2022 rxn4chemistry - Mark Martori Lopez

import json

"""def storeResults(d,filename,outputdir):
    # # Json output file:
    outputfile = filename.replace(".png",".json")
    # Txt output file:
    with open(outputdir+outputfile,"w") as f:
        json.dump(d,f)"""


def storeResults(d, filename, outputdir):
    outputfile = filename.replace(".png", ".json")
    with open(outputdir+outputfile, "w") as f:
        json.dump(d, f, indent=4)  # Add indent=4 for pretty printing