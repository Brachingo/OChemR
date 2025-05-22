import json

"""
The following functions are very simple:

    1. Order_Reaction(): Used to retrive the dictionary of the reaction from the initial inference 
                         and reorganize it in question of the Arrow IDs from the reaction to ensure that the
                         order of the reaction is correct.
                         
    2. Annotate_Reaction(): Once the order of the reaction is correct, this function will annotate the reaction with the
                            correct information. The function will take the dictionary of the ordered reaction and 
                            annotate it properly in a JSON file.

"""

def Order_Reaction(data):
    ordered = dict(
            sorted(
                data.items(),
                key=lambda pair: int(pair[0].replace('arrow',''))
                )
            )
    keys = list(ordered.keys())
    for i in range(1, len(keys)):
        prev_key = keys[i - 1]
        curr_key = keys[i]
        ordered[curr_key]['prev_mol'] = ordered[prev_key]['post_mol']
        
    return ordered

def Annotate_Reaction(data, outputdir, filename):
    with open(outputdir+filename.replace(".png",".json"), "w") as f:
        json.dump(data, f, indent=4)
        
        