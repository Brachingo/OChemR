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

        
    return ordered

def Annotate_Reaction(data, outputdir, filename):
    with open(outputdir+filename.replace(".png",".json"), "w") as f:
        json.dump(data, f, indent=4)
        
        