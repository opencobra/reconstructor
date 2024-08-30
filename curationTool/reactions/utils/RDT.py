# This Python file contains functionality to process a chemical reaction represented by SMILES strings
# using the Reaction Decoder Tool (RDT). It involves calling Java-based RDT to generate atom-atom mapping (AAM)
# and visualizing the reaction, then replacing placeholders with actual values from a dictionary,
# and finally moving the generated image to a specific directory.

import os
import shutil
import subprocess

def RDT(rxn_file_path, destination_path_png='ECBLAST_smiles_AAM.png', destination_path_rxn='ECBLAST_smiles_AAM.rxn'):
    """
    Processes a chemical reaction using RDT to generate atom-atom mapping and visualization.
    It also replaces placeholders in the reaction file with labels.

    Input:
    - reaction_smiles (str): SMILES string representing the reaction.
    - replacement_dict (dict): Dictionary mapping placeholders to actual values for labelling the png.

    Output:
    - (dict): Dictionary containing the path to the reaction visualization image and a flag indicating if it's found in VMH.
    """
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    rdt_path = os.path.join(os.path.dirname(base_dir), 'rdt-2.4.1-jar-with-dependencies.jar')
    cmd = ['java', '-jar', rdt_path, '-Q', 'RXN', '-q', rxn_file_path, '-g', '-j', 'AAM', '-f', 'XML']
    subprocess.Popen(cmd).wait()
    # destination_path_png = os.path.join('reactions/static/reactions', destination_path_png)        
    # destination_path_rxn = os.path.join('reactions/static/reactions', destination_path_rxn)
    img_path = 'ECBLAST_temp_AAM.png'
    rxn_file_path = 'ECBLAST_temp_AAM.rxn'
    shutil.move(img_path, destination_path_png)
    shutil.move(rxn_file_path, destination_path_rxn)

    # Return the constructed formula and image path
    img_path = destination_path_png.replace('media/', '')
    return {
        "visualizations": [img_path],
        "vmh_found": False
    }
