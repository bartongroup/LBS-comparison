# Clean PDB Script
# ================
# 
# 
# PDB File Issues
# ---------------
# 
# - Multiple models ** DONE **
# - Multiple occupancies/alternate locations  ** DONE **
#   - Pick highest occupancy, remove alternate locations  ** DONE **
#   - Set occupancies to 1.00  ** DONE **
# - Missing atoms, residues ** NOT DEALT WITH **
# - Chain breaks ** DONE **
#   - Output to file: or CA or N or C position ** DONE **
# - Selenomets to Mets ** DONE **
# - Nonstandard res to ATOM records ** DONE **

# IMPORTS
import argparse
import logging
import operator
import os
import sys
import traceback
from functools import reduce

from collections import OrderedDict

from Bio.PDB import PDBParser
from Bio.PDB import Select
from Bio.PDB.Polypeptide import PPBuilder

# CONSTANTS
PDB_LINE_TEMPLATE = '{record: <6}{serial: >5} {atom_name: ^4}{altloc: ^1}{resname: ^3} {chain_id: ^1}{resnum: >4}{icode: ^1}   {x: >8.3f}{y: >8.3f}{z: >8.3f}{occ: >6.2f}{tfac: >6.2f}          {element: >2}{charge: >2}'

# MAIN
if __name__ == '__main__':
    
    # ARGUMENT PARSING
    parser = argparse.ArgumentParser(description='''

#############
# CLEAN PDB #
#############

A program for cleaning PDB files.

Dependencies:
- Python (v2.7)
- BioPython (>= v1.60)

''', formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('pdb', type=str, help='Path to the PDB file to be cleaned.')
    parser.add_argument('-rmw', '--remove-waters', action='store_true', help='Remove waters.')
    parser.add_argument('-kh', '--keep-hydrogens', action='store_true', help='Keep hydrogens.')
    parser.add_argument('-if', '--informative_filenames', action='store_true', help='Keep a record of the flags used for cleaning the output filename.')
    
    args = parser.parse_args()
    
    pdb_path = args.pdb
    pdb_noext, pdb_ext = os.path.splitext(pdb_path)
    pdb_ext = pdb_ext.replace('.', '')
    
    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure(os.path.split(os.path.splitext(pdb_path)[0])[1], pdb_path)
    
    # OUTPUT LABEL
    output_label = 'clean'
    
    if args.informative_filenames:
        if args.remove_waters:
            output_label = output_label + '_dry'
        
        if args.keep_hydrogens:
            output_label = output_label + '_kh'
    
    # REMOVE MULTIPLE MODELS
    # BY TAKING THE FIRST MODEL
    model = structure[0]
    
    # RAISE AN ERROR FOR TOO MANY ATOMS
    if len(list(model.get_atoms())) > 99999:
        try:
            raise ValueError('More than 99999 atoms in the PDB model!')
        except:
            traceback.print_exc(file=sys.stdout)
            exit(9)
    
    # DETERMINE POLYPEPTIDES AND CHAIN BREAKS
    ppb = PPBuilder()
    polypeptides = ppb.build_peptides(model, aa_only=False)
    
    # MAKE DATA STRUCTURES FOR CHAIN POLYPEPTIDES
    chain_ids = set([x.id for x in model.child_list])
    chain_pieces = OrderedDict()
    chain_polypeptides = OrderedDict()
    chain_break_residues = OrderedDict()
    chain_sequences = OrderedDict()
    
    for chain_id in chain_ids:
        chain_pieces[chain_id] = 0
        chain_break_residues[chain_id] = []
        chain_polypeptides[chain_id] = []
    
    # GET ALL POLYPEPTIDE RESIDUES IN THE MODEL
    polypeptide_residues = []
    
    for pp in polypeptides:
        for res in pp:
            polypeptide_residues.append(res)
    
    # GET THE CHAIN_ID(S) ASSOCIATED WITH EACH POLYPEPTIDE
    polypeptide_chain_id_sets = [set([k.get_parent().id for k in x]) for x in polypeptides]
    
    for e, polypeptide_chain_id_set in enumerate(polypeptide_chain_id_sets):
    
        # WARN IF NOT JUST ONE CHAIN ID ASSOCIATED WITH THE POLYPEPTIDE
        if len(polypeptide_chain_id_set) != 1:
            logging.warn('A polypeptide had {} chains associated with it: {}'.format(len(polypeptide_chain_id_set),
                                                                                   polypeptide_chain_id_set))
    
        for polypeptide_chain_id in polypeptide_chain_id_set:
            chain_pieces[polypeptide_chain_id] = chain_pieces[polypeptide_chain_id] + 1
    
            # ADD FIRST AND LAST RESIDUE TO THE CHAIN BREAK RESIDUES (POLYPEPTIDE TERMINAL RESIDUES)
            chain_break_residues[polypeptide_chain_id] = chain_break_residues[polypeptide_chain_id] + [polypeptides[e][0], polypeptides[e][-1]]
            chain_polypeptides[polypeptide_chain_id] = chain_polypeptides[polypeptide_chain_id] + [polypeptides[e]]
            
    # POP OUT THE FIRST AND LAST RESIDUES FROM THE CHAIN BREAK RESIDUES
    # TO REMOVE THE GENUINE TERMINI
    for chain_id in chain_break_residues:
        chain_break_residues[chain_id] = chain_break_residues[chain_id][1:-1]
    
    all_chain_break_residues = reduce(operator.add, chain_break_residues.values())
    
    # MAKE THE CHAIN SEQUENCES FROM THE CHAIN POLYPEPTIDE PIECES
    for chain_id in chain_polypeptides:
    
        pp_seqs = [str(x.get_sequence()) for x in chain_polypeptides[chain_id]]
        
        if pp_seqs:
            chain_sequences[chain_id] = reduce(operator.add, pp_seqs)
    
    # WRITE OUT CLEANED PDB
    # MANY OF THE ISSUES ARE SOLVED DURING THE WRITING OUT
    with open('.'.join((pdb_noext, output_label, pdb_ext)), 'wb') as fo:
        
        atom_serial = 1
        
        for residue in model.get_residues():
        
            # REMOVE WATERS IF FLAG SET
            if args.remove_waters:
                if residue.get_full_id()[3][0] == 'W':
                    continue
            
            record = 'ATOM'
            
            # SET HETATM RECORD IF IT WAS ORIGINALLY A HETATM OR WATER
            if residue.get_full_id()[3][0] == 'W' or residue.get_full_id()[3][0].startswith('H_'):
                record = 'HETATM'
            
            # SET ATOM RECORD IF THE RESIDUE IS IN A POLYPEPETIDE
            if residue in polypeptide_residues:
                record = 'ATOM'
            
            # LOOP THROUGH ATOMS TO OUTPUT
            for atom in residue.child_list:
                
                # DEAL WITH DISORDERED ATOMS
                if atom.is_disordered():
                    atom = atom.disordered_get()
                
                # REMOVE HYDROGENS
                if not args.keep_hydrogens:
                    if atom.element.strip() == 'H':
                        continue
		    
            if atom.element.strip() == 'D': # these are bloody DEUTERIUM, we also want them out. JSU dixit.
                continue
                
                # CONVERT SELENOMETHIONINES TO METHIONINES
                if residue in polypeptide_residues and (residue.resname == 'MSE' or residue.resname == 'MET'):
                    
                    residue.resname = 'MET'
                    
                    if atom.name == 'SE' and atom.element == 'SE':
                        atom.name = 'SD'
                        atom.element = 'S'
                
                # FIX ATOM NAME BUG
                if len(atom.name) == 3:
                    atom.name = ' ' + atom.name
                
                # PDB OUTPUT
                # ATOM SERIALS ARE RENUMBERED FROM 1
                # ALTLOCS ARE ALWAYS BLANK
                # CHARGES ARE ALWAYS BLANK(?)
                # OCCUPANCIES ARE ALWAYS 1.00
                output_line = PDB_LINE_TEMPLATE.format(record=record,
                                                       serial=atom_serial,
                                                       atom_name=atom.name,
                                                       altloc=' ',
                                                       resname=residue.resname,
                                                       chain_id=residue.get_parent().id,
                                                       resnum=residue.get_id()[1],
                                                       icode=residue.get_id()[2],
                                                       x=float(atom.coord[0]),
                                                       y=float(atom.coord[1]),
                                                       z=float(atom.coord[2]),
                                                       occ=1.00,
                                                       tfac=atom.bfactor,
                                                       element=atom.element,
                                                       charge='')
                
                fo.write('{}\n'.format(output_line))
                
                atom_serial += 1
                
                # RAISE AN ERROR IF WE'VE NOW GOT TOO MANY ATOMS
                if (atom_serial - 1) > 99999:
                    
                    try:
                        raise ValueError('More than 99999 atoms in the PDB when renumbered!')
                    except:
                        traceback.print_exc(file=sys.stdout)
                        exit(9)
    
    # WRITE OUT COORDINATES FOR CHAIN BREAKS FOUND WITH THE PDB FILE
    with open('.'.join((pdb_noext, pdb_ext, 'breaks')), 'wb') as fo, \
         open('.'.join((pdb_noext, pdb_ext, 'break_residues')), 'wb') as fo2:
        
        for chain_break_residue in all_chain_break_residues:
            
            fo2.write('{},{}{}`{}\n'.format(chain_break_residue.get_parent().id,
                                             chain_break_residue.get_id()[1], # RESIDUE NUMBER
                                             chain_break_residue.get_id()[2].strip(), # INSERTION CODE
                                             chain_break_residue.resname.strip()))
            
            if 'CA' in chain_break_residue.child_dict:
                break_coord = list(chain_break_residue.child_dict['CA'].coord)
                
            elif 'N' in chain_break_residue.child_dict:
                break_coord = list(chain_break_residue.child_dict['N'].coord)
                
            elif 'C' in chain_break_residue.child_dict:
                break_coord = list(chain_break_residue.child_dict['C'].coord)
                
            else:
                raise ValueError('Chain break residue {} had no mainchain atom to extract coordinates from.'.format(chain_break_residue))
            
            fo.write('{}\n'.format(str(break_coord)))
