#!/usr/bin/python

program_description = "Script to show live results of FoldIt competition"
__author__ = "Kyle Barlow"

import os
import sys
import argparse
import time
import subprocess

# Constants
dssp_binary_location = os.path.expanduser('~/bin/mkdssp')

def check_for_new_pdbs(input_directory, previously_discovered_pdbs):
    new_pdbs = set()

    for directory_root, directory, filenames in os.walk(input_directory):
        for filename in filenames:
            if 'pdb' in filename.lower():
                full_filename = os.path.join(directory_root, filename)
                if full_filename.endswith('.gz'):
                    try:
                        subprocess.check_call( ['gunzip', full_filename] )
                    except Exception:
                        print 'Failed unzipping file', full_filename
                else:
                    if full_filename not in previously_discovered_pdbs:
                        new_pdbs.add(full_filename)

    return new_pdbs

def run_dssp(pdb_filepath):
    '''
    Returns a dictionary
      Key: Secondary structure type
      Value: List of residue numbers of that structure
    '''

    try:
        dssp_output = subprocess.check_output([dssp_binary_location, '-i', pdb_filepath])
    except:
        print 'DSSP failed on ' + pdb_filepath + '\nCheck your file if this is yours!'
        return None

    parsing = False
    return_dict = {}
    for line in dssp_output.split('\n'):
        if not parsing and line[2:25] == '#  RESIDUE AA STRUCTURE':
            parsing=True
            continue

        if parsing:
            if line[5:10].strip() != '':
                res_num = int(line[5:10])
                structure = line[16:17]
                if structure not in return_dict.keys():
                    return_dict[structure]=[]
                return_dict[structure].append(res_num)

    return return_dict

def helical_content_from_dict(dssp_dict):
    '''
    Returns total helical content as a percentage
    '''
    # Secondary structure codes as in output explanation from http://swift.cmbi.ru.nl/gv/dssp/
    total_residues = 0
    helix_residues = 0
    for key in dssp_dict:
        total_residues += len( dssp_dict[key] )
        if key == 'H':
            helix_residues += len( dssp_dict[key] )
    
    return float(helix_residues) / float(total_residues)

def longest_helix_from_dict(dssp_dict):
    '''
    Returns length of longest alpha helix in structure
    '''
    longest_helix = 0
    if 'H' in dssp_dict:
        helix_residues = sorted(dssp_dict['H'])
        if len(helix_residues) == 0 or len(helix_residues) == 1:
            return len(helix_residues)

        this_helix = 1
        last_residue = helix_residues[0]
        for residue in helix_residues[1:]:
            if last_residue+1 == residue:
                this_helix += 1
            else:
                if this_helix > longest_helix:
                    longest_helix = this_helix
                this_helix = 1
            last_residue = residue
            
    return longest_helix

def main(input_directory):
    known_pdbs = set()

    while(True):
        new_pdbs = check_for_new_pdbs(input_directory, known_pdbs)

        for new_pdb in new_pdbs:
            print '\nProcessing new pdb: %s' % os.path.basename(new_pdb)
            dssp_dict = run_dssp(new_pdb)
            longest_helix = longest_helix_from_dict(dssp_dict)
            helical_content = helical_content_from_dict(dssp_dict)
            print 'Longest helix: %d residues' % longest_helix
            print 'Total helical content: %.2f%%' % (helical_content*100.0)

        known_pdbs.update(new_pdbs)

        time.sleep(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=program_description)

    parser.add_argument('input_directory',
                        help = 'Name of directory to monitor for live addition of PDB files')

    args = parser.parse_args()

    assert( os.path.isdir(args.input_directory) )

    main(args.input_directory)
