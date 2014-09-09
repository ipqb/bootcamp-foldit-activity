#!/usr/bin/python

program_description = "Script to show live results of FoldIt competition"
__author__ = "Kyle Barlow"

import os
import sys
import argparse
import datetime
import time
import subprocess

# Constants
dssp_binary_location = os.path.expanduser('~/bin/mkdssp')
time_print_interval = datetime.timedelta(minutes=5) # How often to print regular time updates
short_time_print_interval = datetime.timedelta(seconds=15) # How often to print short time updates
shortest_time_print_interval = datetime.timedelta(seconds=1) # How often to print shortest time updates

def ordinal(value):
    """
    Hattip: http://code.activestate.com/recipes/576888-format-a-number-as-an-ordinal/

    Converts zero or a *postive* integer (or their string 
    representations) to an ordinal value.

    >>> for i in range(1,13):
    ...     ordinal(i)
    ...     
    u'1st'
    u'2nd'
    u'3rd'
    u'4th'
    u'5th'
    u'6th'
    u'7th'
    u'8th'
    u'9th'
    u'10th'
    u'11th'
    u'12th'

    >>> for i in (100, '111', '112',1011):
    ...     ordinal(i)
    ...     
    u'100th'
    u'111th'
    u'112th'
    u'1011th'

    """
    try:
        value = int(value)
    except ValueError:
        return value

    if value % 100//10 != 1:
        if value % 10 == 1:
            ordval = u"%d%s" % (value, "st")
        elif value % 10 == 2:
            ordval = u"%d%s" % (value, "nd")
        elif value % 10 == 3:
            ordval = u"%d%s" % (value, "rd")
        else:
            ordval = u"%d%s" % (value, "th")
    else:
        ordval = u"%d%s" % (value, "th")

    return ordval

def mktime(timestring):
    return datetime.datetime.strptime(timestring, '%Y-%m-%d %I:%M%p')

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

def print_remaining_time(end_time):
    time_left = end_time - datetime.datetime.now()
    total_seconds = int(time_left.total_seconds())
    hours, remainder = divmod(total_seconds, 60*60)
    minutes, seconds = divmod(remainder, 60)
    if hours > 0:
        print '\n######## Contest ends in: %d hours, %d minutes ########' % (hours, minutes)
    else:
        if minutes > 5:
            print '\n######## Contest ends in: %d minutes ########' % (minutes)
        else:
            if minutes > 0:
                print '\n######## Contest ends in: %d minutes, %d seconds ########' % (minutes, seconds)
            else:
                print '\n######## Contest ends in: %d seconds ########' % (seconds)

def main(input_directory, end_time):
    known_pdbs = set()
    helix_lengths = []
    helical_contents = []

    print_remaining_time(end_time)
    last_time_update = datetime.datetime.now()
    time_to_compare = time_print_interval
    while(end_time > datetime.datetime.now()):
        time_since_update = datetime.datetime.now() - last_time_update
        remaining_time = end_time - datetime.datetime.now()

        if remaining_time < datetime.timedelta(seconds=30):
            time_to_compare = shortest_time_print_interval
        elif remaining_time < datetime.timedelta(minutes=2):
            time_to_compare = short_time_print_interval
        else:
            time_to_compare = time_print_interval

        if time_since_update > time_to_compare:
            print_remaining_time(end_time)
            last_time_update = datetime.datetime.now()

        new_pdbs = check_for_new_pdbs(input_directory, known_pdbs)

        for new_pdb in new_pdbs:
            print '\nFound new pdb: %s' % os.path.basename(new_pdb)
            dssp_dict = run_dssp(new_pdb)

            longest_helix = longest_helix_from_dict(dssp_dict)
            helical_content = helical_content_from_dict(dssp_dict)

            print 'Longest helix: %d residues' % longest_helix
            helix_lengths.append( (longest_helix, new_pdb) )

            print 'Total helical content: %.2f%%' % (helical_content*100.0)
            helical_contents.append( (helical_content, new_pdb) )
            

        known_pdbs.update(new_pdbs)

        time.sleep(1)

    # Contest done, print winners!
    helix_lengths.sort(reverse=True)
    helical_contents.sort(reverse=True)
    print '\nContest over - runner-up results will be shown in 10 seconds (for suspense)\n'
    time.sleep(10)

    print '#### CONTEST RESULTS!!!! ####'
    print 'Helix length category runner-ups:'
    for i, helix_length in enumerate(helix_lengths[1:5]):
        print '%s place: %s, with a longest helix length of: %d' % (ordinal(i+2),
                                                                    os.path.basename(helix_length[1]),
                                                                    helix_length[0])

    print '\nHelical content category runner-ups:'
    for i, helix_length in enumerate(helical_contents[1:5]):
        print '%s place: %s, with a total helical content of: %.4f%%' % (ordinal(i+2),
                                                                         os.path.basename(helix_length[1]),
                                                                         helix_length[0])

    time.sleep(10)
    if len(helix_lengths) > 0:
        print '\nHELIX LENGTH CHAMPION:'
        print helix_lengths[0][1]
        print '%s place: %s, with a longest helix length of: %d' % (ordinal(1),
                                                                         os.path.basename(helix_lengths[0][1]),
                                                                         helix_lengths[0][0])
    if len(helical_contents) > 0:
        print '\nHELIX CONTENT CHAMPION:'
        print helical_contents[0][1]
        print '%s place: %s, with a total helical content of: %.4f%%' % (ordinal(1),
                                                                         os.path.basename(helical_contents[0][1]),
                                                                         helical_contents[0][0])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=program_description)

    parser.add_argument('input_directory',
                        help = 'Name of directory to monitor for live addition of PDB files')
    parser.add_argument('-e', '--end_time',
                        type = mktime,
                        default = None,
                        help = 'Contest end time in format YY-MM-DD HH:MMpm')

    args = parser.parse_args()

    assert( os.path.isdir(args.input_directory) )
    if args.end_time:
        end_time = args.end_time
    else:
        end_time = datetime.datetime.now() + datetime.timedelta(seconds=30)

    assert( end_time > datetime.datetime.now() )

    main(args.input_directory, end_time)
