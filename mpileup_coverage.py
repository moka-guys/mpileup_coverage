######
# This script is used to identify regions in a BED file which are nt covered at the required read depth
# in a mpileup file. Any regions not covered 100% at the required depth are reported to fail.
# A Jones 05/10/2020
######

import argparse
import os
import sys

def cli_arguments(args):
    """Parses command line arguments.
    Args:
        args: A list containing the expected commandline arguments. Example:
            ['mpileup_coverage.py', '-b', 'path/to/bedfile/', '-m',
             'path/to/mpileup/', '-o', '/path/to/output_file', '-c', '600']
    Returns:
        An argparse.parser object with methods named after long-option command-line arguments. Example:
            -b "path/to/bedfile" --> parser.parse_args(args).bedfile
    """
    # Define arguments.
    parser = argparse.ArgumentParser()
    # The runfolder string argument is immediately passed to os.path.expanduser using the *type* argument, and this
    # value is contained as the .runfolder() method in the object returned by parser.parser_args().
    # Os.path.expanduser allows expands tilde signs (~) to a string containing the user home directory.
    parser.add_argument('-b', '--bedfile', required=True, help='path_to_bedfile', type=os.path.expanduser)
    parser.add_argument('-m', '--mpileup', required=True, help='path to mpileup file', type=os.path.expanduser)
    parser.add_argument('-c', '--coverage', required=True, help='The minimum coverage required for an amplicon')
    parser.add_argument('-o', '--output_file', required=True, help='path to output file', type=os.path.expanduser)
    
    # Collect arguments and return
    return parser.parse_args(args)

def read_bedfile(args):
    """
    This function parses an input sambamba bed file and captures all bases to be assessed for coverage.
    A dictionary (region_dict) contains a entry for each region, with the value as a further dictionary containing the chromosome, start, stop and description
    A list (position_list) containing each base in the form of chr:pos used to quickly parse the mpileup file
    Input: 
        command line arguments

    Returns:
        - region_dict (dictionary)
        - position_list(list) 
    """
    position_list = []
    region_dict = {}
    with open(args.bedfile,'r') as bedfile:
        for line in bedfile.readlines():
            # check it's a sambama file, needed for the amplicon description
            assert len(line.split("\t")) == 8 , "has sambamba bed file been supplied???"
            # name relevant columns
            chr, start, stop, combined_genomic_coord, col5, col6, amplicon_description, entrez_geneid = line.split("\t")
            # create dictionary entry using genomic coord, as this will be unique (amplicon_description may not be)
            region_dict[combined_genomic_coord] = {"chr":str(chr),"start":int(start),"stop":int(stop), "description":amplicon_description}
            # add each base to the position_list in form chr:pos
            for base in range(int(start),int(stop)):
                position_list.append(chr+":"+str(base))
    return region_dict,position_list

def parse_mpileup(args,position_list):
    """
    This function parses the input mpilup file and extracts the relevant columns for the bases in the BED file (using position_list)
    A list (mpileup_list) is generated containing a tupe dictionary for each base, with values (chr, pos, depth)
    Input: 
        - command line arguments
        - position_list

    Returns:
        - mpileup_list (list) 
    """
    mpileup_list=[]
    with open(args.mpileup,'r') as mpileup:
        for line in mpileup.readlines():
            chr, pos, ref, depth, base_call_list, qual_list = line.split("\t")
            # check if it's a base in the bed file
            if str(chr)+":"+str(pos) in position_list:
                mpileup_list.append((str(chr),int(pos),int(depth)))
    return mpileup_list

def region_coverage(region_dict, mpileup_list, args):
    """
    This function loops through each region in the region_dict, and parses the mpileup_list, checking if the read depth at that base is
    below the minimum required coverage (provided argument).
    It checks that the base is present in the mpileup file, if not it will mark it as having low coverage (in case the mpileup file was 
    not run with -a)
    A key pair is added to the amplicon dictionary in region_dict, recording a boolean denoting if the amplicon is not completely covered
    at the required depth (low_coverage)
    A count for the number of amplicons completely covered above the expected value, and those not covered completelyis added to region_dict.
    Input: 
        - command line arguments
        - mpileup_list
        - region_dict

    Returns:
        - region_dict (dictionary) 
    """
    low_coverage_count = 0
    ok_coverage_count = 0
    for region in region_dict:
        # set flags which can be changed to reflect an action.
        # low coverage is used to flag the amplicon if a base is found to be below the required cutoff.
        low_coverage = False
        # take into account zero based, open ended BED file coordinates by adding one to start
        for base in range(region_dict[region]["start"]+1, region_dict[region]["stop"]):
            # set flag to check if the base is included in the mpileup file
            seen = False
            for mpileup_base in mpileup_list:
                if region_dict[region]["chr"] == mpileup_base[0] and base == mpileup_base[1]:
                    seen = True
                    if mpileup_base[2] < int(args.coverage):
                        low_coverage = True
            # if the base has not been seen after parsing the mpileup we can't pass the amplicon so mark as low_coverage
            if not seen:
                low_coverage = True

        # add the result to the dictionary        
        region_dict[region]["low_coverage"] = low_coverage
        
        # add to counts of regions with low or ok coverage
        if low_coverage:
            low_coverage_count += 1
        else:
            ok_coverage_count += 1
    
    # check all regions have been counted
    assert low_coverage_count+ ok_coverage_count == len(region_dict)
    
    # add to dictionary as a seperate field
    region_dict["low_coverage_count"] = low_coverage_count
    region_dict["ok_coverage_count"] = ok_coverage_count
    return region_dict

def report_low_covered_regions(region_dict,args):
    """
    This function parses the region_dict and writes the coverage report. 
    Any amplicons not covered sufficiently are listed, and a count of number of passing amplicons is also stated.
    Input: 
        - command line argument
        - region_dict

    Returns:
        - none
    """
    with open(args.output_file,'w') as output_file:
        output_file.write("The listed amplicons were not completely covered at the required coverage (%sX)\n" % args.coverage)
        for region in region_dict:
            # want to skip these entries in the dictionary as they don't relate to amplicons
            if region != "low_coverage_count" and region != "ok_coverage_count":
                if region_dict[region]["low_coverage"] == True:
                    output_file.write("\t".join([region_dict[region]["chr"],str(region_dict[region]["start"]),str(region_dict[region]["stop"]),region_dict[region]["description"]+"\n"]))
        output_file.write("The remaining %s amplicons were covered above the required coverage (%sX)\n" % (str(region_dict["ok_coverage_count"]),args.coverage))
        

def main(args):
    # Get command line arguments
    parsed_args = cli_arguments(args)
    # read bed file to get list of amplicons/bases to calculate coverage for
    region_dict, position_list = read_bedfile(parsed_args)
    # extract relevant lines and columsn from mpileup file
    mpileup_list = parse_mpileup(parsed_args,position_list)
    # flag any amplicons with insufficient coverage
    region_dict = region_coverage(region_dict, mpileup_list, parsed_args)
    # write coevrage report
    report_low_covered_regions(region_dict,parsed_args)

if __name__ =="__main__":
    main(sys.argv[1:])