# mpileup_coverage v1.0
The input to Varscan2 is a mpileup file. Where mate-pairs overlap there is a difference in how mpileup and sambamba select reads, meaning the coverage reports produced by Sambamba does not represent the coverage available to varscan.

Therefore, in order to ensure an amplicon has sufficient depth to pass the QC requirements, an alternative coverage report is required to sambamba, performed by this script.

The script has 4 inputs.
1. -b /path/to/bedfile. This must be a sambamba BED file, in order to include the more informative amplicon description.
2. -m /path/to/mpileup_file. This should be run with the same parameters applied by varscan (eg min BQ). 
3. -c 500. Required depth for every base in the amplicon for the amplicon to pass QC
4. -o /path/to/output.txt. path to create coverage report

The script parses through each region in the BED file, and then assesses the reported coverage for each base in the mpileup file.
If a single base has coverage below the required level, or is missing from the mpileup file the amplicon will be marked as a fail.

We have checked that softclipped bases are not counted by mpileup, so coverage will not be inflated by the primers of overlapping amplicons.

