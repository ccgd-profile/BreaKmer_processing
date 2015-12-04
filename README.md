# BreaKmer_processing
Scripts to process breakmer output and format it to report to collaborators

process (in development):

run link_breakmer_files.rcprime.py on the cohort of interest
command: ./python link_breakmer_files.rcprime.py <projectID> <subprojectID>
where these refer to the naming conventions used in the paths on the isilon

this should produce a list file of the aggregates being considered. Split this list into tumor_[subprojectID_sample_list].txt and normal_[subprojectID_sample_list].txt

If you don't want to filter any breakmer calls (right now that should be the case) touch a filters.txt file in that folder.

The process_cohort_results.py should point to this folder. Run it
command: ./python process_cohort_results.py <projectID> <tumorsamplelist> <normalsamplelist>

This produces the outputs

filtered events section in the Summary tab:
only mentioned events that are requested to be filtered out.
Always mentioned Any events mapped to non-primary chromosomes (1-22,X,Y)
