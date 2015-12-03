# BreaKmer_processing
Scripts to process breakmer output and format it to report to collaborators

process:

run link_breakmer_files.rcprime.py on the cohort of interest
command: ./python link_breakmer_files.rcprime.py <projectID> <subprojectID>
where these refer to the naming conventions used in the paths on the isilon

this should produce a list file of the aggregates being considered. Split this list into tumor_[subprojectID_sample_list].txt and normal_[subprojectID_sample_list].txt

If you don't want to filter any breakmer calls (right now that should be the case) touch a filters.txt file in that folder.

The process_cohort_results.py should point to this folder. Run it
command: ./python process_cohort_results.py <projectID> <tumorsamplelist> <normalsamplelist>

This produces the outputs

For using the summary R script, I typically copy the first section and change the projectId variable,
start R from the directory containing the processed output files, paste in the code and it will generate:
1. breakmer.target_summary.txt
2. breakmer.sample_summary.txt

```
projectId <- 'SCB0002_TSC1-TSC2'

svs <- read.table(paste(projectId, "breakmer.df.txt", sep='_'), header=T, sep="\t", as.is=T)                                                                                                                                                                                                                                            
tFreqs <- table(svs$sampleid, svs$target_name)

targetGeneSummary <- data.frame("Nevents" = colSums(tFreqs), "Nsamples_with_event"=colSums(tFreqs>0))
write.table(targetGeneSummary, quote=F, row.names=T, col.names=T, sep="\t", file=paste(projectId, "breakmer.target_summary.txt", sep="_"))

tFreqs2 <- table(svs$sampleid, svs$rearr_type)
sampleSummary <- data.frame("Nevents"=rowSums(tFreqs), "Nindels"=tFreqs2[,1], "Nrearr"=tFreqs2[,2]) #, "Ninversions"=tFreqs2[,3])
write.table(sampleSummary, quote=F, row.names=T, col.names=T, sep="\t", file=paste(projectId, "breakmer.sample_summary.txt", sep="_"))
```
