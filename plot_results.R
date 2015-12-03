projectId <- 'SCB0002_TSC1-TSC2'

svs <- read.table(paste(projectId, "breakmer.df.txt", sep='_'), header=T, sep="\t", as.is=T)                                                                                                                                                                                                                                            
tFreqs <- table(svs$sampleid, svs$target_name)

targetGeneSummary <- data.frame("Nevents" = colSums(tFreqs), "Nsamples_with_event"=colSums(tFreqs>0))
write.table(targetGeneSummary, quote=F, row.names=T, col.names=T, sep="\t", file=paste(projectId, "breakmer.target_summary.txt", sep="_"))

tFreqs2 <- table(svs$sampleid, svs$rearr_type)
sampleSummary <- data.frame("Nevents"=rowSums(tFreqs), "Nindels"=tFreqs2[,1], "Nrearr"=tFreqs2[,2]) #, "Ninversions"=tFreqs2[,3])
write.table(sampleSummary, quote=F, row.names=T, col.names=T, sep="\t", file=paste(projectId, "breakmer.sample_summary.txt", sep="_"))

library("ggplot2")
#------------------------------------------------------------------------------
results_by_target_bar <- function(data) {

  rearr_idx <- which(data$rearr_type != 'indel')
  indel_idx <- which(data$rearr_type == 'indel')

  rearr_target_freq <- table(data$sampleid[rearr_idx], data$target_name[rearr_idx])
  indel_target_freq <- table(data$sampleid[indel_idx], data$target_name[indel_idx])

  rearr_target_freq <- colSums(rearr_target_freq > 0)
  indel_target_freq <- colSums(indel_target_freq>0)

  rearr_df <- data.frame("target_gene" = factor(names(rearr_target_freq), levels=names(rearr_target_freq)[order(rearr_target_freq, decreasing=F)], ordered=T), "rearr_freq" = as.numeric(rearr_target_freq))
  xl <- "Target captured gene name"
  yl <- "Number of samples with rearr event"
  gp <- ggplot(rearr_df[which(rearr_df[,2] > 5),], aes(x=target_gene, y=rearr_freq) ) +
               geom_bar(stat="identity") +
               labs(x=xl, y=yl, title="Target capture gene rearrangement freq (> 5 samples)") +
               theme_bw() +
               guides(fill=guide_legend(title=NULL)) +
               theme(text = element_text(size=18)) +
               theme(legend.text=element_text(size=18)) +
               theme(strip.text=element_text(size=18)) +
               coord_flip() +
               theme(axis.text=element_text(size=14), axis.title=element_text(size=16))
  ggsave(filename="target_rearr_freq_bar.pdf", plot=gp, width=14, height=12)

  indel_df <- data.frame("target_gene" = factor(names(indel_target_freq), levels=names(indel_target_freq)[order(indel_target_freq, decreasing=F)], ordered=T), "indel_freq" = as.numeric(indel_target_freq))
  xl <- "Target captured gene name"
  yl <- "Number of samples with indel event"
  gp <- ggplot(indel_df[which(indel_df[,2] > 5),], aes(x=target_gene, y=indel_freq) ) +
               geom_bar(stat="identity") +
               labs(x=xl, y=yl, title="Target gene indel freq (> 5 samples)") +
               theme_bw() +
               guides(fill=guide_legend(title=NULL)) +
               theme(text = element_text(size=18)) +
               theme(legend.text=element_text(size=18)) +
               theme(strip.text=element_text(size=18)) +
               coord_flip() +
               theme(axis.text=element_text(size=16), axis.title=element_text(size=16))
  ggsave(filename="target_indel_freq_bar.pdf", plot=gp, width=14, height=10)

  rearr_pair_freq <- table(data$sampleid[rearr_idx], data$gene_list[rearr_idx])
  rearr_pair_freq <- colSums(rearr_pair_freq > 0)
  rearr_pair_df <- data.frame("rearr_genes" = factor(names(rearr_pair_freq), levels=names(rearr_pair_freq)[order(rearr_pair_freq, decreasing=F)], ordered=T), "rearr_freq" = as.numeric(rearr_pair_freq))
  xl <- "Rearranged genes"
  yl <- "Number of samples with rearr pairs"
  gp <- ggplot(rearr_pair_df[which(rearr_pair_df[,2] > 2),], aes(x=rearr_genes, y=rearr_freq) ) +
               geom_bar(stat="identity") +
               labs(x=xl, y=yl, title="Gene rearrangements freq (> 5 samples)") +
               theme_bw() +
               guides(fill=guide_legend(title=NULL)) +
               theme(text = element_text(size=18)) +
               theme(legend.text=element_text(size=18)) +
               theme(strip.text=element_text(size=18)) +
               coord_flip() +
               theme(axis.text=element_text(size=14), axis.title=element_text(size=16))
  ggsave(filename="rearr_pair_freq_bar.pdf", plot=gp, width=14, height=12)

  rearrs <- data[rearr_idx, ]
  indels <- data[indel_idx, ]

  xl <- "Split read counts"
  yl <- "Discordant read counts"
  gp <- ggplot(rearrs, aes(x=log2(max_splitread_counts+1), y=log2(max_discread_counts+1), color=max_repeat_overlap) ) +
               geom_point(aes(size=1/min_uniq_align), alpha=0.7) +
               labs(x=xl, y=yl, title="Rearrangement evidence") +
               theme_bw() +
               guides(fill=guide_legend(title=NULL)) +
               theme(text = element_text(size=18)) +
               theme(legend.text=element_text(size=18)) +
               theme(strip.text=element_text(size=18)) +
               coord_flip() +
               theme(axis.text=element_text(size=14), axis.title=element_text(size=16))
  ggsave(filename="rearr_counts_scatter.pdf", plot=gp, width=14, height=12)

  xl <- "log2(Read counts)"
  yl <- "Contig length"
  gp <- ggplot(rearrs, aes(x=log2(max_splitread_counts+max_discread_counts), y=len_contig, color=max_repeat_overlap) ) +
               geom_point(aes(size=1/min_uniq_align), alpha=0.7) +
               labs(x=xl, y=yl, title="Rearrangement evidence vs. contig length") +
               theme_bw() +
               guides(fill=guide_legend(title=NULL)) +
               theme(text = element_text(size=18)) +
               theme(legend.text=element_text(size=18)) +
               theme(strip.text=element_text(size=18)) +
               coord_flip() +
               theme(axis.text=element_text(size=14), axis.title=element_text(size=16))
  ggsave(filename="rearr_counts_vs_len_scatter.pdf", plot=gp, width=14, height=12)

  xl <- "log2(Read counts)"
  yl <- "Contig length"
  gp <- ggplot(indels, aes(x=log2(max_splitread_counts), y=len_contig, color=subprojectid) ) +
               geom_point(size=4, alpha=0.7) +
               labs(x=xl, y=yl, title="Indel evidence vs. contig length") +
               theme_bw() +
               guides(fill=guide_legend(title=NULL)) +
               theme(text = element_text(size=18)) +
               theme(legend.text=element_text(size=18)) +
               theme(strip.text=element_text(size=18)) +
               coord_flip() +
               theme(axis.text=element_text(size=14), axis.title=element_text(size=16))
  ggsave(filename="indel_counts_vs_len_scatter.pdf", plot=gp, width=14, height=12)

  xl <- "Target genes"
  yl <- "Read count support"
  gp <- ggplot(rearrs, aes(x=target_name, y=log2(max_splitread_counts+max_discread_counts)) ) +
               geom_boxplot() +
               labs(x=xl, y=yl, title="Indel evidence vs. contig length") +
               theme_bw() +
               guides(fill=guide_legend(title=NULL)) +
               theme(text = element_text(size=18)) +
               theme(legend.text=element_text(size=18)) +
               theme(strip.text=element_text(size=18)) +
               coord_flip() +
               theme(axis.text=element_text(size=7), axis.title=element_text(size=16))
  ggsave(filename="rearr_target_counts_box.pdf", plot=gp, width=14, height=20)
}
#------------------------------------------------------------------------------
