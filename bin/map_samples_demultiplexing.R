# assign demultiplexing IDs to samples

# library
    library(data.table)
    library(ggplot2)
    library(stringr)
    library(ggpubr)

# input data
    sample_pos = fread('/project/holstegelab/Share/pacbio/data_processed/HLA_project/standard_pipeline/demultiplex_plate3/assign_samples/sample_positions_plate4.txt', h=T, stringsAsFactors=F)
    barcod_pos = fread('/project/holstegelab/Share/pacbio/data_processed/HLA_project/standard_pipeline/demultiplex_plate3/assign_samples/barcode_positions.txt', h=F, stringsAsFactors=F)
    counts_info = fread('/project/holstegelab/Share/pacbio/data_processed/HLA_project/standard_pipeline/demultiplex_plate3/m64346e_221010_124738.ccs.primrose.hifi.demultiplexed.lima.counts', h=T, stringsAsFactors=F)

# main
# assign column names and rownames to barcode_pos
    rownames(barcod_pos) = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H')
    colnames(barcod_pos) = as.character(seq(1, 12))
# reshape barcodes
    barcod_df = data.frame()
    for (i in 1:nrow(barcod_pos)){
        for (j in 1:ncol(barcod_pos)){
            tmp = data.frame(barcode = as.character(barcod_pos[i, ..j]), position = paste(rownames(barcod_pos)[i], colnames(barcod_pos)[j], sep=""))
            barcod_df = rbind(barcod_df, tmp)
        }
    }
# map samples
    table(sample_pos$Slot %in% barcod_df$position)
    #colnames(sample_pos)[22] = 'Volume2'
    samples_mapped = merge(sample_pos, barcod_df, by.x = 'Slot', by.y = 'position')
    head(samples_mapped)
# also add counts
    table(samples_mapped$barcode %in% counts_info$IdxFirstNamed)            # H10, 11 and 12 are missing
    samples_mapped_with_counts = merge(samples_mapped, counts_info, by.x = 'barcode', by.y = 'IdxCombinedNamed', all.x = T)
    write.table(samples_mapped_with_counts, '/project/holstegelab/Share/pacbio/data_processed/HLA_project/standard_pipeline/demultiplex_plate3/assign_samples/mapping_samples_plate.txt', quote=F, row.names=F, col.names=F, sep = "\t")
# plots
    p0 = ggplot(samples_mapped_with_counts, aes(x = Slot, y = Counts, fill = gender)) + geom_bar(stat = 'identity') + facet_grid(cols = vars(pheno), scale = 'free') + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
    p1 = ggplot(samples_mapped_with_counts, aes(x = MeanScore, y = Counts, color = pheno)) + geom_jitter(size = 6, alpha = 0.7, width = 0.25)
    p2 = ggplot(samples_mapped_with_counts, aes(x = MeanScore, y = Counts, color = pheno)) + geom_jitter(size = 6, alpha = 0.7, width = 0.25) + ylim(0, 3000)
    combined = ggarrange(plotlist=list(p1, p2), nrow = 1, ncol = 2)
    pdf('/project/holstegelab/Share/pacbio/data_processed/HLA_project/standard_pipeline/demultiplex_plate3/assign_samples/qc_samples.pdf', height=8, width = 16, onefile=TRUE)
    p0
    combined
    dev.off()
