library(karyoploteR)
library(regioneR)
library(GenomicRanges)
library(IRanges)
library(GenomeInfoDb)
library(rtracklayer)
library(dplyr)
library(tidyr)

# 0) force chr/start/end

as_chr_start_end <- function(df, chr_col = 1, start_col = 2, end_col = 3) {
  df <- as.data.frame(df)
  colnames(df)[chr_col] <- "chr"
  colnames(df)[start_col] <- "start"
  colnames(df)[end_col] <- "end"
  df$start <- as.numeric(df$start)
  df$end   <- as.numeric(df$end)
  df$chr   <- as.character(df$chr)
  df <- df[!is.na(df$chr) & !is.na(df$start) & !is.na(df$end), ]
  df
}

# 1) Load genome from .fai

genome_df <- read.table("PshNP85002.v4.final.fasta.fai", sep = "\t", header = FALSE,
                        stringsAsFactors = FALSE)[, 1:2]
colnames(genome_df) <- c("chr", "end")
genome_df$start <- 0
genome_df$chr <- as.character(genome_df$chr)
genome_df$end <- as.numeric(genome_df$end)
genome_df <- genome_df[, c("chr", "start", "end")]


# 2) Load genes (GFF3 -> df)

genes_gff <- rtracklayer::import("PshNP85002_hapA.v4.final.genes.gff3")
genes_df <- as.data.frame(genes_gff) %>%
  dplyr::filter(type == "gene") %>%
  dplyr::transmute(chr = as.character(seqnames),
                   start = as.numeric(start),
                   end = as.numeric(end))

genes_gff_b <- rtracklayer::import("PshNP85002_hapB.v4.final.genes.gff3")
genes_df_b <- as.data.frame(genes_gff_b) %>%
  dplyr::filter(type == "gene") %>%
  dplyr::transmute(chr = as.character(seqnames),
                   start = as.numeric(start),
                   end = as.numeric(end))

genes_df <- dplyr::bind_rows(genes_df, genes_df_b)

# 3) Load TEs 

te_a <- read.table("PshNP85002_hapA.REPET_TEannot.match.overlap_TEs_merged.bed",
                   sep = "\t", header = FALSE, stringsAsFactors = FALSE)
te_b <- read.table("PshNP85002_hapB.REPET_TEannot.match.overlap_TEs_merged.bed",
                   sep = "\t", header = FALSE, stringsAsFactors = FALSE)

# BED is usually 0-based start; convert to 1-based for GRanges
te_a <- data.frame(chr = te_a[,1], start = te_a[,2] + 1, end = te_a[,3])
te_b <- data.frame(chr = te_b[,1], start = te_b[,2] + 1, end = te_b[,3])
te_df <- dplyr::bind_rows(te_a, te_b)
te_df$chr <- as.character(te_df$chr)
te_df$start <- as.numeric(te_df$start)
te_df$end <- as.numeric(te_df$end)

# 4) Load centromeres + gaps

centro_a <- read.table("PshNP85002_hapA.centromeres_pos.filt.bed",
                       sep = "\t", header = FALSE, stringsAsFactors = FALSE)
centro_b <- read.table("PshNP85002_hapB.centromeres_pos.filt.bed",
                       sep = "\t", header = FALSE, stringsAsFactors = FALSE)

centro_a <- data.frame(chr = centro_a[,1], start = centro_a[,2] + 1, end = centro_a[,3])
centro_b <- data.frame(chr = centro_b[,1], start = centro_b[,2] + 1, end = centro_b[,3])
centro_df <- dplyr::bind_rows(centro_a, centro_b)
centro_df$chr <- as.character(centro_df$chr)
centro_df$start <- as.numeric(centro_df$start)
centro_df$end <- as.numeric(centro_df$end)

gaps_bed <- read.table("assembly_gaps (1).bed",
                       sep = "\t", header = FALSE, stringsAsFactors = FALSE)
gaps_df <- data.frame(chr = gaps_bed[,1], start = gaps_bed[,2] + 1, end = gaps_bed[,3])
gaps_df$chr <- as.character(gaps_df$chr)
gaps_df$start <- as.numeric(gaps_df$start)
gaps_df$end <- as.numeric(gaps_df$end)

# 5) Plot function 

plot_aligned_4track <- function(outfile = "PshNP85002_4track_aligned.png") {
  
  message("Creating aligned 4-track per chromosome...")
  
  # chromosome order
  chr_numbers <- 1:18
  ordered_chrs <- c()
  for (chr_num in chr_numbers) {
    chr_A <- paste0("chr", chr_num, "A")
    chr_B <- paste0("chr", chr_num, "B")
    if (chr_A %in% genome_df$chr && chr_B %in% genome_df$chr) {
      ordered_chrs <- c(ordered_chrs, chr_A, chr_B)
    }
  }
  
  existing_chrs <- intersect(ordered_chrs, genome_df$chr)
  genome_filtered_df <- genome_df[genome_df$chr %in% existing_chrs, ]
  genome_filtered_gr <- toGRanges(genome_filtered_df)
  
  # output PDF
  out_pdf <- sub("\\.png$", ".pdf", outfile)
  cairo_pdf(out_pdf, width = 8, height = 20)
  
  # plot params
  plot.params <- getDefaultPlotParams(plot.type = 1)
  plot.params$ideogramheight <- 1
  plot.params$trackheight <- 12
  plot.params$interchromosomicdist <- 0
  plot.params$topmargin <- 5
  plot.params$bottommargin <- 5
  
  kp <- plotKaryotype(
    genome = genome_filtered_gr,
    plot.type = 1,
    chromosomes = existing_chrs,
    ideogram.plotter = NULL,
    labels.plotter = NULL,
    plot.params = plot.params
  )
  
  # chromosome names: only label hap A
  final_chr_labels <- character(length(existing_chrs))
  names(final_chr_labels) <- existing_chrs
  for (chr_name in existing_chrs) {
    chr_num <- gsub("chr([0-9]+)[AB]", "\\1", chr_name)
    haplotype <- substr(chr_name, nchar(chr_name), nchar(chr_name))
    final_chr_labels[chr_name] <- ifelse(haplotype == "A", paste0("Chr", chr_num), "")
  }
  kpAddChromosomeNames(kp, cex = 0.8, font = 2, chr.names = final_chr_labels)
  
  # tracks
  track_TE_A   <- autotrack(current.track = 1, total.tracks = 4, margin = -0.01)
  track_GENE_A <- autotrack(current.track = 2, total.tracks = 4, margin = -0.01)
  track_GENE_B <- autotrack(current.track = 3, total.tracks = 4, margin = -0.01)
  track_TE_B   <- autotrack(current.track = 4, total.tracks = 4, margin = -0.01)
  
  # params
  window_size <- 10000
  gene_colors <- c("#FFE0B2", "#FF9800", "#E65100", "#BF360C")
  te_fill_color <- "#0277BD"
  
  
    # TE bp fraction per window 
  calc_te_frac <- function(chr_name, chr_length, windows, window_ends) {
    chr_tes <- te_df[te_df$chr == chr_name, ]
    if (nrow(chr_tes) == 0) return(rep(0, length(windows)))
    
    chr_tes <- chr_tes[!is.na(chr_tes$start) & !is.na(chr_tes$end), ]
    if (nrow(chr_tes) == 0) return(rep(0, length(windows)))
    
    chr_tes$start <- pmax(1, chr_tes$start)
    chr_tes$end   <- pmin(chr_length, chr_tes$end)
    chr_tes <- chr_tes[chr_tes$start <= chr_tes$end, ]
    if (nrow(chr_tes) == 0) return(rep(0, length(windows)))
    
    te_gr <- GRanges(seqnames = chr_name, ranges = IRanges(start = chr_tes$start, end = chr_tes$end))
    seqlengths(te_gr) <- setNames(chr_length, chr_name)
    
    cov01 <- pmin(coverage(te_gr)[[chr_name]], 1)
    
    vapply(seq_along(windows), function(i) {
      s <- windows[i]; e <- window_ends[i]
      as.numeric(sum(cov01[s:e])) / (e - s + 1)
    }, numeric(1))
  }
  
  # draw per chromosome
  for (chr_name in existing_chrs) {
    
    chr_length <- genome_df$end[match(chr_name, genome_df$chr)]
    if (is.na(chr_length)) next
    chr_length <- as.numeric(chr_length)
    
    windows <- seq(1, chr_length, by = window_size)
    window_ends <- pmin(windows + window_size - 1, chr_length)
    midpoints <- (windows + window_ends) / 2
    
    haplotype <- substr(chr_name, nchar(chr_name), nchar(chr_name))
    
    # TE coverage
    te_frac <- calc_te_frac(chr_name, chr_length, windows, window_ends)
    te_plot <- pmin(te_frac, cap) / cap
    
    te_y <- te_plot
    if (break_area_on_zero) te_y[te_y == 0] <- NA
    
    if (haplotype == "A") {
      kpArea(kp, chr = chr_name, x = midpoints, y = te_y,
             r0 = track_TE_A$r0, r1 = track_TE_A$r1,
             col = te_fill_color, border = NA)
    } else {
      kpArea(kp, chr = chr_name, x = midpoints, y = te_y,
             r0 = track_TE_B$r0, r1 = track_TE_B$r1,
             col = te_fill_color, border = NA)
    }
    
    # genes (counts per window, normalized per chromosome)
    chr_genes <- genes_df[genes_df$chr == chr_name, ]
    gene_counts <- vapply(seq_along(windows), function(i) {
      s <- windows[i]; e <- window_ends[i]
      sum(chr_genes$start >= s & chr_genes$start <= e)
    }, numeric(1))
    if (max(gene_counts) > 0) gene_counts <- gene_counts / max(gene_counts)
    
    gene_df_plot <- data.frame(chr = chr_name, start = windows, end = window_ends, density = gene_counts)
    gene_gr_plot <- toGRanges(gene_df_plot)
    gene_gr_plot$genedensity <- gene_df_plot$density
    
    if (haplotype == "A") {
      kpHeatmap(kp, data = gene_gr_plot, y = gene_gr_plot$genedensity,
                col = gene_colors, r0 = track_GENE_A$r0, r1 = track_GENE_A$r1)
    } else {
      kpHeatmap(kp, data = gene_gr_plot, y = gene_gr_plot$genedensity,
                col = gene_colors, r0 = track_GENE_B$r0, r1 = track_GENE_B$r1)
    }
  }
  
  # Centromeres: drawn as one block per haplotype
  # (spans TE+Gene for A; spans Gene+TE for B)
 
  centro_gr <- toGRanges(centro_df)
  
  if (length(centro_gr) > 0) {
    
    all_merged_centros <- GRanges()
    
    for (chr_name in unique(as.character(seqnames(centro_gr)))) {
      chr_centros <- centro_gr[seqnames(centro_gr) == chr_name]
      if (length(chr_centros) == 0) next
      
      if (chr_name == "chr9A") {
        merged <- reduce(chr_centros, min.gapwidth = 5000)
      } else if (chr_name %in% c("chr2B", "chr9B")) {
        merged <- reduce(chr_centros, min.gapwidth = 50000)
      } else {
        merged <- reduce(chr_centros, min.gapwidth = 200000)
      }
      
      all_merged_centros <- c(all_merged_centros, merged)
    }
    
    # Filter tiny merged segments
    w <- width(all_merged_centros)
    chr <- as.character(seqnames(all_merged_centros))
    
    keep <- logical(length(all_merged_centros))
    for (i in seq_along(all_merged_centros)) {
      if (chr[i] == "chr9A") {
        keep[i] <- w[i] >= 5000
      } else if (chr[i] %in% c("chr2B", "chr9B")) {
        keep[i] <- w[i] >= 50000
      } else {
        keep[i] <- w[i] >= 200000
      }
    }
    centro_final <- all_merged_centros[keep]
    
    # full-height inside the combined block
    cent_y0 <- 0.00
    cent_y1 <- 1.00
    
    # style
    cent_fill   <- "#00000027"
    cent_border <- "#00000066"  
    
    # these r0/r1 span TWO tracks -> one continuous band
    hapA_r0 <- track_TE_A$r0
    hapA_r1 <- track_GENE_A$r1
    hapB_r0 <- track_GENE_B$r0
    hapB_r1 <- track_TE_B$r1
    
    if (length(centro_final) > 0) {
      for (chr_name in unique(as.character(seqnames(centro_final)))) {
        chr_cent <- centro_final[seqnames(centro_final) == chr_name]
        hap <- substr(chr_name, nchar(chr_name), nchar(chr_name))
        
        if (hap == "A") {
          kpRect(kp, data = chr_cent, y0 = cent_y0, y1 = cent_y1,
                 r0 = hapA_r0, r1 = hapA_r1,
                 border = cent_border, lwd = 0.8, col = cent_fill)
        } else {
          kpRect(kp, data = chr_cent, y0 = cent_y0, y1 = cent_y1,
                 r0 = hapB_r0, r1 = hapB_r1,
                 border = cent_border, lwd = 0.8, col = cent_fill)
        }
      }
    }
  }
  
  # gaps markers (on gene tracks)
  gaps_gr <- toGRanges(gaps_df)
  if (length(gaps_gr) > 0) {
    for (chr_name in unique(as.character(seqnames(gaps_gr)))) {
      chr_gaps <- gaps_gr[seqnames(gaps_gr) == chr_name]
      gap_mid <- (start(chr_gaps) + end(chr_gaps)) / 2
      hap <- substr(chr_name, nchar(chr_name), nchar(chr_name))
      if (hap == "A") {
        kpPoints(kp, chr = chr_name, x = gap_mid, y = 0.5,
                 r0 = track_GENE_A$r0, r1 = track_GENE_A$r1,
                 pch = 21, cex = 1.0, bg = "#F44336", col = "white")
      } else {
        kpPoints(kp, chr = chr_name, x = gap_mid, y = 0.5,
                 r0 = track_GENE_B$r0, r1 = track_GENE_B$r1,
                 pch = 21, cex = 1.0, bg = "#F44336", col = "white")
      }
    }
  }
  
  # telomeres
  missing_5prime_telomeres <- c("chr2B", "chr15A", "chr15B")
  for (chr_name in existing_chrs) {
    chr_len <- genome_df$end[match(chr_name, genome_df$chr)]
    if (is.na(chr_len)) next
    small_offset <- chr_len * 0.01
    hap <- substr(chr_name, nchar(chr_name), nchar(chr_name))
    
    if (hap == "A") {
      if (!chr_name %in% missing_5prime_telomeres) {
        kpText(kp, chr = chr_name, x = 1 - small_offset, y = 0.5,
               r0 = track_GENE_A$r0, r1 = track_GENE_A$r1, labels = "►", cex = 0.8, col = "black")
      }
      kpText(kp, chr = chr_name, x = chr_len + small_offset, y = 0.5,
             r0 = track_GENE_A$r0, r1 = track_GENE_A$r1, labels = "◄", cex = 0.8, col = "black")
    } else {
      if (!chr_name %in% missing_5prime_telomeres) {
        kpText(kp, chr = chr_name, x = 1 - small_offset, y = 0.5,
               r0 = track_GENE_B$r0, r1 = track_GENE_B$r1, labels = "►", cex = 0.8, col = "black")
      }
      kpText(kp, chr = chr_name, x = chr_len + small_offset, y = 0.5,
             r0 = track_GENE_B$r0, r1 = track_GENE_B$r1, labels = "◄", cex = 0.8, col = "black")
    }
  }
  
  dev.off()
  message("Done: ", out_pdf)
}

# run
plot_aligned_4track("PshNP85002_4track_aligned.png")

