process SPLIT_GENOME {
    label 'process_low'

    module 'miniconda3'
    conda 'conda-forge::r-base'

    input:
    path fai 

    output:
    path "segment_*.txt", emit: segments

    script:
    """
    cat <<'RSCRIPT' > split_genome.R
    split_genome <- function(genome_df, num_segments) {
      total_genome_length <- sum(genome_df\$chrom_len)
      target_length <- total_genome_length / num_segments
      segments <- list()
      current_segment <- c()
      current_length <- 0
      seg_index <- 1

      for (j in 1:nrow(genome_df)) {
        chrom_name <- genome_df\$chrom_name[j]
        chrom_len  <- genome_df\$chrom_len[j]
        pos <- 1

        while (pos <= chrom_len) {
          remaining_scaffold <- chrom_len - pos + 1
          remaining_target   <- target_length - current_length

          if (remaining_scaffold <= remaining_target) {
            current_segment <- c(current_segment, paste0(chrom_name, ":", pos, "-", chrom_len))
            current_length  <- current_length + remaining_scaffold
            pos <- chrom_len + 1

            # Always close the segment at the end of a scaffold, even if it
            # hasn't reached target_length. This keeps segments from spanning
            # multiple scaffolds, since downstream HMM position calculations
            # (dupHMM.R) assume positions increase monotonically within a
            # segment and break at scaffold boundaries otherwise.
            segments[[seg_index]] <- data.frame(region = current_segment)
            seg_index       <- seg_index + 1
            current_segment <- c()
            current_length  <- 0
          } else {
            chunk_end <- pos + remaining_target - 1
            current_segment <- c(current_segment, paste0(chrom_name, ":", pos, "-", as.integer(chunk_end)))
            segments[[seg_index]] <- data.frame(region = current_segment)
            seg_index       <- seg_index + 1
            current_segment <- c()
            current_length  <- 0
            pos <- as.integer(chunk_end) + 1
          }
        }
      }

      if (length(current_segment) > 0) {
        segments[[seg_index]] <- data.frame(region = current_segment)
      }

      return(segments)
    }

    fai <- read.csv("${fai}", sep = "\\t", header = FALSE)
    colnames(fai) <- c("chrom_name", "chrom_len", "offset", "line_bases", "line_width")

    data <- fai[, c("chrom_name", "chrom_len")]
    data <- data[-grep("ctg", data\$chrom_name), ]

    no_chunks <- ${params.num_genome_segments}

    for (i in 1:nrow(data)) {
      data\$order[i] <- sum(data\$chrom_len[1:i])
    }
    data\$chrom_start <- c(1, data\$order + 1)[1:nrow(data)]
    data <- data[, c("chrom_name", "chrom_len", "chrom_start", "order")]

    genome_segments <- split_genome(data, no_chunks)

    for (i in seq_along(genome_segments)) {
      write.table(genome_segments[[i]],
                  file = paste0("segment_", i, ".txt"),
                  row.names = FALSE, sep = "\\t",
                  col.names = FALSE, quote = FALSE,
                  fileEncoding = "UTF-8")
    }
RSCRIPT

    Rscript split_genome.R
    """
}