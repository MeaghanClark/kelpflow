process DEPTH_HIST {
    label 'process_medium'

    module 'r/4.4.1'
    
    input:
    path coverage_file

    output:
    path "depth_histogram.png", emit: histogram
    path "depth_cutoff.txt",    emit: cutoff

    script:
    """
    Rscript - <<'EOF'
    args <- commandArgs(trailingOnly=TRUE)

    # Read coverage data
    data <- read.table("${coverage_file}", header=FALSE, sep="\\t", col.names=c("chr", "pos", "depth"))
    depths <- data\$depth

    # Calculate median and cutoff
    median_depth <- median(depths)
    depth_cutoff <- median_depth * ${params.depth_cutoff}

    # Create PNG histogram
    png("depth_histogram.png", width=800, height=600)

    # Full range histogram
    hist(depths,
         main="Depth Distribution Across Genome",
         xlab="Depth",
         ylab="Frequency",
         breaks=100,
         col="steelblue",
         border="white")

    # Add median line
    abline(v=median_depth, col="black", lwd=2, lty=2)

    # Add cutoff line
    abline(v=depth_cutoff, col="red", lwd=2, lty=1)

    # Add legend
    legend("topright",
           legend=c(paste("Median:", round(median_depth, 2)),
                   paste("Cutoff (median ×", ${params.depth_cutoff}, "):", round(depth_cutoff, 2))),
           col=c("black", "red"),
           lty=c(2, 1),
           lwd=2)

    dev.off()

    # Save cutoff value to file
    cat(paste("Depth Cutoff:", round(depth_cutoff, 2)), file="depth_cutoff.txt")
    EOF
    """
}
