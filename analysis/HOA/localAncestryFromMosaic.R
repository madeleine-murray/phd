#!/usr/bin/env Rscript

# Import local ancestry
require(MOSAIC) # import package
library(ggplot2)
library(dplyr)
load("MOSAIC_RESULTS6_allContinents/localanc_Mayan_3way_1-18_1-22_3884_60_0.99_100.RData") # import local ancestry results
fpath="/st1/hdd/pg/human_data_processing/data_merges/hs37d5/snlab/HN00147677_HN00147680/HOA/mosaic_input/" # directory for input (to get SNP position)
chrnos=1:22 # chromosome variable
local_pos=grid_to_pos(localanc, fpath, g.loc, chrnos) # convert gridpoints (on recombination distances) to SNP position
posGen <- read.table("HOA_07122022.pileup.bed", header = FALSE, colClasses = "numeric") # import genomic positions

# Output coordinates for ancestry for each individual
find_ancestry <- function(individual, outputFile) {
	
	changes <- data.frame(Chr = numeric(0), Start = numeric(0), End = numeric(0), Ancestry = character(0))
	for (chromosome in 1:22) {

		# Combine ancestries and name rows with genomic positions from the appropriate chromosome
		ancestryData <- as.data.frame(cbind(local_pos[[chromosome]][1,individual,], local_pos[[chromosome]][2,individual,], local_pos[[chromosome]][3,individual,]))
		rownames(ancestryData) <- posGen$V2[posGen$V1 == chromosome]
	
		# Find the max ancestry for each position
		classified <- apply(ancestryData, 1, function(row) {
		max_col <- which.max(row[1:3])
		if (max_col == 1) { "American"
		} else if (max_col == 2) { "European"
		} else { "African"}
		})
		ancestryData$Ancestry <- classified

		current_ancestry <- ancestryData$Ancestry[1]
		current_rowname <- rownames(ancestryData)[1]

		# Iterate though Ancestry and record rowname (genomic position) whenever the Ancestry changes
		for (i in 2:nrow(ancestryData)) {

			# Record changes in ancestry
			if (ancestryData$Ancestry[i] != current_ancestry) {
				changes <- rbind(changes, data.frame(Chr = chromosome, Start = current_rowname, End = rownames(ancestryData)[i-1], Ancestry = current_ancestry))
				current_ancestry <- ancestryData$Ancestry[i]
				current_rowname <- rownames(ancestryData)[i]
			}

			# Record ancestry even when there's no change
			if (i == nrow(ancestryData)) {
				changes <- rbind(changes, data.frame(Chr = chromosome, Start = current_rowname, End = rownames(ancestryData)[i], Ancestry = current_ancestry))
			}
		}
	}
  
	# Create output and save in file: 4 columns (Chr, Start and End positions, ancestry affiliation)
	changes <- rbind(changes, data.frame(Chr = chromosome, Start = current_rowname, End = rownames(ancestryData)[nrow(ancestryData)], Ancestry = current_ancestry))
	write.table(changes, file = outputFile, sep = "\t", row.names = FALSE, quote = FALSE)
}

# Do it for all 18 individuals (but haploids so go to 36 eg: ind1 is hap 1&2, ind2 is hap3&4...)
for (hap in 1:36) {
  output_file <- paste("ancestry_hap-", hap, ".tsv", sep = "")
  find_ancestry(hap, output_file)
}

# Plot a boxplot summarising ancestry proportions for each individual for each chromosome
allData <- list()
for (hap in 1:36) {
	
	# Read data for the current individual
	ancestryData <- read.table(paste("ancestry_hap-", hap, ".tsv", sep=""), head=TRUE)
	
	# Calculate the cumulate segment length (eg: if American appears several times, add the segments length)
	chromosomeCumulative <- ancestryData %>%
		group_by(Chr, Ancestry) %>%
		summarize(segment_length = sum(End - Start))

	# Calculate proportion for each ancestry on each chromosome
	chromosomeProportions <- chromosomeCumulative %>%
		group_by(Chr) %>%
		mutate(chromosome_length = sum(segment_length),
		Proportion = segment_length / chromosome_length)
	
	# Highlight SGDP individuals (others are HOA only)
	chromosomeProportions$SGDP <- ifelse(hap %in% c(3,4,7,8), "yes", "no")

	# Store the processed data in the list
	allData[[hap]] <- chromosomeProportions
}
plotData <- bind_rows(allData) # combine data for all individuals
plotData$Ancestry <- factor(plotData$Ancestry, levels = c("American","European","African")) # reorder the ancestries
plotChr <- ggplot(data = plotData) +
	geom_boxplot(aes(x = factor(Chr), y = Proportion, fill = Ancestry),
		lwd = 5, outlier.shape = 16, outlier.size = 8) + # make boxplots bigger
	geom_point(data = plotData %>% filter(SGDP == 'yes'),
		aes(x = factor(Chr), y = Proportion),
		shape = 24, size = 18, fill = 'gold', stroke = 5) + # highlight SGDP (others are HOA only)
	facet_wrap(~ Ancestry, ncol = 1) + # one plot of three rows for each ancestry
	labs(title = "Ancestry proportions along chromosomes",
		x = "", y = "Proportion (as haploid)") +
	scale_fill_manual(values=c("American"="red","European"="deepskyblue","African"="green")) + # set ancestry colours
	guides(fill = 'none') + # remove the legend
	theme_linedraw() +
	theme(plot.title = element_text(size = 90, face = 'bold'), strip.text = element_text(size = 70, face = "bold"),
		axis.title = element_text(size = 80), axis.text.y = element_text(size = 50), axis.text.x = element_text(size = 60))

# Plot a boxplot summarising ancestry proportions for each individual genome-wide
allData <- list()
for (hap in 1:36) {
	
	# Read data for the current individual
	ancestryData <- read.table(paste("ancestry_hap-", hap, ".tsv", sep=""), head=TRUE)

	# Calculate the sum for each Ancestry genome-wide
	ancestryLength <- ancestryData %>%
		group_by(Ancestry) %>%
		summarize(segmentLength = sum(End - Start))

	# Calculate the length of the genome
	genomeLength <- sum(ancestryData$End - ancestryData$Start)

	# Calculate the proportion for each Ancestry compared the the whole-genome
	genomeProportions <- ancestryLength %>%
		mutate(Proportion = segmentLength / genomeLength)

	# Highlight SGDP individuals (others are HOA only)
	genomeProportions$SGDP <- ifelse(hap %in% c(3,4,7,8), "yes", "no")

	# Store the processed data in the list
	allData[[hap]] <- genomeProportions
}
plotData <- bind_rows(allData) # combine data for all individuals
plotData$Ancestry <- factor(plotData$Ancestry, levels = c("American","European","African")) # reorder the ancestries
plotGw <- ggplot(data = plotData) +
	geom_boxplot(aes(x = Ancestry, y = Proportion, fill = Ancestry),
		lwd = 5, outlier.shape = 16, outlier.size = 8) + # make boxplots bigger
	geom_point(data = plotData %>% filter(SGDP == 'yes'),
		aes(x = Ancestry, y = Proportion),
		shape = 24, size = 18, fill = 'gold', stroke = 5) + # highlight SGDP (others are HOA only)
	labs(title = "Ancestry proportions genome-wide",
		x = "", y = "Proportion (as haploid)") +
	scale_fill_manual(values = c("American" = "red", "European" = "deepskyblue", "African" = "green")) +
	guides(fill = 'none') +
	theme_linedraw() +
	theme(plot.title = element_text(size = 90, face = 'bold'),
		axis.title = element_text(size = 80), axis.text.y = element_text(size = 50), axis.text.x = element_text(size = 70))
#print(as_tibble(plotData), n = 108) 
library(gridExtra)
png(filename = 'AncestryProportions_chr.png', width=3840, height=3840)
plotChr
dev.off()
png(filename = 'AncestryProportions_gw.png', width=2160, height=3840)
plotGw
dev.off()
png(filename = 'AncestryProportions.png', width=6000, height=3840)
grid.arrange(plotGw, plotChr, ncol = 2, widths = c(1, 2))
dev.off()

# Output American coordinates for each diploid (putting together both haploids)
find_nonAmerican <- function(individual, outputFile) {

	changes <- data.frame(Chr = numeric(0), Start = numeric(0), End = numeric(0), American = character(0))
	for (chromosome in 1:22) {
		
		# Combine ancestries and name rows with genomic positions from the appropriate chromosome
		diploidData <- as.data.frame(cbind(local_pos[[chromosome]][1,individual,], local_pos[[chromosome]][2,individual,], local_pos[[chromosome]][3,individual,],
						local_pos[[chromosome]][1,individual+1,], local_pos[[chromosome]][2,individual+1,], local_pos[[chromosome]][3,individual+1,]))
		rownames(diploidData) <- posGen$V2[posGen$V1 == chromosome]

		# If the max ancestry is American in both chromosomes, then assign American, else nonAmerican
		classified <- apply(diploidData, 1, function(row) {
			if (which.max(row[1:3]) == 1 && which.max(row[4:6]) == 1) { "yes" } else { "no" }
		})
		diploidData$American <- classified
		
		current_ancestry <- diploidData$American[1]
		current_rowname <- rownames(diploidData)[1]

		# Iterate through American and record rename (genomic position) whenever the non-American ancestry appears
		for (i in 2:nrow(diploidData)) {
			
			# Record changes in ancestry
			if (diploidData$American[i] != current_ancestry) {
				changes <- rbind(changes, data.frame(Chr = chromosome, Start = current_rowname, End = rownames(diploidData)[i-1], American = current_ancestry))
				current_ancestry <- diploidData$American[i]
				current_rowname <- rownames(diploidData)[i]
			}
			
			# Record ancestry even when there's no change
			if (i == nrow(diploidData)) {
				changes <- rbind(changes, data.frame(Chr = chromosome, Start = current_rowname, End = rownames(diploidData)[i], American = current_ancestry))
			}
		}
	}

	# Create output and save in file: 4 columns (Chr, Start and End positions, "yes" or "no" for fully American) and exclude the last because gets repeated for some reason
	changes <- rbind(changes, data.frame(Chr = chromosome, Start = current_rowname, End = rownames(diploidData)[nrow(diploidData)], American = current_ancestry))
	write.table(changes[-nrow(changes),], file = outputFile, sep = "\t", row.names = FALSE, quote = FALSE)
}

# Do it for all 18 individuals (as diploids)
for (hap in seq(1,36,by=2)) {
	output_file <- paste0("masking/american_ind-", ceiling(hap/2), ".tsv")
	find_nonAmerican(hap, output_file)
}

q()

# Plot ancestry for each individual (karyograms) --> sanity check
ancestry <- function(individual, chromosome) {
	pos <- posGen$V2[posGen$V1 == chromosome] # get the appropriate positions for the chromomome
	american <- as.data.frame(cbind(pos, local_pos[[chromosome]][2,individual,]))
	european <- as.data.frame(cbind(pos, local_pos[[chromosome]][3,individual,]))
	african  <- as.data.frame(cbind(pos, local_pos[[chromosome]][1,individual,]))
	ggplot() +
		geom_line(data=american, linewidth=6, aes(x=pos, y=V2, colour="American")) +
		geom_line(data=european, linewidth=6, aes(x=pos, y=V2, colour="European")) +
		geom_line(data=african,  linewidth=6, aes(x=pos, y=V2, colour="African")) +
		ggtitle(paste("Ancestry for haploid", individual, "on chromosome", chromosome, sep=" ")) +
		labs(x="Genomic position", y="Ancestry proportion") +
		scale_colour_manual(name = "Ancestry",
			values = c("American" = "red", "European" = "blue", "African" = "green")) +
		theme_classic() +
		theme(axis.title = element_text(size = 65), axis.text = element_text(size = 50),
			plot.title = element_text(size = 80, face = "bold"),
			legend.title = element_text(size = 65), legend.text = element_text(size = 50), legend.key.size = unit(5, "lines"))

}
png(filename = 'karyogram_hap8_chr16.png', width=3840, height=2160)
ancestry(8, 16) # arguments: individual & chromosome
dev.off()
