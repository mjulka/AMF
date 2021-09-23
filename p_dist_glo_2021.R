#2021 September
#p distances calculations and diagrams 
getwd()
setwd("D:/Yulia/AMF/2021/september")
library("ape")
library('phangorn')

# function for making histogram for fasta alignment 


p.histogram <- function(fasta) {
  alignment<-read.FASTA(file=fasta, type = "DNA")
  p_distances <-dist.dna(alignment, model="raw", pairwise.deletion = TRUE)
  taxa_number <- length(alignment)
  #hist_breaks <- taxa_number/2
  #picture_name <- paste (fasta, "_", taxa_number, ".jpg", sep = "")
  picture_name <- paste (fasta, "_", taxa_number, "_140", ".jpg", sep = "")
  print (picture_name)
  jpeg(picture_name,width = 540, height = 545, quality = 95)
  #hist(p_distances, hist_breaks, main= fasta, xlab= "p-distances")
  #hist(p_distances, taxa_number, main= fasta, xlab= "p-distances") #for small taxa number
  hist(p_distances, 140, main= fasta, xlab= "p-distances")
  dev.off()
}
#data from Petr
p.histogram ("ref_ITS1_alined.1.fas")
p.histogram ("ref_ITS2_aligned.1.fas")

#by genera

p.distances.genus <- function(fasta) {
  alignment<-read.FASTA(file=fasta, type = "DNA")
  p_distances <-dist.dna(alignment, model="raw", pairwise.deletion = TRUE)
  average<- mean(p_distances)
  med<- median(p_distances)
  minimum <- min(p_distances)
  maximum<-max(p_distances)
  taxa_number <- length(alignment)
  file_name <- paste (fasta, "_p-distance_genus_", taxa_number, ".txt", sep = "")
  sink(file_name)
  cat("p-distances for", fasta, taxa_number, "samples")
  cat("\n")
  cat('average', average)
  cat("\n")
  cat('mediana', med)
  cat("\n")
  cat('minimum', minimum)
  cat("\n")
  cat('maximum', maximum)
  sink()
}
#ITS2
#input files - aligned by muscle (via Aliview) fastas, each species represented only by 1 seq
#only genera with >1 species analysed
setwd("D:/Yulia/AMF/2021/september/ITS2_Genera")
getwd()
p.distances.genus("Acaulospora_ITS2.fa")
p.distances.genus("Ambispora_ITS2.fa")
p.distances.genus("Archaeospora_ITS2.fa")
p.distances.genus("Cetraspora_ITS2.fa")
p.distances.genus("Claroideoglomus_ITS2.fa")
p.distances.genus("Corymbiglomus_ITS2.fa")
p.distances.genus("Dentiscutata_ITS2.fa")
p.distances.genus("Diversispora_ITS2.fa")
p.distances.genus("Dominikia_ITS2.fa")
p.distances.genus("Entrophospora_ITS2.fa")
p.distances.genus("Funneliformis_ITS2.fa")
p.distances.genus("Gigaspora_ITS2.fa")
p.distances.genus("Glomus_ITS2.fa")
p.distances.genus("Kamienskia_ITS2.fa")
p.distances.genus("Paraglomus_ITS2.fa")
p.distances.genus("Racocetra_ITS2.fa")
p.distances.genus("Redeckera_ITS2.fa")
p.distances.genus("Rhizoglomus_ITS2.fa")
p.distances.genus("Rhizophagus_ITS2.fa")
p.distances.genus("Sacculospora_ITS2.fa")
p.distances.genus("Scutellospora_ITS2.fa")
p.distances.genus("Septoglomus_ITS2.fa")
#input file with samples genera by 1 species (type where it possible)
p.distances.genus("between_genera_ITS2.fa")
#ITS1
#input files - aligned by muscle (via Aliview) fastas, each species represented only by 1 seq
#only genera with >1 species analysed
setwd("D:/Yulia/AMF/2021/september/ITS1_Genera")
getwd()
p.distances.genus("Acaulospora_ITS1.fa")
p.distances.genus("Ambispora_ITS1.fa")
p.distances.genus("Archaeospora_ITS1.fa")
p.distances.genus("Cetraspora_ITS1.fa")
p.distances.genus("Claroideoglomus_ITS1.fa")
p.distances.genus("Corymbiglomus_ITS1.fa")
p.distances.genus("Dentiscutata_ITS1.fa")
p.distances.genus("Diversispora_ITS1.fa")
p.distances.genus("Dominikia_ITS1.fa")
p.distances.genus("Entrophospora_ITS1.fa")
p.distances.genus("Funneliformis_ITS1.fa")
p.distances.genus("Gigaspora_ITS1.fa")
p.distances.genus("Glomus_ITS1.fa")
p.distances.genus("Kamienskia_ITS1.fa")
p.distances.genus("Paraglomus_ITS1.fa")
p.distances.genus("Racocetra_ITS1.fa")
p.distances.genus("Redeckera_ITS1.fa")
p.distances.genus("Rhizoglomus_ITS1.fa")
p.distances.genus("Rhizophagus_ITS1.fa")
p.distances.genus("Sacculospora_ITS1.fa")
p.distances.genus("Scutellospora_ITS1.fa")
p.distances.genus("Septoglomus_ITS1.fa")
#input file with samples genera by 1 species (type where it possible)
p.distances.genus("between_genera_ITS1.fa")
