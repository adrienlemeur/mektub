{
  rm(list=ls())
  graphics.off()
  gc()

  for(lib in c("seqinr", "jackalope", "optparse", "Biostrings")){
    suppressPackageStartupMessages(library(lib, character.only = T))
  }
  cat(paste("Packages loaded...", "\n"))
} #start

{
  args = commandArgs(trailingOnly=TRUE)
  option_list = list(
    make_option(c("-r", "--reference"), type="character", default=NULL, 
      help="reference sequence to mutate (fasta)", metavar="character"),
    make_option(c("-t", "--transposon"), type="character", default=NULL, 
                help="transposon positions in the reference sequence (bed)", metavar="character"),
    make_option(c("-c", "--haplotypecount"), type='numeric', default=1,
      help="number of haplotype to compute (>=2)"),
    make_option(c("-n", "--name"), type='character', default=2,
                help="artificial genome name prefix")
    );
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
  opt$haplotypecount = opt$haplotypecount+1

  if (is.null(opt$reference)){
    print_help(opt_parser)
    stop("A reference sequence must be supplied", call.=FALSE)
  } else if (is.null(opt$transposon)){
    print_help(opt_parser)
    stop("A bed file with transposons positions must be provided", call. = FALSE)
  } else if (is.null(opt$name)){
    print_help(opt_parser)
    stop("An output name for artificial genome must be provided", call. = FALSE)
  }
} #input parameters

{
  reference_sequence <- list()
  #reference_sequence$sequence <- read.fasta(file = "BLrg.fasta", set.attributes = F)[[1]]
  reference_sequence$sequence <- read.fasta(file = opt$reference, set.attributes = F)[[1]]
  reference_sequence$position = c(1:length(reference_sequence$sequence))

  neo_sequence <- list()
  neo_sequence$sequence <- reference_sequence$reference ; neo_sequence$position = reference_sequence$position

  #transposons_loci <- apply(read.table("BLrg_IS6110.bed", header = F, sep = "\t", col.names = c("chrom", "start", "stop")), 1, function(x) as.list(sapply(x, trimws) ))
  transposons_loci <- apply(read.table(opt$transposon, header = F, sep = "\t", col.names = c("chrom", "start", "stop")), 1, function(x) as.list(sapply(x, trimws) ))

  for(i in 1:length(transposons_loci)){
    transposons_loci[[i]]$start <- as.numeric(transposons_loci[[i]]$start)+1
    transposons_loci[[i]]$stop <- as.numeric(transposons_loci[[i]]$stop)
  }
} #initialisation

{
  seqCut <- function(sequence, bed){
  output = list()
  output$sequence <- sequence$sequence[-c(bed$start:bed$stop)]
  output$position <- sequence$position[-c(bed$start:bed$stop)]
  return(output)
}
  seqPaste <- function(reference_sequence, neo_sequence, bed){
  output = list()
  sequence <- reference_sequence$sequence[bed$start:bed$stop]
  position <- reference_sequence$position[bed$start:bed$stop]

  insert_pos <- sample(1:length(neo_sequence$sequence), 1)
  #new_length = length(neo_sequence$sequence)
  #a = length(neo_sequence$sequence[1:insert_pos])
  #b = length(neo_sequence$sequence[insert_pos+1:new_length])
  #print(length(reference_sequence$sequence) - a)
  #print(length(reference_sequence$sequence) - b)
  #print(length(reference_sequence$sequence) - (a + b))
  output$sequence = c(neo_sequence$sequence[1:insert_pos-1], sequence, neo_sequence$sequence[insert_pos:length(neo_sequence$sequence)])
  output$position = c(neo_sequence$position[1:insert_pos-1], position, neo_sequence$position[insert_pos:length(neo_sequence$position)])
  return(output)
}

  neo_sequence = list() ; neo_sequence$sequence <- reference_sequence$sequence ; neo_sequence$position = reference_sequence$position

  for(i in 1:length(transposons_loci)){
  neo_sequence <- seqCut(neo_sequence, transposons_loci[[i]])
}

  for(i in 1:length(transposons_loci)){
  neo_sequence <- seqPaste(reference_sequence, neo_sequence, transposons_loci[[i]])
}

  random_deletion <- floor(abs(rnorm(3, 900, 2150)))
  deletion_position <- sample(1:length(neo_sequence$sequence), 3, replace = F)

  for(i in 1:length(deletion_position)){
    neo_sequence$sequence <- neo_sequence$sequence[-c(deletion_position[i]:(deletion_position[i]+random_deletion[i]))]
    neo_sequence$position <- neo_sequence$position[-c(deletion_position[i]:(deletion_position[i]+random_deletion[i]))]
  }
  equivalence = as.data.frame(cbind(1:length(neo_sequence$position), neo_sequence$position))
  equivalence = equivalence[equivalence$V1 != equivalence$V2,]

  cat(paste("Structural variants added, writing equivalence table...", "\n"))
  write.table(equivalence, file = paste(opt$name, "_position_equivalence.txt", sep = ""), quote = F, sep = "\t", row.names = F, col.names = F)
} # structural variation

#deletion_size <- unlist(read.csv("DEL_SIZE.txt", F))
#library(EnvStats)
#epois(deletion_size, ci = TRUE, conf.level = 0.9)
#enorm(deletion_size, ci = TRUE, conf.level = 0.9)

#trans_colors = c(rgb(1,0,0,0.5), rgb(1,1,0,0.5))
#hist(deletion_size, breaks= 50, freq = F, col = "firebrick4")
#hist(abs(rnorm(1000, 900, 2150)), freq = F, add = T, breaks = 25, col = trans_colors[2])

write.fasta(neo_sequence$sequence, file.out = "tmp.fa", names = opt$name)
ref <- read_fasta("tmp.fa", cut_names = TRUE)

unlink("tmp.fa", recursive = FALSE, force = FALSE)

{
  #GTR = list()
  #GTR$nuc$A = 0.172 ; GTR$nuc$T = 0.172 ; GTR$nuc$C = 0.327 ; GTR$nuc$G = 0.329
  #GTR$mut$A_C = 1.01590 ; GTR$mut$A_G = 3.56960 ; GTR$mut$A_T = 0.43746 ; GTR$mut$C_G = 0.91385 ; GTR$mut$C_T = 3.56960 ; GTR$mut$G_T = 1.00000 ; GTR$invariants = 0.957
  #transver <- c(GTR$mut$C_T, GTR$mut$A_T, GTR$mut$G_T, GTR$mut$A_C, GTR$mut$C_G, GTR$mut$A_G)
  #tuberculosis_alignement <- readDNAStringSet("MSA_25_strains.fasta", use.names = T)
  #sub_matrix <- oligonucleotideTransitions(tuberculosis_alignement, left=1, right=1, as.prob=FALSE)[c(4:1), c(4:1)]/105876744
  #sub_GTR_parameters <- c(sub_matrix[1,2], sub_matrix[1,3], sub_matrix[1,4], sub_matrix[2,3], sub_matrix[2,4], sub_matrix[3,4])
  nuc <- c(0.172, 0.327, 0.172, 0.329)
  sub_GTR_parameters <- c(0.06286848, 0.06086589, 0.01691312, 0.11487831, 0.06114962, 0.06270075)
  GTR_tuberculosis <- sub_GTR(pi_tcag = nuc, abcdef = sub_GTR_parameters, invariant = 0.957)
  haplotypes <- create_haplotypes(ref, haps_theta(theta = 5*10^3*2*1.23*10^-7, n_haps = opt$haplotypecount), sub = GTR_tuberculosis, ins = indels(2.685*10^-9), del = indels(2.685*10^-9), n_threads = 10)
  haplotypes$set_names(paste("H", 1:opt$haplotypecount, sep = ""))
  #Θ=2*Ne*μ
  haplotypes$rm_haps(paste("H", opt$haplotypecount, sep = ""))

  cat(paste("Sequence mutated, writing fasta and vcf...", "\n"))
  write_fasta(haplotypes, opt$name, compress = 1, n_threads = 10, overwrite = T, text_width = 60)
  write_vcf(haplotypes, opt$name, compress = 1, overwrite = T)

} #model defintion
cat(paste("All done :)", "\n"))