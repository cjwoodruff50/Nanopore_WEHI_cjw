# Rscript to access Sereika M etal ENA data and process it to extract nanopore
# reads covering the 16S and (modified) 23S subunits, as well as the Z2 and Z5 
# sub-regions of the 23S and the rrn of each bacterial genome.
#
# Revision built on version v02b  of 3Feb20224.
#
# source("/stornext/Bioinf/data/lab_speed/cjw/microbiome/scripts/Rscripts/extract_Sereika_16S23S_v02b.R")

args = commandArgs(trailingOnly=TRUE)

#######################################################################################
##########################                               ##############################
##########################         INITIALISATIONS       ##############################
##########################                               ##############################
#######################################################################################
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
packages <- c("stringr","seqinr","tictoc","latex2exp","ShortRead", "DECIPHER","bfsl","dplyr",
              "nnet", "MASS","stringdist","phyloseq","vegan","ape","collapse")    # Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  BiocManager::install(packages[!installed_packages])
  #  BiocManager::install(packages[!installed_packages],repos="https://cran.ms.unimelb.edu.au/")
}# Packages loading
invisible(lapply(packages, library, character.only = TRUE))


# Alignment parameters
mat=nucleotideSubstitutionMatrix(match = 1, mismatch = -2, baseOnly = TRUE)
substitute_mat1 = mat
substitute_mat2 = mat
align_params = list(align_type ="global", gapOpen = 2, gapExtend = 2)
align_parameters1 = list(align_type="global-local", gapOpen=2, gapExtend=2)
align_parameters2 = list(align_type="overlap", gapOpen=1, gapExtend=0.5)
gapOpen = 2
gapExtend = 0

out_dateString = "03022024"
numcores = 56
#######################################################################################
##########################                               ##############################
##########################           FUNCTIONS           ##############################
##########################                               ##############################
#######################################################################################
source("/stornext/Bioinf/data/lab_speed/cjw/microbiome/scripts/Rscripts/rrn_operon_in_silico_16S_ITS_23S_primers_04Aug2021.R")
nprimers = length(primers)
# Note that primer_sets 1 to 4 give the 16S 5', 16S 3', 23S 5' and 23S 3' primers, respectively.
source("/stornext/Bioinf/data/lab_speed/cjw/microbiome/scripts/Rscripts/rrn_16S_ITS_23S_boundaries.R")

primerAlign = function(i,frags){
  # Alignment of a single primer, identified by the list member i, against a list
  # of fragment sequences frags[jf] where jf is in 1:length(frags)
  # Returns 3 alignment MOPs (alscore, PID, Levdist) and start_loc of the primer
  # alignment for each fragment sequence.
  # This function can be easily called by mclapply.
  primer = i
  skip = FALSE
  # First check to see if this is a primer set associated with a boundary that has already been determined.
  # This can be done by parsing the primer type of one entry in this primer_set, together with keeping a
  # running tally of which Z region 3' boundaries have already been determined.
  # 
  do.check = FALSE
  if (do.check){
    zstr = s2c(primer$type)
    if (zstr[1]=="Z"){
      this_bdry = as.integer(zstr[2])
      if ((zstr[4]=="3") & (zstr[5]=="p")){ # So this is genuinely a 3' boundary
        zprime3p[this_bdry] = 1
      } else if (this_bdry>1){ 
        if (zprime3p[this_bdry]==1){ 
          # So this is a 5' boundary and the immediately prior 3' boundary has been computed. 
          # Just use the previously computed AA01 matrices for the set of primers in this primer_set.
          skip = TRUE
        }
      }
    }    #   end    conditional    block  for avoiding calculation of Zregion 5' boundaries except for Z1_5'  (zstr[1]=="Z")
  }      #   end    conditional    block  checking to avoid redundant computation of boundary alignments.
  
  if (!skip){
    A = matrix(-1,nrow=length(frags),ncol=4)
    colnames(A) = c("score","PID","Levdist","startpos")
    for (jf in 1:length(frags)){
      fragseq = frags[jf]    # frags_reord[jf]
      algnf = pairwiseAlignment(primer$sequence,fragseq, type="global-local",substitutionMatrix=mat,
                                gapOpening=align_params$gapOpen,gapExtension=align_params$gapExtend)
      algnr = pairwiseAlignment(reverseComplement(DNAString(primer$sequence)),fragseq,
                                type="global-local",substitutionMatrix=mat,
                                gapOpening=align_params$gapOpen,gapExtension=align_params$gapExtend)
      if(score(algnf)>score(algnr)){algn=algnf} else {algn=algnr}
      Levdist = stringDist(c(as.character(pattern(algn))[1],as.character(subject(algn))))
      A[jf,] = c(score(algn), pid(algn), Levdist, start(subject(algn)) )
    }
  }
  out = A
}


select_boundarySet_primers = function(boundarySet){
  # Returns the set of primers for each of the rrn boundaries indexed by boundarySet.
  # Output is a list of integer vectors.
  # 9 March 2023                          [cjw]
  rrnsubuPrimers = primers
  sucount = 0
  for (jb in boundarySet){
    np = length(primer_sets[[jb]])
    for (kp in 1:np){
      sucount = sucount + 1
      rrnsubuPrimers[[sucount]] = primers[[primer_sets[[jb]][kp]]]
    }
  }
  nsu = sucount
  out = rrnsubuPrimers[1:nsu]
}

read_filter_frags = function(inname,fastqpath,thresh,numcores,makeplot=FALSE){
  # Accesses the read fastq sequences and quality scores. Length filter these
  # reads and returns the full set of fragments and an index to those that 
  # meet the filter low -threshold cutoff length.
  # Can be used to give a histogram of lengths of reads that are above the 
  # lower threshold, thresh, and shorter than 20001 bases.
  # 9 March 2023                          [cjw]
  Fq = readFastq(file.path(fastqpath,inname))
  Q = quality(Fq)
  Nr = length(Q)
  Fqlen = function(i){width(Fq[i])}
  rL = mclapply(1:Nr,Fqlen,mc.cores = numcores)
  readLen = unlist(rL)
  ilen = which(readLen>thresh)
  Nrg = length(ilen)
  readLenfilt = readLen[ilen]
  if (makeplot){
    hist(readLenfilt[which(readLenfilt<20001)],main=paste(inname," minlength ",thresh,sep=" "),xlab="Read Length")
  }
  frags = sapply(1:length(ilen),function(j){unlist(sread(Fq[ilen[j]]))}) 
  out = list(ifilt = ilen, frags = frags, Fq = Fq)
}

chrName = function(shead){
  p1 = str_locate(shead,"SN:")
  p2 = str_locate(shead,"contig_")
  out = substr(shead, p1[2]+1,p2[2]+1)
}

chrLen = function(shead){
  p1 = str_locate(shead,"LN:")
  out = as.integer(substr(shead, p1[2]+1,nchar(shead)))
}

reorient16S = function(fr16, refOrient16S){
  # Identifies which of the sequences in fr16 have the opposite orientation
  # to refOrient16S, then reverseComplements each of these.
  # Returns the reoriented sequences, together with the identification of 
  # original strandedness (needed to match teh quality data for each sequence.)
  y1 = DNAStringSet(c(refOrient16S,fr16))
  which_ref = 1
  zzor = OrientNucleotides(y1,reference=which_ref,type="orientations")
  refseqor = OrientNucleotides(y1,reference=which_ref,type="sequences")
  zzor = zzor[2:length(zzor)]
  refseqor = refseqor[2:length(refseqor)]
  out = list(seqs=refseqor, strand = zzor)
}

reorient23S = function(fr23, refOrient23S){
  # Identifies which of the sequences in fr23 have the opposite orientation
  # to refOrient23S, then reverseComplements each of these.
  # Returns the reoriented sequences, together with the identification of 
  # original strandedness (needed to match teh quality data for each sequence.)
  y1 = DNAStringSet(c(refOrient23S,fr23))
  which_ref = 1
  zzor = OrientNucleotides(y1,reference=which_ref,type="orientations")
  refseqor = OrientNucleotides(y1,reference=which_ref,type="sequences")
  zzor = zzor[2:length(zzor)]
  refseqor = refseqor[2:length(refseqor)]
  out = list(seqs=refseqor, strand = zzor)
}

reorientZ = function(frZ, refOrientZ){
  # Identifies which of the sequences in frZ have the opposite orientation
  # to refOrientZ, then reverseComplements each of these.
  # Returns the reoriented sequences, together with the identification of 
  # original strandedness (needed to match teh quality data for each sequence.)
  y1 = DNAStringSet(c(refOrientZ,frZ))
  which_ref = 1
  zzor = OrientNucleotides(y1,reference=which_ref,type="orientations")
  refseqor = OrientNucleotides(y1,reference=which_ref,type="sequences")
  zzor = zzor[2:length(zzor)]
  refseqor = refseqor[2:length(refseqor)]
  out = list(seqs=refseqor, strand = zzor)
}
#######################################################################################
############################     RUN INITIALISATION      ##############################
#######################################################################################
#basepath = "/vast/projects/rrn/microbiome/Sereika"
#fastqpath = file.path(basepath,"fastq");  inpath = fastqpath
#basepath = "/vast/projects/rrn/microbiome/paper"
#basepath = "/vast/scratch/users/woodruff.c"
#basepath = "/vast/scratch/users/woodruff.c/output"
basepath = "/vast/scratch/users/woodruff.c/Zymo_hmw_r104"
basepath = "/vast/projects/rrn/microbiome/papercheck"
fastqpath = file.path(basepath,"insplitfastq") 
inpath = fastqpath
outpath = file.path(basepath,"output")
datasetID = "D6322B"    
# Assuming we are working only with the R10.4 MAG Zymo D6322 reference genomes the 
# reference genome path (not used, 30 March 2023) is as below:-
#reffastapath = file.path(basepath,"import","ZYMO_assemblies/r104")
refOrient16S= "AGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAACAGCTTGCTGTTTCGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTA"
refOrient23S_full= "CCGAACTGTCTCACGACGTTCTGAACCCAGCTCGCGTGCCGCTTTAATGGGCAAACAGCCCAACCCTTGGGACCGACTACAGCCCCAGGATGCGACGAGCCGACATCGAGGTGCCAAACCTCCCCGTCGATGTGGACTCTTGGGGGAGATAAGCCTGTTATCCCCAGGGTAGCTTTTATCCGTTGAGCGATGGCCCTTCCATGCGGAACCACCGGATCACTAAGCCCGACTTTCGTCCCTGCTCGACTTGTAGGTCTCGCAGTCAAGCTCCCTTCTGCCTTTACACTCTTCGAAATGATTGAACCATTCTGAGGGAACCTTTGGGCGCCTCCGTTACCTTTTAGGAGGCGACCGCCCCAGTCAAACTGCCCATCTGACACTGTCTCCCACCACCAATTAGTGGTGCGGGTTAGAGGGTTCATAACACAAGGGTAGTATCCCACCAGCGCCTCCTCCGAAACTAGCGTTCCGGTCTCATCGGCTCCTACCTATCCTGTACATGTGGTACAAACACTCAATATCAAACTACAGTAAAGCTCCATGGTATTTCCGTCCTGTCGCGGGTAACCTGCATCTTCACAGGTACTAAAATTTCACCGAGTCTCTCGTTGAGACAGTGCCCAAATCGTTACGCCTTTCGTGCGGGTCAGGAACTTACCCGACAAGGAATTTCACTTGTTCGGACCGTTATAGTTACGGCCGCCGTTTACTGGGGCTTCAATTCTGAGCTTCGCCGAAGCTAACCCATCCTCTTAACCTTCCAGCACCGGGCAGGCGTCAGCCCCTATACTTCATCTTACGATTTTGCAGAGACCTGTGTTTTTGATAAACAGTCGCTTGGGCCTATTCACTGCGGCTGACCGAAGTCAGCACCCCTTCTCCCGAAGTTACGGGGTCATTTTTGCCGAGTTCCTTAACGAGAGTTCGCTCGCTCACCTTAGGATACTCTCCTCGACTACACTGTCGGTTTACGGTACGGGCAGTTGTTTTCTCACTAGAAGCTTTTCTTGGCAGTGTGACATCAGGAACTTCGCTACTATTATTTCGCTCCCCGTCACAACTTGTCCTTAAAGAGTAAAGCATTTAACTCTACTAAACTTGGTTTCTGCTTGGACATGCACTTCCAATCGCATGCATTCCTTAGCCTACTGCGTCCCTCCATTGCTCAAACAAAAACAACTGGTACAGGAATATCAACCTGTTGTCCATCGCCTACGCCTATCGGCCTCGGCTTAGGTCCCGACTAACCCTGGGCGGACGAGCCTTCCCCAGGAAACCTTAGTCATACGGTGGACAGGATTCTCACCTGTCTTTCGCTACTCATACCGGCATTCTCACTTCTAAGCGCTCCAGCCGTCCTCACGATCGACCTTCAACGCCCTTAGAACGCTCTCCTACCACTACACCTAATGGTGTAGTCCACAGCTTCGGTAATATGTTTAGCCCCGGTACATTTTCGGCGCAGGGTCACTCGACTAGTGAGCTATTACGCACTCTTTAAATGGTGGCTGCTTCTAAGCCGACGTCCTAGTTGTCTGTGCAACCCCACATCCTTTTCCACTTAACATATATTTTGGGACCTTAGCTGGTGGTCTGGGCTGTTTCCCTTTCGACTACGGATCTTATCACTCGCAGTCTGACTCCCGGATATAAATGAATGGCATTCGGAGTTTATCTGAATTCGGTAACCCGAGATGGGCCCCTAGTCCAAACAGTGCTCTACCTCCATCATTCTCAATTCCGAGGCTAGCCCTAAAGCTATTTCGGAGAGAACCAGCTATCTCCAAGTTCGTTTGGAATTTCTTCCCCCGCTACCCACACCTCATCCCCGCACTTTTCAACGTACGTGGGTTCGGTCCTCCAGTGCGTTTTACCGCACCTTCAACCTGGACATGGGTAGATCACATGGTTTCGGGTCTACGACTACATACTCATTCGCCCTATTCAGACTCGCTTTCGCTGCGGCTCCGTCTCTTCGACTTAACCTCGCATGCAATCGTAACTCGCCGGTTCATTCTACAAAAGGCACGCCATCACTCATTAACGAGCTTTGACTTGTTGTAGGCACACGGTTTCAGGATCTATTTCACTCCCCTTTTGCAGAGATCTCTCACCTTTCCCTCACGGTACTGGTTCACTATCGGTCACTAGGGAGTATTTAGCCTTGCGGGATGGTCCCCGCGGATTCCGACGGAATTTCTCGTGTTCCGCCGTACTCAGGATCCTCCTAGGTGTTGTCAGCATTTCGTCTACGGGGCTTTCACCCTCTTTAGCGGAATTTTCCAAATCCTTCAACTATACTAACAGACTACCATATTGGAGTCCTACAACCCCAACAAGCAAGCTTGTTGGTTTGGGCTCTTCCCGTTTCGCTCGCCGCTACTCAGGGAATCGAATTTTCTTTCTCTTCCTGCAGGTACTAAGATGTTTCAGTTCTCTGCGTCTACCTCTAATCAGCTATGTATTCACTGAAAAGTAATATCCTATAAAAGATATTGG"


#######################################################################################
#######################################################################################
##########################                               ##############################
##########################           MAIN BODY           ##############################
##########################                               ##############################
#######################################################################################
#######################################################################################
slurm = TRUE
# If slurm==TRUE it is likely that a batch run is being undertaken with the read .fastq file to be 
# processed being an input parameter of the slurm call.
# If slurm==FALSE it is probably the case that the .fastq file of reads is imported from ENA or 
# some such and has the extension .fastq
# Note that, in both cases, the fastq files to be processed are in the same directory, file.path(basepath,"fastq").

cat("\nRunning extract_Sereika_16S23S_v10.R code.  \n\n")
lothresh16S = 1250;  hithresh16S = 1700
lothresh23S = 2000;  hithresh23S = 2650
lothreshrrn = 3800;  hithreshrrn = 5600
lothreshZ2 = 150;   hithreshZ2 = 350
lothreshZ5 = 500;   hithreshZ5 = 750

boundarySet = c(1,2,3,24,26,30,32,4)   # Configured for 16S, 23S, Z2, Z5 (not 23Sm)
mapjb = rbind(c(1,2),c(3,8),c(4,5),c(6,7))
rownames(mapjb) = c("16S","23S","Z2","Z5");  colnames(mapjb) = c("5'","3'")
# Note that, through the definition of the boundarySet values for 23S (for example) the
# nominal "23S" boundaries can be changed - e.g. if boundarySet[c(3,8)] = c(22,38) the
# nominal "23S" boundaries are c( Z1_5', Z8_3').
npr = rep(0,length(boundarySet)) 
for (k in 1:length(boundarySet)){npr[k] = length(primer_sets[[boundarySet[k]]])}
ibstart = c(1,cumsum(npr)[1:(length(npr)-1)]+1)
rrnsubuPrimers = select_boundarySet_primers(boundarySet)
nsu = length(rrnsubuPrimers)

# PART 1: 
# Process each of the files that contains the results of aligning the set of primers
# to a block of fastq read sequences.
iseven=function(j){(j%%2)==0}
innameStem0 = "split_zymo_hmw_r104_"        #  "split_zymo_hmw_r104"        # "split_R104-202001-"
# innameStem = "rrn_subunits_Z2Z5_align_indices_names_03042023_split_zymo_hmw_r104"        
# "rrn_subunits_align_indices_names_03042023split_zymo_hmw_r104"
boundSets = list(c(ibstart[1]:(ibstart[1]+npr[1]-1)),c(ibstart[2]:(ibstart[2]+npr[2]-1)),
                   c(ibstart[3]:(ibstart[3]+npr[3]-1)),c(ibstart[4]:(ibstart[4]+npr[4]-1)),
                     c(ibstart[5]:(ibstart[5]+npr[5]-1)),c(ibstart[6]:(ibstart[6]+npr[6]-1)),
                       c(ibstart[7]:(ibstart[7]+npr[7]-1)),c(ibstart[8]:(ibstart[8]+npr[8]-1)) )
kblock =  0
# Depending on how the Sereika et al. zymo_hmw_r104 data is split there will be variation in how 
# many split files to process.  The unix split function uses lower case alphabetic indexing.
# Here we assume the letter indexing runs from aa to ck.
llb =  c(1,1,1)       
lub = c(26,26,24)
BNDS = vector("list",sum(lub))
Names16 = vector("list",sum(lub))
Names23 = vector("list",sum(lub))
NamesZ2 = vector("list",sum(lub))
NamesZ5 = vector("list",sum(lub))
for (k1 in 1:length(llb)){
  j1 = letters[k1]
  for (k2 in llb[k1]:lub[k1]){
    j2 = letters[k2]
    kblock = kblock + 1
    ib =  vector("list",nsu)
    inname0 = paste(innameStem0,paste(j1,j2,sep=""),sep="")
    cat("Processing file ",inname0,"\n")
    F = read_filter_frags(inname0,fastqpath,lothresh16S,numcores,makeplot=FALSE)  
    ilen = F$ifilt
    frags = F$frags
    Fq = F$Fq
    tic() 
    suAlignPrimer = mclapply(rrnsubuPrimers,primerAlign,frags,mc.cores = numcores)
    toc()
    sucount = 0
    # Identify which reads have a low Levdist, and either a high score or PID
    for (jb in boundarySet){
      np = length(primer_sets[[jb]])
      for (kp in 1:np){
        sucount = sucount + 1
        test1 = suAlignPrimer[[sucount]]
        plen = nchar(primers[[primer_sets[[jb]][kp]]]$sequence)
        i11 = which(test1[,"score"]>(plen-3*2))
        i12 = which(test1[,"PID"]>100*((plen-2)/plen))
        i13 = which(test1[,"Levdist"]<ifelse(plen>15,3,2))
        ib[[sucount]] = list(intersectAll=intersect(i13,union(i11,i12)),i11=i11,i12=i12,i13=i13)
        #        print(c(jb, np, kp,length(i11),length(i12),length(i13),length(ib[[sucount]]$intersectAll)))
      }    #   end   kp   loop
    }      #   end   jb   loop
    #     rm(F)
    fragsQ = quality(Fq[ilen])
    
    zz=as.character(subseq(Fq[ilen]@id, start=1, end=36, width=NA))
    fragNames = unlist(zz)
    # We will save this data to facilitate downstream development without the need for boundary re-computation.
    outname = paste("rrn_subunits_Z2Z5_align_indices_names_",out_dateString,"_",inname0,".RData",sep="")
    save(file=file.path(outpath,outname),ilen,suAlignPrimer,ib,frags,fragsQ,fragNames)
    toc()
    
    # Now identify alignments for each of the 4 boundary types that are of
    # sufficient quality to be accepted as probably correct.
    # 16S 5' and 3' - both boundaries have 2 primers
    ibSelectMethod = "union"
    sucount = 0
    jbm = mapjb["16S","5'"]
    jb = boundarySet[jbm]
    np = length(primer_sets[[jb]])
    if (ibSelectMethod=="intersect"){
      ind16S5p = intersect(ib[[sucount+1]]$intersectAll,ib[[sucount+2]]$intersectAll)
    } else {
      ind16S5p = union(ib[[sucount+1]]$intersectAll,ib[[sucount+2]]$intersectAll)
    }
    sucount = sucount + np
    jbm = mapjb["16S","3'"]
    jb = boundarySet[jbm]
    np = length(primer_sets[[jb]])
    if (ibSelectMethod=="intersect"){
      ind16S3p = intersect(ib[[sucount+1]]$intersectAll,ib[[sucount+2]]$intersectAll)
    } else {
      ind16S3p = union(ib[[sucount+1]]$intersectAll,ib[[sucount+2]]$intersectAll)
    }
    sucount = sucount + np
    ind16S = intersect(ind16S5p,ind16S3p)
    ind16Sx = union(ind16S5p,ind16S3p)
    
    # 23S  5' and 3' (possibly modified - handled by boundarySet)
    ibSelectMethod = "union"
    jbm = mapjb["23S","5'"]
    jb = boundarySet[jbm]
    np = npr[jbm]
    sucount = ibstart[jbm] - 1
    if (ibSelectMethod=="intersect"){
      ind23S5p = intersect(ib[[sucount+1]]$intersectAll,ib[[sucount+2]]$intersectAll)
    } else {
      ind23S5p = union(ib[[sucount+1]]$intersectAll,ib[[sucount+2]]$intersectAll)
    }
    if (np>2){
      for (kp in (3:np)){
        if (ibSelectMethod=="intersect"){
          ind23S5p = intersect(ind23S5p,ib[[sucount+kp]]$intersectAll)
        } else {
          ind23S5p = union(ind23S5p,ib[[sucount+kp]]$intersectAll)
        }
      }
    }
    jbm = mapjb["23S","3'"]
    jb = boundarySet[jbm]
    np = npr[jbm]
    sucount = ibstart[jbm] - 1
    if (ibSelectMethod=="intersect"){
      ind23S3p = intersect(ib[[sucount+1]]$intersectAll,ib[[sucount+2]]$intersectAll)
    } else {
      ind23S3p = union(ib[[sucount+1]]$intersectAll,ib[[sucount+2]]$intersectAll)
    }
    if (np>2){
      for (kp in (3:np)){
        if (ibSelectMethod=="intersect"){
          ind23S3p = intersect(ind23S3p,ib[[sucount+kp]]$intersectAll)
        } else {
          ind23S3p = union(ind23S3p,ib[[sucount+kp]]$intersectAll)
        }
      }
    }
    ind23S = intersect(ind23S5p,ind23S3p)
    
    # Z2 of 23S  5' and 3'
    ibSelectMethod = "union"
    jbm = mapjb["Z2","5'"]
    jb = boundarySet[jbm]
    np = npr[jbm]
    sucount = ibstart[jbm] - 1
    if (ibSelectMethod=="intersect"){
      indZ25p = intersect(ib[[sucount+1]]$intersectAll,ib[[sucount+2]]$intersectAll)
    } else {
      indZ25p = union(ib[[sucount+1]]$intersectAll,ib[[sucount+2]]$intersectAll)
    }
    if (np>2){
      for (kp in (3:np)){
        if (ibSelectMethod=="intersect"){
          indZ25p = intersect(indZ25p,ib[[sucount+kp]]$intersectAll)
        } else {
          indZ25p = union(indZ25p,ib[[sucount+kp]]$intersectAll)
        }
      }
    }
    jbm = mapjb["Z2","3'"]
    jb = boundarySet[jbm]
    np = npr[jbm]
    sucount = ibstart[jbm] - 1
    if (ibSelectMethod=="intersect"){
      indZ23p = intersect(ib[[sucount+1]]$intersectAll,ib[[sucount+2]]$intersectAll)
    } else {
      indZ23p = union(ib[[sucount+1]]$intersectAll,ib[[sucount+2]]$intersectAll)
    }
    if (np>2){
      for (kp in (3:np)){
        if (ibSelectMethod=="intersect"){
          indZ23p = intersect(indZ23p,ib[[sucount+kp]]$intersectAll)
        } else {
          indZ23p = union(indZ23p,ib[[sucount+kp]]$intersectAll)
        }
      }
    }
    indZ2 = intersect(indZ25p,indZ23p)
    
    # Z5 of 23S  5' and 3'
    ibSelectMethod = "union"
    jbm = mapjb["Z5","5'"]
    jb = boundarySet[jbm]
    np = npr[jbm]
    sucount = ibstart[jbm] - 1
    if (ibSelectMethod=="intersect"){
      indZ55p = intersect(ib[[sucount+1]]$intersectAll,ib[[sucount+2]]$intersectAll)
    } else {
      indZ55p = union(ib[[sucount+1]]$intersectAll,ib[[sucount+2]]$intersectAll)
    }
    if (np>2){
      for (kp in (3:np)){
        if (ibSelectMethod=="intersect"){
          indZ55p = intersect(indZ55p,ib[[sucount+kp]]$intersectAll)
        } else {
          indZ55p = union(indZ55p,ib[[sucount+kp]]$intersectAll)
        }
      }
    }
    jbm = mapjb["Z5","3'"]
    jb = boundarySet[jbm]
    np = npr[jbm]
    sucount = ibstart[jbm] - 1
    if (ibSelectMethod=="intersect"){
      indZ53p = intersect(ib[[sucount+1]]$intersectAll,ib[[sucount+2]]$intersectAll)
    } else {
      indZ53p = union(ib[[sucount+1]]$intersectAll,ib[[sucount+2]]$intersectAll)
    }
    if (np>2){
      for (kp in (3:np)){
        if (ibSelectMethod=="intersect"){
          indZ53p = intersect(indZ53p,ib[[sucount+kp]]$intersectAll)
        } else {
          indZ53p = union(indZ53p,ib[[sucount+kp]]$intersectAll)
        }
      }
    }
    indZ5 = intersect(indZ55p,indZ53p)
    print(c(length(ind16S),length(ind23S),length(indZ2),length(indZ5)))

    suAP = suAlignPrimer
    n16S = length(ind16S);   n23S = length(ind23S);  nZ2 = length(indZ2);    nZ5 = length(indZ5) 
    
    # Create the dataframe whose columns correspond to the different primers and give
    # the start position of the alignment of that primer, together with a final column
    # which provides a median-based best start position.
    Df = vector("list",length(npr))
    for (jbm in 1:length(npr)){
      np = npr[jbm]
      cnames = paste("p",1:np,sep="")
      if (jbm %in% mapjb["16S",]){indS=ind16S} else if (jbm %in% mapjb["23S",]) {indS=ind23S}  else if (jbm %in% mapjb["Z2",]) {indS=indZ2} else if (jbm %in% mapjb["Z5",]) {indS=indZ5}
      df = data.frame(p1 = suAP[[boundSets[[jbm]][1]]][indS,"startpos"])
      for (k in 2:np){
        z = suAP[[boundSets[[jbm]][k]]][indS,"startpos"]
        df = mutate(df,new=z)
        colnames(df) = cnames[1:k]
      }  #    end    k       loop
      bestcolval = rowMedians(as.matrix(df))
      ix = which(rowVars(as.matrix(df))>0)
      if (length(ix)>0){
        xr = df[ix,]
        mxr = as.matrix(xr)
        bestcolval[ix] = sapply(1:nrow(mxr),function(j){fmode(mxr[j,],ties="first")})
      }
      df = mutate(df,pos=bestcolval)
      Df[[jbm]] = df
    }    #    end    jb      loop
    # Now identify those pairs of boundaries that are consistent with being 16S, or 23S, or Z2, or Z5.
    suLen16 = sapply(1:n16S,function(j){x=Df[[mapjb["16S","3'"]]][j,"pos"] - Df[[mapjb["16S","5'"]]][j,"pos"]})
    suStrand16 = sign(suLen16)
    i16g = which(sapply(1:n16S,function(j){out = abs(suLen16[j])>lothresh16S && abs(suLen16[j])<hithresh16S}))
    
    suLen23 = sapply(1:n23S,function(j){x=Df[[mapjb["23S","3'"]]][j,"pos"] - Df[[mapjb["23S","5'"]]][j,"pos"]})
    suStrand23 = sign(suLen23)
    i23g = which(sapply(1:n23S,function(j){out = abs(suLen23[j])>lothresh23S && abs(suLen23[j])<hithresh23S}))
    
    suLenZ2 = sapply(1:nZ2,function(j){x=Df[[mapjb["Z2","3'"]]][j,"pos"] - Df[[mapjb["Z2","5'"]]][j,"pos"]})
    suStrandZ2 = sign(suLenZ2)
    iZ2g = which(sapply(1:nZ2,function(j){out = abs(suLenZ2[j])>lothreshZ2 && abs(suLenZ2[j])<hithreshZ2}))
    
    suLenZ5 = sapply(1:nZ5,function(j){x=Df[[mapjb["Z5","3'"]]][j,"pos"] - Df[[mapjb["Z5","5'"]]][j,"pos"]})
    suStrandZ5 = sign(suLenZ5)
    iZ5g = which(sapply(1:nZ5,function(j){out = abs(suLenZ5[j])>lothreshZ5 && abs(suLenZ5[j])<hithreshZ5}))
    
    # Now identify matching pairs of 16S 5' and 23S 3' so as to identify possible rrn sequences.
    indrrn = intersect(ind16S,ind23S);   nrrn = length(indrrn)
    irrn16 = unlist(sapply(1:length(indrrn), function(j){which(ind16S==indrrn[j])}))
    irrn23 = unlist(sapply(1:length(indrrn),function(j){which(ind23S==indrrn[j])}))
    suLenrrn = sapply(1:nrrn,function(j){x=Df[[mapjb["23S","3'"]]][irrn23[j],"pos"] - Df[[mapjb["16S","5'"]]][irrn16[j],"pos"]})
    suStrandrrn = sign(suLenrrn)
    irrng = which(sapply(1:nrrn,function(j){out = abs(suLenrrn[j])>lothreshrrn && abs(suLenrrn[j])<hithreshrrn}))
    cat("\n Subunit Identification Performance: 16S ",length(i16g), n16S, length(i16g)/n16S*100,
        "   23S",length(i23g), n23S, length(i23g)/n23S*100,
        "   rrn",length(irrng), nrrn, length(irrng)/nrrn*100,"\n")
    # Now create fr16, qual16, and fr16Names where these are for those readsof the ONT-sequenced data that
    # have a 16S sequence contained within them.  Analogously for the 23S, Z2 and Z5 regions.
    fr16 = rep("",length(i16g))
    n16 = length(i16g)
    fr16Names = rep("",n16)
    qual16 = vector("list",n16)
    bounds16 = matrix(0,nrow=n16,ncol=2)
    for (j in 1:n16){
      p1 = Df[[1]][i16g[j],"pos"];    p2 = Df[[2]][i16g[j],"pos"]
      fwd = p2>p1
      lb = min(p1,p2);  ub = max(p1,p2) + ifelse(fwd,nchar(primers[[3]]$sequence),nchar(primers[[1]]$sequence)) - 1
      bounds16[j,] = c(lb,ub)
      fr16[j] = substr(toString(frags[ind16S[i16g[j]]][[1]]),start=bounds16[j,1],stop=bounds16[j,2])
      fr16Names[j] = fragNames[ind16S[i16g[j]]]
      qual16[[j]] = BString(fragsQ[[ind16S[i16g[j]]]],start=lb,nchar=(ub-lb+1))
      
    }
    Qual16 = BStringSet(qual16)
    
    fr23 = rep("",length(i23g))
    n23 = length(i23g)
    fr23Names = rep("",n23)
    qual23 = vector("list",length(i23g))
    bounds23 = matrix(0,nrow=length(i23g),ncol=2)
    for (j in 1:n23){
      p3 = Df[[mapjb["23S","5'"]]][i23g[j],"pos"];    p4 = Df[[mapjb["23S","3'"]]][i23g[j],"pos"]
      lb = min(p3,p4);  ub = max(p3,p4)
      bounds23[j,] = c(lb,ub)
      fr23[j] = substr(toString(frags[ind23S[i23g[j]]][[1]]),start=lb,stop=ub)
      fr23Names[j] = fragNames[ind23S[i23g[j]]]
      qual23[[j]] = BString(fragsQ[[ind23S[i23g[j]]]],start=lb,nchar=(ub-lb+1))
    }
    Qual23 = BStringSet(qual23)
    
    frZ2 = rep("",length(iZ2g))
    nZ2 = length(iZ2g)
    frZ2Names = rep("",nZ2)
    qualZ2 = vector("list",length(iZ2g))
    boundsZ2 = matrix(0,nrow=length(iZ2g),ncol=2)
    for (j in 1:nZ2){
      p5 = Df[[mapjb["Z2","5'"]]][iZ2g[j],"pos"];    p6 = Df[[mapjb["Z2","3'"]]][iZ2g[j],"pos"]
      lb = min(p5,p6);  ub = max(p5,p6)
      boundsZ2[j,] = c(lb,ub)
      frZ2[j] = substr(toString(frags[indZ2[iZ2g[j]]][[1]]),start=lb,stop=ub)
      frZ2Names[j] = fragNames[indZ2[iZ2g[j]]]
      qualZ2[[j]] = BString(fragsQ[[indZ2[iZ2g[j]]]],start=lb,nchar=(ub-lb+1))
    }
    QualZ2 = BStringSet(qualZ2)
    
    frZ5 = rep("",length(iZ5g))
    nZ5 = length(iZ5g)
    frZ5Names = rep("",nZ5)
    qualZ5 = vector("list",length(iZ5g))
    boundsZ5 = matrix(0,nrow=length(iZ5g),ncol=2)
    for (j in 1:nZ5){
      p7 = Df[[mapjb["Z5","5'"]]][iZ5g[j],"pos"];    p8 = Df[[mapjb["Z5","3'"]]][iZ5g[j],"pos"]
      lb = min(p7,p8);  ub = max(p7,p8)
      boundsZ5[j,] = c(lb,ub)
      frZ5[j] = substr(toString(frags[indZ5[iZ5g[j]]][[1]]),start=lb,stop=ub)
      frZ5Names[j] = fragNames[indZ5[iZ5g[j]]]
      qualZ5[[j]] = BString(fragsQ[[indZ5[iZ5g[j]]]],start=lb,nchar=(ub-lb+1))
    }
    QualZ5 = BStringSet(qualZ5)
    
    do.rrn = FALSE
    if (do.rrn){
      frrrn = rep("",length(irrng))
      qualrrn = vector("list",length(irrng))
      boundsrrn = matrix(0,nrow=length(irrng),ncol=2)
      for (j in 1:length(irrng)){
        pr1 = Df[[mapjb["16S","5'"]]][indrrn[irrng[j]],"pos"];    pr3 = Df[[mapjb["23S","5'"]]][indrrn[irrng[j]],"pos"]
        pr2 = Df[[mapjb["16S","3'"]]][indrrn[irrng[j]],"pos"];    pr4 = Df[[mapjb["23S","3'"]]][indrrn[irrng[j]],"pos"]
        lb = min(pr1,pr2,pr3,pr4);  ub = max(pr1,pr2,pr3,pr4)
        boundsrrn[j,] = c(lb,ub)
        frrrn[j] = substr(toString(frags[indrrn[irrng[j]]][[1]]),start=lb,stop=ub)
        qualrrn[[j]] = BString(fragsQ[[indrrn[irrng[j]]]],start=lb,nchar=(ub-lb+1))
      }
      Qualrrn = BStringSet(qualrrn)
    }     #     end     do.rrn    conditional loop
    
    
    # Re-orientation 
    fror = reorient16S(fr16, refOrient16S)
    fr16 = fror$seq
    ireor16 = which(fror$strand=="rc")
    qual16[ireor16] = sapply(ireor16,function(j){out=reverse(qual16[[j]])})
    
#   (For 23Sm) lsuBnds = c(64,2139)  # See AdHoc code below (~line 675) for determination of 23Sm, Z2 and Z3 boundary locations.
    lsuBnds = c(1,nchar(refOrient23S_full))  # See AdHoc code below (~line 928) for determination of 23Sm, Z2 and Z3 boundary locations.
    refOrient23S = substr(refOrient23S_full,start=lsuBnds[1],stop=lsuBnds[2])    
    fror = reorient23S(fr23, refOrient23S)
    fr23 = fror$seq
    ireor23 = which(fror$strand=="rc")
    qual23[ireor23] = sapply(ireor23,function(j){out=reverse(qual23[[j]])})
    
    Z2Bnds = c(116,369)
    refOrientZ2 = substr(refOrient23S_full,start=Z2Bnds[1],stop=Z2Bnds[2])
    fror = reorientZ(frZ2, refOrientZ2)
    frZ2 = fror$seq
    ireorZ2 = which(fror$strand=="rc")
    qualZ2[ireorZ2] = sapply(ireorZ2,function(j){out=reverse(qualZ2[[j]])})
    
    Z5Bnds = c(980,1586)
    refOrientZ5 = substr(refOrient23S_full,start=Z5Bnds[1],stop=Z5Bnds[2])
    fror = reorientZ(frZ5, refOrientZ5)
    frZ5 = fror$seq
    ireorZ5 = which(fror$strand=="rc")
    qualZ5[ireorZ5] = sapply(ireorZ5,function(j){out=reverse(qualZ5[[j]])})
    
   
    # Now create fastq files for each of the set of 16S sequences, the 
    # set of 23S sequences, and the set of rrn sequences.
    # Note that writeFastq(), by default, will output a gzip-compressed fastq file. This can be 
    # changed using   compress=FALSE.  Otherwise the output file would be better named
    # by having the name being <filename_stem.fastq.gz> - to be fixed later, below.
    #
    n16 = length(fr16)
    frD16 = DNAStringSet(x=fr16,start=rep(1,n16), end=sapply(1:n16,function(j){nchar(fr16[j])}), width=NA)
    SRQ16 = ShortReadQ(sread=frD16,quality=BStringSet(qual16), id= BStringSet(fr16Names))
    outfastqname = paste("amplicon16S_",datasetID,"_",j1,j2,"_trimmed_",out_dateString,".fastq",sep="")
    writeFastq(SRQ16,file=file.path(outpath,outfastqname), mode="w",compress=FALSE)
    
    n23 = length(fr23)
    frD23 = DNAStringSet(x=fr23,start=rep(1,n23), end=sapply(1:n23,function(j){nchar(fr23[j])}), width=NA)
    SRQ23 = ShortReadQ(sread=frD23,quality=BStringSet(qual23), id= BStringSet(fr23Names))
    outfastqname = paste("amplicon23S_",datasetID,"_",j1,j2,"_trimmed_",out_dateString,".fastq",sep="")
    writeFastq(SRQ23,file=file.path(outpath,outfastqname), mode="w",compress=FALSE)
    
    nZ2 = length(frZ2)
    frDZ2 = DNAStringSet(x=frZ2,start=rep(1,nZ2), end=sapply(1:nZ2,function(j){nchar(frZ2[j])}), width=NA)
    SRQZ2 = ShortReadQ(sread=frDZ2,quality=BStringSet(qualZ2), id= BStringSet(frZ2Names))
    outfastqname = paste("ampliconZ2_",datasetID,"_",j1,j2,"_trimmed_",out_dateString,".fastq",sep="")
    writeFastq(SRQZ2,file=file.path(outpath,outfastqname), mode="w",compress=FALSE)
    
    nZ5 = length(frZ5)
    frDZ5 = DNAStringSet(x=frZ5,start=rep(1,nZ5), end=sapply(1:nZ5,function(j){nchar(frZ5[j])}), width=NA)
    SRQZ5 = ShortReadQ(sread=frDZ5,quality=BStringSet(qualZ5), id= BStringSet(frZ5Names))
    outfastqname = paste("ampliconZ5_",datasetID,"_",j1,j2,"_trimmed_",out_dateString,".fastq",sep="")
    writeFastq(SRQZ5,file=file.path(outpath,outfastqname), mode="w",compress=FALSE)
    
    if (do.rrn){
      nrrn = length(frrrn)
      frDrrn = DNAStringSet(x=frrrn,start=rep(1,nrrn), end=sapply(1:nrrn,function(j){nchar(frrrn[j])}), width=NA) 
      SRQrrn = ShortReadQ(sread=frDrrn,quality=Qualrrn)
      outfastqname = paste("ampliconrrn",j1,j2,"_",out_dateString,".fastq",sep="")
      writeFastq(SRQrrn,file=file.path(outpath,outfastqname), mode="w",compress=FALSE)
    }     #     end     do.rrn    conditional loop
    BNDS[[kblock]] = list(bnds16 = bounds16, bnds23 = bounds23, bndsZ2 = boundsZ2, bndsZ5 = boundsZ5)
    Names16[[kblock]] = fr16Names
    Names23[[kblock]] = fr23Names
    NamesZ2[[kblock]] = frZ2Names
    NamesZ5[[kblock]] = frZ5Names
    toc() 
  }        #    end    k2     loop 
}          #    end    k1      loop
outname = paste("Sereika_frag_amplicon_bounds_",out_dateString,"_", k1,".RData",sep="")
save(file=file.path(outpath,outname),BNDS,fr16,qual16,Names16,
                   fr23,qual23,Names23,frZ2,qualZ2,NamesZ2,frZ5,qualZ5,NamesZ5,
                     ireor16,ireor23,ireorZ2,ireorZ5)
cat("Completed Run  \n")
  
if (slurm){quit(save = "no", status = 0, runLast = TRUE)}
quit(save = "no", status = 0, runLast = TRUE)

###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
# Find where on refOrient23S  the various Zx boundaries of interest are.

ss = DNAString(refOrient23S_full)
for (jb in c(2:4)){
  for (je in 1:2){
    bS = boundarySet[mapjb[jb,je]]
    sp = DNAString(primers[[primer_sets[[bS]][1]]]$sequence)
    a1 = pairwiseAlignment(sp,ss,type="local",substitutionMatrix = mat, gapOpening = gapOpen, gapExtension = gapExtend)
    a1r = pairwiseAlignment(sp,reverseComplement(ss),type="local",substitutionMatrix = mat, gapOpening = gapOpen, gapExtension = gapExtend)
    if (score(a1)>score(a1r)){a0=a1} else {a0 = a1r}
    Levdist = stringDist(c(as.character(pattern(a0)),as.character(subject(a0))))
    cat("mapjb[",jb,",",je,"] ",score(a0),pid(a0),Levdist,"   start:",start(subject(a0)) ,"\n")
  }
}

# Gave output:-

# mapjb[ 2 , 1 ]  16 94.73684 1    start: 64 
# mapjb[ 2 , 2 ]  18 100 0         start: 2139 
# mapjb[ 3 , 1 ]  16 100 0         start: 116 
# mapjb[ 3 , 2 ]  12 93.33333 1    start: 369 
# mapjb[ 4 , 1 ]  10 92.30769 1    start: 980 
# mapjb[ 4 , 2 ]  10 92.30769 1    start: 1586 

