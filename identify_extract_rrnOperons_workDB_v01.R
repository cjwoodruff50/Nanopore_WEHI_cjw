args = commandArgs(trailingOnly=TRUE)
# Rscript to select regions of all genomes from a specified database of bacterial genomes  
# (e.g. Refseq) that would approximately be the 16S and 23S genes.  This code is specifically
# written to deal with species such as Helicobacter pylori in which the 16S and 23S are not
# embedded in an rrn operon. This code is built on  identify_extract_rrnOperons_v21.R of 
# January 2023.
# Revision of identify_extract_rrnOperons_toyDB_v01.R of 29April2023 to remove a lot of redundant
# code, mainly in to process the function  extract_16S23S_singleFolder().
#
# Revision of identify_extract_rrnOperons_toyDB_v02.R of 30 April 2023.
#   Final DB for strain resolution work using Sereika D6322 dataset and Murrell group's RAD
#   Uses full23S.
#
# 19 July 2023                                                   [cjw]
#
#  source("/stornext/Bioinf/data/lab_speed/cjw/microbiome/scripts/Rscripts/identify_extract_rrnOperons_workDB_v01.R")
#
#
#########################################################################
###################                               #######################
###################         INITIALISATIONS       #######################
###################                               #######################
#########################################################################
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
packages <- c("BiocManager","Biostrings","stringr","seqinr","tictoc","latex2exp","ShortRead","DECIPHER",
              "dplyr", "nnet", "MASS","stringdist")    # Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  BiocManager::install(packages[!installed_packages],repos="https://cran.csiro.au/")
}# Packages loading
invisible(lapply(packages, library, character.only = TRUE))
# seqinr        provides c2s() and s2c() functions 
# ShortRead     provides readFasta, readFastq functions
# DECIPHER      provides the OrientNucleotides function 
# nnet          provides the call  which.is.max  which breaks ties at random. 
# MASS          provides rlm function.

out_dateString = "19072023" 

# Need to provide the number of rrn operons for each species.  That would normally require 
# extracting the rrn operons for all species, so that extraction has to be performed before 
# requiring num_ops_per_species.

# Alignment parameters
mat=nucleotideSubstitutionMatrix(match = 1, mismatch = -2, baseOnly = TRUE)
align_params = list(align_type ="global", gapOpen = 2, gapExtend = 2)
gapOpen = 2;     gapExtend = 2


segment = c("16S","23S","V1","V2","V3","V4","V5","V6","V7","V9","Z0","Z1","Z2","Z3","Z4","Z5","Z6","Z7","Z8","Z9","ITS")
nseg = length(segment)
boundary_names = c("16S_5'","16S_3'","23S_5'","23S_3'","V1_5'","V1_3'","V2_5'","V2_3'","V3_5'",
                   "V3_3'","V4_5'","V4_3'","V5_5'","V5_3'","V6_5'","V6_3'","V7_5'","V7_3'","V9_5'",
                   "V9_3'","Z0_5p","Z0_3p","Z1_5p","Z1_3p","Z2_5p","Z2_3p","Z3_5p","Z3_3p","Z4_5p","Z4_3p","Z5_5p",
                   "Z5_3p","Z6_5p","Z6_3p","Z7_5p","Z7_3p","Z8_5p","Z8_3p","Z9_5p","Z9_3p")



#########################################################################
###################                               #######################
###################           FUNCTIONS           #######################
###################                               #######################
#########################################################################
source("/stornext/Bioinf/data/lab_speed/cjw/microbiome/scripts/Rscripts/rrn_operon_in_silico_16S_ITS_23S_primers_30May2023.R")
nprimers = length(primers)

extract_16S23S_singleFolder = function(genomeFilename){
  # Identify all 16S and 23S genes on genome with filename  genomeFilename (e.g. Dg[jg]),
  # and extract each such gene and save 16S genes separately from 23S genes.
  # The strategy is to independently extract the 16S and 23S genes, then check to see whether
  # there are pairs of 16S/23S genes that are sufficiently close to be considered to lie
  # in an rrn operon.  If so the operon itself will be saved also.
  # Note that different sizes for the stepped window and step size are now required
  # compared to rrn operon search.  Also 16S parameters are smaller value than 23S parameters.
  #
  # Calls directly functions find_start_16S_23S_multi(..), extract_return_16S23S_set(..)
  # Calls function  check_for_boundaries(..) through find_start_16S_23S_multi()
  # 
  # Required parameters are:-  refGenomespath, small_primer_sets
  
  # 30 April 2023                                                      [cjw]
  jg = which(Dg == genomeFilename)
  gF = genomeFilename
  genomespath = refGenomespath
  cat("\n \n Commence processing genome ",jg,gF, " to get details from testing for selected boundaries within \n")
  cat("all overlapping windows along the genome, using  small_primer_sets  to select boundaries. \n")
  #  E.coli K-12 MG1655 locs with my code are (27,1492,120,357,538,798,922,1080)
  opwinlen16=3000; ovlp16=1600; opwinlen23=4500; ovlp23=3000
  # These value for opwinlen16 etc. (above) must be consistent with the declared values for
  # function  find_start_16S_23S_multi(...) - declared, at present, as below:- 
  #    find_start_16S_23S_multi(refGenomespath, gF, small16S_primer_sets, small23S_primer_sets,
  #                   opwinlen16=opwinlen16, ovlp16=ovlp16, opwinlen23=opwinlen23, ovlp23=ovlp23)
  tic()
  # The following function call takes about 22 minutes to run.  It cannot be parallelised as it
  # is within a function that is itself normally run from mclapply().
  # Note that the last window, for both 16S and 23S, is generally not shifted by as much as
  # all the other windows, to allow for the non-integral number of step lengths in a genome.
  B1623 = find_start_16S_23S_multi(refGenomespath, gF, small16S_primer_sets, small23S_primer_sets)
  toc()
  B16 = B1623$boundaries16
  B23 = B1623$boundaries23
  chrom.Len = B1623$chromLength
  shortSteps = B1623$finalStep    # gives length of the final truncated windows for 16S and 23S step scanning.
  cat("   Completed processing genome ",jg,gF,"for possible boundaries.  \n\n")
  
  # Gather the stats for each of the overlapping windows.
  nsteps16 = length(B16)
  align_stats16 = matrix(0,nrow=nsteps16,ncol=4*nb16+2)
  for (js in 1:(nsteps16)){
    align_stats16[js,1:nb16] = B16[[js]]$seg_bounds
    align_stats16[js,(nb16+1):(2*nb16)] = B16[[js]]$align_scores
    align_stats16[js,(2*nb16+1):(3*nb16)] = B16[[js]]$percent_id
    align_stats16[js,(3*nb16+1):(4*nb16)] = B16[[js]]$Lev
    align_stats16[js,(4*nb16+1):(4*nb16+2)] = B16[[js]]$segord
  }
  
  nsteps23 = length(B23)
  align_stats23 = matrix(0,nrow=nsteps23,ncol=4*nb23+2)
  for (js in 1:(nsteps23)){
    align_stats23[js,1:nb23] = B23[[js]]$seg_bounds
    align_stats23[js,(nb23+1):(2*nb23)] = B23[[js]]$align_scores
    align_stats23[js,(2*nb23+1):(3*nb23)] = B23[[js]]$percent_id
    align_stats23[js,(3*nb23+1):(4*nb23)] = B23[[js]]$Lev
    align_stats23[js,(4*nb23+1):(4*nb23+2)] = B23[[js]]$segord
  }
  cat("Genome ",jg,gF," nsteps16, nsteps23 ",nsteps16, nsteps23,"\n")
  # Now determine whether there are any 16S genes and, if so, store details of the 16S 
  # boundaries identified and whether the 16S is reversed.
  nwin16 = nsteps16
  lbf16 = length(bseqfwd16)
  diffs = matrix(0,nrow=nwin16,ncol=(lbf16-1))
  for (jw in 1:nwin16){
    diffs[jw,] = sapply(2:lbf16,function(js){align_stats16[jw,bseqfwd16[js]]-align_stats16[jw,bseqfwd16[js-1]]})
  }
  # The first expression for i16Swin below requires all subregions being used to be of consistent sign.  This
  # is too demanding for broad application, so this expression is listed but not currently used.  The second
  # expression is robust to one bad boundary.
  # i16Swin = which( (sapply(1:nwin16,function(jw){abs(sum(sign(diffs[jw,1:(lbf16-1)])))>=(lbf16-2)})) && 
  #                           (sapply(1:nwin16,function(jw){t1 = abs(sum(diffs[jw,])); out=t1>1100 && t1<2200})) &&
  #                               (sapply(1:nwin16,function(jw){sum(align_stats16[jw,nb16+bseqfwd16])>(lbf16*10)})) )
  i16Swin = which(sapply(1:nwin16,function(jw){
                                    b1 = median(align_stats16[jw,(nb16+1):(2*nb16)])>8
                                    b2 = sum(align_stats16[jw,nb16+bseqfwd16])>(lbf16*10)
                                    out = b1 && b2
                                  }))
  if (length(i16Swin)==0){ # The case where no operon has yet been identified.  Loosen the criteria
    i16Swin = which(sapply(1:length(align_stats16[,1]),function(jw){median(align_stats16[jw,9:16])>8}))
  }
  if (length(i16Swin)>0){
    nops = length(i16Swin)
    delw16 = (opwinlen16-ovlp16)
    if (nops>1){
      bounds16 = matrix(0,nrow=nops,ncol=length(bseqfwd16 ))
      bounds16all = matrix(0,nrow=nops,ncol=nb16)
      for (jop in 1:nops){
        jw = i16Swin[jop]
        if (jw<nwin16){
          bounds16[jop,] =sapply(1:length(bseqfwd16),function(k){(jw-1)*delw16 + align_stats16[jw,bseqfwd16[k]]})
          bounds16all[jop,] =sapply(1:nb16,function(k){(jw-1)*delw16 + align_stats16[jw,k]})
        } else { # The end effect dealt with
          bounds16[jop,] =sapply(1:length(bseqfwd16),function(k){(jw-2)*delw16 + B1623$finalStep[1] + align_stats16[jw,bseqfwd16[k]]})
          bounds16all[jop,] =sapply(1:nb16,function(k){(jw-2)*delw16 + B1623$finalStep[1] + align_stats16[jw,k]})
        }
      }
    } else {
      jw = i16Swin[1]
      if (jw<nwin16){
        bounds16 = sapply(1:length(bseqfwd16),function(k){(jw-1)*delw16 + align_stats16[jw,bseqfwd16[k]]})
        bounds16all = sapply(1:nb16, function(k){(jw-1)*delw16 + align_stats16[jw,k]})
      } else {
        bounds16 = sapply(1:length(bseqfwd16),function(k){(jw-2)*delw16 + align_stats16[jw,bseqfwd16[k]]})
        bounds16all = sapply(1:nb16, function(k){(jw-2)*delw16 + B1623$finalStep[1] + align_stats16[jw,k]})
      }
    }
    
    # Now refine the bounds
    cat("Now refining 16S bounds with regression and residuals analysis followed by correction if necessary.",jg,gF,"\n")
    # Don't need refinement if there is no variance in the absolute value of largest difference.
    check = sapply(1:nrow(bounds16),function(j){out=abs(bounds16[j,ncol(bounds16)]-bounds16[j,1])})
    if (var(check)<1){ # No refinement necessary
      newBnds16 = list(bounds16=bounds16, bounds16all=bounds16all) 
      cat("No refinement of 16S boundaries required.\n")
    } else {
      newBnds16 = linfitBnds16(bounds16all,i16Swin)
      bounds16 = newBnds16$bounds16;    bounds16all = newBnds16$bounds16all
      cat("Checking 16S regression success. \n");  print(str(newBnds16))
    }
  } else {  #   !(length(i16Swin)>0)
    print(paste("No 16S found for genome ",gF,sep=""))
    genes16 = list(bounds16=rep(-1,length(bseqfwd16)),is_reversed=-2)
  }
  nops = nrow(bounds16all)
  if (nops==1){
    is_reversed = bounds16[length(bseqfwd16)]<bounds16[1]
    genes16 = list(bounds16=bounds16,is_reversed=is_reversed)
  } else if (nops>1){
    is_reversed = sapply(1:length(bounds16[,1]),function(j){bounds16[j,length(bseqfwd16)]<bounds16[j,1]})
    genes16 = list(bounds16=bounds16,is_reversed=is_reversed)
  }  
  cat("Have constructed the matrix of 16S sequence end points.",jg,gF, "\n")
  
  n16S = 0     # Default value
  alllocs16S = matrix(-9,2,2)     # Default value
  if (length(i16Swin)>0){
    n16S = length(genes16$bounds16)%/%length(bseqfwd16)
    if (n16S>1){ # There are multiple 16S entries
      which_nonzero16S = rep(TRUE,n16S) 
      alllocs16S = matrix(-1,nrow=n16S,ncol=2) 
      colnames(alllocs16S) = c("start16S","stop16S")
      alllocs16S = bounds16[,c(1,lbf16)]
    } else if (n16S==1){
      alllocs16S[1,] = bounds16[c(1,lbf16)]
    }
  } else {
    # Case !(length(i16Swin)>0)
    cat("\n\n HAZARD:     No 16S found Genome",jg,gF,"\n\n")
  }

  cat("Computed genes16 jg=",jg,"\n")
  print(genes16)
  
  # Now determine whether there are any 23S genes and, if so, store details of the 23S 
  # boundaries identified and whether the 23S is reversed.
  nwin23 = nrow(align_stats23)
  lbf23 = length(bseqfwd23)
  diffs = matrix(0,nrow=nwin23,ncol=(length(bseqfwd23)-1))
  for (jw in 1:nwin23){
    diffs[jw,] = sapply(2:length(bseqfwd23),function(js){align_stats23[jw,bseqfwd23[js]]-align_stats23[jw,bseqfwd23[js-1]]})
  }
  # i23Swin = which( sapply(1:nwin23,function(jw){(abs(sum(sign(diffs[jw,1:(length(bseqfwd23)-1)])))==(length(bseqfwd23 )-1)) && 
  #                    (abs(sum(diffs[jw,]))>1600 && abs(sum(diffs[jw,]))<3500) &&
  #                        (sum(align_stats23[jw,nb23+bseqfwd23])>(length(bseqfwd23)*10))}) )
  i23Swin = which(sapply(1:nwin23,function(jw){
    b1 = median(align_stats23[jw,9:16])>8
    b2 = sum(align_stats23[jw,nb23+bseqfwd23])>(lbf23*8)
    out = b1 && b2
  }))
  if (length(i23Swin)==0){ # The case where no operon has yet been identified.  Loosen the criteria
    i23Swin = which(sapply(1:length(align_stats23[,1]),function(jw){median(align_stats23[jw,9:16])>8}))
  }
  if (length(i23Swin)>0){
    delw23 = (opwinlen23 - ovlp23)
    nops = length(i23Swin)
    if (nops>1){
      bounds23 = matrix(0,nrow=nops,ncol=length(bseqfwd23 ))
      bounds23all = matrix(0,nrow=nops,ncol=nb23)
      for (jop in 1:nops){
        jw = i23Swin[jop]
        if (jw<nwin23){
          bounds23[jop,] =sapply(1:length(bseqfwd23),function(k){(jw-1)*delw23 + align_stats23[jw,bseqfwd23[k]]})
          bounds23all[jop,] =sapply(1:nb23,function(k){(jw-1)*delw23 + align_stats23[jw,k]})
        } else { # The end effect dealt with
          bounds23[jop,] =sapply(1:length(bseqfwd23),function(k){(jw-2)*delw23 + B1623$finalStep[2] + align_stats23[jw,bseqfwd23[k]]})
          bounds23all[jop,] =sapply(1:nb23,function(k){(jw-2)*delw23 + B1623$finalStep[2] + align_stats23[jw,k]})
        }
      }
    } else {
      jw = i23Swin[1]
      if (jw<nwin23){
        bounds23 = sapply(1:length(bseqfwd23),function(k){(jw-1)*delw23 + align_stats23[jw,bseqfwd23[k]]})
        bounds23all = sapply(1:nb23, function(k){(jw-1)*delw23 + align_stats23[jw,k]})
      } else {
        bounds23 = sapply(1:length(bseqfwd23),function(k){(jw-2)*delw23 + align_stats23[jw,bseqfwd23[k]]})
        bounds23all = sapply(1:nb23, function(k){(jw-2)*delw23 + B1623$finalStep[1] + align_stats23[jw,k]})
      }
    }
    
    # Now refine the bounds
    cat("Now refining 23S bounds with regression and residuals analysis followed by correction if necessary.",jg,gF,"\n")
    # Don't need refinement if there is no variance in the absolute value of largest difference.
    check = sapply(1:nrow(bounds23),function(j){out=abs(bounds23[j,ncol(bounds23)]-bounds23[j,1])})
    if (var(check)<1){ # No refinement necessary
      newBnds23 = list(bounds23=bounds23, bounds23all=bounds23all) 
      cat("No refinement of 23S boundaries required.\n")
    } else {
      newBnds23 = linfitBnds23(bounds23all,i23Swin)
      bounds23 = newBnds23$bounds23;    bounds23all = newBnds23$bounds23all
      cat("Checking 23S regression success. \n");  print(str(newBnds23))
    }
  } else {  #   !(length(i123Swin)>0)
    print(paste("No 23S found for genome ",gF,sep=""))
    genes23 = list(bounds23=rep(-1,length(bseqfwd23)),is_reversed=-2)
  }
  nops = length(bounds23all) %/% nb23
  if (nops==1){
    is_reversed = bounds23[length(bseqfwd23)]<bounds23[1]
    genes23 = list(bounds23=bounds23,is_reversed=is_reversed)
  } else if (nops>1){
    is_reversed = sapply(1:length(bounds23[,1]),function(j){bounds23[j,length(bseqfwd23)]<bounds23[j,1]})
    genes23 = list(bounds23=bounds23,is_reversed=is_reversed)
  }  
  
  cat("Computed genes23",jg,"\n")
  print(genes23)
  
  n23S = 0     # Default value
  alllocs23S = matrix(-9,2,2)     # Default value
  if (length(i23Swin)>0){
    n23S = length(genes23$bounds23)%/%length(bseqfwd23)
    if (n23S>1){ # There are multiple 23S entries
      which_nonzero236S = rep(TRUE,n23S) 
      alllocs23S = matrix(-1,nrow=n23S,ncol=2) 
      colnames(alllocs23S) = c("start23S","stop23S")
      alllocs23S = bounds23[,c(1,lbf23)]
    } else if (n23S==1){
      alllocs23S[1,] = bounds23[c(1,lbf23)]
    }
  } else {
    # Case !(length(i23Swin)>0)
    cat("\n\n HAZARD:     No 23S found for Genome",jg,gF,"\n\n")
  }
  
  cat("Have constructed the matrix of 23 sequence end points.",jg,gF, "\n")

  
  # Now extract the flanked 16S and 23S sequences.
  gseqs = matrix("",nrow=max(2,n16S,n23S),ncol=3)
  locs16S23S = matrix(-1,nrow=n16S+n23S,ncol=3)
  colnames(locs16S23S) = c("regID","start","end")
  extracted = extract_return_16S23S_set(genomespath, gF, alllocs16S, alllocs23S)
  if (n16S>0){
    ig16 = which(nchar(extracted$seq16S)>1000)
    gseqs[1:length(ig16),2] = extracted$seq16S[ig16]
  } 
  n16S = length(ig16)
  if (n23S>0){
    ig23 = which(nchar(extracted$seq23S)>1700)
    gseqs[1:length(ig23),3] = extracted$seq23S[ig23]
  }
  n23S = length(ig23)
  # Now consolidate the end-boundary locations.
  if ((n16S+n23S) > 0){
    locs16S23S = matrix(-1,nrow=n16S+n23S,ncol=3)
    colnames(locs16S23S) = c("regID","start","end")
    if (n16S>0){
      for (js in 1:length(ig16)){
        locs16S23S[js,2:3] = alllocs16S[js,1:2]
        locs16S23S[js,1] = "16S"
      }
    }
    if (n23S>0){
      for (js in (n16S+1):(n16S+n23S)){
        locs16S23S[js,2:3] = alllocs23S[(js-n16S),1:2]
        locs16S23S[js,1] = "23S"
      }
    }
  } else {  # Dummy matrix with empty string in the "regID" column
    locs16S23S = matrix(-9,nrow=2,ncol=3)
    locs16S23S[1:2,1] = c("","")
  }
  
  out = list(segSeqs=gseqs, stats16 = align_stats16, stats23 = align_stats23, locs=locs16S23S,
               i16Swin=i16Swin, i23Swin=i23Swin, allbounds16=bounds16all,allbounds23=bounds23all)
}

find_start_16S_23S_multi = function(refGenomespath, gF, small16S_primer_sets, small23S_primer_sets,
                                    opwinlen16=3000, ovlp16=1600, opwinlen23=4500, ovlp23=3000){
  # Returns details from each of the step-scanned windows used in searching for possible locations
  # of 16S and 23S sub-region (V and Z) boundaries. 
  # Input parameters:   
  #    The genome of interest is specified by  refGenomespath, gF (gF is the filename)
  #    small16S_primer_sets, small23S_primer_sets  identify the primer sets for the boundaries of interest.
  #           - indices into the list object  primer_sets. These primers target features
  #                           such as the 16S_5', 16S_3', V3_5' etc boundaries
  #    opwinlen16, ovlp16  give the stepped window length and the overlap of consecutive windows, respectively, 
  #                      for the 16S.  Likewise for 23S;
  # This code checks to see whether the .fna file has multiple contigs that are likely
  # to be chromosomes (e.g. Ceriebacter sphaeroides has 2 chromosomes, 2-5 plasmids).  If so
  # it processes all such contigs, and pools the rrn operons found.
  # Update (23 April 2023):  Incorporated a treament of the partial window step needed for the 
  #                          3' end of the genome, for both 16S and 23S.  Includes making these
  #                          two step sizes avavilable external to this function.
  # 23 April 2023                                                         [cjw]
  
  # Read the fasta file specified, check whether it has more than one sequence identified as a chromosome
  # in it. If so concatenate these sequences into gf.
  # Also, ensure that only ACGT entries are present - replacing non-canonical bases  
  # by a random choice between "A","CX","G","T". (although, if IUPAC ambiguity codes it would 
  # be better to constrain the choice according to this code.)  
  # The output from this section is  bacseq.
  jg = which(Dg==gF)
  gen_fa = read.fasta(file=file.path(refGenomespath,gF), seqtype="DNA",forceDNAtolower = FALSE)
  chrom = which(unlist(sapply(1:length(gen_fa),function(j){grep("chromosome",attributes(gen_fa[[j]])$Annot)==1})))
  if (length(chrom)==0){chrom=1}
  nchr = max(1,length(chrom))
  if (nchr==1){
    gen_fa[[1]] = gen_fa[[chrom]]
    chrom.Len = length(gen_fa[[1]])
  } else if (nchr>1) {
    cat("\n Genome ",gF, " has ",nchr, "chromosomes identified.\n")
    chrom.Len = rep(0,nchr)
    gen_fa[[1]] = gen_fa[[chrom[1]]]
    chrom.Len = sapply(1:nchr,function(jc){out = length(gen_fa[[chrom[jc]]])})
#    gen_fa[[1]] = sapply(2:nchr,function(jc){t1 = c(gen_fa[[1]],gen_fa[[chrom[jc]]])})
    gf = gen_fa[[1]]
    for (jj in 2:nchr){gf = append(gf,gen_fa[[jj]])}
  }
  # Some of the genomes have a few "N" entries.Replace non-canonical base by "A" 
  # and advise number of replacements.
  if (nchr==1){
    noncanon = which(!(gen_fa[[1]] %in% c("A","G","C","T") ))
    nc_codes = gen_fa[[1]][noncanon]
    if (length(noncanon)>0){
      gen_fa[[1]][noncanon] = rep(sample(c("A","C","G","T"),1),length(noncanon),replace=TRUE)
    }
    bacseq = c2s(gen_fa[[1]])
  } 
  if (nchr>1){
    noncanon = which(!(gf %in% c("A","G","C","T") ))
    nc_codes = gf[noncanon]
    if (length(noncanon)>0){
      gf[noncanon] = rep(sample(c("A","C","G","T"),1),length(noncanon),replace=TRUE)
    }
    bacseq = c2s(gf)
  }
  
  # Set up for the 16S stepping window, including the alignment process for identifying the boundaries 
  # based on the in-silico primer use.  
  delw16 = (opwinlen16-ovlp16)
  nsteps16 = (nchar(bacseq)-ovlp16)%/%delw16
  shortstep16 = nchar(bacseq)-(nsteps16-1)*delw16 -opwinlen16
  B16 = vector("list",nsteps16+1)
  #  cat("    Checking for genome ",jg," boundaries in ", nsteps, " windows. \n")
  for (jw in 1:nsteps16){
    winseq = substr(bacseq,start=(jw-1)*(delw16) + 1, stop = (jw-1)*delw16 + opwinlen16) 
    B16[[jw]] = check_for_boundaries(winseq,small16S_primer_sets,primer_sets,primers)
  } 
  winseq = substr(bacseq,start=nchar(bacseq)-opwinlen16+1, stop = nchar(bacseq))
  B16[[nsteps16+1]] = check_for_boundaries(winseq,small16S_primer_sets,primer_sets,primers)
  
  delw23 = (opwinlen23-ovlp23)
  nsteps23 = (nchar(bacseq)-ovlp23)%/%delw23
  shortstep23 = nchar(bacseq)-(nsteps23-1)*delw23 - opwinlen23
  B23 = vector("list",nsteps23+1)
  for (jw in 1:nsteps23){
    winseq = substr(bacseq,start=(jw-1)*(delw23) + 1, stop = (jw-1)*delw23 + opwinlen23) 
    B23[[jw]] = check_for_boundaries(winseq,small23S_primer_sets,primer_sets,primers)
  } 
  winseq = substr(bacseq,start=nchar(bacseq)-opwinlen23+1, stop = nchar(bacseq))
  B23[[nsteps23+1]] = check_for_boundaries(winseq,small23S_primer_sets,primer_sets,primers)
  out = list(boundaries16=B16,boundaries23=B23, finalStep = c(shortstep16,shortstep23),
                non_canonical = cbind(noncanon,nc_codes),chromLength=chrom.Len)
}

primerAlign = function(i){
  # Alignment of a single primer, identified by the list member i, against a list
  # of fragment sequences, indexed by indf, frags[jf] where jf is in indf.
  # Returns 3 alignment MOPs (alscore, PID, Levdist) and start_loc of the primer
  # alignment for each fragment sequence.
  # This function can be easily called by mclapply.
  primer = i
  skip = FALSE
  # First check to see if this is a primer set associated with a boundary that has already been determined.
  # This can be done by parsing the primer type of one entry in this primer_set, together with keeping a
  # running tally of which Z region 3' boundaries have already been determined.
  # 
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
  
  if (!skip){
    A = matrix(-1,nrow=nfrags,ncol=4)
    colnames(A) = c("score","PID","Levdist","startpos")
    for (jf in indf){
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

check_for_boundaries = function(sequence,primer_set_index,primer_sets,primers){
  # A modification of   identify_ITS_segment()  to use a specified subset of primers 
  # and work on a single sequence.  Returns the same list of data as   identify_ITS_segment()
  # but just for single sequence.
  # 8 July 2019                                             [cjw]
  substmat = nucleotideSubstitutionMatrix(match = 1, mismatch = -2, baseOnly = TRUE)
  align_parameters = list(align_type="global-local", gapOpen=2, gapExtend=2)
  nseg = length(primer_set_index)
  ascore = rep(0,nseg)
  apid = rep(0,nseg)
  aL = rep(0,nseg)
  aLev = rep(0,nseg)
  D2align = vector("list",nseg)
  
  if (class(sequence)=="SeqFastadna"){
    s2 = DNAString(c2s(sequence))
  } else if (class(sequence)=="character"){
    s2 = DNAString(sequence)
  } else {print("Not suitable fragment type"); stop()}
#  cat("    Computing boundary primer alignments to sequence of length ",length(s2),"\n")
  
  for (jseg in 1:nseg){ # jsb = j_seg_boundary
    jsb = primer_set_index[jseg]
    nt = length(primer_sets[jsb])
    D1align = vector("list",nt)
    scores = rep(0,nt)
    for (jnt in 1:nt){
      reversed = FALSE
      i_prim = primer_sets[[jsb]][jnt]
#      s1 = DNAString(primers[[i_prim]]$sequence)
#      s1r = reverseComplement(s1)
      s1 = primers[[i_prim]]$sequence
      s1r = reverseComplement(DNAString(s1))
      algnf = pairwiseAlignment(s1,s2,type=as.character(align_parameters$align_type),substitutionMatrix=substmat,
                                gapOpening=align_parameters$gapOpen,gapExtension=align_parameters$gapExtend)
      reversed=TRUE
      algnr = pairwiseAlignment(s1r,s2,type=as.character(align_parameters$align_type),substitutionMatrix=substmat,
                                  gapOpening=align_parameters$gapOpen,gapExtension=align_parameters$gapExtend)
      if (score(algnf)>score(algnr)){algn = algnf} else {algn= algnr}
      scores[jnt] = score(algn)
      D1align[[jnt]] = algn
    }
    jnt_max = which.max(scores)
    align = D1align[[jnt_max]]
    s1 = align@pattern
    D2align[[jseg]] = align
    ascore[jseg] = scores[jnt_max]
    apid[jseg] = pid(align)
    aL[jseg] = start(Views(align,start=0))
    s2start = align@subject@range@start
    s2stop = s2start + align@subject@range@width-1
    s3=substr(s2,start=s2start, stop=s2stop)
    aLev[jseg] = stringDist(c(as.character(s1),as.character(s3)))
  }
  # Now check to see if the segment boundaries are consistently ordered.
  diffsign = sapply(2:nseg,function(jseg){sign(aL[jseg]-aL[jseg-1])})
  is_properly_ordered = (abs(sum(diffsign))==(nseg-1))
  is_reversed = (is_properly_ordered & diffsign[1]<1)
  segord = c(is_properly_ordered,is_reversed)
  out = list(seg_bounds = aL, align_scores = ascore, percent_id = apid, Lev = aLev, 
             segord = cbind(is_properly_ordered,is_reversed))
#  list(seg_bounds = aL, align_scores = ascore, percent_id = apid, Lev = aLev, 
#             segord = cbind(is_properly_ordered,is_reversed), alignments = D2align)
}

extract_return_16S23S_set = function(genomespath, bac_fasta_filename, alllocs16S, alllocs23S){
  # Given a set of single-genome locations at which 16S and 23S start and end, returns
  # the set of 16S and 23S sequences with a 40 base flanking region at each end.
  # Input parameters:   genomespath, bac_fasta_filename, to specify the genome, 
  #                     alllocs16S is a matrix of start and end genomic locations for
  #                       each 16S found;  alllocs23S is the 23S corresponding data.
  # Derived from extract_return_16S23Srrn() of 21 September 2022 which could process 
  #    genomes with more than 1 chromosome, and also dealt separately with 16S and 23S.
  # 19 October 2022
  
  # Read the fasta file specified, check whether it has more than one sequence in it, and if so select the 
  # the longest of these sequences.  Also, ensure that only ACGT entries are present - replacing 
  # non-canonical bases by a random sample of the canonical bases.   The output from this section is  bacseq.
  jg = which(Dg==bac_fasta_filename)
  gF = bac_fasta_filename
  n16S = ifelse(alllocs16S[1,2]==(-9),0,nrow(alllocs16S))
  n23S = ifelse(alllocs23S[1,2]==(-9),0,nrow(alllocs23S))
  gen_fa = read.fasta(file=file.path(genomespath,bac_fasta_filename), seqtype="DNA",forceDNAtolower = FALSE)
  chrom = which(unlist(sapply(1:length(gen_fa),function(j){grep("chromosome",attributes(gen_fa[[j]])$Annot)==1})))
  if (length(chrom)==0){chrom=1}
  nchr = max(1,length(chrom))
  if (nchr==1){
    gen_fa[[1]] = gen_fa[[chrom]]
    chrom.Len = length(gen_fa[[1]])
  } else if (nchr>0) {
    cat("\n Genome ",gF, " has ",nchr, "chromosomes identified.\n")
    chrom.Len = rep(0,nchr)
    gen_fa[[1]] = gen_fa[[chrom[1]]]
    chrom.Len = sapply(1:nchr,function(jc){out = length(gen_fa[[chrom[jc]]])})
    #    gen_fa[[1]] = sapply(2:nchr,function(jc){t1 = c(gen_fa[[1]],gen_fa[[chrom[jc]]])})
    gf = sapply(2:nchr,function(jc){t1 = c(gen_fa[[1]],gen_fa[[chrom[jc]]])})
  }
  # Some of the genomes have a few "N" entries.Replace non-canonical base by "A" 
  # and advise number of replacements.
  if (nchr==1){
    noncanon = which(!(gen_fa[[1]] %in% c("A","G","C","T") ))
    nc_codes = gen_fa[[1]][noncanon]
    if (length(noncanon)>0){
      gen_fa[[1]][noncanon] = rep(sample(c("A","C","G","T"),1),length(noncanon),replace=TRUE)
    }
    bacseq = c2s(gen_fa[[1]])
  } 
  if (nchr>1){
    noncanon = which(!(gf %in% c("A","G","C","T") ))
    nc_codes = gf[noncanon]
    if (length(noncanon)>0){
      gf[noncanon] = rep(sample(c("A","C","G","T"),1),length(noncanon),replace=TRUE)
    }
    bacseq = c2s(gf)
  }
  extend = 40
  # Process the 16S. Note that we need to deal with start and end points that are too
  # close to the genome ends to have a full flank padding - i.e. extend may be too large.
  seq16S = rep("",max(1,n16S))
  if (n16S>0){
    for (js in 1:n16S){
      lb16S = max(1,min(alllocs16S[js,])-extend)
      ub16S = min(max(alllocs16S[js,])+extend,sum(chrom.Len))
      seq16S[js] = substr(bacseq,start=lb16S, stop=ub16S)
    }
  }
  # Process the 23S
  seq23S = rep("",max(1,n23S))
  if (n23S>0){
    for (js in 1:n23S){
      lb23S = max(1,min(alllocs23S[js,])-extend)
      ub23S = min(max(alllocs23S[js,])+extend,sum(chrom.Len))
      seq23S[js] = substr(bacseq,start=lb23S, stop=ub23S)
    }
  }
  out = list(seq16S=seq16S, seq23S=seq23S)
}

linfitBnds16 = function(bounds16all,i16Swin){
  di16win = sapply(2:length(i16Swin),function(j){i16Swin[j] - i16Swin[j-1]})
  i16win2 = which(di16win>1)  # Find the separators of clustered windows.
  i16win2 = append(i16win2,length(i16Swin))
  ntb16 = length(i16win2)
  newbounds16all = matrix(0,nrow=ntb16,ncol=ncol(bounds16all))
  jnc = 0
  for (k in 1:ntb16){
    ncwin = ifelse(k==1,i16win2[1], i16win2[k]-i16win2[k-1] )
    if (ncwin>1){ # Need to deal with 2 or 3 adjacent windows that might have a 16S within them
      df = data.frame(x=EcoliReflocs16,y1=bounds16all[i16win2[k]+1-ncwin,],y2=bounds16all[i16win2[k]+2-ncwin,])
      if (ncwin==3){df =mutate(df,y3=bounds16all[i16win2[k],])}  # It is assumed that ncwin is no greater than 3
      # Compare ordinary least squares (lsq) fit to the data with a robust lsq fit. If the robust version
      # is clearly better then use it. 
      ols = vector("list",ncwin);   robust = ols
      ols[[1]] = lm(y1~x, data=df);   ols[[2]] = lm(y2~x, data=df)   
      robust[[1]] = rlm(y1~x, data=df);   robust[[2]] = rlm(y2~x, data=df)
      fitsum = matrix(0,nrow=ncwin,ncol=2)
      fitsum[1,] = c(summary(ols[[1]])$sigma,summary(robust[[1]])$sigma)
      fitsum[2,] = c(summary(ols[[2]])$sigma,summary(robust[[2]])$sigma) 
      if (ncwin==3){
        ols[[3]] = lm(y3~x, data=df);   robust[[3]] = rlm(y3~x, data=df)
        fitsum[3,] = c(summary(ols[[3]])$sigma,summary(robust[[3]])$sigma)
      }
      jbest = which.min(fitsum);  jcol = 1+(jbest-1)%/%ncwin;  jrow =  jbest - (ncwin*(jcol-1))
      if (jcol==1){
        newbounds16all[k,] = ols[[jrow]]$fitted.values
      } else {
        newbounds16all[k,] = robust[[jrow]]$fitted.values
      }
      jnc = jnc + ncwin
    } else {
      jnc = jnc + 1
      newbounds16all[k,] = bounds16all[jnc,]
    }
  }        #    end    k     loop
  newbounds16all = round(newbounds16all)
  out = list(bounds16=newbounds16all[,bseqfwd16], bounds16all = newbounds16all)
}


linfitBnds23 = function(bounds23all,i23Swin){
  findbnds = function(df){
    # Processes a set of up to 3 rows of 23S bounds to identify a complete 23S.
    # 5 July 2023             [cjw]
    ncwin = ncol(df)-1
    if (ncwin>3){cat("Too many rows in dataframe given to findbnds(). Stopping. \n"); stop}
    if (ncwin<2){cat("Too few rows in dataframe given to findbnds(). Stopping. \n"); stop}
    outbounds = rep(0,nrow(df))
    ols = vector("list",ncwin);   robust = ols
    ols[[1]] = lm(y1~x, data=df);   ols[[2]] = lm(y2~x, data=df)   
    robust[[1]] = rlm(y1~x, data=df);   robust[[2]] = rlm(y2~x, data=df)
    fitsum = matrix(0,nrow=ncwin,ncol=2)
    fitsum[1,] = c(summary(ols[[1]])$sigma,summary(robust[[1]])$sigma)
    fitsum[2,] = c(summary(ols[[2]])$sigma,summary(robust[[2]])$sigma) 
    if (ncwin==3){
      ols[[3]] = lm(y3~x, data=df);   robust[[3]] = rlm(y3~x, data=df)
      fitsum[3,] = c(summary(ols[[3]])$sigma,summary(robust[[3]])$sigma)
    }
    jbest = which.min(fitsum);  jcol = 1+(jbest-1)%/%ncwin;  jrow =  jbest - (ncwin*(jcol-1))
    if (jcol==1){
      bounds = ols[[jrow]]$fitted.values
    } else {
      bounds = robust[[jrow]]$fitted.values
    }
    out = bounds
  }
  
  di23win = sapply(2:length(i23Swin),function(j){i23Swin[j] - i23Swin[j-1]})
  i23win2 = which(di23win>1)  # Find the separators of clustered windows.
  i23win2 = append(i23win2,length(i23Swin))
  ntb23 = length(i23win2)
  newbounds23all = matrix(0,nrow=ntb23+5,ncol=ncol(bounds23all))
  jnc = 0
  allkdone = FALSE
  k = 0;  kw = 0
  while(!allkdone){
    k = k+1; kw = kw+1
    double23S = FALSE
    ncwin = ifelse(k==1,i23win2[1], i23win2[k]-i23win2[k-1] )
    if (ncwin>3){
      double23S = TRUE
      # Need to split the rows of bounds in to two groups such that the two groups
      # probably mainly give data on different 23S sequences. I will use the internal
      # boundaries 3,4 and 5 for this.
      ncwh = ncwin %/%2
      m1 = mean(bounds23all[(i23win2[k]+1-ncwin):(i23win2[k]+ncwh-ncwin),c(3,4,5)])
      m2 = mean(bounds23all[(i23win2[k]+ncwh+1-ncwin):(i23win2[k]),c(3,4,5)])
      if (abs(m2-m1)>4000){ # Clear evidence of two separate 23S sequences. 
        df1 = data.frame(x=EcoliReflocs23,y1=bounds23all[i23win2[k]+1-ncwin,],y2=bounds23all[i23win2[k]+2-ncwin,])
        if (ncwh==3){df =mutate(df,y3=bounds23all[i23win2[k]+ncwh-ncwin,])} 
        newbounds23all[kw,] = findbnds(df1)
        kw = kw+1
        df2 = data.frame(x=EcoliReflocs23,y1=bounds23all[i23win2[k]+ncwh+1-ncwin,],y2=bounds23all[i23win2[k]+ncwh+2-ncwin,])
        if ((ncwh+2) < ncwin){df =mutate(df,y3=bounds23all[i23win2[k],])} 
        newbounds23all[kw,] = findbnds(df2)
      } else{ # Assume only one 23S sequence
        df3 = data.frame(x=EcoliReflocs23,y1=bounds23all[i23win2[k]+2-ncwin,],y2=bounds23all[i23win2[k]+4-ncwin,])
        newbounds23all[kw,] = findbnds(df1)
      }
      jnc = jnc + ncwin
    } else if (ncwin>1) { # Need to deal with 2 or 3 adjacent windows that might have a 23S within them
      df = data.frame(x=EcoliReflocs23,y1=bounds23all[i23win2[k]+1-ncwin,],y2=bounds23all[i23win2[k]+2-ncwin,])
      if (ncwin==3){df =mutate(df,y3=bounds23all[i23win2[k],])}  # It is true that ncwin is no greater than 3
      # Compare ordinary least squares (lsq) fit to the data with a robust lsq fit. If the robust version
      # is clearly better then use it. 
      newbounds23all[kw,] = findbnds(df)
      jnc = jnc + ncwin
    } else {
      jnc = jnc + 1
      newbounds23all[kw,] = bounds23all[jnc,]
    }
    allkdone = (k==length(i23win2))
  }        #    end    k     loop
  newbounds23all = round(newbounds23all[1:kw,])
  if (kw==1){
    zz = list(bounds23=newbounds23all[bseqfwd23], bounds23all = newbounds23all)
  } else {
    zz = list(bounds23=newbounds23all[,bseqfwd23], bounds23all = newbounds23all)
  }
  out = zz
}


refineRR = function(RR,bndsSubunits){
  # Removes rows which have 16S  or 23S lengths outside specified bounds.
  # These bounds are given by bndsSubunits which is a 2x2 matrix giving lower
  # and upper bounds for (row 1) 16S, and (row2) 23S.
  # Removal must be done for both segSeqs and locs, but segSeqs modification
  # is not by row but by single subunit.
  # 19 April 2023                                              [cjw]
  
  #    gseqs = matrix("",nrow=max(2,n16S,n23S),ncol=3)
  for (jg in 1:length(RR)){
    rr = RR[[jg]];  rr0 = rr
    i16S = which(rr$locs[,1]=="16S");    n16S = length(i16S)
    i23S = which(rr$locs[,1]=="23S");    n23S = length(i23S)  
    gseqs = matrix("",nrow=max(2,n16S,n23S),ncol=3)
    gseqs = rr$segSeqs
    id16 = which(sapply(i16S,function(j){
      t1 = abs(as.numeric(rr$locs[j,3]) - as.numeric(rr$locs[j,2]))
      out=(t1>bndsSubunits[1,1]) && (t1<bndsSubunits[1,2])  }  )  )
    noCorrection16 = length(i16S)==length(id16)
    nn16S = length(id16);  i16Sn = 1:nn16S
    if (!noCorrection16){
      cat("Genome ",jg," requires correction of entry ",setdiff(i16S,id16),"\n")
      rr$locs = rr$locs[-c(setdiff(i16S,id16)),]
      rownames(rr$locs) = 1:nrow(rr$locs)
      rr$segSeqs[i16Sn,2] = rr$segSeqs[id16,2]
      rr$segSeqs[(nn16S+1):n16S,2] = "" 
      cat(" 16S correction for genome ",jg,"  Correcting ",setdiff(i16S,id16),"\n")
    }
    # Note that RR[[jg]]$locs might already have had rows removed from the 16S processing just above.
    i23S = which(rr$locs[,1]=="23S");    n23S = length(i23S)  
    id23 = which(sapply(i23S,function(j){
      t1 = abs(as.numeric(rr$locs[j,3]) - as.numeric(rr$locs[j,2]))
      out=(t1>bndsSubunits[2,1]) && (t1<bndsSubunits[2,2])  }  )  )
    noCorrection23 = length(i23S)==length(id23)
    nn23S = length(id23)
    if (!noCorrection23){
      cat("Genome ",jg," requires correction of entry ",setdiff(i23S,i23S[id23]),"\n")
      rr$locs = rr$locs[-c(setdiff(i23S,i23S[id23])),]
      rownames(rr$locs) = 1:nrow(rr$locs)
      id23a = i23S[id23] - nn16S;  
      rr$segSeqs[1:nn23S,3] = rr$segSeqs[id23a,3]
      rr$segSeqs[(nn23S+1):n23S,3] = ""
      cat(" 23S correction for genome ",jg,"  Correcting row ",setdiff(i23S,i23S[id23]),"\n")
    } else {
      id23a = i23S[id23] - nn16S;  
    }
    rr$segSeqs = rr$segSeqs[1:(max(2,nn16S,nn23S)),]
    # It occasionally happens that a pair of 23S locations have the same start location but they are
    # oppositely directed.  It can happen that both are retained by the processing just above, but 
    # one should be rejected.  Here we identify which 16S, if any, is within an ITS range of this
    # start location. Then we choose to retain that 23S that has the same orientation as this 16S, and
    # reject the other 23S.  The locs table has rows ordered by start value, hence duplicated start
    # points are in consecutive rows.
    i23Sn = which(rr$locs[,1]=="23S");    n23Sn = length(i23Sn)
    idup = which(duplicated(rr$locs[i23Sn,2]))
    nu23 = length(unique(rr$locs[i23Sn,2]))
    npairs = n23Sn - nu23
    jp = 0
    ig = rep(-1,npairs)
    while (jp<npairs){
      jp = jp + 1
      x1 = as.numeric(rr$locs[i23Sn[idup[jp]],2])
      diffs = sapply(1:nn16S,function(j){out = min(abs(x1-as.numeric(rr$locs[i16Sn[j],2])),abs(x1-as.numeric(rr$locs[i16Sn[j],3])))})
      ib = which(diffs<1200)
      if (length(ib>1)){
        # Weird!  Just take ib[1]
        ib = ib[1]
      }
      sg16 = sign(as.numeric(rr$locs[i16Sn[ib],3]) - as.numeric(rr$locs[i16Sn[ib],2]))
      sg23_1 =  sign(as.numeric(rr$locs[i23Sn[idup[jp]]-1,3]) - as.numeric(rr$locs[i23Sn[idup[jp]]-1,2]))
      sg23_2 =  sign(as.numeric(rr$locs[i23Sn[idup[jp]],3]) - as.numeric(rr$locs[i23Sn[idup[jp]],2]))
      ig[jp] = which(sg23_1*sg16>0,sg23_2*sg16>0)
      # Note that ig[jp] is either 1 or 2. It gives the index of the entry to be retained - the
      # one to be deleted is 3-ig[jp]
    }
    # So we remove rows ir = i23Sn[idup[jp]]-2+(3-ig[jp]) from rr$locs AND sequences rr$segSeqs[ir-n16Sn,3]
    if (npairs>0){
      for (jp in 1:npairs){
        ir = i23Sn[idup[jp]]-2+ig[jp]
        id = i23Sn[idup[jp]]-2 +(3-ig[jp])   # The locs row to be deleted
        cat(" Further 23S refinement.  Rows ",id,"\n")
        rr$locs = rr$locs[-id,]
        rownames(rr$locs) = 1:nrow(rr$locs)
        id23b = setdiff(i23Sn,id) - nn16S;    nn23S = length(id23b)
        rr$segSeqs[1:length(id23b),3] = rr$segSeqs[id23b,3]
        rr$segSeqs[nn23S+1,3] = ""
      }
      nn23S = nn23S - npairs
      rr$segSeqs = rr$segSeqs[1:(max(2,nn16S,nn23S)),]
    }
    RR[[jg]] = rr
  }
  out = RR
}



get_rrn_into_RR = function(RR,bac_fasta_name,bndsSubunits){
  # From the 16S and 23S data already in RR, determine whether there are rrn operons
  # implied by proximity of 16S and 23S sub-units.  Embed rrn information into
  # RR via the RR[[jg]]$locs component, which becomes a data frame.  Also, save 
  # the rrn sequences in RR[[jg]]$segSeqs[,1].
  #
  # Calls  extract_return_rrn()
  # 21 April 2023                                 [cjw]
  rrnBounds = vector("list",length(RR))
  for (jg in 1:length(RR)){
    if (length(RR[[jg]])==4){
      cat("In  get_rrn_into_RR jg ",jg,"\n")
      n16S = length(which(RR[[jg]]$locs[,1]=="16S"));  cat("    n16S ",n16S,"\n")
      n23S = length(which(RR[[jg]]$locs[,1]=="23S"));  cat("    n23S ",n23S,"\n")
      if ((n16S>0) && (n23S>0)){
        nmaxsops = min(n16S,n23S)
        bnds = matrix(-1,nrow=max(2,nmaxsops),ncol=2)
        seps = matrix(-1,nrow=max(2,n16S),ncol=n23S)
        for (j1 in 1:n16S){
          seps[j1,] = sapply(1:n23S,function(j2){
            df = abs(as.numeric(RR[[jg]]$locs[j1,3]) - as.numeric(RR[[jg]]$locs[n16S+j2,2]))
            dr = abs(as.numeric(RR[[jg]]$locs[j1,2]) - as.numeric(RR[[jg]]$locs[n16S+j2,3]))
            out = min(df,dr)
          })
        }
        if (length(which(rowMins(seps)<1200))<nrow(seps)){ # Value 1200 is used as an upper bound on length(ITS)
          ir = which(rowMins(seps)>=1200)
          cat("Do not consider 16S row ",ir," for forming an rrn.","\n")
          i1 = which(rowMins(seps)<1200);  ni1 = length(i1)
          for (j1 in 1:length(i1)){
            j2 = which.min(seps[i1[j1],])
            j2bnds = c(as.numeric(RR[[jg]]$locs[i1[j1],2:3]),as.numeric(RR[[jg]]$locs[n16S+j2,2:3]))
            ub = max(j2bnds);  lb = min(j2bnds)
            bnds[j1,1] = lb
            bnds[j1,2] = ub
          }
          nmaxsops = min(nmaxsops,nrow(bnds))
        } else {
          cat("Do not consider 23S row ",n16S+which(colMins(seps)>=1200)," for forming an rrn.  This is being handled correctly at present.","\n")
          isrrn = rep(FALSE,n16S)
          nmaxsops = min(n16S,n23S)
          bnds = matrix(-1,nrow=max(2,nmaxsops),ncol=2)
          for (j1 in 1:nmaxsops){
            seps = sapply(1:n23S,function(j2){
              df = abs(as.numeric(RR[[jg]]$locs[j1,3]) - as.numeric(RR[[jg]]$locs[n16S+j2,2]))
              dr = abs(as.numeric(RR[[jg]]$locs[j1,2]) - as.numeric(RR[[jg]]$locs[n16S+j2,3]))
              out = min(df,dr)
            })
            isrrn[j1] = min(seps)<1200 
            if (isrrn[j1]){
              j2 = which.min(seps)
              j2bnds = c(as.numeric(RR[[jg]]$locs[j1,2:3]),as.numeric(RR[[jg]]$locs[n16S+j2,2:3]))
              ub = max(j2bnds);  lb = min(j2bnds)
              bnds[j1,1] = lb
              bnds[j1,2] = ub
            }
          }    #   end    j1    loop
          rrnBounds[[jg]] = bnds
          rrnLen = sapply(1:nrow(rrnBounds[[jg]]),function(jop){abs(rrnBounds[[jg]][jop,2]-rrnBounds[[jg]][jop,1])})
          temprrn = extract_return_rrn(refGenomespath,Dg[jg],bnds)
          kc = 0;   bndsG = matrix(-1,nrow=max(2,nmaxsops),ncol=2)
          for (jjr in 1:length(temprrn)){
            if ((rrnLen[jjr]>bndsSubunits[3,1]) && (rrnLen[jjr]<bndsSubunits[3,2])){
              kc = kc + 1;   bndsG[kc,] = bnds[jjr,]
              RR[[jg]]$segSeqs[jjr,1] =  temprrn[jjr]  
            } else {
              RR[[jg]]$segSeqs[jjr,1] =  ""
            }
          }
          bndsG = bndsG[1:kc,]    
          if (kc>0){
            regID=c(rep("16S",n16S),rep("23S",n23S),rep("rrn",kc))
            if (kc>1){
              start=c(as.numeric(RR[[jg]]$locs[1:(n16S+n23S),2]),bndsG[,1])
              stop=c(as.numeric(RR[[jg]]$locs[1:(n16S+n23S),3]),bndsG[,2])
            } else {
              start=c(as.numeric(RR[[jg]]$locs[1:(n16S+n23S),2]),bndsG[1])
              stop=c(as.numeric(RR[[jg]]$locs[1:(n16S+n23S),3]),bndsG[2])
            }
            newlocs = data.frame(regID=regID,start=start,stop=stop)
            RR[[jg]]$locs = newlocs
          } else {
            cat("\n No rrn found for jg = ",jg, "  genomes ",bac_fasta_name[jg],"\n")
          }
        }
      } else { # No basis for inserting rrn as the 16S/23S part is not complete
        cat("Genome ",jg,bac_fasta_name[jg]," has no identified rrn \n")
      }      #   end   (n16S>0) && (n23S>0) conditional    block
    }        #   end   length(RR[[jg]])==4  conditional    block
  }          #   end    jg    loop
  out = RR
}

make_bac_filenames = function(bacnames_filepath, bacfilename="bac_names.txt"){
  # Generates a text file with name "bac_names_modified.txt" derived from the file "bac_names.txt"
  # in the specified path.  The lines of the new file are valid unix-style filenames constructed from the
  # corresponding lines in the input file.  The purposes of this function is to support generation of
  # individual organism 16S gene, and rrn operon, fasta files for each genome in some database, such as Refseq.
  # It is called,for instance, in   identify_extract_rrnOperon_v01.R  .
  # 7 September                                                     [cjw]
  outfilename = "bac_names_modified.txt"
  L = readLines(con=file.path(bacnames_filepath, bacfilename))
  newL = rep("",length(L))
  replace_set = c(" ","/",":",">",",")
  nbac = length(L)
  for (j in 1:nbac){
    w1 = s2c(L[j])
    for (jc in 2:length(w1)){w1[jc] = ifelse(w1[jc] %in% replace_set,"_",w1[jc])}
    w2 = w1[2:length(w1)]
    w3 = rep(w2[1],length(w2))
    w3[2:length(w3)] = sapply(2:length(w3),function(j){double = (w2[j]=="_") && (w2[j-1]=="_"); ifelse(double,"",w2[j])})
    newL[j] = c2s(w3)
  }
  writeLines(newL,con=file.path(bacnames_filepath, outfilename))
  out = newL
}

pAlign = function(pat,subj,align_type="global-local"){
  # Not used in the code of identify_extract_rrnOperons_v21.R (6 Oct 2022).
  # Computes both forward and reverse-complement pairwiseAlignments of pat and subj,
  # and returns the higher score of the two, together with a logical indicating
  # whether the forward (TRUE) or reverse-complement (FALSE) gave that higher score.
  # 20 June 2022                               [cjw]
  s1 = DNAString(pat);   s2 = DNAString(subj)
  submat = nucleotideSubstitutionMatrix(match = 2, mismatch = 0, baseOnly = FALSE, type = "DNA")
  gapOpening=1; gapExtension=1
  align1 = pairwiseAlignment(s1,s2, type=align_type,substitutionMatrix=submat,
                             gapOpening=1, gapExtension=1)
  sc1 = score(align1)
  align2 = pairwiseAlignment(s1,reverseComplement(s2), type=align_type,substitutionMatrix=submat,
                             gapOpening=1, gapExtension=1)
  sc2 = score(align2)
  orient = sc1>sc2
  if (sc1>sc2){best_align=align1} else {best_align=align2}
  out = list(align=best_align,fwd=orient)
}

fna_Refseq_header_cleanup = function(bac_fasta_name){
  badChar = c(">","|",",","`","'")
  replaceChar = cbind(c(" ",".","/",":"),c("_","p","s","c"))
  new_bac_fasta_name = rep("",length(bac_fasta_name))
  for (j in 1:length(bac_fasta_name)){
    if (!is.na(bac_fasta_name[j])){
      t1 = unlist(strsplit(bac_fasta_name[j],split=" "))
      nt = length(t1)
      for (k in 1:nt){
        x = s2c(t1[k])
        x = sapply(1:length(x),function(k1){out = ifelse(x[k1] %in% badChar,"",x[k1])})
        x = sapply(1:length(x),function(k2){xx=x[k2];
        for (k21 in 1:nrow(replaceChar)){
          if(xx == replaceChar[k21,1]){xx = replaceChar[k21,2]} }
        out = xx })
        t1[k] = c2s(x)
      }
      new_bac_fasta_name[j] =paste(t1,sep="_",collapse="_")
    }
  }
  out = new_bac_fasta_name
}  


extract_return_rrn = function(genomespath,bac_fasta_filename,rrnbnds){
  # NOTE: Called by  get_rrn_into_RR(...) (19 Oct 2022)
  # For genome file.path(genomespath,bac_fasta_filename) with lower and upper bounds of its
  # operons at rrnbnds[1:nops,1:2] returns the sequences of those operons, after replacing 
  # non-canonical bases by randomly sampled canonical bases.
  # Returns a set of genomic locations at which 16S, 23S and rrn operons are found to start, together
  # with the rrn sequences for those starting locations.
  # Input parameters:   genomespath,bac_fasta_filename, rrnbnds)
  # Derived from extract_return_16S23Srrn() of 29 September 2022 which could process 
  # genomes with more than 1 chromosome, and dealt separately with 16S and 23S.
  
  # Read the fasta file specified, check whether it has more than one sequence in it, and if so select the 
  # the longest of these sequences.  Also, ensure that only ACGT entries are present - replacing 
  # non-canonical bases by a random sample of the canonical bases.   The output from this section is  bacseq.
  locsRegs = rep(0,6)
  jg = which(Dg==bac_fasta_filename)
  gF = Dg[jg]
  gen_fa = read.fasta(file=file.path(genomespath,gF), seqtype="DNA",forceDNAtolower = FALSE)
  chrom = which(unlist(sapply(1:length(gen_fa),function(j){grep("chromosome",attributes(gen_fa[[j]])$Annot)==1})))
  if (length(chrom)==0){chrom=1}
  nchr = max(1,length(chrom))
  if (nchr==1){
    gen_fa[[1]] = gen_fa[[chrom]]
    chrom.Len = length(gen_fa[[1]])
  } else if (nchr>0) {
    cat("\n Genome ",gF, " has ",nchr, "chromosomes identified.\n")
    chrom.Len = rep(0,nchr)
    gen_fa[[1]] = gen_fa[[chrom[1]]]
    chrom.Len = sapply(1:nchr,function(jc){out = length(gen_fa[[chrom[jc]]])})
    #    gen_fa[[1]] = sapply(2:nchr,function(jc){t1 = c(gen_fa[[1]],gen_fa[[chrom[jc]]])})
    gf = sapply(2:nchr,function(jc){t1 = c(gen_fa[[1]],gen_fa[[chrom[jc]]])})
  }
  # Some of the genomes have a few "N" entries.Replace non-canonical base by "A" 
  # and advise number of replacements.
  if (nchr==1){
    noncanon = which(!(gen_fa[[1]] %in% c("A","G","C","T") ))
    nc_codes = gen_fa[[1]][noncanon]
    if (length(noncanon)>0){
      gen_fa[[1]][noncanon] = rep(sample(c("A","C","G","T"),1),length(noncanon),replace=TRUE)
    }
    bacseq = c2s(gen_fa[[1]])
  } 
  if (nchr>1){
    noncanon = which(!(gf %in% c("A","G","C","T") ))
    nc_codes = gf[noncanon]
    if (length(noncanon)>0){
      gf[noncanon] = rep(sample(c("A","C","G","T"),1),length(noncanon),replace=TRUE)
    }
    bacseq = c2s(gf)
  }
  extend = 40
  nops = nrow(rrnbnds)
  rrnSeqs = rep("",nops)
  for (jop in 1:nops){
    rrnSeqs[jop] = substr(bacseq,start=max(1,rrnbnds[jop,1]-extend), stop=min(rrnbnds[jop,2]+extend,sum(chrom.Len)))
  }
  out = rrnSeqs
}


row_subset_correction = function(bnds){
  # Given a matrix, bnds, whose rows are locations of a common set of operon boundaries
  # this code identifies whether any boundary location is incongruent and, if so,
  # adjusts based on a consensus of such boundaries. Called in extract_16S23S_singleFolder()
  # for both 16S and 23S (22 April 2023).
  # 12 October 2022                                                      [cjw]
  make_pre_plot = FALSE
  nrb = nrow(bnds);   ncb = ncol(bnds)
  absdiffs = t(sapply(1:nrb,function(jg){
    sapply(1:(ncb-1),function(jc){
      abs(bnds[jg,jc+1]-bnds[jg,jc])
    })
  }))
  medabsdiffs = sapply(1:(ncb-1),function(jc){median(absdiffs[,jc])})
  # Now identify boundary locations that might not be reasonable. 
  # The approach is to use the fact that absdiffs values that are close 
  # to the median for the separation of consecutive boundaries are almost
  # certainly good locations.  So identify good boundaries.  Then assume 
  # any other boundaries are not good.  These are the ones that are then 
  # to be adjusted.
  goodBnds = matrix(-1,nrow=nrb,ncol=ncb)
  for (jg in 1:nrb){
    t1=-1
    for (jc in 1:(ncb-1)){
      if (abs(absdiffs[jg,jc]-medabsdiffs[jc])<100){t1 = append(t1,c(jc,jc+1))}
    }
    t1 = unique(t1[2:length(t1)]) 
    goodBnds[jg,1:length(t1)] = t1
  }
  badBnds = matrix(goodBnds<0,nrow=nrow(goodBnds),ncol=ncol(goodBnds))
  goodBnds = ifelse(badBnds,0,1)
  strand = rep(0,nrb)
  for (jg in 1:nrb){
    ig = which(goodBnds[jg,]==1)
    nig = length(ig)
    strand[jg] = sign(median(bnds[jg,ig[(nig-1)]:ig[nig]]) - median(bnds[jg,ig[1:2]]))
  }
  # Now identify a suitable centering column.
  # Create a modified set of bounds to process 
  orbnds = t(sapply(1:nrb,function(jg){out=strand[jg]*bnds[jg,]}))
  cSds = rep(0, ncb)
  csum = rep(-1,ncb)
  for (jbc in 1:ncb){
    bndsc = orbnds
    for (jg in 1:nrb){
      bndsc[jg,] = sapply(1:ncb,function(jb){orbnds[jg,jb]-orbnds[jg,jbc]})
    }
    cSds = colSds(bndsc)
    cat("Centering on column ",jbc,"\n")
    print(cSds)
    csum[jbc] = sum(cSds)
    cat("Sum of column standard deviations: ",csum[jbc],"\n")
  }
  jbc = which.min(csum)
  cat("Centring using column  ",jbc,"\n")
  # Now centre each genome's boundaries.
  bndsc = orbnds
  cadj = rep(0,nrb)
  for (jg in 1:nrb){
    bndsc[jg,] = sapply(1:ncb,function(jb){orbnds[jg,jb]-orbnds[jg,jbc]})
  }
  
  refbnds = sapply(1:ncb,function(jb){median(bndsc[,jb])})
  # Now find adjustments for those boundaries that are not identified as good.
  br = which(sapply(1:nrb,function(jg){sum(goodBnds[jg,])<ncol(goodBnds)}))
  newbndsc = bndsc
  for (jg in br){
    ig = which(goodBnds[jg,]==1)
    nig = length(ig)
    fit = lm(bndsc[jg,ig] ~ refbnds[ig])
    ib = setdiff(1:ncb,ig)
    for (jc in ib){
      newbndsc[jg,jc] = fit$coefficients[2]*refbnds[jc] + fit$coefficients[1]
    }
  }
  # Now decentre all boundary locations.
  newbnds = bnds
  print(dim(newbndsc))
  print(dim(orbnds))
  for (jg in 1:nrb){
    newbnds[jg,] = strand[jg]*(newbndsc[jg,] + orbnds[jg,jbc])
  }
  out = newbnds  
}


#######################################################################################
#######################################################################################
##########################                               ##############################
##########################           MAIN BODY           ##############################
##########################                               ##############################
#######################################################################################
#######################################################################################
slurm =   TRUE   
if (slurm){
  # Deal with input parameters.
  # Test if there is exactly 1 argument: if not, return an error; if so extract key parameters
  if (!(length(args)==1)) {
    stop("One parameter is required,  numcores.", call.=FALSE)
  } else { 
    numcores = as.numeric(args[1])
  } 
} else {
  numcores = 10
}

# Select database for which rrn regions are to be extracted.
DBid = "workDB1"
message1 = paste("Commenced run of identify_extract_rrnOperons_workDB_v01.R for ",
                   DBid, ", with",numcores,"cores.   Date: ", date(),sep=" ")
print(message1)

# Various  path  definitions
#  /vast/projects/rrn/microbiome/blastdb/workDB/genomes
#basepath = "/stornext/Bioinf/data/lab_speed/cjw/microbiome/workDB1/D6322.refseq/Genomes"
basepath = "/stornext/Bioinf/data/lab_speed/cjw/microbiome"
#refGenomespath = file.path(basepath,DBid,"D6322.refseq/Genomes")
vastbasepath = "/vast/projects/rrn/microbiome/blastdb"
refGenomespath = file.path(vastbasepath,DBid,"genomes")
outpath = file.path(vastbasepath,DBid,"output")
refopspath = file.path(outpath,"rrn")
outRDatapath = file.path(outpath,"RData")
plotpath = file.path(outpath,"plots")
bounds_path = file.path(outpath,"boundaries")

m = as.numeric(system2("ls",args=paste(" -lrt ",refGenomespath," | wc -l",sep=""),stdout=TRUE))
Ng = m-2
jgLow = 1;   jgHi = Ng
print(paste("jgLow: ",jgLow,"   jgHi: ",jgHi,sep=""))
generate_rrn = TRUE
generate_rrn_operons_16S = TRUE
generate_rrn_operons_23S = TRUE
generate_rrn_operons_ITS = FALSE
create_bac_names = TRUE
full23S = TRUE
small16S_primer_sets = c(1,2,6,8,10,12,14,16) # 16S_5p, 16S_3p, V1_3p,  V2_3p, V3_3p, V4_3p, V5_3p, V6_3p (27,1492,120,357,538,798,922,1080)
bseqfwd16 = c(1,4,5,6,8,2)   # NB: First bound must be 16S_5', and last bound 16S_3'.   so (16S_5',V2_3',V3_3',V4_3',V6_3',16S_3')
EcoliReflocs16 = c(27,1492,120,357,538,798,922,1080)  #  Corresponding to small16S_primer_sets above.
if (full23S){
  small23S_primer_sets = c(3,4,23,25,27,29,33,35)  # 23S_5', 23S_3', Z1_5', Z2_5', Z3_5', Z4_5', Z6_5', Z7_5' (Note: Z<n>_5'==Z<n-1>_3')
  bseqfwd23 = c(1,3,5,7,8,2)
  bndsSubunits = rbind(c(1200,1600),c(2200,2800),c(4000,6000))
  EcoliReflocs23 = c(10,2488,84,136,385,728,1591,1725)  #  Corresponding to small23S_primer_sets above.
} else {
  small23S_primer_sets = c(23,25,27,29,33,35,37,39)   # Z1_5', Z2_5', Z3_5', Z4_5', Z6_5', Z7_5', Z8_5', Z9_5' (Note: Z<n>_5'==Z<n-1>_3')
  bseqfwd23 = c(1,3,4,5,7,8)
  bndsSubunits = rbind(c(1200,1600),c(1800,2400),c(3650,5650))
  EcoliReflocs23 = c(84,136,385,728,1591,1725,1821,2145)  #  Corresponding to small23S_primer_sets above.
}
colnames(bndsSubunits) = c("lowb","uppb");   rownames(bndsSubunits) = c("16S","23S","rrn")
nb16 = length(small16S_primer_sets)
nb23 = length(small23S_primer_sets)
pattern_str =  ".fna"   #  "fasta"     
Dg0 = dir(refGenomespath,pattern=pattern_str)
ind_bacteria = which(sapply(1:length(Dg0),function(jg){
  z1 = file.size(file.path(refGenomespath,Dg0[jg])) 
  out = ((z1>130000) & (z1<12000000))
}))
jg_highest = length(ind_bacteria)
Dg = Dg0[ind_bacteria]
jgHi = min(jgHi, length(Dg))
# Can access the fasta file header for each file as follows:-
bac_fasta_name = rep("",length(ind_bacteria))
#setwd(refopspath)
namesMap = matrix("",nrow=length(Dg),ncol=2)
# namesMap links names of the .fna files to the names of the sequences in the file.
# Note that these sequence names become the entries in bac_names.txt, and these 
# are subsequently cleaned up - see fna_Refseq_header_cleanup(), used below -  to 
# become names of individual 16S, 23S and rrn .fasta files that are saved.  
if (create_bac_names){
  system2("rm",args=file.path(outpath,"text","bac_names.txt"))
  for (jb in 1:length(ind_bacteria)){    # length(Dg)
    argstr = paste(" -1 ",file.path(refGenomespath,Dg[jb]), " >> ",
                   file.path(outpath,"text","bac_names.txt"),sep="")
    system2("head",args= argstr)
    argstr1 = paste(" -1 ",file.path(refGenomespath,Dg[jb]),sep="")
    zz = system2("head",args= argstr1,stdout=TRUE)  # Get the first line of the .fna file
    namesMap[jb,1] = Dg[jb]
    namesMap[jb,2] = substr(zz,start=2,stop=nchar(zz))
  }
}
bac_fasta_name = readLines(con=file.path(outpath,"text","bac_names.txt"))
new_bac_fasta_name  = fna_Refseq_header_cleanup(bac_fasta_name)
writeLines(new_bac_fasta_name, con=file.path(outpath,"text","bac_names_cleaned.txt"))
bac_fasta_name = new_bac_fasta_name[1:(jgHi-jgLow+1)]

tempExclude = FALSE
if (!tempExclude){
#  RR = vector("list",length(ind_bacteria))
  #for (jg in 1:length(ind_bacteria)){   # length(ind_bacteria)
  cat("Commencing identification of the 16S and 23S subunits \n\n")
#  for (jg in 1:length(Dg)){  
#    RR[[jg]] = extract_16S23S_singlespecies(Dg[jg])
#  }
  RR = mclapply(Dg,extract_16S23S_singleFolder, mc.cores = numcores)
#  RR = lapply(Dg,extract_16S23S_singleFolder)
  cat("Completed identification of the 16S and 23S subunits \n\n")
  outname  = paste("ribosomal_subunits_",DBid,".RData",sep="")
  RR0 = RR
  # The following is an interim save, because of the computational load of the 
  # computation of RR - which would be lost in a subsequent code failure.
  outname2  = paste("ribosomal_subunits_interim_",DBid,".RData",sep="")
  save(file=file.path(outRDatapath,outname2),RR0, bndsSubunits)
#  RR = refineRR(RR,bndsSubunits)     # I 
  # Note that details on rrn operons need to be generated based on the 16S and 23S 
  # data returned for each genome.  The following code achieves this.
  # Need to compare start/end coordinates of each 16S with start/end coordinates
  # of all 23S to see if any pairs are sufficiently close to be from a single rrn.
  RR = get_rrn_into_RR(RR,Dg,bndsSubunits)
  outname  = paste("ribosomal_subunits_",DBid,".RData",sep="")
  save(file=file.path(outRDatapath,outname),RR, bndsSubunits,small16S_primer_sets,bseqfwd16,
               small23S_primer_sets,bseqfwd23)
} else {
  outname  = paste("ribosomal_subunits_",DBid,".RData",sep="")
  load(file=file.path(outRDatapath,outname))
  cat("\n Completed tempExclude block.\n")
}
## iibad = which(sapply(1:length(RR),function(j){length(RR[[j]])<2}))
# For 16S, 23S and rrn, identify the indices of RR for which such regions
# are found on the genome.
cat("\n\n Now processing to use the RR data to identify and extract 16S and 23S sequences. \n")
Ng = length(RR)
indp = rep(0,Ng)
indpOp = rep(0,Ng)
indp16S = rep(0,Ng)
indp23S = rep(0,Ng)
ip = 0;  ipOp = 0;  ip16S = 0;   ip23S = 0
for (jg in 1:Ng){
  if (length(RR[[jg]])>1){
    n16S = length(which((RR[[jg]]$locs[,1]=="16S")))
    n23S = length(which((RR[[jg]]$locs[,1]=="23S")))
    nops = length(which((RR[[jg]]$locs[,1]=="rrn")))
    present_rrn = nops>0
    present_16S = n16S>0
    present_23S = n23S>0
    if (!present_rrn){
      cat("No rrn operon present in genome ",jg,"\n")
      cat("   - number of 16S is ",n16S,", of 23S ",n23S,"\n")
    } 
    present = (present_rrn || (present_16S || present_23S))
    if (present){
      ip = ip+1
      indp[ip] = jg
      if (present_rrn){
        ipOp = ipOp+1
        indpOp[ipOp] = jg
      }
      if (present_16S){
        ip16S = ip16S+1
        indp16S[ip16S] = jg
      }
      if (present_23S){
        ip23S = ip23S+1
        indp23S[ip23S] = jg
      }
    }
  }
}
indp = indp[1:ip]
indpOp = indpOp[1:ipOp]
indp16S = indp16S[1:ip16S]
indp23S = indp23S[1:ip23S]
cat("\n ip, ipOp, ip16S, ip23S: ",ip,ipOp,ip16S,ip23S,"\n")

# Now write all rrn operon sequences from one genome to a single file, all 16S to a separate single 
# file, and all 23S sequences to another separate file.
# Also, assign values to variable  num_ops_per_species, num_16S_per_species, num_23S_per_species  .
# Note that indp indexes any genomes that have at least one of an rrn operon or a 16S or a 23S.
num_ops_per_species = rep(0,Ng)
num_16S_per_species = rep(0,Ng)
num_23S_per_species = rep(0,Ng)
for (jjg in 1:length(indp)){
  open_typeOp = "w"
  open_type16S = "w"
  open_type23S = "w"
  jg = indp[jjg]
  base_bac_name = unlist(strsplit(new_bac_fasta_name[jg],split=","))[1]
  base_bac_name4 = paste(unlist(strsplit(base_bac_name,split="_"))[1:4],sep="_",collapse="_")
  if (length(RR[[jg]])>1){
    if (length(RR[[jg]]$segSeqs)==3){ # Is failure in RR[[jg]] computation, or a dummy entry, or genuine single operon genome
      if (nchar(RR[[jg]]$segSeqs[1])<3000){
        cat("Genome index ",jg," has no operons identified \n")
        num_ops_per_species[jg] = 0
      } else {
        if (!RR[[jg]]$segSeqs[1]==""){ # There is a single rrn operon
          num_ops_per_species[jg] = 1
          ops_fasta_name = paste(base_bac_name4,"_rrn_operons_1_",out_dateString,".fa",sep="")
          write.fasta(RR[[jg]]$segSeqs[1],names=paste(base_bac_name,"rrn_operon_1",sep=" "),
                      file.out=file.path(refopspath,ops_fasta_name),open=open_typeOp)
          open_typeOp = "a"
        }
      }
      # There may be 16S or 23S sequences.  These need to be identified and processed here.
      if (!RR[[jg]]$segSeqs[2]==""){
        num_16S_per_species[jg] = 1
        ssu_fasta_name = paste(base_bac_name4,"_16S_genes_1_",out_dateString, ".fa",sep="")
        write.fasta(RR[[jg]]$segSeqs[2],names=paste(base_bac_name,"16S_gene_1",sep=" "),
                    file.out=file.path(refopspath,ssu_fasta_name),open=open_type16S)
        open_type16S = "a"
      }
      if (!RR[[jg]]$segSeqs[3]==""){
        num_23S_per_species[jg] = 1
        lsu_fasta_name = paste(base_bac_name4,"_23S_genes_1_",out_dateString, ".fa",sep="")
        write.fasta(RR[[jg]]$segSeqs[3],names=paste(base_bac_name,"23S_gene_1",sep=" "),
                    file.out=file.path(refopspath,lsu_fasta_name),open=open_type23S)
        open_type23S = "a"
      }
    } else { # There are 2 or more operons for this genome
      n16S = length(which((RR[[jg]]$locs[,1]=="16S")))
      n23S = length(which((RR[[jg]]$locs[,1]=="23S")))
      nops = length(which((RR[[jg]]$locs[,1]=="rrn")))
      num_ops_per_species[jg] = 0
      num_16S_per_species[jg] = 0
      num_23S_per_species[jg] = 0
      if (nops>0){
        for (jop in 1:nops){
          len.rrn = nchar(RR[[jg]]$segSeqs[jop,1])
          if ((len.rrn>bndsSubunits[3,1]) && (len.rrn<bndsSubunits[3,2])){
            ops_fasta_name = paste(base_bac_name4,"_rrn_operons_",nops,"_",out_dateString, ".fa",sep="")
            write.fasta(RR[[jg]]$segSeqs[jop,1],names=paste(base_bac_name,"rrn_operon_",jop,sep=""),
                        file.out=file.path(refopspath,ops_fasta_name),open=open_typeOp)
            num_ops_per_species[jg] =  num_ops_per_species[jg] + 1
            open_typeOp = "a"
          }
        }
      }
      if (n16S>0){
        for (j16 in 1:n16S){
          len.16S = nchar(RR[[jg]]$segSeqs[j16,2])
          if ((len.16S>bndsSubunits[1,1]) && (len.16S<bndsSubunits[1,2])){
            ssu_fasta_name = paste(base_bac_name4,"_16S_genes_",n16S,"_",out_dateString, ".fa",sep="")
            write.fasta(RR[[jg]]$segSeqs[j16,2],names=paste(base_bac_name,"16S_gene_",j16,sep=""),
                        file.out=file.path(refopspath,ssu_fasta_name),open=open_type16S)
            num_16S_per_species[jg] =  num_16S_per_species[jg] + 1
            open_type16S = "a"
          }
        }
      }
      if (n23S>0){
        for (j23 in 1:n23S){
          len.23S = nchar(RR[[jg]]$segSeqs[j23,3])
          if ((len.23S>bndsSubunits[2,1]) && (len.23S<bndsSubunits[2,2])){
            lsu_fasta_name = paste(base_bac_name4,"_23S_genes_",n23S,"_",out_dateString, ".fa",sep="")
            write.fasta(RR[[jg]]$segSeqs[j23,3],names=paste(base_bac_name,"23S_gene_",j23,sep=""),
                        file.out=file.path(refopspath,lsu_fasta_name),open=open_type23S)
            num_23S_per_species[jg] =  num_23S_per_species[jg] + 1
            open_type23S = "a"
          }
        }
      }
    }          #       end    conditional on length(RR[[jg]]$segSeqs)     block
  } else {     #       end    conditional on length(RR[[jg]])>1           block
    num_ops_per_species[jg] = 0
    num_16S_per_species[jg] = 0
    num_23S_per_species[jg] = 0
  }            #       end    conditional on length(RR[[jg]])>1           block
}
lbo = c(1,1+cumsum(num_ops_per_species)[1:(length(num_ops_per_species)-1)])  # Lower index bounds for operons associated with each genome.
ubo = cumsum(num_ops_per_species)        # Upper index bounds for operons associated with each genome.
xbo = 1:ubo[length(ubo)]
IP = list(indp=indp, indpOp=indpOp, indp16S=indp16S, indp23S=indp23S)
NperReg = list(nops=num_ops_per_species, n16S=num_16S_per_species, n23S=num_23S_per_species)
if (full23S){
  outname = paste("extraction_data_rrn_operon_",DBid,"_genomes_",out_dateString,".RData",sep="")
} else {
  outname = paste("extraction_data_rrn_operon_",DBid,"_genomes_modified23S_",out_dateString,".RData",sep="")
}
save(file=file.path(outRDatapath,outname),RR,NperReg,lbo,ubo,xbo,IP,new_bac_fasta_name,namesMap)
cat("Primary output to ",file.path(outRDatapath,outname),"\n")

# If only wanting to reload data from a former run execute the following:-
#  inname = paste("extraction_data_rrn_operon_",microbiomeID,"_genomes_",in_dateString,".RData",sep="")
#  load(file=file.path(outRDatapath,inname))
#  OR  load(file=file.path(outRDatapath,outname))


if (slurm){quit(save = "no", status = 0, runLast = TRUE)}
quit(save = "no", status = 0, runLast = TRUE)

# 26 April 2023:  The following code is investigative and developmental.  It is part of a close
#                 analysis of why I appear to be getting some unsound 16S bounds results.
#                 The resolution being attempted is to use scanning-window adjacency and robust
#                 linear regression to identify the best sets of 16S bounds.
#                 The code should eventually sit inside function extract_16S23S_singleFolder(),
#                 but external to   find_start_16S_23S_multi()
print(cbind(NperReg[[2]],NperReg[[3]],NperReg[[1]]))
i12 = which(sapply(1:length(NperReg[[2]]),function(j){out = !(NperReg[[2]][j]==NperReg[[3]][j])}))
i13 = which(sapply(1:length(NperReg[[2]]),function(j){out = !(NperReg[[2]][j]==NperReg[[1]][j])}))
i23 = which(sapply(1:length(NperReg[[2]]),function(j){out = !(NperReg[[3]][j]==NperReg[[1]][j])}))
print(i12)
print(i13)
print(i23)

opwinlen16=3000; ovlp16=1600; opwinlen23=4500; ovlp23=3000
jg = 1
startloc0 = 4035042   # 210422  
jw0 = ((startloc0-ovlp16)/(opwinlen16-ovlp16) %/% 1) -2
print(cbind(RR[[jg]]$stats16[jw0,],RR[[jg]]$stats16[jw0+1,],RR[[jg]]$stats16[jw0+2,],RR[[jg]]$stats16[jw0+3,],
                RR[[jg]]$stats16[jw0+4,],RR[[jg]]$stats16[jw0+5,]))
jw01 = 2130
print(cbind(align_stats16[jw01,],align_stats16[jw01+1,],align_stats16[jw01+2,],align_stats16[jw01+3,],align_stats16[jw01+4,]))
startloc1 = 216087
jw1 = ((startloc1-ovlp23)/(opwinlen23-ovlp23) %/% 1) - 2
print(cbind(RR[[jg]]$stats23[jw1,],RR[[jg]]$stats23[jw1+1,],RR[[jg]]$stats23[jw1+2,],RR[[jg]]$stats23[jw1+3,],
            RR[[jg]]$stats23[jw1+4,],RR[[jg]]$stats23[jw1+5,]))

plot(EcoliReflocs,bounds16all[1,]-min(bounds16all[1,]),ylab="Relative Location",xlab="Boundary indices",pch=1,col=1,
        ylim = c(-100,3000))
for (j in 2:3){
  points(EcoliReflocs,bounds16all[j,]-min(bounds16all[j,]),pch=j,col=j)
}
legend("left",legend=c("1","2","3"),col=1:3,pch=1:3)

# Will now try to use robust linear regression, preceded by identification of 
# adjacent windows from i16Swin, to identify unique 16S regions.
linfitBnds16 = function(bounds16all,i16Swin){
  di16win = sapply(2:length(i16Swin),function(j){i16Swin[j] - i16Swin[j-1]})
  i16win2 = which(di16win>1)  # Find the separators of clustered windows.
  i16win2 = append(i16win2,length(i16Swin))
  ntb16 = length(i16win2)
  newbounds16all = matrix(0,nrow=ntb16,ncol=ncol(bounds16all))
  jnc = 0
  for (k in 1:ntb16){
    if (k==ntb16){
      ncwin = length(i16Swin) - i16win2[k-1]
    } else {
      ncwin = ifelse(k==1,i16win2[1], i16win2[k]-i16win2[k-1] +1 )
    }
    if (ncwin>1){ # Need to deal with 2 or 3 adjacent windows that might have a 16S within them
      df = data.frame(x=EcoliReflocs,y1=bounds16all[i16win2[k]+1-ncwin,],y2=bounds16all[i16win2[k]+2-ncwin,])
      if (ncwin==3){df =mutate(df,y3=bounds16all[i16win2[k],])}  # It is assumed that ncwin is no greater than 3
      # Compare ordinary least squares (lsq) fit to the data with a robust lsq fit. If the robust version
      # is clearly better then use it. 
      ols = vector("list",ncwin);   robust = ols
      ols[[1]] = lm(y1~x, data=df);   ols[[2]] = lm(y2~x, data=df)   
      robust[[1]] = rlm(y1~x, data=df);   robust[[2]] = rlm(y2~x, data=df)
      fitsum = matrix(0,nrow=ncwin,ncol=2)
      fitsum[1,] = c(summary(ols[[1]])$sigma,summary(robust[[1]])$sigma)
      fitsum[2,] = c(summary(ols[[2]])$sigma,summary(robust[[2]])$sigma) 
      if (ncwin==3){
        ols[[3]] = lm(y3~x, data=df);   robust[[3]] = rlm(y3~x, data=df)
        fitsum[3,] = c(summary(ols[[3]])$sigma,summary(robust[[3]])$sigma)
      }
      jbest = which.min(fitsum);  jcol = 1+(jbest-1)%/%ncwin;  jrow =  1+(jbest-1)%/%2
      if (jcol==1){
        newbounds16all[k,] = ols[[jrow]]$fitted.values
      } else {
        newbounds16all[k,] = robust[[jrow]]$fitted.values
      }
      jnc = jnc + ncwin
    } else {
      jnc = jnc + 1
      newbounds16all[k,] = bounds16all[jnc,]
    }
  }        #    end    k     loop
  newbounds16all = round(newbounds16all)
  out = list(bounds16=newbounds16all[bseqfwd16], bounds16all = newbounds16all)
}
  
df = data.frame(x=EcoliReflocs,y1=bounds16all[1,],y2=bounds16all[2,],y3=bounds16all[3,])
ols1 = lm(y1~x, data=df);   ols2 = lm(y2~x, data=df);   ols3 = lm(y3~x, data=df)  
robust1 = rlm(y1~x, data=df);   robust2 = rlm(y2~x, data=df);   robust3 = rlm(y3~x, data=df)
print(rbind(c(summary(ols1)$sigma,summary(robust1)$sigma),c(summary(ols2)$sigma,summary(robust2)$sigma),
          c(summary(ols3)$sigma,summary(robust3)$sigma)) )


plot(df$y3, ols3$residuals, ylab='Standardized Residuals', xlab='y', pch="o",
        ylim=c(min(ols3$residuals,robust3$residuals),max(ols3$residuals,robust3$residuals))) 
abline(h=0)
points(df$y3, robust3$residuals, pch="+")
################################################################################################
################################################################################################
################################################################################################
################################################################################################
# TASK 1: Create a stand-alone function that takes a matrix of boundary locations, where
#         each row is one genome, and iteratively corrects a subset of the rows based on the 
#         complementary subset of rows.  function iterative_row_subset_correction()
#
# Test data from "Actinobacter baumannii", derived from RR[[1]].
align_stats23 = RR[[2]]$stats23
nwin23 = nrow(align_stats23)
lbf23 = length(bseqfwd23)
diffs = matrix(0,nrow=nwin23,ncol=(lbf23-1))
for (jw in 1:nwin23){
   diffs[jw,] = sapply(2:lbf23,function(js){align_stats23[jw,bseqfwd23[js]]-align_stats23[jw,bseqfwd23[js-1]]})
}
i23Swin = which(sapply(1:nwin23,function(jw){
                              b1 = median(align_stats23[jw,9:16])>8
                              b2 = sum(align_stats23[jw,nb23+bseqfwd23])>(lbf23*8)
                              out = b1 && b2
                           }))
print(diffs[i23Swin,])  
# This suggests that rows c(1,3,5,6,8,9,10,11) be used for initial development purposes.

diffst = diffs[i23Swin[c(1,3,5,6,8,9,10,11)],]
bnds23t = align_stats23[i23Swin[c(1,3,5,6,8,9,10,11)],bseqfwd23]
# Bring each genome into an approximate alignment by centering using the 3rd boundary.
# (Note: Several boundaries should be tried as the centering one and that which gives 
#        the least "sum(column_variance)" then used.)

row_subset_correction = function(bnds){
  # Given a matrix, bnds, whose rows are locations of a common set of operon boundaries
  # this code identifies whether any boundary location is incongruent and, if so,
  # adjusts based on a consensus of such boundaries.
  # 11 October 2022                                                      [cjw]
  make_pre_plot = FALSE
  nrb = nrow(bnds);   ncb = ncol(bnds)
  # Must first ensure all bound sets correspond to the same orientation.
  # Use sign of difference in medians of last 3 and first 3 boundaries.
  strand = sapply(1:nrb,function(jg){
    out = sign(median(bnds[jg,(ncb-1):ncb]) - median(bnds[jg,1:2]))
  })
  # Create a modified set of bounds to process 
  orbnds = t(sapply(1:nrb,function(jg){out=strand[jg]*bnds[jg,]}))
  cSds = matrix(0,nrow=nrb,ncol=ncb)
  csum = rep(-1,ncb)
  for (jbc in 1:ncb){
    bndsn = orbnds
    for (jg in 1:nrb){
      bndsn[jg,] = sapply(1:ncb,function(jb){orbnds[jg,jb]-orbnds[jg,jbc]})
    }
    cSds[jbc,] = colSds(bndsn)
    cat("Centering on column ",jbc,"\n")
    print(cSds[jbc,])
    csum[jbc] = sum(cSds[jbc,])
    cat("Sum of column standard deviations: ",csum[jbc],"\n")
  }
  jbc = which.min(csum)
  bndsn = orbnds
  for (jg in 1:nrb){
    bndsn[jg,] = sapply(1:ncb,function(jb){orbnds[jg,jb]-orbnds[jg,jbc]})
  }
  
  refbnds = sapply(1:ncb,function(jb){median(bndsn[,jb])})
  if (make_pre_plot){
    plot(refbnds,bndsn[1,],pch=1,col=1,xlab="Ref Loc",ylab="Actual Loc",main=paste(microbiomeID," jg=1 23S",sep=""),
         ylim=c(min(bndsn),max(bndsn)))
    for (jg in 2:nrow(bndsn)){points(refbnds,bndsn[jg,],pch=jg,col=jg)}
  }
  
  newbnds = orbnds
  adjust = matrix(0,nrow=nrb,ncol=ncb)
  for (jg in 1:nrb){
    rfit = rlm(bndsn[jg,] ~ refbnds)
    ib = as.numeric(which(abs(rfit$residuals)>100))
    if (length(ib)>0){
      adjust[jg,ib] = as.numeric(rfit$residuals[ib])
      newbnds[jg,ib] = strand[jg]*orbnds[jg,ib] - adjust[jg,ib]
    }
  }
  out = newbnds
}

bnds23t.adjusted = row_subset_correction(bnds23t)
# Check the adjustment via the following plot.
bnds23tn = bnds23t.adjusted
for (jg in 1:nrow(bnds23t.adjusted)){
  bnds23tn[jg,] = sapply(1:ncol(bnds23t.adjusted),function(jb){bnds23t.adjusted[jg,jb]-bnds23t.adjusted[jg,jbc]})
}

plot(refbnds,bnds23tn[1,],pch=1,col=1,xlab="Ref Loc",ylab="Actual Loc",main=paste(microbiomeID," jg=1 23S",sep=""),
     ylim=c(min(bnds23tn),max(bnds23tn)))
for (jg in 2:nrow(bnds23tn)){points(refbnds,bnds23tn[jg,],pch=jg,col=jg)}


# Another attempt!
row_subset_correction = function(bnds){
  # Given a matrix, bnds, whose rows are locations of a common set of operon boundaries
  # this code identifies whether any boundary location is incongruent and, if so,
  # adjusts based on a consensus of such boundaries.
  # 12 October 2022                                                      [cjw]
  make_pre_plot = FALSE
  nrb = nrow(bnds);   ncb = ncol(bnds)
  absdiffs = t(sapply(1:nrb,function(jg){
    sapply(1:(ncb-1),function(jc){
      abs(bnds[jg,jc+1]-bnds[jg,jc])
    })
  }))
  medabsdiffs = sapply(1:(ncb-1),function(jc){median(absdiffs[,jc])})
  # Now identify boundary locations that might not be reasonable. 
  # The approach is to use the fact that absdiffs values that are close 
  # to the median for the separation of consecutive boundaries are almost
  # certainly good locations.  So identify good boundaries.  Then assume 
  # any other boundaries are not good.  These are the ones that are then 
  # to be adjusted.
  goodBnds = matrix(-1,nrow=nrb,ncol=ncb)
  for (jg in 1:nrb){
    t1=-1
    for (jc in 1:(ncb-1)){
      if (abs(absdiffs[jg,jc]-medabsdiffs[jc])<100){t1 = append(t1,c(jc,jc+1))}
    }
    t1 = unique(t1[2:length(t1)]) 
    goodBnds[jg,1:length(t1)] = t1
  }
  badBnds = matrix(goodBnds<0,nrow=nrow(goodBnds),ncol=ncol(goodBnds))
  goodBnds = ifelse(badBnds,0,1)
  strand = rep(0,nrb)
  for (jg in 1:nrb){
    ig = which(goodBnds[jg,]==1)
    nig = length(ig)
    strand[jg] = sign(median(bnds[jg,ig[(nig-1)]:ig[nig]]) - median(bnds[jg,ig[1:2]]))
  }
  # Now identify a suitable centering column.
  # Create a modified set of bounds to process 
  orbnds = t(sapply(1:nrb,function(jg){out=strand[jg]*bnds[jg,]}))
  cSds = matrix(0,nrow=nrb,ncol=ncb)
  csum = rep(-1,ncb)
  for (jbc in 1:ncb){
    bndsc = orbnds
    for (jg in 1:nrb){
      bndsc[jg,] = sapply(1:ncb,function(jb){orbnds[jg,jb]-orbnds[jg,jbc]})
    }
    cSds[jbc,] = colSds(bndsc)
    cat("Centering on column ",jbc,"\n")
    print(cSds[jbc,])
    csum[jbc] = sum(cSds[jbc,])
    cat("Sum of column standard deviations: ",csum[jbc],"\n")
  }
  jbc = which.min(csum)
  cat("Centering columns is ",jbc,"\n")
  # Now centre each genome's boundaries.
  bndsc = orbnds
  cadj = rep(0,nrb)
  for (jg in 1:nrb){
    bndsc[jg,] = sapply(1:ncb,function(jb){orbnds[jg,jb]-orbnds[jg,jbc]})
  }
  
  refbnds = sapply(1:ncb,function(jb){median(bndsc[,jb])})
  # Now find adjustments for those boundaries that are not identified as good.
  br = which(sapply(1:nrb,function(jg){sum(goodBnds[jg,])<ncol(goodBnds)}))
  newbndsc = bndsc
  for (jg in br){
    ig = which(goodBnds[jg,]==1)
    nig = length(ig)
    fit = lm(bndsc[jg,ig] ~ refbnds[ig])
    ib = setdiff(1:ncb,ig)
    for (jc in ib){
      newbndsc[jg,jc] = fit$coefficients[2]*refbnds[jc] + fit$coefficients[1]
    }
  }
  # Now decentre all boundary locations.
  newbnds = bnds
  print(dim(newbndsc))
  print(dim(orbnds))
  for (jg in 1:nrb){
    newbnds[jg,] = strand[jg]*(newbndsc[jg,] + orbnds[jg,jbc])
  }
  out = newbnds  
}


adjBnds = row_subset_correction(bnds)

  ##################################################
shortrrn = rep(FALSE,6*Ng);   shortrrn2 = rep(FALSE,6*Ng)
jc1=0;  jc2 = 0

for (jjg in 1:length(IP$indpOp)){
  jg = IP$indpOp[jjg]
  if (length(RR[[jg]]$segSeqs)==3){
    nops=1
    jc1 = jc1 + 1
    shortrrn[jc1] = !(nchar(RR[[jg]][[1]][1])>4000)
    jc2 = jc2 + 1
    shortrrn2[jc2] = shortrrn[jc1]
  }
  if (length(RR[[jg]]$segSeqs)>3){
    nops = length(which((RR[[jg]]$locs[,1]=="rrn")))
    for (jop in 1:nops){
      jc1 = jc1 + 1
      shortrrn[jc1] = !(nchar(RR[[jg]][[1]][jop,1])>4000)
    }
  }
}
shortrrn = shortrrn[1:jc1]
shortrrn2 = shortrrn2[1:jc2]
  
#######################################################################  
#######################################################################  
#######################################################################
# Code to identify those RR elements that have an rrn locs entry that is
# a dummy value (i.e. set to -1), and which row(s) of RR[[j]]$locs is (are)
# dummy entries.  This will probably be needed for rrn_subregions_Gmatrix..R code.
xx=sapply(1:length(RR),function(k){
          out=which(sapply(1:nrow(RR[[k]]$locs),function(kk){
                            RR[[k]]$locs[kk,1]=="rrn" && RR[[k]]$locs[kk,2]<0
                          }))
           })
jgBadrrn = which(sapply(1:length(xx),function(j){length(xx[[j]])>0}))
iLocRowBad = unlist(xx)
print(cbind(jgBadrrn,iLocRowBad))
# Print the locs table for a specific one of the errant genomes
print( RR[[jgBadrrn[23]]]$locs)

#######################################################################  
#######################################################################  
#######################################################################
# Create code to process the 16S and 23S entries from RR[[jg]]$locs to  
# remove any invalid length sub-units. Shuffle up to fill in removed lines.
refineRR = function(RR,bnds){
  # Removes rows which have 16S  or 23S lengths outside specified bounds.
  # These bounds are given by bnds which is a 2x2 matrix giving lower
  # and upper bounds for (row 1) 16S, and (row2) 23S.
  # Removal must be done for both segSeqs and locs, but segSeqs modification
  # is not by row but by single subunit.
  # 19 April 2023                                              [cjw]
  
  #    gseqs = matrix("",nrow=max(2,n16S,n23S),ncol=3)
  for (jg in 1:length(RR)){
    rr = RR[[jg]];  rr0 = rr
    i16S = which(rr$locs[,1]=="16S");    n16S = length(i16S)
    gseqs = matrix("",nrow=max(2,n16S,n23S),ncol=3)
    gseqs = rr$segSeqs
    id16 = which(sapply(i16S,function(j){
      t1 = abs(as.numeric(rr$locs[j,3]) - as.numeric(rr$locs[j,2]))
      out=(t1>bnds[1,1]) && (t1<bnds[1,2])  }  )  )
    noCorrection16 = length(i16S)==length(id16)
    nn16S = length(id16)
    if (!noCorrection16){
      cat("Genome ",jg," requires correction of entry ",setdiff(i16S,id16),"\n")
      rr$locs = rr$locs[-c(setdiff(i16S,id16)),]
      rr$segSeqs[1:nn16S,2] = rr$segSeqs[id16,2]
      rr$segSeqs[(nn16S+1):n16S,2] = "" 
      cat(" 16S correction for genome ",jg,"  Correcting ",setdiff(i16S,id16),"\n")
    }
    # Note that RR[[jg]]$locs might already have had rows removed from the 16S processing just above.
    i23S = which(rr$locs[,1]=="23S");    n23S = length(i23S)  
    id23 = which(sapply(i23S,function(j){
      t1 = abs(as.numeric(rr$locs[j,3]) - as.numeric(rr$locs[j,2]))
      out=(t1>bnds[2,1]) && (t1<bnds[2,2])  }  )  )
    noCorrection23 = length(i23S)==length(id23)
    nn23S = length(id23)
    if (!noCorrection23){
      cat("Genome ",jg," requires correction of entry ",setdiff(i23S,i23S[id23]),"\n")
      rr$locs = rr$locs[-c(setdiff(i23S,i23S[id23])),]
      id23a = i23S[id23] - nn16S;  
      rr$segSeqs[1:nn23S,3] = rr$segSeqs[id23a,3]
      rr$segSeqs[(nn23S+1):n23S,3] = ""
      cat(" 23S correction for genome ",jg,"  Correcting row ",setdiff(i23S,i23S[id23]),"\n")
    }
    rr$segSeqs = rr$segSeqs[1:(max(2,nn16S,nn23S)),]
    RR[[jg]] = rr
  }
  out = RR
}


# load(file.path(outRDatapath,outname))




######################################################################################################### 
#########################################################################################################  
####################            Probably Obsolete Code (22 April 2023)             ######################
#########################################################################################################  
#########################################################################################################

extract_return_16S23Srrn = function(genomespath,bac_fasta_filename, start_loc16S, stop_loc16S,
                                    start_loc23S, stop_loc23S,chromLength){
  # NOTE: Probably of no further value (6 Oct 2022)
  # Returns a set of genomic locations at which 16S, 23S and rrn operons are found to start, together
  # with the rrn sequences for those starting locations.
  # Input parameters:   genomespath,bac_fasta_filename, rrn_start_loc, rrn_16S_stop_loc,is_reversed
  # Derived from extract_return_rrn() of 29 August 2022 which could process genomes with more than 1 chromosome.
  # Revised 21 September 2022 to deal separately with 16S and 23S.
  
  # Read the fasta file specified, check whether it has more than one sequence in it, and if so select the 
  # the longest of these sequences.  Also, ensure that only ACGT entries are present - replacing 
  # non-canonical bases by a random sample of the canonical bases.   The output from this section is  bacseq.
  locsRegs = rep(0,6)
  jg = which(Dg==bac_fasta_filename)
  gen_fa = read.fasta(file=file.path(genomespath,bac_fasta_filename), seqtype="DNA",forceDNAtolower = FALSE)
  chrom = which(unlist(sapply(1:length(gen_fa),function(j){grep("chromosome",attributes(gen_fa[[j]])$Annot)==1})))
  nchr = length(chrom)
  if (nchr==1){
    gen_fa[[1]] = gen_fa[[chrom]]
    chrom.Len = length(gen_fa[[1]])
  } else if (nchr>0) {
    cat("\n Genome ",bac_fasta_filename, " has ",nchr, "chromosomes identified.\n")
    chrom.Len = rep(0,nchr)
    gen_fa[[1]] = gen_fa[[chrom[1]]]
    chrom.Len = sapply(1:nchr,function(jc){out = length(gen_fa[[chrom[jc]]])})
    #    gen_fa[[1]] = sapply(2:nchr,function(jc){t1 = c(gen_fa[[1]],gen_fa[[chrom[jc]]])})
    gf = sapply(2:nchr,function(jc){t1 = c(gen_fa[[1]],gen_fa[[chrom[jc]]])})
  } # Some of the genomes have a few "N" entries.Replace non-canonical base by "A" 
  # and advise number of replacements.
  noncanon = which(!(gen_fa[[1]] %in% c("A","G","C","T") ))
  nc_codes = gen_fa[[1]][noncanon]
  if (length(noncanon)>0){ # Replace non-canonical base by "A" and advise number of replacements 
    print(paste("Genome ",jg," has ",length(noncanon)," non_canonical entries.",sep=""))
    #    print(noncanon)
    gen_fa[[1]][noncanon] = rep(sample(c("A","G","C","T"),1,replace=TRUE),length(noncanon))
  }
  extend = 40
  bacseq = c2s(gen_fa[[1]])
  
  lb16S = min(start_loc16S,stop_loc16S)-extend;   ub16S = max(start_loc16S,stop_loc16S)+extend
  seq16S = substr(bacseq,start=lb16S, stop=ub16S)
  lb23S = min(start_loc23S,stop_loc23S)-extend;   ub23S = max(start_loc23S,stop_loc23S)+extend
  seq23S = substr(bacseq,start=lb23S, stop=ub23S)
  locsRegs[1:4] = c(lb16S,ub16S,lb23S,ub23S)
  out = list(seq16S=seq16S, seq23S=seq23S,locs=locsRegs)
}



#############################################
for (jg in 1:length(Dg)){
  argstr = paste("'^>' ",file.path(refGenomespath,Dg[jg]),sep=" ")
  cat("\n",jg,"\n")
  system2("grep",arg = argstr)
}



