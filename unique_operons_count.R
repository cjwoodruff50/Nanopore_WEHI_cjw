# Script to identify the unique 16S/23S/rrn operons of a genome and how many
# of each unique sequence are present.
#
# source("/stornext/Bioinf/data/lab_speed/cjw/microbiome/scripts/Rscripts/unique_operons_count.R")
#
# (10 July 2023                                                        [cjw])
# 22 January 2024                                                      [cjw]  
# Minor tidy-up 3 February 2024                                        [cjw]  

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

# Alignment parameters
mat=nucleotideSubstitutionMatrix(match = 1, mismatch = -2, baseOnly = TRUE)
substitute_mat1 = mat
substitute_mat2 = mat
align_params = list(align_type ="global", gapOpen = 2, gapExtend = 2)
align_parameters1 = list(align_type="global-local", gapOpen=2, gapExtend=2)
align_parameters2 = list(align_type="overlap", gapOpen=1, gapExtend=0.5)
gapOpen = 2
gapExtend = 2


out_dateString = "03022024"         # "10072023" 
baseset = c("A","C","G","T")

#########################################################################
###################                               #######################
###################           FUNCTIONS           #######################
###################                               #######################
#########################################################################

allPIDS1strain = function(bacseq){
  # Returns a matrix of PID and Score values for all pairwise alignments of the sequences in bacseq
  # 20 June 2023                                        [cjw]
  Nb = length(bacseq)
  mopsP = matrix(0,nrow=Nb,ncol=Nb)
  mopsS = matrix(0,nrow=Nb,ncol=Nb)
  for (j1 in 1:Nb){
    s2 = bacseq[j1]
    for (j2 in setdiff(1:Nb,j1)){
      s1 = bacseq[j2]
      algnf = pairwiseAlignment(s1,s2,type="global",substitutionMatrix=mat,
                                gapOpening=gapOpen,gapExtension=gapExtend)
      algnr = pairwiseAlignment(s1,reverseComplement(DNAString(s2)),type="global",substitutionMatrix=mat,
                                gapOpening=gapOpen,gapExtension=gapExtend)
      Al = if (score(algnf)>score(algnr)){Al=algnf} else {Al=algnr}
      #    writePairwiseAlignments(Al[[k]])
      mopsP[j1,j2] = pid(Al)
      mopsS[j1,j2] = score(Al)
    }
  }
  out = list(PID=mopsP, Score=mopsS)
}



uniqueByAlign = function(specid,seqtyp,seqPath){
  # Returns a list(uniques, PIDmat,Scoremat),where uniques gives the set of identical "operon"
  # for each strain of the species, specid, from the files in seqPath having the identifier
  # seqtyp in their name. For example, specid="coli", seqtyp="16S", 
  # seqPath= "/stornext/Bioinf/data/lab_speed/cjw/genomes/D6322.refseq/ssrRNAs"
  # 
  # 10 July 2023                                       [cjw]
  pat = paste(specid,sep="_")
#  Dspec = dir(file.path(seqPath,seqtyp),pattern=pat)
  Dspec = dir(seqPath,pattern=pat)
  is.fasta = which(sapply(1:length(Dspec),function(k){Len=nchar(Dspec[k]); out=substr(Dspec[k],start=Len-2,stop=Len)==".fa"}))
  D16S = Dspec[is.fasta]
  cat("Accessing directory ",file.path(seqPath,seqtyp),"\n Searching for pattern ",pat,"\n","  Returned ",length(D16S)," files. \n")
  print(D16S)
  
  EqPID = vector("list",length(D16S))
  outline = rep("",length(D16S))
  M = vector("list",length(D16S))
  S = vector("list",length(D16S))
  for (jg in 1:length(D16S)){
    cat("jg ",jg," file: ",D16S[jg],"\n")
#    gen_fa = read.fasta(file=file.path(seqPath,seqtyp,D16S[jg]), seqtype="DNA",forceDNAtolower = FALSE)
    gen_fa = read.fasta(file=file.path(seqPath,D16S[jg]), seqtype="DNA",forceDNAtolower = FALSE)
    n16S = length(gen_fa)
    if (n16S<2){
      EqPID[[jg]] = 1
      s2 = paste(as.character(EqPID[[jg]][1]),sep="    ",collapse="    ")
      t1 = unlist(strsplit(D16S[jg],split="_"))
      isp = which(sapply(1:length(t1),function(j){t1[j]==specid}))
      if (is.na(str_locate(D16S[jg],"D6322")[1])){
        outline[jg] = paste(t1[isp-1],specid,t1[isp-2],s2,sep="   ")
      } else {
        outline[jg] = paste(t1[isp-1],specid,"D6322",s2,sep="   ")
      }
  #    outline[jg] = "Only one operon."
      M[[jg]] = 100
      S[[jg]] = length(gen_fa[[1]])
    } else {
      bacseq=rep("",n16S)
      nameStem=rep("",n16S)
      for (jop in 1:n16S){
        bacseq[jop] = c2s(getSequence(gen_fa[[jop]]))
        nameStem[jop] = attributes(gen_fa[[jop]])$name
      }
      mops = allPIDS1strain(bacseq)
      nops = nrow(mops$PID)
      # Identify groups of sequences that have a PID of 100. These are 
      # considered to be identical sequences - though they might have a 
      # gapped alignment with no mis-matches.
      rset = 1:n16S
      eqp = vector("list",n16S)
      jc = 0;  jcc = 0
      while (length(rset)>0){
        jc = jc + 1
        izz = which(mops$PID[jc,rset]>99.99)
        if (length(izz)>0){
          jcc = jcc + 1
          eqp[[jcc]] = c(jc,rset[izz])
        } else if (jc %in% rset){
          jcc = jcc + 1
          eqp[[jcc]] = jc
        }
        #        cat("jc ",jc,"  jcc ",jcc,"  rset:",rset,"   eqp[[jcc]]:",eqp[[jcc]],"\n")
        rset = setdiff(rset,c(jc,eqp[[jcc]]))
      }
      EqPID[[jg]] = eqp[1:jcc]
      nuniq=jcc
      t1 = unlist(strsplit(D16S[jg],split="_"))
      isp = which(sapply(1:length(t1),function(j){t1[j]==specid}))
      s1 = as.character(EqPID[[jg]][1])
      if (jcc>1){
        for (k in 2:jcc){s1 = append(s1,as.character(EqPID[[jg]][k]))}
      }
      s2 = paste(s1,sep="    ",collapse="    ")
      if (is.na(str_locate(D16S[jg],"D6322")[1])){
        outline[jg] = paste(t1[isp-1],specid,t1[isp-2],s2,sep="   ")
      } else {
        outline[jg] = paste(t1[isp-1],specid,"D6322",s2,sep="   ")
      }
      M[[jg]] = mops$PID
      S[[jg]] = mops$Score
    }
    cat("jg ",jg," uniques: ",outline[jg],"\n")
  }
  out = list(uniques=outline, PIDmat=M,Scoremat=S)
}



#######################################################################################
#######################################################################################
##########################                               ##############################
##########################           MAIN BODY           ##############################
##########################                               ##############################
#######################################################################################
#######################################################################################
# This code uses the call of uniqueByAlign().  More relevant might be alignment and checking 
# whether PID=100%.
speciesSet = c("aeruginosa","aureus","beijerinki","coli","difficile","enterica",
               "faecalis","fermentum","monocytogenes","subtilis")
seqtypeSet = c("16S","23S","rrn")
pathworkDB1 = "/stornext/Bioinf/data/lab_speed/cjw/microbiome/workDB1"   #  "/vast/projects/rrn/microbiome/blastdb/workDB"                #"/stornext/Bioinf/data/lab_speed/cjw/microbiome/toyDB2/output/rrn" 
basepath = "/vast/projects/rrn/microbiome/papercheck"
pathworkDB1 = file.path(basepath,"workDB1")
outpath = file.path(basepath,"output")

# 24 June:  Want to get uniques of strains in toyDB216S_17MayAdj located in folder
#            /vast/projects/rrn/microbiome/blastdb/toyDB2/toy16SDB2_17May
#           Folder /stornext/Bioinf/data/lab_speed/cjw/microbiome/toyDB2/output/rrn
#           has sub-folder 16S and 23S which havethe individual .fa files of sub-unit 
#           sequences for each strain in toyDB2.
#               - includes NRRL-B-1109 with its CP ID as part of name in conventional form.


# Which of the species?
indspec =  c(10,7,4,9,1,6,2)  # This gives results in alphabetical genus order.
indtyp = c(1,2);   nindtyp = length(indtyp) # Note indtype MUST be either 1 or c(1,2) or c(1,2,3)
cat("Set of species chosen:- \n")
print(speciesSet[indspec])
IU = vector("list",length(indspec)*nindtyp)   
for (iseqtyp in seq_along(indtyp)){
  seqtype = seqtypeSet[iseqtyp]
  seqPath = file.path(pathworkDB1,seqtype)   # Assumes separate sub-directories for 16S, 23S, rrn sequences.
  for (isp in seq_along(indspec)){
    ispecid = indspec[isp]
    uind = ifelse(seqtype=="16S",nindtyp*(isp-1)+1,ifelse(seqtype=="23S",nindtyp*(isp-1)+2,nindtyp*(isp-1)+3))
    specid = speciesSet[ispecid]
    # out = list(uniques=outline, PIDmat=M, Scoremat=S)
    cat("Progress:  Started ",seqtype, speciesSet[ispecid],"\n")
    UbA = uniqueByAlign(specid,seqtype,seqPath)
    IU[[uind]] = UbA
#    outname = paste("strain_operon_multiplicity_",seqtype,"_",specid,"_",out_dateString,".RData",sep="")
#    save(file=file.path(seqPath,seqtype,outname),UbA)
#    save(file=file.path(seqPath,outname),UbA)
  }  #      end    specid    loop
}    #      end    seqtyp    loop

# There appears to be a greater number of unique operons for 16S sequences than
# for 23S sequences across the various species.  The following seeks to quantify
# the variation in number of unique sequences for each sub-unit.
maxStrains = 15   # Assuming there are no more than 15 strains per species
numUniqs = array(0,dim=c(length(indtyp),length(indspec),maxStrains))
for (iseqtyp in seq_along(indtyp)){
  seqtype = seqtypeSet[iseqtyp]
  seqPath = pathworkDB1    # Assumes 16S, 23S, rrn are in the same folder - needs attention.
  for (isp in seq_along(indspec)){
    ispecid = indspec[isp]
    uind = ifelse(seqtype=="16S",2*(isp-1)+1,2*isp)
    nops = length(IU[[uind]]$PID)
    for (jop in 1:nops){
      s1 = unlist(strsplit(IU[[uind]]$uniques[jop],split="  "))
      nu = length(which(sapply(4:length(s1),function(j){out=!(s1[j]=="")})))
      numUniqs[iseqtyp,isp,jop] = nu
    }
  }
}
nu16S = rowSums(numUniqs[1,,])
nu23S = rowSums(numUniqs[2,,])
print(cbind(nu16S,nu23S))
Numuniqs = data.frame(Nuniq16S=nu16S,Nuniq23S=nu23S)
rownames(Numuniqs) = speciesSet[indspec]
print(Numuniqs)

outname = paste("strain_operon_multiplicity_All_",out_dateString,".RData",sep="")
save(file=file.path(outpath,outname),IU,Numuniqs,numUniqs)





    
    
    
    
    