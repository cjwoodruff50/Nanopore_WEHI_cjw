# Rscript to sub-sample from selected species of the Sereika D6322 in-silico amplicon datasets.
# Method outline in diary notes of 27 July 2023.
#
# Major revision 07/01/2024 to include sub-sampling as done initially - that is, 4 sub-samplings
# with manually selected sub-sampling factors.  This revision was on version of 21 November 2023.
#
# Revising 4 Feb 2024 to build code around a simplified file directory structure as outlined in 
# diary of 2-3 February.
# Simplified verion of  subsample_D6322_fast1.R   of 04022024 for pub lication.
# 
# Note that the 16S and 23 primary D6322 fastq files as used as input to RAD must be in a 
# stand-alone directory for accessing by readFastq(). 
# Create a sub-directory  in  under basePath  and set inD6322fastqPath = file.path(basePath"in")
#  - see line 392 (as of 04022024).
#
# 19 February 2024                                                                       [cjw]

if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
packages <- c("stringr","seqinr","tictoc","latex2exp","ShortRead", "DECIPHER","bfsl","dplyr",
              "nnet", "MASS","stringdist","phyloseq","vegan","ape","dada2")    # Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  BiocManager::install(packages[!installed_packages])
  #  BiocManager::install(packages[!installed_packages],repos="https://cran.ms.unimelb.edu.au/")
}# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

speciesSet = c("subtilis","faecalis","coli","monocytogenes","aeruginosa","enterica","aureus")
Ns = length(speciesSet)

#######################################################################################
##########################                               ##############################
##########################           FUNCTIONS           ##############################
##########################                               ##############################
#######################################################################################

makeDesignSet = function(Full, Nv,wtsSet=c(0.01,0.02,0.05,0.1),same="TRUE"){
  # Full is a 2 x Nv matrix with row 1 being the count of 16S reads, row2 the count of 
  # 23S reads. Subsampling is performed independently on both rows.
  # Returns choose(Nv,2) of length Nv vectors of sub-sampling weights, subject to a set of 
  # constraints. These subsets are returned as the columns of a matrix.
  # Constraints are as follows:-
  #    1. The minimum number of reads from a particular species in the sub-sample is minCount.
  #    2. Each sub-sampled dataset has exactly 2 species subsampled.
  #    3. Sub-sampling must not reduce the total full dataset count by more than 50%.
  #    4. All pairs of species are included in the sample of datasets.
  #    5. If desired (set by input same), the same down-sampling is used for pairs of 16S and 23S  
  #
  # 23 November 2023                                                          [cjw]
  minCount = 6
  nw = length(wtsSet)
  Des16 = matrix(1,nrow=Nv,ncol=choose(Nv,2))
  Des23 = matrix(1,nrow=Nv,ncol=choose(Nv,2))
  Dsub16 = matrix(1,nrow=Nv,ncol=choose(Nv,2))
  Dsub23 = matrix(1,nrow=Nv,ncol=choose(Nv,2))
  jcol = 0
  for (jsp1 in 1:(Nv-1)){
    for (jsp2 in (jsp1+1):Nv){
      jcol = jcol+1
      # Randomly choose a pair of weights, check to see whether constrains 1 or 3 are contravened,
      # and resample weights if so.
      OK = FALSE
      OK16 = FALSE
      OK23 = FALSE
      while (!OK) {
        vec16 = Full[1,]
        vec23 = Full[2,]
        iw = sample(1:nw,2)
        vec16[jsp1] = Full[1,jsp1]*wtsSet[iw[1]]
        vec16[jsp2] = Full[1,jsp2]*wtsSet[iw[2]]
        b1 = (vec16[jsp1]>(minCount-1)) &&  (vec16[jsp2]>(minCount-1))  # Constraint 1
        b2 = (sum(vec16)>0.5*sum(Full[1,])) &&  (sum(vec23)>0.6*sum(Full[2,]))
        OK16 = b1 && b2
        if (same){
          vec23 = Full[2,]
          vec23[jsp1] = Full[2,jsp1]*wtsSet[iw[1]]
          vec23[jsp2] = Full[2,jsp2]*wtsSet[iw[2]]
          b1 = (vec23[jsp1]>(minCount-1)) &&  (vec23[jsp2]>(minCount-1))  # Constraint 1
          b2 = sum(vec23)>0.5*sum(Full[2,]) 
          OK23 = b1 && b2
          OK = OK16 && OK23
        }
        else {
          OK = OK16
        }
      }
      Des16[c(jsp1,jsp2),jcol] = wtsSet[iw]
      if (same){
        Des23[c(jsp1,jsp2),jcol] = wtsSet[iw]
      } else {
        OK = FALSE
        while (!OK) {
          vec23 = Full[2,]
          iw = sample(1:nw,2)
          vec23[jsp1] = Full[2,jsp1]*wtsSet[iw[1]]
          vec23[jsp2] = Full[2,jsp2]*wtsSet[iw[2]]
          b1 = (vec23[jsp1]>(minCount-1)) &&  (vec23[jsp2]>(minCount-1))  # Constraint 1
          b2 = sum(vec23)>0.5*sum(Full[2,]) 
          OK = b1 && b2
        }
        Des23[c(jsp1,jsp2),jcol] = wtsSet[iw]
      }
      Dsub16[,jcol] = round(vec16);   Dsub23[,jcol] = round(vec23)
    }
  }
  out = list(D16=Dsub16, D23=Dsub23, Design16=Des16, Design23=Des23)
}


#######################################################################################
############################     RUN INITIALISATION      ##############################
#######################################################################################
whichD6322Set = 2  
whichDesign = 5*(whichD6322Set-1)

extracting = TRUE

out_dateString = ifelse(whichD6322Set==1,"18012024", "19022024")

#######################################################################################
#######################################################################################
##########################                               ##############################
##########################           MAIN BODY           ##############################
##########################                               ##############################
#######################################################################################
#######################################################################################
# Note that ordering of species in giving counts immediately below is (Bs, Ef, Ec, Lm, Pa, Se, Sa)
if (whichD6322Set==1){   
  D6322full16.specCounts = c(558, 938, 478, 798, 223, 979, 529)
  D6322full23.specCounts = c(374, 739, 454, 501, 192, 837, 615)   
  stem11 = "amplicon_D6322_"
  stem12 =  "_05012024_"        #     "_trimmed_05012024_filtered_denoise_"
  stem21 = "amplicon_D6322_"
  stem22 = "_05012024_D6322_05012024_filtered_denoise"
  stem31 = "amplicon_D6322_"
  stem32 = "_trimmed_05012024"
  basePath = "/stornext/Bioinf/data/lab_speed/cjw/microbiome/paper_microbiome_ASV_check"
  outRDatapath = file.path(basePath,"output/RData")
  out_fastqPath = file.path(basePath,"Sub")
} else if (whichD6322Set==2){ # 
  D6322fullcheckA16.specCounts = c(812,1343,658,1103,320,1407,744)       # This vector is ordered by genus alphabetically
  D6322fullcheckA23.specCounts = c(539,1053,633,713,265,1210,859) 
  
  D6322_denoise_dateString = "05022024" 
  denoise_dateString = "06022024"   # Only a fill-in date, pending a new set of primary D6322 dataset denosing being completed.
  D6322full16.specCounts = c(812,1343,658,1103,320,1407,744)       # This vector is ordered by genus alphabetically
  D6322full23.specCounts = c(539,1053,633,713,265,1210,859)        # This vector is ordered by genus alphabetically   
  stem11 = "amplicon_D6322_"
  stem12 =  paste("_",D6322_denoise_dateString,"_filtered_denoise",sep="")        #     "_trimmed_05012024_filtered_denoise_"
  stem21 = "amplicon_D6322_"
  stem22 = paste("_",D6322_denoise_dateString,"_filtered_denoise",sep="")
  stem31 = "amplicon_D6322_"
  stem32 = paste("_",D6322_denoise_dateString,"_filtered",sep="")
  basePath = "/vast/projects/rrn/microbiome/papercheck2"
  outRDatapath = file.path(basePath,"output/RData")
  out_fastqPath = basePath
} else {
  cat("\n INVALID whichD6322Set value. \n"); break
}
Ns = length(speciesSet)

seqtype = "16S"   #   or  "23S"
for (seqtype in c("16S","23S")){
  totASVset = NULL
  whichDesign = 5*(whichD6322Set-1)    # +ifelse(seqtype=="16S",1,2)
  inname1 = paste(stem11,seqtype,stem12, "_names","_out.txt",sep="")
  inname2 = paste(stem21,seqtype,stem22, "_indices","_out.txt",sep="")
  infastqname = ifelse(seqtype=="16S",paste(stem31,seqtype,stem32,".fastq",sep=""),
                       paste(stem31,seqtype,stem32,".fastq",sep=""))
  cat(seqtype,infastqname,"\n", inname1,"\n")  

  # From RAD text output in the "names" file, parse each name to extract the index into the fastq
  # sequences. Store in the data.frame A.df.
  A = readLines(con=file.path(basePath,inname1))
  Na = length(A)
  A.ind = rep(0,Na)
  A.ee = rep(0,Na)
  for (j in 1:Na){
    t1 = unlist(strsplit(A[j],split="[|]"))
    A.ind[j] = as.numeric(substr(t1[1],start=4,stop=nchar(t1[1]))) 
    A.ee[j] = as.numeric(unlist(strsplit(t1[2],"="))[2])
  }
  A.df = data.frame(index=A.ind, EstErr = A.ee, seqname = A)
  
  # Using the primary output data from the code blastn_output_analysis_v02.R, for a single species
  # from each strain identify all ASVs that optimally align to that species.  Then use the RAD
  # output text file, "indices", which indexes into A.df$index, to get indices into the original 
  # fastq file fed to RAD for each (species,strain) set.
  # Now we need to select a D6322 species and strain (if there is ambiguity for that species) 
  spec_ind = vector("list",Ns)
  for (jspec in 1:length(speciesSet)){
    species = speciesSet[jspec] #   strain = speciesStrain[jspec]
    # Identify all ASVs that optimally align to that species/strain.
    which_subunit = seqtype
    strainAnalysis_dateString =  ifelse(whichD6322Set==3,"05012024","06022024")       # "15072023"
    whichDN = ifelse(seqtype=="16S",2*whichDesign+1,2*whichDesign+2)
    inname = paste("strainAnalysis_",which_subunit,"_",whichDN,"_D6322_",strainAnalysis_dateString,".RData",sep="")
    load(file=file.path(outRDatapath,inname))   #   Loads gss16/23,strainsASVset16/23.df,Allstrains.df) 
    if (seqtype=="16S"){strainsASVset.df = strainsASVset16.df} else {strainsASVset.df  = strainsASVset23.df}
    k = which(sapply(1:nrow(strainsASVset.df),function(j){
      t1 = unlist(strsplit(strainsASVset.df[j,"GenusSpecies"],split="_")) 
      out=t1[2]==species}))
    klen = length(unique(strainsASVset.df[k,"strainASVs"]))
    if (klen>1){ # This species has multiple strains having at least 1 non-identical ASVset.
      # Unlist the strainASVs strings for each strain and pool the indices of all ASVs. Then form the
      # vector of unique ASV indices, calling it ASVset
      lenk = length(k)
      asv = NULL
      for (jk in 1:lenk){
        zz=as.numeric(unlist(strsplit(strainsASVset.df[k[jk],"strainASVs"],split="_")))
        asv = append(asv,zz)
      }
      ASVset = unique(asv)
    } else {
      k = k[1]
      ASVset = as.numeric(unlist(strsplit(strainsASVset.df[k,"strainASVs"],split="_")))
    }
    #    print(ASVset)
    totASVset = append(totASVset,ASVset)
    
    # Now get the index of read sequences associated with each of these ASVs into a single integer vector.
    B = readLines(con=file.path(basePath,inname2))
    
    Nb = length(B)
    seqsVec = rep(-1,Na)
    jc = 0
    for (j in 1:length(ASVset)){
      jj = ASVset[j]
      t1 = as.numeric(unlist(strsplit(B[jj],split="\t")))
      seqsVec[(jc+1):(jc+length(t1))] = t1
      jc = jc+length(t1)
    }
    seqsVec = seqsVec[1:jc]
    seqsVecord = sort(seqsVec)
    
    # A.ind[seqsVecord] should be giving the indices of the fastq sequences that probably are from
    # the currently-selected species and that were processed by RAD.
    spec_ind[[jspec]] = A.ind[seqsVecord]
  }
  #  for (j in 1:Ns){cat(speciesSet[j],length(spec_ind[[j]]),"\n")}
  totreads = sum(sapply(1:length(speciesSet),function(j){length(spec_ind[[j]])}));  cat("Total reads ",totreads,"\n")
  if (seqtype == "16S") {
    spec_ind16 = spec_ind
    totASVset16 = totASVset
    ObsCounts16 = sapply(1:7,function(j){length(spec_ind[[j]])})
  } else {
    spec_ind23 = spec_ind
    totASVset23 = totASVset
    ObsCounts23 = sapply(1:7,function(j){length(spec_ind[[j]])})
  }
}     #   end    seqtype   loop
print(c(max(totASVset16),length(totASVset16),max(totASVset23),length(totASVset23)))
max.Exp = max(D6322fullcheckA16.specCounts,D6322fullcheckA23.specCounts)
min.Exp = min(D6322fullcheckA16.specCounts,D6322fullcheckA23.specCounts)
max.Obs = max(ObsCounts16,ObsCounts23)
min.Obs = min(ObsCounts16,ObsCounts23)
plot(D6322fullcheckA16.specCounts, ObsCounts16,pch=1, col=1, xlim=c(min.Exp,max.Exp), ylim=c(min.Obs,max.Obs),
       xlab="Expected", ylab = "Observed", main="D6322full Species Read Counts - RAD-derived")
points(D6322fullcheckA23.specCounts, ObsCounts23,pch=2, col=2)
cat("Observed vs. Expected 16S \n")
print(cbind(D6322fullcheckA16.specCounts,ObsCounts16))
cat(" \n Observed vs. Expected 16S \n")
print(cbind(D6322fullcheckA23.specCounts,ObsCounts23))

# Now start the downsampling process.
Nsub = 4
Dsub = matrix(1,nrow=Ns,ncol=Nsub)
Dsub[c(1,6,7),1] = c(0.1,0.5,0.2)
Dsub[c(1,2,3),2] = c(0.5, 0.01,0.1)
Dsub[c(2,7),3] = c(0.02,0.02)
Dsub[c(3,7),4] = c(0.1,0.02)
rownames(Dsub) = c("Bs","Ef","Ec","Lm","Pa","Se","Sa")
colnames(Dsub) = c("Sub21","Sub22","Sub23","Sub24")
subIDSet = paste("Sub",1:Nsub,sep="")
specCounts16 = matrix(-1, nrow= length(subIDSet), ncol = Ns);  specCounts23 = specCounts16
rownames(specCounts16) = subIDSet
colnames(specCounts16) = speciesSet
rownames(specCounts23) = subIDSet
colnames(specCounts23) = speciesSet

DDind = vector("list",14)
for (whichSub in 1:Nsub){
  for (seqtype in c("16S","23S")){
    if (seqtype=="16S"){
      subvec = Dsub[,whichSub]
      spec_ind = spec_ind16
    } else {
      subvec=Dsub[,whichSub]
      spec_ind = spec_ind23
    }
    SubsampleTable = data.frame(species = speciesSet, fraction=subvec)
    subID = paste("Sub",whichD6322Set,whichSub,sep="")
    # Select the sets of indices for each species, implementing the sub-sampling in this process.
    indset1 = sapply(1:Ns,function(j){
      Nind = length(spec_ind[[j]])
      nind = round(SubsampleTable[j,"fraction"]*Nind)
      isamp = sample(1:Nind,nind)
      out = isamp
    })
    Dindset = NULL
    for (j in 1:Ns){Dindset = append(Dindset,spec_ind[[j]][indset1[[j]]])}
    #for (j in 1:Ns){cat("Subsampled Read Counts ",speciesSet[j],length(spec_ind[[j]][indset1[[j]]]),"\n")}
    totSubsampleReads = sum(sapply(1:length(speciesSet),function(j){length(spec_ind[[j]][indset1[[j]]])}))  
    cat("Total reads ",totSubsampleReads,"\n")
    if (seqtype=="16S"){
      specCounts16[whichSub,] = sapply(1:length(speciesSet),function(j){length(spec_ind[[j]][indset1[[j]]])})
    } else {
      specCounts23[whichSub,] = sapply(1:length(speciesSet),function(j){length(spec_ind[[j]][indset1[[j]]])})
    }
    out_fastqName = paste("amplicon_D6322_",seqtype,stem22,"_",subID,"_",out_dateString,".fastq",sep="") 
    inD6322fastqPath = file.path(basePath,"in")
    gen_fq = readFastq(dirPath=inD6322fastqPath,pattern=seqtype)
    sub_gen_fq = gen_fq[Dindset]
    if (extracting){
      # We now extract these sequences.
      writeFastq(sub_gen_fq,file=file.path(out_fastqPath,out_fastqName), mode="w", withIds=TRUE, compress=FALSE)
    }
  }     #   end    seqtype   loop
}       #    end    whichDesign   loop

cat("16S rRNA gene:- \n");    print(specCounts16)
cat("23S rRNA gene:- \n");    print(specCounts23)






