# Rscript  species_strain_analysis_residualAmbiguity_v02.R .
# Rscript to construct a blastn call given a multi-sequence fasta file of queries and
# a blastn reference database, and then to process taxonomic data generated from local 
# blastn runs. 
# This version is intended as a minimalist version suitable for public access via GitHub,
# written in a smarter way than previous version.
# It is a revision of species_strain_analysis_residualAmbiguity_v01.R of 7 February 2024 .
#
#  source("/stornext/Bioinf/data/lab_speed/cjw/microbiome/scripts/Rscripts/species_strain_analysis_residualAmbiguity_v02.R")
#
#
# 21 March 2024                                                  [cjw]

#######################################################################################
##########################                               ##############################
##########################         INITIALISATIONS       ##############################
##########################                               ##############################
#######################################################################################

args = commandArgs(trailingOnly=TRUE)
slurm = FALSE

if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
packages <- c("stringr","seqinr","tictoc","latex2exp","ShortRead", "DECIPHER","bfsl","dplyr",
              "nnet", "MASS","stringdist","phyloseq","vegan","ape","dada2","compositions","scales")    # Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  BiocManager::install(packages[!installed_packages])
  #  BiocManager::install(packages[!installed_packages],repos="https://cran.ms.unimelb.edu.au/")
}# Packages loading
invisible(lapply(packages, library, character.only = TRUE))
numcores = 4    # Upper limit for use of vc7-shared
verbose = FALSE



#######################################################################################
##########################                               ##############################
##########################           FUNCTIONS           ##############################
##########################                               ##############################
#######################################################################################

getStats2 = function(diff,nhits,A,kstartA,B,kstartB,maxTargetSeqsA=50){
  # 
  # 13 February 2024                                            [cjw]
  if (diff>(min(10,maxTargetSeqsA)+2)){
    nhitsA[k] = as.numeric(unlist(strsplit(A[kstartA[k]+3],split=" "))[2])
    numh = rep(0,nhitsA[k]);   chrh = rep("",nhitsA[k])
    stats = data.frame(ID=chrh,accver=chrh,PID=numh,Alen=numh,MM=numh,Gap=numh,Eval=numh,bitscore=numh,score=numh)
    for (k2 in 1:nhitsA[k]){
      line1 = unlist(strsplit(A[kstartA[k]+3+k2],"\t"))
      stats[k2, 2:ncol(stats)] = line1[c(2,3,4,5,6,11,12,13)]
      line2 = unlist(strsplit(B[kstartB[k]+5+k2],"  ")); i2 = which(nchar(line2)>1)
      stats[k2,1] = line2[i2[1]]
    }
  } else { # Probably have zero hits
    nhitsA[k] = 0
    numh = 0;   chrh = "--"
    stats = data.frame(ID=chrh,accver=chrh,PID=numh,Alen=numh,MM=numh,Gap=numh,Eval=numh,bitscore=numh,score=numh)
  }
  out = stats
}

identify_poorly_aligned_ASV = function(stats2,minAlen){
  # Insert a filtering here that identifies unsatisfactorily aligned ASVs.
  out = !(as.numeric(stats2[1,"Alen"]) > minAlen)
}

numHits = function(nASVs,kstartB,kstart){
  # Returns, for each ASV, the number of alignments returned.  Generally
  # this will be maxTargetSeqsB, but not always.
  # 13 February 2024                                           [cjw]
  nhits = rep(0,nASVs)
  for (k in 1:nASVs){
    ib = which.min(sapply(1:length(kstartC),function(j){d = kstartC[j]-kstartB[k]; out=ifelse(d<0,10000,d)}))
    delta = kstartC[ib]-kstartB[k]
    nhits[k] = delta-8
  }
  out = nhits
}

restoreLong = function(D0){
  # Restores full-length genus and species names of entries in the dataframes D16, D23
  # using the tables possibleGenera and possibleSpecies.
  # 13 July 2023  
  D = D0
  for (j1 in (1:nrow(D))){
    m1 = which(possibleGenera[,2]==D[j1,"genus"])
    if (length(m1)>0){
      D[j1,"genus"] = possibleGenera[m1,1]
    }
    m2 = which(possibleSpecies[,2]==D[j1,"species"])
    if (length(m2)>0){
      D[j1,"species"] = possibleSpecies[m2,1]
    }
  }
  out = D
}

extractShortStrainName = function(longName){
  t1 = unlist(strsplit(longName,split="_"))
  if (t1[1] %in% c("CP","NC","NZ")){
    nom=t1[2];  shortName = substr(nom,start=1,stop=min(8,nchar(nom)))
  }else{
    m = which(sapply(1:length(t1),function(j){t1[j]=="D6322"}))
    if (length(m)>0){
      shortName = "D6322"
    } else {
      if (length(t1)==1){
        shortName = t1
      } else {
        nom = ifelse(nchar(t1[1])>nchar(t1[2]),t1[1],t1[2])
        shortName = substr(nom,start=1,stop=min(8,nchar(nom)))
      }
    }
  }
  out = shortName
}


getStrainCellCounts = function(gss,strainsASVset.df, abund){
  # For all species identified as listed in gss, returns the cellular abundances for each strain
  # of a species that has an ASV associated with in gss and strainsASset.df.
  # No attempt at disambiguation is included in this function. 
  # Also returned are the ASVs contributing to each strain's abundance
  #
  # 14 September 2023                                                                 [cjw]
  specSet = unique(gss$Species);    Ns = length(specSet);
  strainSet = sapply(1:nrow(gss),function(j){out=extractShortStrainName(gss[j,"Strain"])})
  Nstr = length(strainSet)
  strainParams = matrix(0,nrow = Nstr, ncol=5)
  colnames(strainParams) = c("Species","StrainID","RawCount","Cellabund","Numops")
  possSpecs = sapply(1:nrow(gss),function(j){out=unlist(strsplit(gss[j,"Species"],split="_"))[2]})
  strainASVs = vector("list",Nstr)
  for (k in 1:Nstr){ # Process a single line of gss and the corresponding line of strainsASVset.df
    ASVs = as.numeric(unlist(strsplit(strainsASVset.df[k,"strainASVs"],split="_")))
    nops = as.numeric(gss[k,"No.Ops of Strain"])
    strainParams[k,1] = unlist(strsplit(gss[k,"Species"],split="_"))[2]
    strainParams[k,2] = strainSet[k]
    strainParams[k,3] = sum(abund[ASVs])
    strainParams[k,5] = nops
    strainParams[k,4] = round(as.numeric(strainParams[k,"RawCount"])/nops)
    strainASVs[[k]] = list(ASVs=ASVs,ASVCpts=strainsASVset.df[k,"ASVsetsCpts"])
  }
  out = list(StCounts = strainParams, strainASVs = strainASVs)
} 

solveSearch = function(B,resids,allsums1=FALSE){
  #
  K = nrow(B)
  solved = FALSE
  p = 1
  while((!solved) && (p<K)){
    p = p+1
    # Check to see if there is an n-strain solution
    # Form a vector that is the sum of n rows. Check if any elements are not equal to 1.
    # If so reject that set of n rows.  Repeat for all possible n-vectors of rows, recording the successful n-vectors.
    if (!solved){
      goodnvecs = matrix(0,nrow=max(2,choose(K,K-2)),ncol=ncol(B))
      jg = 0
      resids = vector("list",K);  for (m in 1:K){resids[[m]]=1:K}
      allCases = recursive.strainSet2(p,K,resids,duplicates=FALSE)
      ncases = length(allCases) %/% p
      cat("Dimension allCases ",dim(allCases),"   ncases ",ncases,"\n")
      jcs = 0
      while (jcs<ncases && !solved){
        #for (jcs in 1:ncases){
        jcs = jcs+1
        if (ncases>1){ 
          vec = colSums(B[allCases[jcs,],])
          # solved if allsums1==TRUE, requires both internal conditions be FALSE to give solved==TRUE
          solved = ifelse(allsums1,!(any(vec<1) || any(vec>1)),!any(vec<1) ) 
          if (solved){jg = jg+1;  goodnvecs[jg:(jg+ncol(allCases)-1),] = B[allCases[jcs,],]; jg = jg+ncol(allCases)-1}
        } else { 
          vec = colSums(B[allCases,])
          # solved if allsums1==TRUE, requires both internal conditions be FALSE to give solved==TRUE
          solved = ifelse(allsums1,!(any(vec<1) || any(vec>1)),!any(vec<1) ) 
          if (solved){jg = jg+1;  goodnvecs[jg,] = B[allCases[1],]; jg = jg+1}
        }
      }
      if (jg>0){
        if (jg==1){goodnvecs = goodnvecs[1,]} 
        else {
          goodnvecs = goodnvecs[1:jg,]
        }
        if (verbose) print(goodnvecs)
      }
      solved = jg>0
      if (solved){
        if (verbose) cat(p,"strain solution(s) found. \n")
        nK=2
      } else {
        if (verbose) cat("No ",p,"- strain solution found. \n")
        goodnvecs = NULL
      }
    }
  }     #   end    while     loop
  if (is.matrix(allCases)){Bind = allCases[jcs,]} else {Bind = allCases}
  out = list(goodnvecs=goodnvecs,is.solved=solved,Bindices=Bind)
}


minBset2 = function(B){
  # The TASK is to find a smallest set of rows of B such that the sum of these rows 
  # is rep(1,ncol(B)), where B is a moderately sparse binary matrix.
  # This function attempts to solve this task.
  # 17 October 2023                                       [cjw]
  K = nrow(B);  P = ncol(B)
  goodnvecs = matrix(0,nrow=choose(K,K-2),ncol=ncol(B))
  goodMono = NULL;  goodPairs = NULL;   goodTriples = NULL;  goodQuads = NULL
  nK = 0
  singleStrain = FALSE
  # Check if there is a single strain solution
  m = which(rowSums(B)==P)
  solved = length(m)>0
  if (solved){
    if (verbose) {
      if (m==1){
        cat("There is one single strain solution. Check row(s) ",m,"\n")
      } else{
        cat("There are multiple single strain solutions. Check row(s) ",m," of B. \n")
      }
    }
    goodnvecs[m,] = B[m,]
    nK=1
    singleStrain = TRUE
    soln = list(goodnvecs=goodnvecs[m,],is.solved=solved,Bindices=m, isSingle = singleStrain)
  } else {
    if (verbose) cat("No single strain solution found. \n")
  }
  # Now try to use recursive.strainSet(n,resids) to check for solutions.
  if ((!solved)){
    resids = vector("list",K);   for (m in 1:K){resids[[m]] = 1:K}
    result1 = solveSearch(B,resids,allsums1=TRUE)
    solved = result1$is.solved
    if (!solved){ # No solution available that does not split an ASV set contribution to strains.
      # Assume that this function is operating on the reduced B matrix that arises from removal
      # of the jcom strain(s) and the columns that they have non-zero entries for. So no rows
      # are privileged.
      # Steps: 1. If any columns have a single non-zero entry identify the strains involved and 
      #           remove them - and also other columns that have non-zero entries for these strains.
      #        2. Check the new reduced, reduced B matrix to see if a minimal set of strains can be 
      #           determined by minBset2 for this matrix. 
      #                a. If so accept that and set  solved = TRUE
      #                b. If not run a minBset2() algorithm that does not require all column sums to be 1.
      kc1 = which(colSums(B)==1)
      if (length(kc1)>0){
        # In each column that sums to 1 identify in which row of B that 1 occurs.
        jr1 = unique(sapply(1:length(kc1), function(j){which(B[,kc1[j]]==1)}))
        if (length(jr1)==1){b1=TRUE; b11=all(B[jr1,]>0)} else {b1=FALSE;  b11=all(colSums(B[jr1,])>0)}
        if (b1 && b11){
          solved = TRUE
          soln = list(goodnvecs=B[jr1,],is.solved=solved,Bindices=jr1, isSingle = singleStrain)
        } else if (!b1 && b11){
          t1 = jr1
          solved = TRUE
          soln = list(goodnvecs=B[jr1,],is.solved=solved,Bindices=jr1, isSingle = singleStrain)
        } else {
          jvec = setdiff(1:nrow(B),jr1)  
          if (length(jr1)==1){
            kvec = setdiff(1:ncol(B),which(B[jr1,]>0))
          } else {
            kvec = setdiff(1:ncol(B),which(colSums(B[jr1,])>0))
          }
          B1 = B[jvec,kvec]
          if (!is.matrix(B1)){ # Exactly 1 of jvec, kvec has length 1
            if (length(jvec)==1){
              t1 = unique(c(jr1,jvec))
              solved = TRUE
              soln = list(goodnvecs=B[t1,],is.solved=solved,Bindices=t1, isSingle = singleStrain)
            } else { # There are length(jvec) strains that could solve the problem. Choose the one with minimal ASVUs present.
              t1 = rowSums(B)
              jt1 = which.min(t1)
              t2 = unique(c(jr1,jvec[jt1[1]]))
              solved = TRUE
              soln = list(goodnvecs=B[t2,],is.solved=solved,Bindices=t2, isSingle = singleStrain)
            }
          } else if (all(as.integer(colSums(B1))==1)){
            t1 = unique(c(jvec,jr1))
            solved = TRUE
            soln = list(goodnvecs=B[t1,],is.solved=solved,Bindices=t1, isSingle = singleStrain)
          } else {# Case 2.b above
            K1 = nrow(B1)
            result2 = solveSearch(B1,resids,allsums1=FALSE)             
            if (!result2$is.solved){
              cat("No result found.\n")   
              solved = FALSE
              soln = list(goodnvecs=NULL,is.solved=solved,Bindices=-1, isSingle = singleStrain)
            } else {
              Bind = c(jr1,jvec[result2$Bindices])
              soln = list(goodnvecs=B[Bind,],is.solved=TRUE,Bindices=Bind, isSingle = singleStrain)
            }
          }
        }
      } else {
        K = nrow(B)
        result3 = solveSearch(B,resids,allsums1=FALSE)             
        if (!result3$is.solved){cat("No result found.\n")}
        soln = result3
      }
    } else {soln = result1}
  }
  out = soln
}


make_plot_legend_label2 = function(thisSpecies){
  # Generate plot legend label given the set of species - e.g. Given "coli" generate "E.coli"
  # Output is a matrix with species in column 1 and label in column 2.
  # 18 January 2024                                                  [cjw]
  Ns = length(thisSpecies)
  labels = sapply(1:Ns,function(j){
    t1 = substr((possibleGenera[thisSpecies])[j],start=1,stop=1)
    out1 = paste(t1,(possibleSpecies[thisSpecies])[j],sep=".")})
  out = cbind(possibleSpecies[thisSpecies],labels)
}

minBsetSearch = function(B){
  nr = nrow(B);    nc = ncol(B)
  foundSoln = FALSE
  solnSets = NULL
  Bsub = NULL
  kch = 0
  while ((!foundSoln) && (kch<nr)){
    kch = kch+1
    if (kch==nr){ # Can directly check if using all possible rows is a valid solution or not.
      validAll = identical(colSums(B),rep(1,nc))
      if (validAll){
        foundSoln = TRUE
        out = list(foundSoln=foundSoln,  solutionSets=1:nr, task=B, Bcompromise=Bsub)
      }
    } else {
      nsets = choose(nr,kch)
      resids = vector("list",kch)
      for (j in 1:kch){
        resids[[j]] = j:(nr-kch+j)
      }
      indicesSets = recursive.strainSet3(kch,nr,resids,duplicates=FALSE)
      if (ncol(indicesSets$goodnvecs)==1){jch=1}
      if (kch==1){
        allowed = which(rowSums(B)==nc)
        foundSoln = length(allowed)>0
        if (foundSoln){
          nallowed = length(allowed)
          solnSets = vector("list", nallowed)
          for (ks in 1:nallowed){
            solnSets[[ks]] = allowed[ks]
          }
        }
      } else {
        check = rep(FALSE,nsets)
        for (js in 1:nsets){
          check[js] = identical(colSums(B[indicesSets$goodnvecs[js,],]),rep(1,nc))
        }
        allowed = which(check)
        foundSoln = length(allowed)>0
        if (foundSoln){
          nallowed = length(allowed)
          solnSets = vector("list", nallowed)
          for (ks in 1:nallowed){
            solnSets[[ks]] = indicesSets$goodnvecs[allowed[ks],]
          }
        }
      }
    }
  }
  
  if (!foundSoln){
    cat("No ideal solution found. \n Compromise solution being generated.\n")
    kch = 1
    while ((!foundSoln) && (kch<nr)){
      kch = kch+1
      if (kch==nr){ # Can directly check if using all possible rows is a valid solution or not.
        Bsub = B
      } else {
        # Form all kch-size sets of rows.  For each set form the colSums(B[<kch-set>,])
        # vector, cvec. If all elements of cvec>0 then a solution will be created for 
        # this value of kch, as follows. Identify any set that has the minimum of cvec
        # components equal to 2, and randomly choose one of these sets. 
        nsets = choose(nr,kch)
        resids = vector("list",kch)
        for (j in 1:kch){
          resids[[j]] = j:(nr-kch+j)
        }
        indicesSets = recursive.strainSet3(kch,nr,resids,duplicates=FALSE)
        if (ncol(indicesSets$goodnvecs)==1){jch=1}
        check = rep(FALSE,nsets)
        allposSums = NULL;  possibleSets = matrix(0,nrow=100,ncol=kch);  nposSets = 0
        for (js in 1:nsets){
          t1 = colSums(B[indicesSets$goodnvecs[js,],])
          if (length(which(t1>0))==nc){
            allposSums = append(allposSums,js)
            nposSets = nposSets + 1
            possibleSets[nposSets,] = indicesSets$goodnvecs[js,]
          }
        }
        if (!is.null(allposSums)){foundSoln = TRUE}
      }
    }
    
    cat("Solution can be derived for",kch, " strains. \n")
    if (kch<nr){
      possibleSets = possibleSets[1:nposSets]
      # Randomly choose one of these index sets.
      jch = ifelse(nposSets==1,allposSums,sample(allposSums,1))
      Bsub = B[indicesSets$goodnvecs[jch,],]
    }
    cvec = colSums(Bsub)
    ctrim = which(cvec>1)
    BsubScoreMat = matrix(0,nrow=nrow(Bsub),ncol=ncol(Bsub)) 
    # Priority ordering of ASVsets, based on smaller numbered SVsets tending to carry more associated reads.
    for (jj in 1:nrow(BsubScoreMat)){
      BsubScoreMat[jj,] = sapply(1:ncol(Bsub),function(j){Bsub[jj,j]*(ncol(Bsub)+1-j)})
    }
    for (jctr in ctrim){
      if ((jctr==1) && (length(which(Bsub[,1]>0))>0)){
        # Need to resolve ambiguity.  Choose the row of BsubScoreMat that gives the maximum cumsum.
        it0 = which(Bsub[,1]==1)
        im = which.max(sapply(it0,function(j){cumsum(BsubScoreMat[j,])[ncol(Bsub)]}))
      } else{
        it0 = which(Bsub[,jctr]==1)
        im = which.max(sapply(it0,function(j){cumsum(BsubScoreMat[j,])[jctr-1]}))
      }
      it1 = setdiff(it0,it0[im])
      Bsub[it1,jctr] = rep(0,length(it1))
      for (jj in 1:nrow(BsubScoreMat)){
        BsubScoreMat[jj,] = sapply(1:ncol(Bsub),function(j){Bsub[jj,j]*(ncol(Bsub)-j)})
      }
    }
    solnSets = indicesSets$goodnvecs[jch,]
    foundSoln = identical(as.integer(colSums(Bsub)),as.integer(rep(1,ncol(Bsub))))
  }
  out = list(foundSoln=foundSoln,  solutionSets=solnSets, task=B, Bcompromise=Bsub)
}

recursive.strainSet3 = function(p,m,resids,duplicates=FALSE){
  # Given a list, resids, of length >=m whose components are integer vectors, constructs
  # all possible combinations of n integers where each component of resids provides one
  # integer, component j providing the j^th vector component. If duplicates=FALSE
  # any combination accepted for output has no duplicate elements.
  # The sets of integers that constitute resids must be ordered 
  # 6 November 2023                                                    [cjw]
  if (p==1){
    x = array(1:m,dim=c(m,1))
    soln = list(goodnvecs=x, is.solved=TRUE,Bindices=x)
    colnames(soln$goodnvecs) = "C1"
  } else {
    AA = recursive.strainSet3(p-1, m, resids, duplicates)
    A0 = AA$goodnvecs
    nr = dim(A0)[1]
    if (length(resids[[p]])==1){
      vec = sapply(1:nr,function(j){sample(setdiff(resids[[j]],A0[j,]),1)})
      A1 = cbind(A0,vec)
    } else {
      A1 = A0
      for (j in 2:length(resids[[p]])){A1 = rbind(A1,A0)}
      
      vec0 = NULL 
      for(jj in 1:length(resids[[p]])){
        vec1 = rep(resids[[p]][jj],nr)
        vec0=append(vec0,vec1)
      }
      A1 = cbind(A1,vec0)
    }
    allCases = unique(t(sapply(1:nrow(A1),function(j){sort(A1[j,])})))
    if (!duplicates){
      t0 = which(sapply(1:nrow(allCases),
                        function(j){length(unique(allCases[j,]))==ncol(allCases)} ) )
      t1 = allCases[t0, ]
      soln = list(goodnvecs=t1,is.solved=TRUE,Bindices=t0)
    } else {soln = A1}
    colnames(soln$goodnvecs) = paste("C",1:p,sep="")
  }
  out = soln
}

recursive.strainSet2 = function(p,m,resids,duplicates=FALSE){
  # Given a list, resids, of length >=m whose components are integer vectors, constructs
  # all possible combinations of n integers where each component of resids provides one
  # integer, component j providing the j^th vector component. If duplicates=FALSE
  # any combination accepted for output has no duplicate elements.
  # The sets of integers that constitute resids must be ordered 
  # 5 November 2023                                                    [cjw]
  if (p==1){
    x = array(1:m,dim=c(m,1))
    soln = list(goodnvecs=x, is.solved=TRUE,Bindices=x)
    colnames(soln$goodnvecs) = "C1"
  } else {
    AA = recursive.strainSet2(p-1, m, resids, duplicates)
    A0 = AA$goodnvecs
    nr = dim(A0)[1]
    if (length(resids[[p]])==1){
      vec = sapply(1:nr,function(j){sample(setdiff(resids[[j]],A0[j,]),1)})
      A1 = cbind(A0,vec)
    } else {
      A1 = A0
      for (j in 2:length(resids[[p]])){A1 = rbind(A1,A0)}
      
      vec0 = NULL 
      for(jj in 1:length(resids[[p]])){
        vec1 = rep(resids[[p]][jj],nr)
        vec0=append(vec0,vec1)
      }
      A1 = cbind(A1,vec0)
    }
    allCases = unique(t(sapply(1:nrow(A1),function(j){sort(A1[j,])})))
    if (!duplicates){
      t0 = which(sapply(1:nrow(allCases),
               function(j){length(unique(allCases[j,]))==ncol(allCases)} ) )
      t1 = allCases[t0, ]
      soln = list(goodnvecs=t1,is.solved=TRUE,Bindices=t0)
    } else {soln = A1}
    colnames(soln$goodnvecs) = paste("C",1:p,sep="")
  }
  out = soln
}


makeBmatrix = function(jcommon,iSu,specASVus){
  # Create the array whose columns are ASV sets associated with a particular species, and rows
  # are the strains of that species.  The matrix is binary. The ASV sets are sorted numerically
  # ascending.
  # 31 October 2023                                         [cjw]
  ASVus = sort(setdiff(unique(as.vector(specASVus[iSu,])),0))
  B = array(0,dim=c(length(iSu),length(ASVus)))
  colnames(B) = as.character(sort(as.numeric(ASVus)))
  niSu = length(iSu)
  for (jj in jcommon){
    row2 = rep(0,length(ASVus))
    row1 = specASVus[iSu[jj],]
    ic = which(row1>0)
    for(k in 1:length(ic)){
      b1=which(colnames(B)==as.character(row1[ic[k]]))
      row2[b1] = rep(1,length(b1))
    }
    B[jj,] = row2
  }
  for (jj in setdiff(1:niSu,jcommon)){
    row2 = rep(0,length(ASVus))
    row1 = specASVus[iSu[jj],]
    ic = which(row1>0)
    for(k in 1:length(ic)){
      b1=which(colnames(B)==as.character(row1[ic[k]]))
      row2[b1] = rep(1,length(b1))
    }
    if (jj>nrow(B)){
      B = rbind(B,row2)
    } else {
      B[jj,] = row2
    }
  }
  # Check to see whether all ASVus have a strain that is associated with them
  #  if (any(colSums(B)==0)){
  #    im = which(colSums(B)==0)
  #    if (length(im)){
  #      jmiss = sapply(1:length(im),function(k){which(unlist(unique(gss16[i16[setdiff(1:length(i16),jcommon16)],"ASVsetsCpts"])) == as.numeric(ASVus16[im[k]]))})
  #      cat("ASVus ",unlist(unique(gss16[i16[setdiff(1:length(i16),jcommon16)],"ASVsetsCpts"]))[jmiss]," not accounted for.\n")
  #    }
  #  }
  out = list(B=B, nASVu = length(ASVus), ASVs = ASVus)
}


makeReducedBmatrix2 = function(B,jcom){
  # Construct the reduced B16 matrix in which the jcom strain(s), and any strain
  # with identical ASVus, have been removed. Also remove the columns corresponding to
  # these ASVus. 
  # 30 October 2023                                                            [cjw]
  
  nr = nrow(B);  nc = ncol(B)
  # First check whether the jcom row or rows solve the problem.
  solved = FALSE
  if (length(jcom)==1){
    solved = (sum(B[jcom,] == ncol(B)))
  } else if (nc==1){
    solved = sum(B)>0
  } else {
    solved = all(colSums(B[jcom,])>0)
  }
  # The following line is intended to identify through jrid any rows of B that are not in  jcom  
  # but are identical to one of these rows.
  jrid = NULL
  krm = NULL
  for (k in 1:length(jcom)){
    t1 = setdiff(which(sapply(1:nr,function(j){identical(B[j,],B[jcom[k],])})),jcom[k])
    if (length(t1)>0){jrid = append(jrid,t1)}
    t2  = which(B[jcom[k],]==1)
    if (length(t2)>0){krm = append(krm,t2)}
  }
  jrid = setdiff(unique(jrid),jcom)
  krm = sort(unique(krm))
  kBp = setdiff(1:nc,krm)
  if (solved){ # The reduced matrix has all columns (ASV sets) retained but, usually, fewer rows.
    kincl = 1:nc;   kexcl=NULL
    jresid = setdiff(1:nrow(B),c(jcom,jrid))
    if (length(kBp)==0){Bp=NULL; jret=jcom} else {Bp=B[jresid,kBp];  jret=jcom}
    out = list(Bp = Bp, jretained=jret)
  } else { # There is at least one ASV set that is not associated with the jcom strains.
    kincl = kBp;    kexcl = setdiff(1:nc,kBp)
    jBp = setdiff(1:nr,c(jcom,jrid))
    Bp = B[jBp,kBp]
    # It is possible that some rows of Bp would have no non-zero entries.  These need to be identified and removed. 
    if (is.matrix(Bp)) {j0 = which(rowSums(Bp)==0)} else {j0=which(Bp==0)} 
    if (length(j0)>0){
      jBp = setdiff(jBp,jBp[j0])
    }
    Bp = B[jBp,kBp]
    out = list(Bp = Bp, jretained=jcom, kincl=kincl, kexcl=kexcl)
  }
}


computeSpecCounts = function(Ns, strainCountsA, gssA, IA, BA, abundASVAu){
  # Ns = number of species. suffix "A" can be 16 or 23, depending on which Amplicon.
  # 7 Nov 2023                                                         [cjw]  
  specCounts = rep(0,Ns)
  for (js in 1:Ns){
    nstr = length(strainCountsA[[js]])
    if (nstr>0){
      sabund = rep(0,nstr)
      nstrainOps = gssA[IA[[js]],"No.Ops of Strain"]
      sabund = as.numeric(gssA[IA[[js]],"Total Counts"])/as.numeric(nstrainOps)
      # There can be 
      #     1. a single strain, with one or multiple ASV sets;
      #     2. Multiple strains, all with the same single ASV sets;
      #     3. Multiple strains, with more than 1 ASV set across the set of strains.
      # Case 3 requires determination of a minimal strain set determined such as to meet 
      # the rule set adopted.
      if (BA[[js]]$nASVu==1){
        soln = round(sabund[sample(1:length(BA[[js]]$B),1)])
      } else if (dim(BA[[js]]$B)[1]==1){
        soln = round(sabund)
      } else if (length(unique(gssA[IA[[js]],"ASVsetsCpts"]))==1){
        # Raw counts are the same for all strains as the ASV set is common.  However NumOps 
        # might vary, and hence the cellular abundance.  Randomly choose one.
        iChoice = sample(length(IA[[js]]),1)
        soln = round(sabund[iChoice])
      } else { # minBset may be required.
        if (identical(as.integer(colSums(BA[[js]]$B)),as.integer(rep(1,BA[[js]]$nASVu)))){
          soln = round(sum(sabund))
        } else {
          # Look at the unique rows of BA[[js]]$B to see if any of them is a valid solution.
          Bu = unique(BA[[js]]$B)
          nuniqr = nrow(Bu)
          isSoln = which(sapply(1:nuniqr,function(j){identical(as.integer(Bu[j,]),as.integer(rep(1,BA[[js]]$nASVu)))}))
          if (length(isSoln)>0){
            jsol = which(sapply(1:nrow(BA[[js]]$B),function(j){identical(BA[[js]]$B[j,], Bu[isSoln,])}))
            njsol = length(jsol)
            jsolChoice = sample(1:njsol,1)
            soln = round(sabund[jsol[jsolChoice]])
          } else {
            MB = minBsetSearch(BA[[js]]$B)
            # Returns  list(foundSoln=foundSoln,  solutionSets=solnSets, task=B, Bcompromise=Bsub)
            if (is.null(MB$Bcompromise)){
              nuniqASVCpts = length(unique(gssA[IA[[js]],"ASVsetsCpts"]))
              if (nuniqASVCpts==1){
                soln = round(sabund[unlist(MB$solutionSets)[1]])
              } else {
                kr = sample(1:length(IA[[js]]),1)
                soln = round(sabund[unlist(MB$solutionSets)[kr]])
              }
            } else {
              # If the solution returned is that of Bcompromise then the strain abundances from gssA
              # are not valid - they must be computed based on the re-assignment of ASV sets to strains.
              ciA = MB$solutionSets   #  Identifies the strain indices in i16 that are the strains of the compromise.
              cCounts = sapply(1:length(ciA),function(jk){
                kASVu = which(as.integer(colnames(BA[[js]]$B)) %in% unlist(gssA[IA[[js]][ciA[jk]],"ASVsetsCpts"]))  
                totCounts = MB$Bcompromise[jk,kASVu]*abundASVAu[as.integer(colnames(BA[[js]]$B))[kASVu]]
                numOps = as.integer(gssA[IA[[js]][kASVu],"No.Ops of Strain"])
                sabund = totCounts/numOps  } )
              sumc = 0
              for (jk in 1:length(ciA)){
                sumc = sumc + sum(unlist(cCounts[[jk]]))
              }
              soln = round(sumc)
            }
          }
        }
      }    #  Completion of minBsetSearch use block
      specCounts[js] = soln
    } else {
      specCounts[js] = 0
    }
  }      #    end    js     loop
  out = specCounts 
}

#######################################################################################
############################     RUN INITIALISATION      ##############################
#######################################################################################
possibleGenera = cbind(c("Bacillus","Clostridia","Enterococcus","Escherichia","Lactobacillus",
                   "Listeria", "Pseudomonas","Salmonella","Shigella","Staphylococcus",
                   "Streptococcus","OtherGenus"),
                   c("Bacillus","Clostrid","Enteroco","Escheric","Lactobac",
                   "Listeria", "Pseudomo","Salmonel","Shigella","Staphylo",
                   "Streptoc","OtherGen"))
colnames(possibleGenera) = c("Full","Trunc8")
Ngenera = nrow(possibleGenera)/2

possibleSpecies = cbind(c("subtilis","difficile","faecalis","coli","fermentum",
                    "monocytogenes","aeruginosa","enterica","boydii","aureus",
                    "agalactiae","OtherSpecies"),
                    c("subtilis","difficil","faecalis","coli","fermentu",
                    "monocyto","aerugino","enterica","boydii","aureus",
                    "agalacti","OtherSpe"))
colnames(possibleSpecies) = c("Full","Trunc8")
Nspecies = nrow(possibleSpecies)/2

num_ops_per_species = c(10,14,4,7,5,6,4,7,7,6,6,5)  # Only 1 value per species so not always correct.
# Streptococcus varies markedly with many at 4 (e.g. suis) most common is 6 (e.g. pyogenes, thermophilus)
# Clostridium - most 9,10 but beijerinckii is 14,15,16
Glengths = c(4214810,2542841,2845422,4639221,2011830,2905187,6792191,4755909,4599354,2874302,1852442,3000000)
# Note that species "other" is assigned a genome length of 3000000, so theoretical abundances based on 
# mass proportions of DNA are only approximately correct for this species, if deemed present (with flow-on
# to all other species).

# Alignment quality thresholds (simple) for 16S and 23S rRNA genes. 
# Values c(10,5, 17,8)  have been used for species-level discrimination.
mmthresh16 = 10
gapthresh16 = 5
mmthresh23 = 17
gapthresh23 = 8

basepath4 = "/vast/projects/rrn/microbiome/papercheck"

# Currently working with two possible sets of datasets, where each such set has 2 primary datasets
# and 8 sub-sampled datasets.
# RADnameSet gives the 10 fasta files of ASVs for each set of datasets, 1 fasta file per dataset.
RADnameSet = c("denoised_amplicon_D6322_16S_trimmed_05012024_filtered.fasta",
               "denoised_amplicon_D6322_23S_trimmed_05012024_filtered.fasta",
               "denoised_amplicon_D6322_16S_05012024_Sub31_14012024_filtered.fasta",
               "denoised_amplicon_D6322_23S_05012024_Sub31_14012024_filtered.fasta",
               "denoised_amplicon_D6322_16S_05012024_Sub32_14012024_filtered.fasta",
               "denoised_amplicon_D6322_23S_05012024_Sub32_14012024_filtered.fasta",
               "denoised_amplicon_D6322_16S_05012024_Sub33_14012024_filtered.fasta",
               "denoised_amplicon_D6322_23S_05012024_Sub33_14012024_filtered.fasta",
               "denoised_amplicon_D6322_16S_05012024_Sub34_14012024_filtered.fasta",
               "denoised_amplicon_D6322_23S_05012024_Sub34_14012024_filtered.fasta",
               "denoised_amplicon_D6322_16S_05022024_filtered.fasta",
               "denoised_amplicon_D6322_23S_05022024_filtered.fasta",
               "denoised_amplicon_D6322_16S_05022024_Sub21_06022024_filtered.fasta",
               "denoised_amplicon_D6322_23S_05022024_Sub21_06022024_filtered.fasta",
               "denoised_amplicon_D6322_16S_05022024_Sub22_06022024_filtered.fasta",
               "denoised_amplicon_D6322_23S_05022024_Sub22_06022024_filtered.fasta",
               "denoised_amplicon_D6322_16S_05022024_Sub23_06022024_filtered.fasta",
               "denoised_amplicon_D6322_23S_05022024_Sub23_06022024_filtered.fasta",
               "denoised_amplicon_D6322_16S_05022024_Sub24_06022024_filtered.fasta",
               "denoised_amplicon_D6322_23S_05022024_Sub24_06022024_filtered.fasta")

cat("There are ",length(RADnameSet), " entries for RADnameSet. \n")
verbose = FALSE
thisGenera = c(1,3,4,6,7,8,10)
thisSpecies = c(1,3,4,6,7,8,10)
speciesSet = c("subtilis","faecalis","coli","monocytogenes","aeruginosa","enterica","aureus")
specSet = speciesSet
Ns = length(speciesSet)
# Note that the preferred order of species is alphabetical except that "other" is last.
initalSpeciesOrder = possibleSpecies[thisSpecies,1]
ireqOrd = order(initalSpeciesOrder)
NopsLengths = data.frame(Genus=possibleGenera,Nops=num_ops_per_species, length=Glengths)

in_dateString = "21032024"   # This variable could be removed and replaced by out_dateString everywhere.       
out_dateString = "21032024"     

# Select a species and the fraction of the reads determined to be present in the full
# D63222 dataset that are to be used for this species.
D6322fullDesign = rep(1,7)
sub01Design = c(0.1,1,1,1,1,0.5,0.2)
sub02Design = c(0.5,0.01,0.1,1,1,1,1)   
sub03Design = c(1,0.02,1,1,1,1,0.02)    
sub04Design = c(1,1,0.1,1,1,1,0.02)     

DesignSet = data.frame(D6322fullcheck = D6322fullDesign, sub11=sub01Design,sub12=sub02Design,sub13=sub03Design,sub14=sub04Design,
                       D6322fullcheckA = D6322fullDesign, sub21=sub01Design,sub22=sub02Design,sub23=sub03Design,sub24=sub04Design)
subIDSet = c("D6322","Sub11","Sub12","Sub13","Sub14",
             "D6322A","Sub21","Sub22","Sub23","Sub24")

# Retain D6322full16* as this is the primary dataset on which the paper is based.
# D6322fullcheck16.specCounts (and 23) was the first lot of checking - 5 Jan 2024 - while
# D6322fullcheck16A.specCounts (and 23) was the second lot of checking.  This dataset did analyse
#  the full set of split zymo_hmw_r104 data - that is, 76 split_zymo* file (aa to cx).  It is the 
#  most appropriate dataset - if major revision of the paper is required (e.g. for re-submission) 
#   then all tables should be updated to the outputs from this dataset.
D6322full16.specCounts = c(729 , 1231 , 247 , 992 , 295 , 1273 , 689)   # This vector is ordered by genus alphabetically
D6322full23.specCounts = c(492 , 775 , 579 , 644 , 241 , 934 , 796)     # This vector is ordered by genus alphabetically
D6322fullcheck16.specCounts = c(558, 938, 478, 798, 223, 979, 529)     # This vector is ordered by genus alphabetically    # E.coli count comes up as 198 from sub-sampling code!
D6322fullcheck23.specCounts = c(374, 739, 454, 501, 192, 837, 615)     # This vector is ordered by genus alphabeticall    # E.faecalis, S.enterica come up as 583, 711 from sub-sampling code.
D6322fullcheckA16.specCounts = c(812,1343,658,1103,320,1407,744)       # This vector is ordered by genus alphabetically
D6322fullcheckA23.specCounts = c(548,1026,595,686,265,1258,902)        # This vector is ordered by genus alphabetically
D6322fullcheckB16.specCounts = c(320,744,658,1407,1343,1103,812)       # This vector is ordered by species alphabetically
D6322fullcheckB23.specCounts = c(265,902,595,1258,1026,686,548)        # This vector is ordered by species alphabetically 
# For running sub-sampled datasets without having run the primary datasets within that same run 
# the computed cellular abundances from the primary datasets need to be supplied explicitly. The following
# lines give this data for the RAD outputs of 05/01/2024 and 05/02/2024. These are from the merged 16S and 23S results.
numOpsperSpecies = c(10,4,7,5,4,7,6)   #  Ordered by genus alphabetically
D6322full16M1 = c(0.072, 0.114, 0.088, 0.181, 0.302, 0.172, 0.072)
D6322full23M1 = c(0.075, 0.159, 0.101, 0.187, 0.289, 0.131, 0.058)
D6322fullcheck16M1 = (D6322fullcheck16.specCounts/numOpsperSpecies)/sum((D6322fullcheck16.specCounts/numOpsperSpecies))
D6322fullcheck23M1 = (D6322fullcheck23.specCounts/numOpsperSpecies)/sum((D6322fullcheck23.specCounts/numOpsperSpecies))
D6322fullcheckA16M1 = (D6322fullcheckA16.specCounts/numOpsperSpecies)/sum((D6322fullcheckA16.specCounts/numOpsperSpecies))
D6322fullcheckA23M1 = (D6322fullcheckA23.specCounts/numOpsperSpecies)/sum((D6322fullcheckA23.specCounts/numOpsperSpecies))

logPath = "/stornext/Bioinf/data/lab_speed/cjw/microbiome/scripts/shell"
logName = paste("species_strain_analysis_v10_log_",out_dateString,".out",sep="")
sink(file = file.path(logPath,logName),type="output")
#######################################################################################
  #######################################################################################
  ##########################                               ##############################
  ##########################           MAIN BODY           ##############################
  ##########################                               ##############################
  #######################################################################################
  #######################################################################################
  # PART A: Outputs basic abundance calculation (genus, species, strain) plus ..
  #            outname = paste("strain_operon_multiplicity_All16S23S_",out_dateString,".RData",sep="")
  #            save(file=file.path(outRDatapath,outname),D16,abund16,S16, D23,abund23,S23)
  #         as well as 16S and 23S versions of 
  #            outname = paste("RADdenoiseBM_",paste(unlist(strsplit(inname1,split="_"))[1:6],collapse="_"),"_",
  #                                  out_dateString,".RData",sep="")
  #            save(file=file.path(outRDatapath,outname),D,ab)

for (whichDesign in 5:9){ #  0:4  for 05012024 RAD output from primary D6322
  message1 = paste("\n \n COMMENCED RUN FOR whichDesign =",whichDesign," Date: ", date(),sep="")
  print(message1)
  
  # A set of initialisations follows for the specific 16S and 23S datasets selected by whichDesign.
  # denoisein_dateString  gives the file date of the RAD output which is used as input in 
  # specifying the file names for indices and names that were output from RAD.
  if (whichDesign %in% c(0:4)){
    basepath = ifelse(whichDesign ==0,"/stornext/Bioinf/data/lab_speed/cjw/microbiome/paper_microbiome_strain_check",
                      "/stornext/Bioinf/data/lab_speed/cjw/microbiome/paper_microbiome_strain_check/Sub")
    denoisein_dateString = ifelse(whichDesign ==0,"05012024","14012024") 
  } else if (whichDesign %in% c(5:9)){
    basepath = "/vast/projects/rrn/microbiome/papercheck"
    denoisein_dateString = ifelse(whichDesign ==5,"05022024","06022024")
    # Assumes primary D6322 datasets through RAD on 05022024 and sub-samples through on 06022024  
  }
  
  # Get key information characterising the multiplicity of rRNA operons of the genomes of 
  # the bacterial strains in the reference databases. 
  # Need number of operons for any strains of interest. (IU from  unique_operons_count.R )
  inPath = file.path(basepath,"workDB1")
  inname = "strain_operon_multiplicity_All_22012024.RData"
  load(file=file.path(inPath,inname))   # Loads IU, numUniqs
  # IU[[n]] has structure  list(uniques=outline, PIDmat=M,Scoremat=S).
  
  RADFastaPath = basepath
  RADpath = basepath
  outRDatapath = file.path(basepath,"output/RData")
  plotpath = file.path(basepath,"output/plots")
  textpath = file.path(basepath,"output/text")
  intextpath = basepath
  
  SubsampleTable = data.frame(species = speciesSet, fraction=DesignSet[,whichDesign+1])
  subID = subIDSet[whichDesign+1]
  currentrunID =  subID
  which_inFasta16 = 2*whichDesign+1
  which_inFasta23 = which_inFasta16 + 1
  
  blastStem = "/vast/projects/rrn/microbiome/papercheck/workDB1"
  
  # Completion of this set of initialisations.
  
  # Approximately 500 lines of code follow which, for each of the 16S and 23S datasets for this
  # microbiome design, 
  #     1. aligns all ASVs to the relevant database
  #     2. For each ASV alignment identifies those strains giving the optimal alignment ot that ASV;
  #     3. Brings together, for each such alignment, data that characterises the strain and its alignment.
  for (which_subunit in c("16S","23S")){  
    # A set of initialisations follows for constructing the blastn alignment calls.
    whichDN = ifelse(which_subunit=="16S",which_inFasta16,which_inFasta23)  
    loThresh = ifelse(which_subunit=="16S",1200,2200) 
    # Specify the reference DB of 16S or 23 rRNA gene sequences.
    blastDBname = paste("workDB1_",which_subunit,"_19012024_nameAdj",sep="") 
    blastDBpath =  blastStem  
    DNpath = basepath;  fastaInpath = basepath;  DNnameSet = RADnameSet; DNid = "RAD"
    abind = 3
    maxTargetSeqsA = 50
    maxTargetSeqsB = 50
    # Completion of rRNA gene and alignment initialisations
    
    # Set up to form the two blastn calls which will be made via system2(). The calls differ in having different output formats
    # which allow production of the required inofrmation for down stream processing.
    blastn_options_str1 = ' -outfmt "7 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore score" '
    blastn_options_str2A = paste('  -max_target_seqs ',maxTargetSeqsA,' -num_threads 4',sep="",collapse="")
    blastn_options_str2B = paste('  -max_target_seqs ',maxTargetSeqsB,' -num_threads 4',sep="",collapse="")
    blastn_call_db = paste(" -db ",file.path(blastDBpath,blastDBname),sep="")
    blastn_call_query = paste(" -query ",file.path(fastaInpath,DNnameSet[whichDN]),sep="")
    outnameA = paste("blastn_",DNid,"_dataset_",whichDN,"_",blastDBname,"_",out_dateString,"A.txt",sep="")
    outnameB = paste("blastn_",DNid,"_dataset_",whichDN,"_",blastDBname,"_",out_dateString,"B.txt",sep="")
    blastn_call_outA = paste(" -out", file.path(DNpath,outnameA))
    blastn_call_outB = paste(" -out", file.path(DNpath,outnameB))
    
    argstrA = paste(blastn_call_db,blastn_call_query,blastn_options_str2A,blastn_options_str1,blastn_call_outA,sep="")
    argstrB = paste(blastn_call_db,blastn_call_query,blastn_options_str2B,blastn_call_outB,sep="")
    system2("/stornext/System/data/apps/ncbi-blast/ncbi-blast-2.11.0+/bin/blastn",args=argstrA)
    system2("/stornext/System/data/apps/ncbi-blast/ncbi-blast-2.11.0+/bin/blastn",args=argstrB)
    # Completion of ASV alignments to the rRNA gene database.
    
    # Parsing the A and B forms of output. B only provides the detailed name of the reference.
    # This also includes filtering for unsatisfactory alignments (with ASVs potentially removed 
    # from further consideration).
    # Set mismatch and gap thresholds that will be used to determine whether an optimal alignment
    # is too poor for any strain in the database being considered a match to that ASV.
    
    # Initialisations for the process of parsing alignments output.
    mmthresh = ifelse(which_subunit=="16S",mmthresh16,mmthresh23)
    gapthresh = ifelse(which_subunit=="16S",gapthresh16,gapthresh23)
    # Construct the names of the text files output by the RAD denoiser which are needed for the following analysis.
    if (whichDesign %in% c(0:4)){
      denoise_indices_name = paste("amplicon_D6322_",which_subunit,"_trimmed_05012024_",
                                   "filtered_denoise_indices_out.txt", sep="")
      denoise_names_name = paste("amplicon_D6322_",which_subunit,"_trimmed_05012024_",
                                 "filtered_denoise_names_out.txt", sep="")
    } else if (whichDesign %in% c(5:9)){
      denoise_indices_name = paste("amplicon_D6322_",which_subunit,"_05022024_",currentrunID,"_",denoisein_dateString,
                                   "_filtered_denoise_indices_out.txt", sep="")
      denoise_names_name = paste("amplicon_D6322_",which_subunit,"_05022024_",currentrunID,"_",denoisein_dateString,
                                 "_filtered_denoise_names_out.txt", sep="")
    } else{
      print("Invalid whichDesign value.")
    }
    minAlen = ifelse(which_subunit=="16S",1100,2100)
    # Completion of initialisation for alignments parsing process.
    
    A = readLines(con=file.path(DNpath,outnameA))
    B = readLines(con=file.path(DNpath,outnameB))
    kstartA = which(sapply(1:length(A),function(j){grep("Query:",A[j])})>0)
    kstartB = which(sapply(1:length(B),function(j){grep("Query=",B[j])})>0)
    kstartC = which(sapply(1:length(B),function(j){grep("^>",B[j])})>0)
    nASVs = length(kstartA)
    # How many hits - derived from fileB.
    nhits = numHits(nASVs,kstartB,kstart)
    
    # Characterise each optimal alignment for each ASV identified. This generates
    # the array D16 (or D23).
    stats1 = data.frame(abund=rep(0,nASVs),ID=rep("",nASVs),accver=rep("",nASVs))
    ab = rep(0,nASVs)
    nhitsA = rep(0,nASVs)
    S = vector("list",nASVs)
    removedASVs = rep(0,10);  nlost = 0
    kc = 0
    for (k in 1:nASVs){
      diff = ifelse(k<nASVs,kstartA[k+1]-kstartA[k],length(A)-kstartA[nASVs]+1)
      ab[k] = as.numeric(unlist(strsplit(A[kstartA[k]],split="_"))[3])
      stats2 = getStats2(diff,nhits,A,kstartA,B,kstartB,maxTargetSeqsA=50)
      check = identify_poorly_aligned_ASV(stats2,minAlen)
      if (check){
        nlost = nlost + 1
        removedASVs[nlost] = k
      } else {
        kc = kc + 1
        S[[kc]] = stats2
      }
    }
    removedASVs = removedASVs[1:nlost]
    iab = setdiff(1:length(ab),removedASVs)
    ab = ab[iab]
    nhits = nhits[iab]
    nhitsA = nhitsA[iab]
    if (which_subunit == "16S"){
      S16 = S[1:kc]
      nASVs16 = nASVs  
    } else if (which_subunit == "23S"){
      S23 = S[1:kc]
      nASVs23 = nASVs 
    } 
    cat("\n Have removed ",removedASVs,which_subunit,"ASVs", "from further consideration (original ASV set indexing). \n")
    nASVs = kc
     
    # The following extended body of code provides a strain-based taxonomic profiling procedure. 
    # It replaces an earlier procedure that simply used the first listed strain to which an ASV aligned.
    # Firstly, the 16S case.
    #         Import text files generated from Julia code call of denoise() followed by 
    #         writing of ASV indices to text file.
    if (which_subunit=="16S"){
      S = S16
    } else if (which_subunit=="23S"){
      S = S23
    }
    
    # Construct filenames of the RAD output text files that are to be used. Then read in the 
    # required data.
    stem1 = ifelse((whichDesign %% 5)==0,"filtered_denoise",paste(denoisein_dateString,"_filtered_denoise",sep=""))
    taskSet = c("indices","names","proj","templates")
    whichdataset = whichDesign+1
    dataset = subIDSet[whichdataset]
    for (whichtask in 1:2){
      task = taskSet[whichtask]
      if (whichtask==1){
        if (whichDesign %in% c(0:4)){
          indname = ifelse((whichDesign %% 5)==0,
                               paste("amplicon_D6322",which_subunit,"05012024",stem1,task,"out.txt",sep="_"),
                                 paste("amplicon_D6322",which_subunit,"05012024",dataset,stem1,task,"out.txt",sep="_"))
        } else { 
          indname = ifelse((whichDesign %% 5)==0,
                           paste("amplicon_D6322",which_subunit,"05022024",stem1,task,"out.txt",sep="_"),
                           paste("amplicon_D6322",which_subunit,"05022024",dataset,stem1,task,"out.txt",sep="_"))
        }
        Jindices = readLines(con=file.path(intextpath,indname))
      } else if (whichtask==2){
        if (whichDesign %in% c(0:4)){
          namesname = ifelse((whichDesign %% 5)==0,
                           paste("amplicon_D6322",which_subunit,"05012024",stem1,task,"out.txt",sep="_"),
                           paste("amplicon_D6322",which_subunit,"05012024",dataset,stem1,task,"out.txt",sep="_"))
        } else { 
          namesname = ifelse((whichDesign %% 5)==0,
                           paste("amplicon_D6322",which_subunit,"05022024",stem1,task,"out.txt",sep="_"),
                           paste("amplicon_D6322",which_subunit,"05022024",dataset,stem1,task,"out.txt",sep="_"))
        }
        Jnames = readLines(con=file.path(intextpath,namesname))
      } else if (whichtask==3){
        if (whichDesign %in% c(0:4)){
          umapname = ifelse((whichDesign %% 5)==0,
                           paste("amplicon_D6322",which_subunit,"05012024",stem1,task,"out.txt",sep="_"),
                           paste("amplicon_D6322",which_subunit,"05012024",dataset,stem1,task,"out.txt",sep="_"))
        } else { 
          umapname = ifelse((whichDesign %% 5)==0,
                           paste("amplicon_D6322",which_subunit,"05022024",stem1,task,"out.txt",sep="_"),
                           paste("amplicon_D6322",which_subunit,"05022024",dataset,stem1,task,"out.txt",sep="_"))
        }
        P = read.table(file=file.path(intextpath,umapname),sep="\t")
      } else {
        if (whichDesign %in% c(0:4)){
          templname = ifelse((whichDesign %% 5)==0,
                           paste("amplicon_D6322",which_subunit,"05012024",stem1,task,"out.txt",sep="_"),
                           paste("amplicon_D6322",which_subunit,"05012024",dataset,stem1,task,"out.txt",sep="_"))
        } else { 
          templname = ifelse((whichDesign %% 5)==0,
                           paste("amplicon_D6322",which_subunit,"05022024",stem1,task,"out.txt",sep="_"),
                           paste("amplicon_D6322",which_subunit,"05022024",dataset,stem1,task,"out.txt",sep="_"))
        }
        Templ = readLines(con=file.path(intextpath,templname))
      }
    }
    # Completed input of RAD text-file data that is required.
    
    # The dataframe, D, is constructed to have data on all optimal alignments for each ASV.
    chrh = rep("",10000)
    numh = rep(0,10000)
    D = data.frame(ID=chrh,ASV=chrh,Abund=numh,maxPID=numh,MM=numh,gaps=numh,
                   genus=chrh, species=chrh, strain=chrh, Operon=numh)
    dcount = 0
    for (k in 1:nASVs){
      maxPID = max(as.numeric(S[[k]]$PID[which(as.numeric(S[[k]]$Alen)>loThresh)]))
      nh = length(S[[k]]$PID)
      joptim = which(sapply(1:nh, function(j){
        out=((as.numeric(S[[k]]$PID[j])==maxPID) && (as.numeric(S[[k]]$Alen[j])>loThresh))
      }))
      noptim = length(joptim)   # The number of optimal alignments for this ASV
      genus = rep("",noptim)
      species = rep("",noptim)
      strain = rep("",noptim)
      ID1 = rep("",noptim)
      op = rep("",noptim)
      for (jo in joptim){
        dcount = dcount + 1
        # Check that the alignment meets quality criteria. If not assign to classification "other".
        is.poor = (as.numeric(S[[k]]$MM[1])>mmthresh) || (as.numeric(S[[k]]$Gap[1])>gapthresh)
        if (is.poor){
          genusNames[k] = "otherG"
          speciesNames[k] = "otherS"
          strainNames[k] = "other"
          strainNamesfull[k] = "other"
          genus[jo] = "otherG"
          species[jo] = "otherS"
          strain[jo] = "other"
          op[jo] = 1
        } else {
          strain0 = S[[k]]$ID[jo]
          # Parse the string strain to identify the genus of the strain and to establish
          # the structure of the naming of this strain.  It is assumed that the genus ID
          # is in one of fields 1 to 4 of the split string, strain.
          s1 = unlist(strsplit(strain0,split="_"))
          found=FALSE;  foundg=FALSE;  foundsp=FALSE
          kjo = 0;  kg=0;  ksp=0
          while ((!found) && (kjo<4)) {
            kjo = kjo+1
            if (!foundg){kg = which(possibleGenera == s1[kjo])[1]; foundg = !(is.na(kg))}
            if (!foundsp){ksp = which(possibleSpecies == s1[kjo])[1]; foundsp = !(is.na(ksp))}
            found = foundg && foundsp
 #           cat(kjo,kg,ksp,kg,ksp,foundg,foundsp,found,"\n")
          }
          if (foundg){genus[jo] = possibleGenera[kg]}
          if (foundsp){species[jo] = possibleSpecies[ksp]}  
          if (species[jo]==""){
            cat("\n ASV ",k," Strain ",strain0, foundg, foundsp,kg,ksp,"\n")
          }
          
          ns1 = length(s1)
          if (kjo>2){
            strain[jo] = paste(s1[1:2],sep="_",collapse="_")
          } else if (kjo>1){
            strain[jo] = paste(s1[(kjo+1):(ns1-1)],sep="_",collapse="_")
          } else {
            strain[jo] = paste(s1[(kjo+1):(ns1-1)],sep="_",collapse="_")
          }
          # Operon indexing may be of the form Op8, Op11 or 8, 11.  The following line deals with this variation.
          op[jo] = ifelse(s2c(s1[ns1])[1]=="O",as.numeric(substr(s1[ns1],3,nchar(s1[ns1]))),as.numeric(substr(s1[ns1],1,nchar(s1[ns1]))))
        }
        ID1[jo] = paste(k,jo,sep="_",collapse="_")
        D[dcount,"ID"] = ID1[jo]
        D[dcount,"ASV"] = k
        D[dcount,"Abund"] = ab[k]
        D[dcount,"maxPID"] = maxPID
        D[dcount,"MM"] = S[[k]]$MM[jo]
        D[dcount,"gaps"] = S[[k]]$Gap[jo]
        D[dcount,"genus"] = genus[jo]
        D[dcount,"species"] = species[jo]
        D[dcount,"strain"] = strain[jo]
        D[dcount,"Operon"] = op[jo]
      }
    }   #    end       k in ASVs    loop
    D = D[1:dcount,]
    uspecies = unique(D[,"species"])
    for (uspec in uspecies){
      iu = which(D[,"species"]==uspec)
      if (verbose){print(D[iu,])}
    }
    
    ustrain = unique(D[,"strain"])
    for (ustr in ustrain){
      iust = which(D[,"strain"]==ustr)
      if (verbose){print(D[iust,])}
    }
    outname = paste("RADdenoiseBM_",paste(unlist(strsplit(denoise_indices_name,split="_"))[1:6],collapse="_"),".RData",sep="",collapse="")
    save(file=file.path(outRDatapath,outname),D,ab)
    
    if (which_subunit == "16S"){
      D16 = D;   abund16 = ab
    }
    if (which_subunit == "23S"){
      D23 = D;  abund23 = ab
    }
    if (which_subunit=="16S"){
      iab16 = iab;  nlost16 = nlost; removedASVs16 = removedASVs
    } else {
      iab23 = iab;  nlost23 = nlost; removedASVs23 = removedASVs
    }
  }     #    end   which_subunit    loop
  
  if ((length(D16)>0) && (length(D23)>0)){
    outname = paste("strain_operon_multiplicity_All16S23S_",currentrunID,"_",out_dateString,".RData",sep="")
    save(file=file.path(outRDatapath,outname),
                 D16,abund16,S16,iab16,nlost16,removedASVs16,
                 D23,abund23,S23,iab23,nlost23,removedASVs23)
         
    in_dateString = out_dateString
  }
  D16 = restoreLong(D16)
  D23 = restoreLong(D23)
  
  
  #PART B:  Develop code to call strains and their proportions.
  #           Key data is  IU from  unique_operons_count.R .
  #           D16, D23 are a summaries of the ASV alignments to the DB. 
  #           IU gives multiplicity information on sequences of "operons" from each strain.
  
  # Identify the complete set of sequences in the DB to which the ASVs optimally align.
  matchesID16 = sapply(1:nrow(D16),function(j){out=paste(D16[j,7:10],sep="_",collapse="_")})
  matchesID16.uniq = unique(matchesID16)
  n16uniq = length(matchesID16.uniq)
  matchesID23 = sapply(1:nrow(D23),function(j){out=paste(D23[j,7:10],sep="_",collapse="_")})
  matchesID23.uniq = unique(matchesID23)
  cat("\n How many optimal alignments of ASVs to DB entries are there? \n")
  cat("Total optimal matches alignments for 16S and 23S ",length(matchesID16),length(matchesID23),"\n",
      "Unique matches",length(matchesID16.uniq),length(matchesID23.uniq),"\n")
  # Index the D16, D23 rows according to the unique sequence to which they align.
  matchesID16.iuniq = sapply(1:n16uniq,function(j){which(matchesID16 == matchesID16.uniq[j])})
  matchesID23.iuniq = sapply(1:nrow(D23),function(j){which(matchesID23.uniq == matchesID23[j])})

  # The following code analyses 16S and 23S data independently.  Later work will integrate 
  # the analysis.
  
  # PART B.16S
  
  which_subunit = "16S";  whichDN = which_inFasta16
  # Now associate with each ASV the indices of the unique sequences - that is, find the
  # set of indices for each ASV.  Later we will look for sets that are identical.
  # The ASVs are ordered numerically.  So extracting the last entry in D16 and parsing to get
  # the ASV number gives the total number of ASVs - see nASV16 below.
  nASV16 = as.numeric(unlist(strsplit(D16[nrow(D16),1],split="_"))[1])
  splitID16 = sapply(1:nrow(D16),function(j){out= as.numeric(unlist(strsplit(D16[j,1],split="_"))[1])})
  iuASV16 = vector("list",nASV16)  
  # iuASV16 is a list with each element having the rows of D16 containing data of each operon of each 
  #  strain aligning to a specific ASV.
  ID16u = vector("list",nASV16)
  for (k in 1:nASV16){
    iuASV16[[k]] = which(splitID16==k)
    n1 = length(iuASV16[[k]])
    ID16u[[k]] = sapply(1:n1,function(j){which(matchesID16.uniq == matchesID16[iuASV16[[k]][j]])})
  }
  # This gives indexing into D16 rows via iuASV16
  #   - e.g.   k = 3;  D16[iuASV16[[k]][1]:(iuASV16[[k+1]][1]-1),]  
  #            returns all D16 data on ASV 3

  # An ASV set is the set of ASVs where each ASV has the same pattern of strains 
  # giving an optimal alignment (or alignments).
  # Generate the list defining these ASV sets.
  rset = 1:nASV16
  jc=0
  ASV16u = vector("list",nASV16)
  while (length(rset)>0){
    jc = jc + 1
    x = which(sapply(rset,function(j){identical(ID16u[[j]],ID16u[[rset[1]]])}))
    ASV16u[[jc]] = rset[x]
    rset = setdiff(rset,ASV16u[[jc]])
  }
  nA16u = jc
  ASV16u = ASV16u[1:nA16u]
  if (verbose){
    for (ja in 1:nA16u){
      str1 = D16[which(splitID16==ASV16u[[ja]][1])[1],c("genus","species")]
      cat("\n ",str1$genus,str1$species,":  ASV set ",ja," has ASVs",ASV16u[[ja]],"\n")
      cat("                Abundances:",abund16[ASV16u[[ja]]],"  Total:",sum(abund16[ASV16u[[ja]]]),"\n") 
    }
  }
  
  # Implementing parsimony principle for 16S.
  # Note that ASV-associated strains as given by D are actually (strain, Operon) entities.
  # Firstly, gather the unique strain data for each ASV16u set. 
  # CAVEAT! The alignments are to sequences of all rRNA genes of every strain in the DB.
  #         Thus multiple ASVs may align to different rRNA gene sequences of the same strain.
  #         This may give rise to multiple ASV sets associating with a particular strain.
  #         Also, multiple strains will often align optimally with some ASV sets.
 
  # Determine how many rRNA operons each strain of each species has.
  numOpsStrain16 = vector("list",Ns)
  numOpsStrain16 = sapply(1:Ns,function(js){
    nstrains = length(IU[[2*(js-1)+1]]$uniques)
    nops = sapply(1:nstrains,function(ks){
      out=ifelse(length(IU[[2*(js-1)+1]]$PID[[ks]])>1, nrow(IU[[2*(js-1)+1]]$PID[[ks]]),1)
    })
  })  
  # Completed - note that this information is the same for 23S rRNA data.
  
  # For each ASV set identify the species with which it is associated, and for each strain of that species
  # that gave an optimal alignment specify the index within the strain of the rRNA gene giving the
  # optimal alignment. This information is stored in  Ustrain.
  Ustrain = vector("list",nA16u)
  for (ja in 1:nA16u){
    # For ASV ja (and others in its set of identically-matching ASVs) we have the list of (strain/Ops)
    strOp = matchesID16.uniq[ID16u[[ASV16u[[ja]][1]]]];   nstrOp = length(strOp)
    # Now extract the strain string and Op index for these strains. Also get genus and species for later.
    opnum = rep(0,nstrOp);    strain = rep("",nstrOp)
    
    for (j in 1: nstrOp){
      t1 = unlist(strsplit(strOp[j], split="_"));  nt1 = length(t1)
      if (j==1){genus = t1[1];  species = t1[2]}
      opnum[j] = as.numeric(t1[nt1])
      strain[j] = paste(t1[3:(nt1-1)],sep="_",collapse="_")
    }
    # And reduce consideration to unique strains and the set of operons for each
    ustrain = unique(strain)
    nustr = length(ustrain)
    opsets = vector("list",nustr)
    opnums = rep("",nustr)
    for (j in 1:nustr){
      m = which(strain==ustrain[j])
      opsets[[j]] =  opnum[m]
      opnums[j] = length(m)
    }
    opsetsStr = sapply(1:nustr,function(k){v=opsets[[k]]; vs = paste(v[1:length(v)],sep="_",collapse="_");out=vs})
    Ustrain[[ja]] = data.frame(ustrain=ustrain, opsSets = opsetsStr, numops = opnums)
    #    print(cbind(ustrain,opsetsStr))
  }
  
  # Now, for each species, get information on each strain that has an ASV optimally aligned to it.
  # Consider ASV sets in 2 categories 
  #    1. Those having only a single strain to which they optimally align;
  #    2. Those having multiple strains to which they optimally align;
  # This second category is further divided into those    
  # 
  
  # Gather key data for each ASV set  having only a single strain optimally aligned. This
  # includes total reads associated with the ASVs of the ASV set concerned, and the implied
  # cellular abundance for the particular strain.
  isingle16 = which(sapply(1:nA16u,function(ja){length(Ustrain[[ja]]$ustrain)==1}))
  cat("\n 16S Singleton strains \n")
  for (j in isingle16){print(Ustrain[[j]]$ustrain)}
  singletons = unique(sapply(isingle16,function(k){Ustrain[[k]]$ustrain}))
  # Access the genus and species information on these strains.
  singSpecies16 = rep("",length(isingle16))
  singStrain16 = rep("",length(isingle16))
  numops = rep(0,length(isingle16))
  cabund1 = rep(0,length(isingle16))
  sabund = rep(-1,length(isingle16))
  for (jj in 1:length(isingle16)){
    j = isingle16[jj]   # The indexing of isingle16 is the same as that of ASV16u
    m = length(ASV16u[[j]])  # Gives the number ASVs in this ASV set
    numops[jj] = length(which(D16[,"ASV"]==ASV16u[[j]][1]))   # Only need to deal with the first as the other ASVs have identical strains(s)
#    for (jm in 1:m){
#      print(D16[which(D16[,"ASV"]==ASV16u[[j]][jm]),c("ASV","Abund","genus","species","strain","Operon")])
#    }
    singSpecies16[jj] = paste(D16[which(D16[,"strain"]==Ustrain[[j]]$ustrain)[1],c("genus","species","strain")],sep="_",collapse="_")
    singStrain16[jj] = D16[which(D16[,"strain"]==Ustrain[[j]]$ustrain)[1],"strain"]
    t1 = unlist(strsplit(singStrain16[jj],split="_"))
    shortsingStrain = ifelse(nchar(t1[1])>2,t1[1],t1[2])
    sabund[jj] = sum(abund16[ASV16u[[j]]])
    # Find the species and strain indices into the list that gives the number of operons for this strain.
    # This number is needed for method 1 of cellular abundance calculation.
    js = which(speciesSet==D16[which(D16[,"strain"]==Ustrain[[j]]$ustrain)[1],"species"])
    ks = which(sapply(1:length(IU[[2*(js-1)+1]]$uniques),function(j){
      t1 = unlist(strsplit(IU[[2*(js-1)+1]]$uniques[j],split="   "))[3]
      if (nchar(t1)>nchar(shortsingStrain)){t1 = substr(t1,start=1,stop=nchar(shortsingStrain))}; out = t1}) 
      == shortsingStrain)
    cabund1[jj] = sabund[jj]/numOpsStrain16[[js]][ks]
  } 
  singletons.df = data.frame(IDStr=singSpecies16,Strain=singStrain16,ReadCounts=sabund,
                             CellAbund1=round(cabund1), NumOps=numops,ASV16uSets=isingle16)
  if (verbose){cat("\n 16 singletons.df \n"); print(singletons.df)}
  # Completed singleton data collection 
  

  # The ASV16u sets not yet accounted for are given by
  imult16 = setdiff(1:nA16u,isingle16)   # The index of ASV sets that are not singletons in respect of strains.
  nmult16 = length(imult16)
  # Which ASVs are first in each of these non-singleton ASVsets?
  kASV16u = sapply(1:nmult16,function(j){out=ASV16u[[imult16[j]]][1]})
  multSpecies16 = sapply(kASV16u,function(k){out = paste(D16[which(D16[,"ASV"]==k),c("genus","species")][1,],sep="_",collapse="_")})
  if (verbose){
    cat("\n ASV sets with ambiguous strains  \n")
    print(data.frame(ASVset=imult16,firstASV=kASV16u,Species=multSpecies16))
  }
  SS = vector("list",nmult16)  
  multASVStrD16 = vector("list",nmult16)
  numops = rep(0,nmult16)
  sabund = rep(0,nmult16)
  for (jj in 1:nmult16){
    # For this ASV set, represented by the lowest indexed ASV of this set, do the following:-
    #   1. Identify the unique strains having optimal alignments;
    #   2. For each unique strain determine the number of operons having optimal alignments;
    #   3. Compute the cellular abundance for each unique strain;
    #   4. Identify the genus and species for this set of unique strains (assuming all have the same genus and species);
    #   5. Construct an analogous dataframe to singletons.df for each ASV set
#    j = unresolved[jj]
    j = imult16[jj]
    multASVStrD16[[jj]] = D16[which(D16[,"ASV"]==kASV16u[jj]),c("Abund","genus","species","strain")]
    uniqunresASVStr16 = unique(multASVStrD16[[jj]][,"strain"])
    cat("\n Species in mult resolved part:",multSpecies16[jj],"\n")
    uniqStr = unique(D16[which(D16[,"ASV"]==kASV16u[jj]),"strain"])
    nuniqStr = length(uniqStr)
    sabund[jj] = sum(abund16[ASV16u[[j]]])
    if (verbose){
      cat("Abundance (fragment counts for ASVset) ",sabund[jj],"\n")
      cat("Possible strains are \n")
      print(uniqStr)
    }
    # Determine cellular abundance of each strain, assuming all raw counts belong to the single strain.
    cabund1 = rep(0,nuniqStr) 
    for (jus in 1:nuniqStr){
      # Find the species and strain indices into the list that gives the number of operons for this strain.
      # This number is needed for method 1 of cellular abundance calculation.
      thisStrain = Ustrain[[imult16[jj]]]$ustrain[jus]
      t1 = unlist(strsplit(thisStrain,split="_"))
      shortsingStrain = ifelse(nchar(t1[1])>2,t1[1],t1[2])
      #    shortStrainName = extractShortStrainName(thisStrain)
      js = which(speciesSet==D16[which(D16[,"strain"]==thisStrain)[1],"species"])
      ks = which(sapply(1:length(IU[[2*(js-1)+1]]$uniques),function(j){
        t1=unlist(strsplit(IU[[2*(js-1)+1]]$uniques[j],split="   "))[3]
        shortName = substr(t1,start=1,stop=min(nchar(shortsingStrain),nchar(t1)))
      })
      == shortsingStrain)
      cabund1[jus] = sabund[jj]/numOpsStrain16[[js]][ks]
      print(c(jus,shortsingStrain,js,ks,sabund[jj],numOpsStrain16[[js]][ks]))
    }
    # Number of operons contributing for each strain
    numops = sapply(1:length(uniqStr),function(k){out = length(which(D16[which(D16[,"ASV"]==kASV16u[jj]),"strain"]==uniqStr[k]))})
    if (verbose){
      cat("Number of contributing operons for each strain: ",numops,"\n")
      cat("Cellular Count for each strain of species \n")
           print(cbind(multSpecies16[jj], round(cabund1)))
    }
    multAbund = rep(sabund[jj],nuniqStr)
    S = data.frame(IDStr=rep(multSpecies16[jj],nuniqStr), Strain=uniqunresASVStr16, ReadCounts=multAbund,
                   CellAbund1=round(cabund1), NumOps=numops, ASV16uSets=rep(imult16[jj],nuniqStr))
    SS[[jj]] = S
  }
  if (verbose){
    for (kk in 1:nmult16){cat("\n For species ",multSpecies16[kk],"\n");  print(SS[[kk]])}
  }
  
  
  # Merge the two sources of taxonomic data to give a single record for each strain. As an
  # interim step create a dataframe with all records, with same strains on consecutive records.
  # Add a new column which is the PID value for the ASV from which the record was extracted.
  zz1 = singletons.df
  iord1 = order(zz1[,"IDStr"])
  zz1ord = zz1[iord1,]
  zz2 = SS[[1]]
  for (ir in 2:nmult16){zz2 = rbind(zz2,SS[[ir]])}
  iord2 = order(zz2[,"IDStr"])
  zz2ord = zz2[iord2,]
  zzord = rbind(zz1ord,zz2ord)
  ibad = which(zzord[,"CellAbund1"]==0)
  igood = setdiff(1:nrow(zzord),ibad)
  Allstrains.df = zzord[igood,]
  PIDmax = sapply(1:nrow(Allstrains.df),function(k){r = Allstrains.df[k,"ASV16uSets"]; out=D16[which(D16[,"ASV"]==ASV16u[[r]][1])[1],"maxPID"]})
  Allstrains.df$maxPID = PIDmax 
  #  print(Allstrains.df)
  # Merging multiple records of a single strain.
  Uallstrains = unique(Allstrains.df[,"Strain"])
  nustrains = length(Uallstrains)
  cellularAbund1 = rep(0,nustrains)
  IUgenus_species = sapply(1:(length(IU)%/%2),function(kk){
    paste(unlist(strsplit(IU[[2*(kk-1)+1]]$uniques[1],split="   "))[c(1,2)],sep="_",collapse="_")})
  # Need number of operons for any strains of interest.
  # Now available through  numOpsStrain16 calculation - see above.
  ASVsetsCpts = vector("list",nustrains)
  gss = matrix("",nrow=nustrains,ncol=5)
  strainASVset16 = rep("",nustrains)  
  gs = rep("",nustrains)
  for (j in 1:nustrains){
    strainindex = which(Allstrains.df[,"Strain"]==Uallstrains[j])
    Asets = Allstrains.df[strainindex,"ASV16uSets"]
    X = NULL
    for (kk in 1:length(Asets)){
      xx = sapply(1:length(ASV16u[[as.numeric(Asets[kk])]]),function(j){which(D16[,"ASV"]==ASV16u[[as.numeric(Asets[kk])]][j])})
      if (class(xx)[1]=="list"){xxv = unlist(xx)} else {xxv = as.vector(xx)}
      y = D16[xxv,"ASV"]
      X = append(X,y)
      if (verbose) print(D16[xxv,])
    }
    Xu = unique(X)
    strainASVset16[j] = paste(Xu,sep="_",collapse="_")
    t0 = unlist(strsplit(strainASVset16[j],split="_"))
    k=0;  sub = NULL
    while (length(t0)>0){
      k = k+1
      t1 = setdiff(t0,ASV16u[[k]])
      if (length(t1)<length(t0)){
        sub = append(sub,k)
        t0 = t1
      }
    }
    ASVsetsCpts[[j]] = sub
    uAsets = unique(Asets)
    mP = Allstrains.df[strainindex,"maxPID"]
    umP = unique(mP)
    abundCounts = sum(Allstrains.df[strainindex,"ReadCounts"])    # Allstrains.df[strainindex,"CellAbund"]*Allstrains.df[strainindex,"NumOps"]
    genus = unlist(strsplit(Allstrains.df$IDStr[strainindex[1]],split="_"))[1]
    species = unlist(strsplit(Allstrains.df$IDStr[strainindex[1]],split="_"))[2]
    genus_species = paste(genus,species,sep="_",collapse="_")
    gs[j] = genus_species
    js = which(speciesSet==species)
    # For some species there is no ambiguity of strain that aligns optimally for any ASV set.
    # For these species the raw count is simply the sum of counts over all ASVs for any such ASV sets.
    # So check which of the ASV sets associated with this strain occur only once in Allstrains.df[strainindex,"ASV16uSets"]
    # CAVEAT: I don't understand what I was on about here!
    ASVulist = unique(Allstrains.df[strainindex,"ASV16uSets"])
    if (length(ASVulist) == length(strainindex)){
      abundCounts = sum(Allstrains.df[strainindex,"ReadCounts"])
    } else{ 
      # There will be a multiplicity of abundance counts for any ASV, but for the strain of 
      # current focus there will be only the count from any ASVs that align to it. Hence
      abundCounts = sum(Allstrains.df[strainindex,"ReadCounts"])
    }
    IUstrains = sapply(1:length(IU[[2*(js-1)+1]]$uniques),function(kk){unlist(strsplit(IU[[2*(js-1)+1]]$uniques[kk],split="   "))[3]})
    nIUstrains = length(IUstrains)
    mIU = which(sapply(1:nIUstrains,function(kk){out=!(is.na(str_locate(Uallstrains[j],unlist(strsplit(IUstrains[kk],split="p"))[1])[1]))}))
    StrainTot = sum(abundCounts)
    thisStrain = Ustrain[[imult16[jj]]]$ustrain[jus]
    shortStrainName = extractShortStrainName(Uallstrains[j])
    js = which(speciesSet==D16[which(D16[,"strain"]==Uallstrains[j])[1],"species"])
    ks = which(sapply(1:length(IU[[2*(js-1)+1]]$uniques),function(j){
      t1=unlist(strsplit(IU[[2*(js-1)+1]]$uniques[j],split="   "))[3]
      shortName = substr(t1,start=1,stop=min(8,nchar(t1)))
    }) == shortStrainName)
    numOps =  numOpsStrain16[[js]][ks]
    cellularAbund1[j] = abundCounts/numOps
    gss[j,1] = genus_species
    gss[j,2] = Uallstrains[j]
    gss[j,3:5] = c(StrainTot,numOps,round(cellularAbund1[j]))
  }
  colnames(gss) = c("Species","Strain","Total Counts","No.Ops of Strain","CellularAbundance1")
  gss16 = as.data.frame(gss)
  gss16$ASVsetsCpts = ASVsetsCpts
#  print(gss16[,c(1:4,6)])
  strainsASVset16.df = data.frame(GenusSpecies=gs,UniqStrains=Uallstrains,strainASVs=strainASVset16)
  strainsASVset16.df$ASVsetsCpts=ASVsetsCpts
#  print(strainsASVset16.df)
  
  C16 = getStrainCellCounts(gss16,strainsASVset16.df, abund16) 
  I16 = sapply(1:Ns,function(j){which(C16$StCounts[,"Species"]==(speciesSet[ireqOrd])[j])})
  speciesSet = c("subtilis","faecalis","coli","monocytogenes","aeruginosa","enterica","aureus")
  
  MD16 = data.frame(species=C16$StCounts[,"Species"], strain=C16$StCounts[,"StrainID"],
                    rawcount=C16$StCounts[,"RawCount"])
  zz = sapply(1:length(C16$strainASVs),function(j){C16$strainASVs[[j]]$ASVCpts})
  maxzzlen = max(sapply(1:length(zz),function(j){length(unlist(zz[[j]]))}))
  specASVus16 = matrix(0,nrow=length(zz),ncol=maxzzlen)
  for (j in 1:length(zz)){t1 = unlist(zz[[j]]); specASVus16[j,1:length(t1)]=t1}
  colnames(specASVus16) = paste("ASVu",1:maxzzlen,sep="")
  MD16 = cbind(MD16,specASVus16)
  zz16 = unique(unlist(as.vector(MD16[,4:ncol(MD16)])))
  usedASVus16 = setdiff(zz16,0)  # Just removing 0 from zz16
  allUsed16 = length(setdiff(1:max(usedASVus16),usedASVus16))==0
  totalRawCounts16 = sum(sapply(1:length(usedASVus16),function(j){sum(abund16[ASV16u[[usedASVus16[[j]]]]])}))
  B16 = vector("list",Ns)
  for (j in 1:Ns){
    B16[[j]] = makeBmatrix(0,I16[[j]],specASVus16)
  }
  abundASV16u = sapply(1:length(ASV16u),function(j){tot = sum(abund16[ASV16u[[j]]])})

  # Some proportions calculations and barplot production.
  plotname = paste("strainAnalysis_barplot_",which_subunit,"_",whichDN,"_",currentrunID,"_",out_dateString,".pdf",sep="")
  pdf(file=file.path(plotpath,plotname),paper="A4r")
  # Note that strainCounts16 and specCounts16 (below) are returned in alphabetic order of species.
  # Also note that the count returned is a cellular abundance count, not raw count.
  strainCounts16 = sapply(1:Ns,function(j){as.numeric(C16$StCounts[unlist(I16[[j]]),"Cellabund"])})
  # Calculation of the species count strictly requires solving the minimum set of strains problem for the species.
  specCounts16 = computeSpecCounts(Ns, strainCounts16, gss16, I16, B16, abundASV16u)
  countsM1_sub_16S = round(specCounts16)   # round(0.5*(as.numeric(specCounts16$CellAbundlwr) + as.numeric(specCounts16$CellAbundupp)) )    
  #    The sub01 vector of counts should be : c(7,307,35,98,74,91,23)
  propsM1_sub_16S = countsM1_sub_16S/sum(countsM1_sub_16S)
  ##   [1] "other"         "faecalis"      "coli"          "monocytogenes" "aeruginosa"    "enterica"      "aureus"
  
  # This change is from 25Sept2023.  It is not the original speciesSet set up at run initialisation.
  # massdesign =  SubsampleTable[,"fraction"]/Glengths[thisSpecies];   massdesign =  massdesign/sum(massdesign)
  design =  SubsampleTable[,"fraction"]/sum(SubsampleTable[,"fraction"])
  design2 = design[ireqOrd] 
  if (whichDesign %in% c(0,5,10)){
    cellAbund16_fullD6322M1 = specCounts16   
    D6322full16M1 = cellAbund16_fullD6322M1/sum(cellAbund16_fullD6322M1)
  } 
  if (whichDesign %in% c(0,5,10)){
    cellAbund = 1/Glengths[(thisSpecies)[ireqOrd]]
    Design16DM1 = cellAbund/sum(cellAbund)
  } else {
    Design16DM1 = design2*D6322full16M1/sum(design2*D6322full16M1)
  }
  cat("\n 16S standalone results.\n")
  print(cbind(Design16DM1,propsM1_sub_16S,D6322full16M1, design2,countsM1_sub_16S),digits=2)
  props16 = cbind(Design16DM1,propsM1_sub_16S)
  MA = acomp(props16[,1:2])
  dMA16  = dist(t(MA))
  m = which(cumsum(1:10)==length(dMA16))
  t1 = rep("",m)
  for (mm in 1:m){t1[mm] = as.character(round(10*dMA16[mm])/10)}
  t1s = paste(t1,sep=",  ",collapse=",  ")
  refstr = ifelse(whichDesign==0,"Theory","DesM1")
  Labels16 = make_plot_legend_label2(thisSpecies)
  barplot(props16, col=c(1:7),cex.names=0.8, xlim=c(0,3),
          names.arg=c(refstr,"subM1"),
          legend.text=c(Labels16[sapply(1:Ns,function(j){which(Labels16==(speciesSet[ireqOrd])[j])}),2]),
          sub=paste("(Aitchison distance of subM1 column to DesDM1: ",t1s,")",sep=""),
          main=paste("Sereika D6322 Mock Microbiome and sub-sample ",subID," - 16S",sep=""))
  print(dMA16,digits=3)
  dev.off()
  Allstrains16.df = Allstrains.df
  cat("\n Allstrains16.df \n")
  print(Allstrains16.df)
  outname = paste("strainAnalysis_",which_subunit,"_",whichDN,"_",currentrunID,"_",out_dateString,".RData",sep="")
  save(file=file.path(outRDatapath,outname),gss16,strainsASVset16.df,Allstrains16.df,ASV16u,specCounts16,props16)
  cat("\n Completed 16S stand-alone section.\n")
  # Need more appropriate design figures - either from full D6322 result or assuming design-modified 
  # equal masses of DNA - use reciprocal genome lengths.
  
  
  # PART B.23S
  
  which_subunit = "23S";  whichDN = which_inFasta23
  # Now associate with each ASV the indices of the unique sequences - that is, find the
  # set of indices for each ASV.  Later we will look for sets that are identical.
  # The ASVs are ordered numerically.  So extracting the last entry in D16 and parsing to get
  # the ASV number gives the total number of ASVs - see nASV23 below.
  nASV23 = as.numeric(unlist(strsplit(D23[nrow(D23),1],split="_"))[1])
  splitID23 = sapply(1:nrow(D23),function(j){out= as.numeric(unlist(strsplit(D23[j,1],split="_"))[1])})
  iuASV23 = vector("list",nASV23)  
  # iuASV23 is a list with each element having the rows of D23 containing data of each operon of each 
  #  strain aligning to a specific ASV.
  ID23u = vector("list",nASV23)
  for (k in 1:nASV23){
    iuASV23[[k]] = which(splitID23==k)
    n1 = length(iuASV23[[k]])
    ID23u[[k]] = sapply(1:n1,function(j){which(matchesID23.uniq == matchesID23[iuASV23[[k]][j]])})
  }
  # This gives indexing into D23 rows via iuASV23
  #   - e.g. ja = ASV23u[[3]][3]; D23[iuASV23[[ja]][1]:(iuASV23[[ja+1]][1]-1),]  
  #            returns all D23 data on ASV 46, since ASV23u[[3]][3]=46
  # Also ID23u[[k]] gives indexing into matchesID23.uniq that pulls out a complete ID string for the 
  # DB entries that ASV set, ASV23u[[k]], contains - e.g. matchesID23.uniq[ID23u[[3]]] gives 14 unique 
  # ID strings of the form <genus>_<species>_<strain>_<Op> ( "Bacillus_subtilis_NZ_CP034943_spizizen_ATCC_10")
  
  # Now construct the list defining the unique ASV sets - i.e. the sets of ASVs that have the same set of 
  # strains to which they optimally align.
  rset = 1:nASV23
  jc=0
  ASV23u = vector("list",nASV23)
  while (length(rset)>0){
    jc = jc + 1
    x = which(sapply(rset,function(j){identical(ID23u[[j]],ID23u[[rset[1]]])}))
    ASV23u[[jc]] = rset[x]
    rset = setdiff(rset,ASV23u[[jc]])
  }
  nA23u = jc
  ASV23u = ASV23u[1:nA23u]
  if (verbose){
    for (ja in 1:nA23u){
      cat("\n ",genus,species,":  ASV set ",ASV23u[[ja]],"\n")
      cat("                Abundances:",abund23[ASV23u[[ja]]],"\n") 
    }
  }
  
  # Implementing parsimony principle for 23S.
  # Note that ASV-associated strains as given by D are actually (strain, Operon) entities.
  # Which ASVs have a single strain giving PID=100 alignment?
  # Firstly, gather the unique strain data for each ASV23u set. CAVEAT! The
  # ASV23u sets are actually (strain,Op) unique sets, and hence some strains 
  # appear in more than 1 ASV23u sets.
  
  # Determine how many rRNA operons each strain of each species has.
  numOpsStrain23 = vector("list",Ns)
  numOpsStrain23 = sapply(1:Ns,function(js){
    nstrains = length(IU[[2*(js-1)+1]]$uniques)
    nops = sapply(1:nstrains,function(ks){
      out=ifelse(length(IU[[2*(js-1)+1]]$PID[[ks]])>1, nrow(IU[[2*(js-1)+1]]$PID[[ks]]),1)
    })
  })  
  
  # For each ASV set identify the species with which it is associated, and for each strain of that species
  # that gave  an optimal alignment specify the index within the strain of the rRNA gene giving the
  # optimal alignment. This information is stored in  Ustrain.
  # Also, 
  Ustrain = vector("list",nA16u)
  for (ja in 1:nA23u){
    # For ASV ja (and others in its set of identically-matching ASVs) we have the list of (strain/Ops)
    strOp = matchesID23.uniq[ID23u[[ASV23u[[ja]][1]]]];   nstrOp = length(strOp)
    # Now extract the strain string and Op index for these strains. Also get genus and species for later.
    opnum = rep(0,nstrOp);    strain = rep("",nstrOp)
    
    for (j in 1: nstrOp){
      t1 = unlist(strsplit(strOp[j], split="_"));  nt1 = length(t1)
      if (j==1){genus = t1[1];  species = t1[2]}
      opnum[j] = as.numeric(t1[nt1])
      strain[j] = paste(t1[3:(nt1-1)],sep="_",collapse="_")
    }
    # And reduce consideration to unique strains and the set of operons for each
    ustrain = unique(strain)
    nustr = length(ustrain)
    opsets = vector("list",nustr)
    opnums = rep("",nustr)
    for (j in 1:nustr){
      m = which(strain==ustrain[j])
      opsets[[j]] =  opnum[m]
      opnums[j] = length(m)
    }
    opsetsStr = sapply(1:nustr,function(k){v=opsets[[k]]; vs = paste(v[1:length(v)],sep="_",collapse="_");out=vs})
    Ustrain[[ja]] = data.frame(ustrain=ustrain, opsSets = opsetsStr, numops = opnums)
  }
  
  # Now, for each species, get information on each strain that has an ASV optimally aligned to it.
  # Consider ASV sets in 2 categories 
  #    1. Those having only a single strain to which they optimally align;
  #    2. Those having multiple strains to which they optimally align.
  
  # Gather key data for each ASV set  having only a single strain optimally aligned. This
  # includes total reads associated with the ASVs of the ASV set concerned, and the implied
  # cellular abundance for the particular strain.
  isingle23 = which(sapply(1:nA23u,function(ja){length(Ustrain[[ja]]$ustrain)==1}))
  cat("\n 23S Singleton strains \n")
  for (j in isingle23){print(Ustrain[[j]]$ustrain)}
  singletons = unique(sapply(isingle23,function(k){Ustrain[[k]]$ustrain}))
  # Access the genus and species information on these strains.
  singSpecies23 = rep("",length(isingle23))
  singStrain23 = rep("",length(isingle23))
  numops = rep(0,length(isingle23))
  cabund1 = rep(0,length(isingle23))
  sabund = rep(-1,length(isingle23))
  for (jj in 1:length(isingle23)){
    j = isingle23[jj]   # The indexing of isingle23 is the same as that of ASV23u
    m = length(ASV23u[[j]])  # Gives the number ASVs in this ASV set
    numops[jj] = length(which(D23[,"ASV"]==ASV23u[[j]][1]))   # Only need to deal with the first as the other ASVs have identical strains(s)
    #    for (jm in 1:m){
    #      print(D23[which(D23[,"ASV"]==ASV23u[[j]][jm]),c("ASV","Abund","genus","species","strain","Operon")])
    #    }
    singSpecies23[jj] = paste(D23[which(D23[,"strain"]==Ustrain[[j]]$ustrain)[1],c("genus","species","strain")],sep="_",collapse="_")
    singStrain23[jj] = D23[which(D23[,"strain"]==Ustrain[[j]]$ustrain)[1],"strain"]
    t1 = unlist(strsplit(singStrain23[jj],split="_"))
    shortsingStrain = ifelse(nchar(t1[1])>2,t1[1],t1[2])
    sabund[jj] = sum(abund23[ASV23u[[j]]])
    # Find the species and strain indices into the list that gives the number of operons for this strain.
    # This number is needed for method 1 of cellular abundance calculation.
    js = which(speciesSet==D23[which(D23[,"strain"]==Ustrain[[j]]$ustrain)[1],"species"])
    ks = which(sapply(1:length(IU[[2*(js-1)+1]]$uniques),function(j){
      t1 = unlist(strsplit(IU[[2*(js-1)+1]]$uniques[j],split="   "))[3]
      if (nchar(t1)>nchar(shortsingStrain)){t1 = substr(t1,start=1,stop=nchar(shortsingStrain))}; out = t1}) 
      == shortsingStrain)
    cabund1[jj] = sabund[jj]/numOpsStrain23[[js]][ks]
  } 
  singletons.df = data.frame(IDStr=singSpecies23,Strain=singStrain23,ReadCounts=sabund,
                             CellAbund1=round(cabund1), NumOps=numops,ASV23uSets=isingle23)
  if (verbose){cat("\n 23S singletons.df \n"); print(singletons.df)}
  # Completed singleton data collection 
  
  
  # The ASV23u sets not yet accounted for are given by
  imult23 = setdiff(1:nA23u,isingle23)   # The index of ASV sets that are not singletons in respect of strains.
  nmult23 = length(imult23)
  # Which ASVs are first in each of these non-singleton ASVsets?
  kASV23u = sapply(1:nmult23,function(j){out=ASV23u[[imult23[j]]][1]})
  multSpecies23 = sapply(kASV23u,function(k){out = paste(D23[which(D23[,"ASV"]==k),c("genus","species")][1,],sep="_",collapse="_")})
  if (verbose){
    cat("\n ASV sets with ambiguous strains  \n")
    print(data.frame(ASVset=imult23,firstASV=kASV23u,Species=multSpecies23))
  }
  SS = vector("list",nmult23)  
  multASVStrD23 = vector("list",nmult23)
  numops = rep(0,nmult23)
  sabund = rep(0,nmult23)
  for (jj in 1:nmult23){
    # For this ASV set, represented by the lowest indexed ASV of this set, do the following:-
    #   1. Identify the unique strains having optimal alignments;
    #   2. For each unique strain determine the number of operons having optimal alignments;
    #   3. Compute the cellular abundance for each unique strain;
    #   4. Identify the genus and species for this set of unique strains (assuming all have the same genus and species);
    #   5. Construct an analogous dataframe to singletons.df for each ASV set
    #    j = unresolved[jj]
    j = imult23[jj]
    multASVStrD23[[jj]] = D23[which(D23[,"ASV"]==kASV23u[jj]),c("Abund","genus","species","strain")]
    uniqunresASVStr23 = unique(multASVStrD23[[jj]][,"strain"])
    cat("\n Species in mult resolved part:",multSpecies23[jj],"\n")
    uniqStr = unique(D23[which(D23[,"ASV"]==kASV23u[jj]),"strain"])
    nuniqStr = length(uniqStr)
    sabund[jj] = sum(abund23[ASV23u[[j]]])
    if (verbose){
      cat("Abundance (fragment counts for ASVset) ",sabund[jj],"\n")
      cat("Possible strains are \n")
      print(uniqStr)
    }
    
    cabund1 = rep(0,nuniqStr) 
    for (jus in 1:nuniqStr){
      # Find the species and strain indices into the list that gives the number of operons for this strain.
      # This number is needed for method 1 of cellular abundance calculation.
      thisStrain = Ustrain[[imult23[jj]]]$ustrain[jus]
      t1 = unlist(strsplit(thisStrain,split="_"))
      shortsingStrain = ifelse(nchar(t1[1])>2,t1[1],t1[2])
      #    shortStrainName = extractShortStrainName(thisStrain)
      js = which(speciesSet==D23[which(D23[,"strain"]==thisStrain)[1],"species"])
      ks = which(sapply(1:length(IU[[2*(js-1)+1]]$uniques),function(j){
        t1=unlist(strsplit(IU[[2*(js-1)+1]]$uniques[j],split="   "))[3]
        shortName = substr(t1,start=1,stop=min(nchar(shortsingStrain),nchar(t1)))
      })
      == shortsingStrain)
      cabund1[jus] = sabund[jj]/numOpsStrain23[[js]][ks]
    }
    # Number of operons contributing for each strain
    numops = sapply(1:length(uniqStr),function(k){out = length(which(D23[which(D23[,"ASV"]==kASV23u[jj]),"strain"]==uniqStr[k]))})
    if (verbose){
      cat("Number of contributing operons for each strain: ",numops,"\n")
      cat("Cellular Count for each strain of species \n")
      print(cbind(multSpecies23[jj], round(cabund1)))
    }
    multAbund = rep(sabund[jj],nuniqStr)
    S = data.frame(IDStr=rep(multSpecies23[jj],nuniqStr), Strain=uniqunresASVStr23, ReadCounts=multAbund,
                   CellAbund1=round(cabund1), NumOps=numops, ASV23uSets=rep(imult23[jj],nuniqStr))
    SS[[jj]] = S
  }
  if (verbose){
    for (kk in 1:nmult23){cat("\n For species ",multSpecies23[kk],"\n");  print(SS[[kk]])}
  }
  
  
  # Merge the two sources of taxonomic data to give a single record for each strain. As an
  # interim step create a dataframe with all records, with same strains on consecutive records.
  # Add a new column which is the PID value for the ASV from which the record was extracted.
  zz1 = singletons.df
  iord1 = order(zz1[,"IDStr"])
  zz1ord = zz1[iord1,]
  zz2 = SS[[1]]
  for (ir in 2:nmult23){zz2 = rbind(zz2,SS[[ir]])}
  iord2 = order(zz2[,"IDStr"])
  zz2ord = zz2[iord2,]
  zzord = rbind(zz1ord,zz2ord)
  ibad = which(zzord[,"CellAbund1"]==0)
  igood = setdiff(1:nrow(zzord),ibad)
  Allstrains.df = zzord[igood,]
  PIDmax = sapply(1:nrow(Allstrains.df),function(k){r = Allstrains.df[k,"ASV23uSets"]; out=D23[which(D23[,"ASV"]==ASV23u[[r]][1])[1],"maxPID"]})
  Allstrains.df$maxPID = PIDmax 
  #  print(Allstrains.df)
  # Merging multiple records of a single strain.
  Uallstrains = unique(Allstrains.df[,"Strain"])
  nustrains = length(Uallstrains)
  cellularAbund1 = rep(0,nustrains)
  cellularAbund2 = rep(0,nustrains)
  IUgenus_species = sapply(1:(length(IU)%/%2),function(kk){
    paste(unlist(strsplit(IU[[2*(kk-1)+1]]$uniques[1],split="   "))[c(1,2)],sep="_",collapse="_")})
  # Need number of operons for any strains of interest.
  # Now available through  numOpsStrain16 calculation - see ~line 740 above.
  ASVsetsCpts = vector("list",nustrains)
  gss = matrix("",nrow=nustrains,ncol=5)
  strainASVset23 = rep("",nustrains)  
  gs = rep("",nustrains)
  for (j in 1:nustrains){
    strainindex = which(Allstrains.df[,"Strain"]==Uallstrains[j])
    Asets = Allstrains.df[strainindex,"ASV23uSets"]
    X = NULL
    for (kk in 1:length(Asets)){
      xx = sapply(1:length(ASV23u[[as.numeric(Asets[kk])]]),function(j){which(D23[,"ASV"]==ASV23u[[as.numeric(Asets[kk])]][j])})
      if (class(xx)[1]=="list"){xxv = unlist(xx)} else {xxv = as.vector(xx)}
      y = D23[xxv,"ASV"]
      X = append(X,y)
      if (verbose) print(D23[xxv,])
    }
    Xu = unique(X)
    strainASVset23[j] = paste(Xu,sep="_",collapse="_")
    t0 = unlist(strsplit(strainASVset23[j],split="_"))
    k=0;  sub = NULL
    while (length(t0)>0){
      k = k+1
      t1 = setdiff(t0,ASV23u[[k]])
      if (length(t1)<length(t0)){
        sub = append(sub,k)
        t0 = t1
      }
    }
    ASVsetsCpts[[j]] = sub
    uAsets = unique(Asets)
    mP = Allstrains.df[strainindex,"maxPID"]
    umP = unique(mP)
    abundCounts = sum(Allstrains.df[strainindex,"ReadCounts"])    # Allstrains.df[strainindex,"CellAbund"]*Allstrains.df[strainindex,"NumOps"]
    genus = unlist(strsplit(Allstrains.df$IDStr[strainindex[1]],split="_"))[1]
    species = unlist(strsplit(Allstrains.df$IDStr[strainindex[1]],split="_"))[2]
    genus_species = paste(genus,species,sep="_",collapse="_")
    gs[j] = genus_species
    js = which(speciesSet==species)
    # For some species there is no ambiguity of strain that aligns optimally for any ASV set.
    # For these species the raw count is simply the sum of counts over all ASVs for any such ASV sets.
    # So check which of the ASV sets associated with this strain occur only once in Allstrains.df[strainindex,"ASV"]
    ASVulist = unique(Allstrains.df[strainindex,"ASV"])
    if (length(ASVulist) == length(strainindex)){
      abundCounts = sum(Allstrains.df[strainindex,"ReadCounts"])
    } else{ 
      # There will be a multiplicity of abundance counts for any ASV, but for the strain of 
      # current focus there will be only the count from any ASVs that align to it. Hence
      abundCounts = sum(Allstrains.df[strainindex,"ReadCounts"])
    }
    IUstrains = sapply(1:length(IU[[2*(js-1)+2]]$uniques),function(kk){unlist(strsplit(IU[[2*(js-1)+2]]$uniques[kk],split="   "))[3]})
    nIUstrains = length(IUstrains)
    mIU = which(sapply(1:nIUstrains,function(kk){out=!(is.na(str_locate(Uallstrains[j],unlist(strsplit(IUstrains[kk],split="p"))[1])[1]))}))
    StrainTot = sum(abundCounts)
    countedOps = Allstrains.df[strainindex,"NumOps"]
    ncO = length(countedOps)
    thisStrain = Ustrain[[imult23[jj]]]$ustrain[jus]
    shortStrainName = extractShortStrainName(Uallstrains[j])
    js = which(speciesSet==D23[which(D23[,"strain"]==Uallstrains[j])[1],"species"])
    ks = which(sapply(1:length(IU[[2*(js-1)+1]]$uniques),function(j){
      t1=unlist(strsplit(IU[[2*(js-1)+1]]$uniques[j],split="   "))[3]
      shortName = substr(t1,start=1,stop=min(8,nchar(t1)))
    }) == shortStrainName)
    numOps =  numOpsStrain23[[js]][ks]
    cellularAbund1[j] = abundCounts/numOps
    gss[j,1] = genus_species
    gss[j,2] = Uallstrains[j]
    gss[j,3:5] = c(StrainTot,numOps,round(cellularAbund1[j]))
    colnames(gss) = c("Species","Strain","Total Counts","No.Ops of Strain","CellularAbundance1")
  }
  gss23 = as.data.frame(gss)
  gss23$ASVsetsCpts = ASVsetsCpts
  #  print(gss23[,c(1:4,6)])
  strainsASVset23.df = data.frame(GenusSpecies=gs,UniqStrains=Uallstrains,strainASVs=strainASVset23)
  strainsASVset23.df$ASVsetsCpts=ASVsetsCpts
  #  print(strainsASVset23.df)
  
  C23 = getStrainCellCounts(gss23,strainsASVset23.df, abund23) 
  I23 = sapply(1:Ns,function(j){which(C23$StCounts[,"Species"]==(speciesSet[ireqOrd])[j])})
  speciesSet = c("subtilis","faecalis","coli","monocytogenes","aeruginosa","enterica","aureus")
  
  MD23 = data.frame(species=C23$StCounts[,"Species"], strain=C23$StCounts[,"StrainID"],
                    rawcount=C23$StCounts[,"RawCount"])
  zz = sapply(1:length(C23$strainASVs),function(j){C23$strainASVs[[j]]$ASVCpts})
  maxzzlen = max(sapply(1:length(zz),function(j){length(unlist(zz[[j]]))}))
  specASVus23 = matrix(0,nrow=length(zz),ncol=maxzzlen)
  for (j in 1:length(zz)){t1 = unlist(zz[[j]]); specASVus23[j,1:length(t1)]=t1}
  colnames(specASVus23) = paste("ASVu",1:maxzzlen,sep="")
  MD23 = cbind(MD23,specASVus23)
  zz23 = unique(unlist(as.vector(MD23[,4:ncol(MD23)])))
  usedASVus23 = setdiff(zz23,0)
  allUsed23 = length(setdiff(1:max(usedASVus23),usedASVus23))==0
  totalRawCounts23 = sum(sapply(1:length(usedASVus23),function(j){sum(abund23[ASV23u[[usedASVus23[[j]]]]])}))
  B23 = vector("list",Ns)
  for (j in 1:Ns){
    B23[[j]] = makeBmatrix(0,I23[[j]],specASVus23)
  }
  abundASV23u = sapply(1:length(ASV23u),function(j){tot = sum(abund23[ASV23u[[j]]])})
  
  # Some proportions calculations and barplot production.
  plotname = paste("strainAnalysis_barplot_",which_subunit,"_",whichDN,"_",currentrunID,"_",out_dateString,".pdf",sep="")
  pdf(file=file.path(plotpath,plotname),paper="A4r")
  # Note that strainCounts23 and specCounts23 (below) are returned in alphabetic order of species.
  # Also note that the count returned is a cellular abundance count, not raw count.
  strainCounts23 = sapply(1:Ns,function(j){as.numeric(C23$StCounts[unlist(I23[[j]]),"Cellabund"])})
  # Calculation of the species count strictly requires solving the minimum set of strains problem for the species.
  specCounts23 = computeSpecCounts(Ns, strainCounts23, gss23, I23, B23, abundASV23u)
  countsM1_sub_23S = round(specCounts23) 
  propsM1_sub_23S = countsM1_sub_23S/sum(countsM1_sub_23S)
  ##   [1] "other"         "faecalis"      "coli"          "monocytogenes" "aeruginosa"    "enterica"      "aureus"
  
  design =  SubsampleTable[,"fraction"]/sum(SubsampleTable[,"fraction"])
  design2 = design[ireqOrd]    
  # Gives the permuted indexing to give alphabetical order of species with "other" being moved to end.
  if (whichDesign %in% c(0,5,10)){
    cellAbund23_fullD6322M1 = specCounts23   
    D6322full16M1 = cellAbund23_fullD6322M1/sum(cellAbund23_fullD6322M1)
  } 
  if (whichDesign %in% c(0,5,10)){
    cellAbund = 1/Glengths[(thisSpecies)[ireqOrd]]
    Design23DM1 = cellAbund/sum(cellAbund)
  } else {
    Design23DM1 = design2*D6322full23M1/sum(design2*D6322full23M1)
  }
  cat("\n 23S standalone results.\n")
  print(cbind(Design23DM1,propsM1_sub_23S,D6322full23M1, design2,countsM1_sub_23S),digits=2)
  props23 = cbind(Design23DM1,propsM1_sub_23S)
  props1623 = cbind()
  MA = acomp(props23[,1:2])
  dMA23  = dist(t(MA))
  MA1623 = acomp(cbind(props16[1,2],props23[,2]))
  dMA1623 = dist(t(MA1623))
  m = which(cumsum(1:10)==length(dMA23))
  t1 = rep("",m)
  for (mm in 1:m){t1[mm] = as.character(round(10*dMA23[mm])/10)}
  t1s = paste(t1,sep=",  ",collapse=",  ")
  refstr = ifelse(whichDesign==0,"Theory","DesM1")
  Labels23 = make_plot_legend_label2(thisSpecies)
  barplot(props23, col=c(1:7),cex.names=0.8, xlim=c(0,3),
          names.arg=c(refstr,"subM1"),
          legend.text=c(Labels23[sapply(1:Ns,function(j){which(Labels23==(speciesSet[ireqOrd])[j])}),2]),
          sub=paste("(Aitchison distance of subM1 column to DesDM1: ",t1s,")",sep=""),
          main=paste("Sereika D6322 Mock Microbiome and sub-sample ",subID," - 23S",sep=""))
  cat("\n Aitchison distances: 16S expected/observed), 23S expected/observed), 16S observed/23S observed. \n")
  print(c(dMA16,dMA23,dMA1623),digits=3)
  dev.off()
  Allstrains23.df = Allstrains.df
  cat("\n Allstrains23.df \n")
  print(Allstrains23.df)
  outname = paste("strainAnalysis_",which_subunit,"_",whichDN,"_",currentrunID,"_",out_dateString,".RData",sep="")
  save(file=file.path(outRDatapath,outname),gss23,strainsASVset23.df,Allstrains23.df,ASV23u,specCounts23,props23)
  cat("\n Completed 23S stand-alone section.\n")
  
}     #    end    whichDesign    loop


# PART C: Merging 16S and 23S

cat("\n \n Commencing PART C  - Merging  \n")
# Can we refine by merging across 16S and 23S to reduce strain ambiguity within a species?
plotname = paste("merged_16S23S_proportions_barplots_",out_dateString, ".pdf",sep="")
outname = paste("merged_16S23S_proportions_distances_",out_dateString, ".RData",sep="")
pdf(file=file.path(plotpath,plotname),paper="a4")
dMA16S23S = vector("list",5)  #  5 because of the set length (primary+ 4 sub-samplings per rRNA gene) of whichDesign sets.
propSets = vector("list",5)
countSets = vector("list",5)
ST = vector("list",5)

for (whichDesign in 5:9){
  cat("\n Design Number ",whichDesign,"\n")
  design = DesignSet[,whichDesign+1]
  design2 = design[c(5,7,3,6,2,4,1)]    
  # Gives the permuted indexing to give alphabetical order of species with "other" being moved to end.)
  currentrunID = subIDSet[whichDesign+1]
  
  inname = paste("strain_operon_multiplicity_All16S23S_",currentrunID,"_",in_dateString,".RData",sep="")
  load(file=file.path(outRDatapath,inname))
  # Need the ASV16u and ASV23u objects.  These are available from the following loads.
  for (which_subunit in c("16S","23S")){# Load the gss16 and gss23 datasubID = subIDSet[whichDesign+1]
    which_inFasta16 = 2*(whichDesign)+1
    which_inFasta23 = which_inFasta16 + 1
    whichDN = ifelse(which_subunit=="16S",which_inFasta16,which_inFasta23) 
    inname = paste("strainAnalysis_",which_subunit,"_",whichDN,"_",currentrunID,"_",in_dateString,".RData",sep="")
    load(file=file.path(outRDatapath,inname))   # Loads gss16/23,strainsASVset16/23.df,Allstrains.df
    cat("Load gss16/23,strainsASVset16/23.df,Allstrains.df \n File loaded is   ") # All strains.df not used in remainder of code.
    print(inname)
  }
  abundASV16u = sapply(1:length(ASV16u),function(j){tot = sum(abund16[ASV16u[[j]]])})
  abundASV23u = sapply(1:length(ASV23u),function(j){tot = sum(abund23[ASV23u[[j]]])})
  C16 = getStrainCellCounts(gss16,strainsASVset16.df, abund16)
  #  specCounts16 = extract_species_from_strains(C16) 
  C23 = getStrainCellCounts(gss23,strainsASVset23.df, abund23)
  #  specCounts23 = extract_species_from_strains(C23) 
  # Need to reset   thisSpecies   to the first assignment (in Run Initialisation)
  speciesSet = c("subtilis","faecalis","coli","monocytogenes","aeruginosa","enterica","aureus")
  
  MD16 = data.frame(species=C16$StCounts[,"Species"], strain=C16$StCounts[,"StrainID"],
                    rawcount=C16$StCounts[,"RawCount"])
  zz = sapply(1:length(C16$strainASVs),function(j){C16$strainASVs[[j]]$ASVCpts})
  maxzzlen = max(sapply(1:length(zz),function(j){length(unlist(zz[[j]]))}))
  specASVus16 = matrix(0,nrow=length(zz),ncol=maxzzlen)
  for (j in 1:length(zz)){t1 = unlist(zz[[j]]); specASVus16[j,1:length(t1)]=t1}
  colnames(specASVus16) = paste("ASVu",1:maxzzlen,sep="")
  MD16 = cbind(MD16,specASVus16)
  zz16 = unique(unlist(as.vector(MD16[,4:ncol(MD16)])))
  usedASVus16 = setdiff(zz16,0)
  allUsed16 = length(setdiff(1:max(usedASVus16),usedASVus16))==0
  totalRawCounts16 = sum(sapply(1:length(usedASVus16),function(j){sum(abund16[ASV16u[[usedASVus16[[j]]]]])}))
  
  MD23 = data.frame(species=C23$StCounts[,"Species"], strain=C23$StCounts[,"StrainID"],
                    rawcount=C23$StCounts[,"RawCount"])
  zz = sapply(1:length(C23$strainASVs),function(j){C23$strainASVs[[j]]$ASVCpts})
  maxzzlen = max(sapply(1:length(zz),function(j){length(unlist(zz[[j]]))}))
  specASVus23 = matrix(0,nrow=length(zz),ncol=maxzzlen)
  for (j in 1:length(zz)){t1 = unlist(zz[[j]]); specASVus23[j,1:length(t1)]=t1}
  colnames(specASVus23) = paste("ASVu",1:maxzzlen,sep="")
  MD23 = cbind(MD23,specASVus23)
  zz23 = unique(unlist(as.vector(MD23[,4:ncol(MD23)])))
  usedASVus23 = setdiff(zz23,0)
  allUsed23 = length(setdiff(1:max(usedASVus23),usedASVus23))==0
  totalRawCounts23 = sum(sapply(1:length(usedASVus23),function(j){sum(abund23[ASV23u[[usedASVus23[[j]]]]])}))
  cat("Total raw counts for 16S, 23S are ",totalRawCounts16,totalRawCounts23,"\n")
  if (allUsed16){cat("All ASV16u sets used in raw counts of species. \n")}
  if (allUsed23){cat("All ASV23u sets used in raw counts of species. \n")}
  
  mergeStrain = vector("list",length(thisSpecies))
  AbundCommon = vector("list",length(thisSpecies))
  # As calculated MD16/23 have the rows ordered alphabetically according to genus.  
  # reorderMD16/23 provide reordering according to species, alphabetically.
  reorderMD16=NULL; for (js in 1:Ns){iz = which(MD16[,"species"]==(speciesSet[ireqOrd])[js]);reorderMD16 = append(reorderMD16,iz)}
  reorderMD23=NULL; for (js in 1:Ns){iz = which(MD23[,"species"]==(speciesSet[ireqOrd])[js]);reorderMD23 = append(reorderMD23,iz)}
  k=0
  totrows = 0
  Nsp = length(possibleSpecies[thisSpecies])
  for (jsp in 1:Nsp){
    #   specid = (speciesSet[ireqOrd])[jsp]
    specid = initalSpeciesOrder[order(initalSpeciesOrder)][jsp]
    str16="";  str23 = ""
    # Identify the strains listed for this species in gss16, gss23
    i16 = which(sapply(1:nrow(gss16),function(j){unlist(strsplit(gss16[j,"Species"],split="_"))[2]==specid}))
    ni16 = length(i16)
    i23 = which(sapply(1:nrow(gss23),function(j){unlist(strsplit(gss23[j,"Species"],split="_"))[2]==specid}))
    ni23 = length(i23)
    specCellabund16 = NULL;  specCellabund23 = NULL
    if ((ni16==0) || (ni23==0)){
      cat("\n PROBLEM:  No strains for at least one of the rRNA genes of species ",specid,"\n")
      if (ni16==0){cat("           16S has no strains identified.\n")}
      if (ni23==0){cat("           23S has no strains identified.\n")}
      cat("No merging possible under this condition for this species.\n")
      cat("AbundCommon[[jsp]] values to be those of the earlier independent 16S and 23S analysis.\n")
      if (ni16==0){
        Str16string = ""
        specCellAbund16 = 0
        StrAbund16string = ""
      } else {
        if (ni16>1){
          Str16string=paste(MD16[i16,"strain"],sep="_",collapse="_")
          strainCellAbund16 = gss16[i16,"CellularAbundance1"]
          StrAbund16string = paste(strainCellAbund16,sep="_",collapse="_")
        } else {
          Str16string=MD16[i16,"strain"]
          StrAbund16string = gss16[i16,"CellularAbundance1"]
        }
        specCellAbund16 = specCounts16[jsp]
      }
      if (ni23==0){
        Str23string = ""
        specCellAbund23 = 0
        StrAbund23string = ""
      } else {
        if (ni23>1){
          Str23string=paste(MD23[i23,"strain"],sep="_",collapse="_")
          strainCellAbund23 = gss23[i23,"CellularAbundance1"]
          StrAbund23string = paste(strainCellAbund23,sep="_",collapse="_")
        } else {
          Str23string=MD23[i23,"strain"]
          StrAbund23string = gss23[i23,"CellularAbundance1"]
        } 
        specCellAbund23 = specCounts23[jsp]
      }
      AbundCommon[[jsp]] = data.frame(Species=specid, Strain16=Str16string, Strain23=Str23string,
                                         specCellAbund16= specCellAbund16,  specCellAbund23= specCellAbund23,
                                             Abund16=StrAbund16string,Abund23=StrAbund23string)
    } else { # (ni16 + ni23)>=2 since neither ni16 or ni23 equals zero
      str16 = gss16[i16,"Strain"]
      str23 = gss23[i23,"Strain"]
      # For each 16S and each 23S identified above extract shortened versions of the strain names.
      str16.trunc = sapply(1:ni16,function(j){paste(unlist(strsplit(str16[j],split="_"))[1:2],sep="_",collapse="_")})
      str23.trunc = sapply(1:ni23,function(j){paste(unlist(strsplit(str23[j],split="_"))[1:2],sep="_",collapse="_")})
      if (whichDesign==2){ # An inappropriate workaround, but being done for now(11Aug2023)
        for (jj in 1:ni23){str23.trunc[jj] = ifelse(str23.trunc[jj]=="CP039753p1_NRRL","D6322_4p84MB",str23.trunc[jj])}
      }
      # If one or other of the rRNA genes has residual ambiguity, identify which of the 23S rRNA gene strains matches each of the
      # ambiguous 16S rRNA genes.
      jcommon16 = which(sapply(1:ni16,function(j){out=length(which(str23.trunc==str16.trunc[j]))>0}))
      if (length(jcommon16)==0){
        jcommon23 = jcommon16
      } else {
        jcommon23 = sapply(1:length(jcommon16),function(j){out=which(str23.trunc==str16.trunc[jcommon16[j]])})
      }
      njc = length(jcommon16)
      #      ksp16 = which(sapply(1:Ns,function(j){specCounts16$strainsforSpecies[[j]][1,"Species"]==specid}))
      #      ksp23 = which(sapply(1:Ns,function(j){specCounts23$strainsforSpecies[[j]][1,"Species"]==specid}))
      ksp16 = which(C16$StCounts[,1]==specid)
      ksp23 = which(C23$StCounts[,1]==specid)
      #      u16A = length(setdiff(unique(specCounts16$strainsforSpecies[[ksp16]][,"ASVu_sets"]),"0"))
      #      u23A = length(setdiff(unique(specCounts23$strainsforSpecies[[ksp23]][,"ASVu_sets"]),"0"))
      #      ASVus16 = unique(unlist(setdiff(unique(specCounts16$strainsforSpecies[[ksp16]][,"ASVu_sets"]),"0")))
      #      ASVus23 = unique(unlist(setdiff(unique(specCounts23$strainsforSpecies[[ksp23]][,"ASVu_sets"]),"0")))
      ASVus16=NULL
      for (j in ksp16){ASVus16 = append(ASVus16,unlist(C16$strainASVs[[j]]$ASVCpts))}
      ASVus16 = sort(unique(ASVus16))
      ASVus23=NULL
      for (j in ksp23){ASVus23 = append(ASVus23,unlist(C23$strainASVs[[j]]$ASVCpts))}
      ASVus23 = sort(unique(ASVus23))
      if (length(unique(gss16[i16,"No.Ops of Strain"]))==1){
        specCellabund16 = round(sum(abundASV16u[setdiff(unique(as.vector(specASVus16[i16,])),0)])/as.numeric(gss16[i16,"No.Ops of Strain"][1]))
      }
      if (length(unique(gss23[i23,"No.Ops of Strain"]))==1){
        specCellabund23 = round(sum(abundASV23u[setdiff(unique(as.vector(specASVus23[i23,])),0)])/as.numeric(gss23[i23,"No.Ops of Strain"][1]))
      }
      # Construct the binary matrix whose columns are the ASV16u (or ASV23u) sets,and rows are the strains.
      
      specStrainCellAbund = function(jcommon,iSu,abundASVSuu,specASVSuus){
        # Returns, for a single species and based on the implied amplicon, the species and strains cellular abundance 
        # values after applying the relevant rules for dealing with ambiguity of strains.  If this is 
        # being used as a central part of the 16s-23S merging process, jcommon may be of non-zero length and value.
        # If not, jcommon=0.
        # 4 November 2023                                                          [2023]
        B = makeBmatrix(jcommon,iSu,specASVus16)
      }
      
      B16mB = makeBmatrix(jcommon16,i16,specASVus16);  B16 = B16mB$B
      B23mB = makeBmatrix(jcommon23,i23,specASVus23);  B23 = B23mB$B
      # Create the matrices, or vectors, corresponding to B16 and B23 but only having the common strain(s), 
      # and also the corresponding objects for the non-common strains.
      njc = length(jcommon16)
      # For 16S
      nr16 = dim(B16)[1];   nc16 = dim(B16)[2]
      B16com = matrix(0,nrow=nrow(B16),ncol=ncol(B16)) 
      B16noncom  = B16com
      R16com = B16com
      if (njc>0){
        B16com[jcommon16,] = B16[jcommon16,]
        R16com = B16com
        if (njc>1){
          k1 = which(colSums(B16com)>1)
          for (kk1 in k1){
            R16com[,kk1] = sapply(1:ni16,function(j){B16com[j,kk1]/sum(B16com[,kk1])})
          }
        }
        jnoncom16 = setdiff(1:ni16,jcommon16)
        if (length(jnoncom16)>0){
          B16noncom[jnoncom16,] = B16[jnoncom16,]
        }
        if (njc>1){
          if (nc16>1){
            kcom16 = which(colSums(B16[jcommon16,])>0)
          } else {
            kcom16 = 1
          }
        } else { # njc == 1
          if (nc16>1){
            kcom16=which(B16[jcommon16,]>0)
          } else {
            kcom16 = 1
          }
        }
        knoncom16 = setdiff(1:ni16,kcom16)
      } else {# njc==0
        jnoncom16 = 1:ni16
        B16noncom[jnoncom16,] = B16[jnoncom16,]
        knoncom16 = 1:ncol(B16);  kcom16 = setdiff(1:ni16,knoncom16)
      }
      if (njc>1){
        if (nc16>1){
          kexcl=which(colSums(B16[jcommon16,])>0)
        } else {kexcl=1}
      } else if (njc==1){
        kexcl=ifelse(nc16>1,which(B16[jcommon16,]>0),1)
      } else {kexcl = NULL}
      kincl = setdiff(1:ni16,kexcl)
      for (k in kexcl){B16noncom[,k] = rep(0,ni16)}
      #      B16red = makeReducedBmatrix2(B16,jcommon16)  # Do I need this - or just kexcl= setdiff(1:ni16,which(colSums(B16[jcommon16,])>0))
      #      for (k in B16red$kexcl){B16noncom[,k] = rep(0,ni16)}
      k1 = which(colSums(B16noncom)>1)
      if (length(k1)>0){ # There needs to be a choice of a single strain
        for (kk1 in k1){
          jk1 = which(B16noncom[,kk1]==1)
          if (length(jk1)>1){
            jk1keep = sample(jk1,1)
            B16noncom[setdiff(jk1,jk1keep),kk1] = rep(0,length(jk1)-1)   
          } else {B16noncom[jk1,kk1] = 1}
        }
      } 
      R16 = R16com+B16noncom
      A16 = matrix(rep(sort(setdiff(unique(specASVus16[i16,]),0)),ni16),ncol=length(setdiff(unique(specASVus16[i16,]),0)),byrow=TRUE)
      strainCounts16Contribs = R16*matrix(abundASV16u[A16],nrow=nrow(A16),ncol=ncol(A16))
      strainCounts16 = rowSums(strainCounts16Contribs)
      strainNumOps16 = sapply(i16,function(j){out=as.numeric(gss16[j,4])})
      strainCellAbund16 = round(strainCounts16/strainNumOps16) 
      specCellAbund16 = round(sum(strainCounts16/strainNumOps16))
      strainNames16 = gss16[i16,1:2]
      
      # For 23S
      nr23 = dim(B23)[1];   nc23 = dim(B23)[2]
      B23com = matrix(0,nrow=nrow(B23),ncol=ncol(B23)) 
      B23noncom  = B23com
      R23com = B23com
      if (njc>0){
        B23com[jcommon23,] = B23[jcommon23,]
        R23com = B23com
        if (njc>1){
          k1 = which(colSums(B23com)>1)
          for (kk1 in k1){
            R23com[,kk1] = sapply(1:ni23,function(j){B23com[j,kk1]/sum(B23com[,kk1])})
          }
        }
        jnoncom23 = setdiff(1:ni23,jcommon23)
        if (length(jnoncom23)>0){
          B23noncom[jnoncom23,] = B23[jnoncom23,]
        }
        if (njc>1){
          if (nc23>1){
            kcom23 = which(colSums(B23[jcommon23,])>0)
          } else {
            kcom23 = 1
          }
        } else { # njc == 1
          if (nc23>1){
            kcom23=which(B23[jcommon23,]>0)
          } else {
            kcom23 = 1
          }
        }
        knoncom23 = setdiff(1:ni23,kcom23)
      } else {# njc==0
        jnoncom23 = 1:ni23
        B23noncom[jnoncom23,] = B23[jnoncom23,]
        knoncom23 = 1:ncol(B23);  kcom23 = setdiff(1:ni23,knoncom23)
      }
      if (njc>1){
        if (nc23>1){
          kexcl=which(colSums(B23[jcommon23,])>0)
        } else {kexcl=1}
      } else if (njc==1){
        if (nc23>1){
          kexcl=which(B23[jcommon23,]>0)
        } else {kexcl=1}
      } else {kexcl = NULL}
      kincl = setdiff(1:ni23,kexcl)
      for (k in kexcl){B23noncom[,k] = rep(0,ni23)}
      k1 = which(colSums(B23noncom)>1)
      if (length(k1)>0){ # There needs to be a choice of a single strain
        for (kk1 in k1){
          jk1 = which(B23noncom[,kk1]==1)
          if (length(jk1)>1){
            jk1keep = sample(jk1,1)
            B23noncom[setdiff(jk1,jk1keep),kk1] = rep(0,length(jk1)-1)   
          } else {B23noncom[jk1,kk1] = 1}
        }
      }
      R23 = R23com+B23noncom
      A23 = matrix(rep(sort(setdiff(unique(specASVus23[i23,]),0)),ni23),ncol=length(setdiff(unique(specASVus23[i23,]),0)),byrow=TRUE)
      strainCounts23Contribs = R23*matrix(abundASV23u[A23],nrow=nrow(A23),ncol=ncol(A23))
      strainCounts23 = rowSums(strainCounts23Contribs)
      strainNumOps23 = sapply(i23,function(j){out=as.numeric(gss23[j,4])})
      strainCellAbund23 = round(strainCounts23/strainNumOps23) 
      specCellAbund23 = round(sum(strainCounts23/strainNumOps23))
      strainNames23 = gss23[i23,1:2]
      
      if (length(jcommon16)>0){
        mergeStrain[[jsp]]$shared = list(shared16=gss16[i16[jcommon16],1:2],shared23=gss23[i23[jcommon23],1:2])
      }
      jnon16 = setdiff(i16[which(rowSums(R16)>0)],i16[jcommon16]);   jnon23 = setdiff(i23[which(rowSums(R23)>0)],i23[jcommon23])
      if (length(jnon16)>0){notsh16 = gss16[jnon16,1:2]} else {notsh16=NULL}
      if (length(jnon23)>0){notsh23 = gss23[jnon23,1:2]} else {notsh23=NULL}
      mergeStrain[[jsp]]$notshared = list(notshared16=notsh16, notshared23=notsh23)
      
      if (ni16>1){
        Str16string=paste(MD16[i16,"strain"],sep="_",collapse="_")
        StrAbund16string = paste(strainCellAbund16,sep="_",collapse="_")
      } else {
        Str16string=MD16[i16,"strain"]
        StrAbund16string = strainCellAbund16
      }
      if (ni23>1){
        Str23string=paste(MD23[i23,"strain"],sep="_",collapse="_")
        StrAbund23string = paste(strainCellAbund23,sep="_",collapse="_")
      } else {
        Str23string=MD23[i23,"strain"]
        StrAbund23string = strainCellAbund23
      }
      AbundCommon[[jsp]] = data.frame(Species=specid, Strain16=Str16string, Strain23=Str23string,
                                      specCellAbund16= specCellAbund16,  specCellAbund23= specCellAbund23,
                                      Abund16=StrAbund16string,Abund23=StrAbund23string)
    }       #     end    (ni16>0  &&  ni23>0)    conditional block      
  }         #     end    jsp (species)         loop
  
  cat("\n AbundCommon \n")
  print(AbundCommon)
  # Now construct proportions data for the species cellular abundances as stored in AbundCommon
  # which has used the minimum for each species (11 Oct 2023).
  if (whichDesign %in% c(0,5,10)){
    cellAbund = 1/Glengths[(thisSpecies)[ireqOrd]]
    theoretical16 = cellAbund/sum(cellAbund)
    theoretical23 = theoretical16
  } else {
    theoretical16 = design2*D6322full16M1/sum(design2*D6322full16M1)
    theoretical23 = design2*D6322full23M1/sum(design2*D6322full23M1)
  } 
  cellCount16 = sapply(1:Nsp,function(j){
    t1 = AbundCommon[[j]]$specCellAbund16
    out = ifelse(length(t1)==1,as.integer(t1),as.integer(t1[1]))})
  proportions16 = cbind(theoretical16,cellCount16/sum(cellCount16))
  cellCount23 = sapply(1:Nsp,function(j){
    t1 = AbundCommon[[j]]$specCellAbund23
    out = ifelse(length(t1)==1,as.integer(t1),as.integer(t1[1]))})
  proportions23 = cbind(theoretical23,cellCount23/sum(cellCount23))
  proportions = cbind(proportions16,proportions23)
  MA16 = acomp(proportions16)
  MA23 = acomp(proportions23)
  MA1623 = acomp(proportions)
  dMA16S23S[[whichDesign+1]]  = c(dist(t(MA16)),dist(t(MA23)), dist(t(MA1623)))
  t1 = c(round(100*dMA16S23S[[whichDesign+1]][1])/100,round(100*dMA16S23S[[whichDesign+1]][2])/100)
  t1s = paste(t1,sep=",  ",collapse=",  ")
  bplotRef16 = ifelse(whichDesign==0,"Theory","Des16")
  bplotRef23 = ifelse(whichDesign==0,"Theory","Des23")
  xuplim = 8               #  ifelse(nci==0,5,ifelse(nci==1,8,ifelse(nci==2,11,ifelse(nci==3,15,5+3*nci))))
  names.str = rep("",4)    #  rep("",4+2*(nci-1))
  names.str[1:2] = c(bplotRef16,"16SD")
  names.str[3:4] = c(bplotRef23,"23SD")
  colnames(proportions) = names.str
  Labels23 = make_plot_legend_label2(thisSpecies)
  barplot(as.matrix(proportions), col=c(1:7),cex.names=0.8, xlim=c(0,xuplim),
          names.arg= names.str,
          legend.text=c(Labels23[sapply(1:Ns,function(j){which(Labels23==(speciesSet[ireqOrd])[j])}),2]),
          sub=paste("(Aitchison distances of other columns to Theoretical: ",t1s,")",sep=""),
          main=paste("Merged 16S and 23S for Dataset ",currentrunID,sep="  "))
  
  # The following is a trimmed down version intended for use in a paper to be submitted (9 Dec 2023).
  
  barplot(as.matrix(proportions), col=c(1:7),cex.names=1.2, xlim=c(0,xuplim), main="",
          names.arg= names.str)
  fullNames = c("Bacillus subtilis","Enterococcus faecalis","Escherichia coli",
                             "Listeria monocytogenes","Pseudomonas aeruginosa", 
                              "Salmonella enterica","Stapholycoccus aureus")
  fullNames.species = c("subtilis","faecalis","coli", "monocytogenes","aeruginosa", "enterica","aureus")
  fNOrder = sapply(1:Ns,function(j){which(fullNames.species==(speciesSet[ireqOrd])[j])})
  legend(x="right",legend=fullNames[fNOrder[Ns:1]], text.font= rep(3,7),col=c(Ns:1),pch=rep(15,7))
  propSets[[whichDesign+1]] = proportions
  countSets[[whichDesign+1]] = cbind(cellCount16,cellCount23)
  cat("\n Proportions (ordered alphabetically on species) \n",speciesSet[ireqOrd],"\n")
  print(propSets[[whichDesign+1]],digits=2)
  cat("\n Aitchison distances \n")
  print(dMA16S23S[[whichDesign+1]],digits=3)  
  outname = paste("merged_final_",currentrunID,"_results.RData",sep="")
  save(file=file.path(outRDatapath,outname),propSets,countSets,dMA16S23S,DesignSet,mergeStrain,AbundCommon,
                         C16,C23,MD16,MD23)
}     #    end     whichDesign    loop

dev.off()

cat("Completed Run  \n")

sink()

######################################################################################################################
######################################################################################################################
######################################################################################################################
# The folllowing code is from species_strain_analysis_v03.R  of 16/01/2024 .

# Now bring together Aitchison metric data from single16S, single 23S and merged 16S and merged 23S analyses
# to allow across-datasets assessment of abundance estimates.
s16S = rep(0,30);  s23S = rep(0,30);  m16S = rep(0,30);  m23S = rep(0,30)
chrID = rep("",30); ambigID = data.frame(Pa=chrID, Sa=chrID, Ec=chrID, Se=chrID, Ef=chrID, Lm=chrID, Bs=chrID)
specAmbig16 = matrix(0,nrow=10,ncol=Ns);   specAmbig23 = matrix(0,nrow=10,ncol=Ns)
colnames(specAmbig16) = c("Pa","Sa","Ec","Se","Ef","Lm","Bs")
colnames(specAmbig23) = c("Pa","Sa","Ec","Se","Ef","Lm","Bs")
rnames = NULL
for (whichDesign in 10:14){
  basepath = ifelse(whichDesign==10,"/stornext/Bioinf/data/lab_speed/cjw/microbiome/paper_microbiome_strain_check",
                    "/stornext/Bioinf/data/lab_speed/cjw/microbiome/paper_microbiome_strain_check/Sub")
  outRDatapath = file.path(basepath,"output/RData")
  # Load the single16S/23S data.  
  currentrunID = subIDSet[whichDesign+1]
  zz0= paste(currentrunID,c("_notmerged","_merged   "),sep="")
  zz1 = rbind(zz0[1],zz0[2])
  rnames = rbind(rnames,zz1)
  inname1 = paste("strainAnalysis_proportions_Ametrics_16S23S_",currentrunID,"_",out_dateString,".RData",sep="")
  load(file=file.path(outRDatapath,inname1))
  # Loads props16,props23,dMA16,dMA23,C16,C23,ASV16u, ASV23u
  s16S[whichDesign+1] = dMA16[1];   s23S[whichDesign+1] = dMA23[1]
  inname2 = paste("merged_final_",currentrunID,"_results.RData",sep="")
  load(file=file.path(outRDatapath,inname2))
  # Loads propSets,countSets,dMA16S23S,DesignSet,mergeStrain)
  m16S[whichDesign+1] = dMA16S23S[[whichDesign+1]][1];   m23S[whichDesign+1] = dMA16S23S[[whichDesign+1]][6]
  for (jsp in 1:Ns){
    specid = initalSpeciesOrder[order(initalSpeciesOrder)][jsp]
    str16="";  str23 = ""
    # Identify the strains listed for this species in C16, C23
    i16 = which(C16$StCounts==specid);   ni16 = length(i16)
    i23 = which(C23$StCounts==specid);   ni23 = length(i23)
  
    wD= whichDesign -10
    specAmbig16[2*(wD-1)+1,jsp] = ni16 
    specAmbig23[2*(wD-1)+1,jsp] = ni23 
    specAmbig16[2*(wD-1)+2,jsp] = as.integer(length(mergeStrain[[jsp]]$shared$shared16[,"Strain"]) + length(mergeStrain[[jsp]]$notshared$notshared16[,"Strain"]))
    specAmbig23[2*(wD-1)+2,jsp] = as.integer(length(mergeStrain[[jsp]]$shared$shared23[,"Strain"]) + length(mergeStrain[[jsp]]$notshared$notshared23[,"Strain"]) )
  }
}
rownames(specAmbig16) = rnames
rownames(specAmbig23) = rnames
# Aitchison distances from barplot figures.
A = data.frame(single16S=s16S, single23S=s23S, merged16S=m16S, merged23S=m23S)
#A = data.frame(single16S=c(1.0,0.9,0.7,0.8,1.5),single23S=c(1.1,0.6,0.2,0.3,1.9),
#                merged16S=c(1.04,0.99,0.70,0.90,0.92),merged23S=c(1.08,0.68,0.21,0.26,0.69))
stats1 = matrix(0,nrow=2,ncol=4)
for (k in 1:4){
  stats1[1,k] = mean(A[,k])
  stats1[2,k] = sd(A[,k])
}
colnames(stats1) = c("single16S","single23S","merged16S","merged16S")
rownames(stats1) = c("Mean","StdDev")
print(stats1,digits=2)
t1 = t.test(A[,1],A[,2],alternative ="two.sided")
print(t1)
t2 = t.test(A[,3],A[,4],alternative ="two.sided")
print(t2)
t12 = t.test(c(A[,1],A[,2]),c(A[,3],A[,4]),alternative ="two.sided")
print(t12)

# Strain disambiguation measure of performance.
disAmbig16 = sapply(1:Ns,function(jsp){sum(specAmbig16[seq(from=1,to=10,by=2),jsp] - specAmbig16[seq(from=2,to=10,by=2),jsp]) })
disAmbig23 = sapply(1:Ns,function(jsp){sum(specAmbig23[seq(from=1,to=10,by=2),jsp] - specAmbig23[seq(from=2,to=10,by=2),jsp]) })
# Can normalise by sum across nonmerged for each species.
disAmbig16norm = disAmbig16/sapply(1:Ns,function(jsp){sum(specAmbig16[seq(from=1,to=10,by=2),jsp])})
disAmbig23norm = disAmbig23/sapply(1:Ns,function(jsp){sum(specAmbig23[seq(from=1,to=10,by=2),jsp])})
DisAmbig = data.frame(disAmbig16=disAmbig16,disAmbig23 = disAmbig23, normeddisambig16=disAmbig16norm, normeddisambig23=disAmbig23norm)
rownames(DisAmbig)= c("P.aerugi","S.aureus","E.coli  ","S.enteri","E.faecal","L.monocy","B.subtil")
print(DisAmbig,digits=2)
ST[[whichDesign+1]] = list(A=A, t.test=list(singlediff=t1, mergeddiff=t2, diffSingleMerged=t12),merged=mergeStrain,AbundCommon)


outname = "overall_summary_D6322_mock_microbiome_16S_23S.RData"
save(file=file.path(outRDatapath,outname),propSets,countSets,DesignSet,ST,dMA16S23S,DisAmbig)


include000 = FALSE
if (include000){
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
  
  
  cellAbund_apportion = function(C){
    # C is an n x (m+1) matrix, such that C[1:n,1:m] is a binary matrix
    # whose rows correspond to equally-rated strains and columns to ASV sets.
    # Column m+1 gives the number of operons of the various strains.
    # This function returns the real valued matrix n x m matrix where some of the non-zero
    # non-zero entries of C have been modified to correspond to equivalued cellular 
    # abundances where ASV sets are shared between strains.  This has application
    # to processing of B_{com} matrices.
    #  28 October 2023                                          [cjw]
    isMat = is.matrix(C)
    if (!isMat){ # Necessarily, C must have only 1 row
      Bmod = C[-length(C)]   # remove the last component
    } else {
      n = nrow(C);   m = ncol(C)-1
      B = matrix(as.integer(C[1:n,1:m]),nrow=n, ncol=m); Nops = C[,m+1]
      Bmod = B
      if (is.matrix(B)){
        kmult = which(colSums(B)>1);   nmult = length(kmult)
        for (kkm in 1:nmult){
          km = kmult[kkm]
          jnz = which(B[,km]>0);  nnz = length(jnz)
          totops = sum(as.integer(C[jnz,m+1]))
          Bmod[jnz,km] = as.integer(C[jnz,m+1])/totops
        }
      } else { # B is either a column vector (so multiple strains but only a single ASV set to consider), or a row vector
        if (m==1){
          jnz = which(B>0);  nnz = length(jnz)
          totops = sum(as.integer(C[jnz,m+1]))
          Bmod[jnz] = as.integer(C[jnz,m+1])/totops
        } else{# n==1
          Bmod = B
        }
      }
    }
    out = Bmod
  }
  
  cellAbund_minBset = function(C){
    # This is the complementary function to cellAbund_apportion(C).
    # It uses the same form of input, and returns the same form of output.
    # It differs in that it retains the minimum number of strains to allow
    # full weight of any ASV sets implied by columns of C.
    # C is an n x (m+1) matrix, such that C[1:n,1:m] is a binary matrix
    # whose rows correspond to equally-rated strains and columns to ASV sets.
    # Column m+1 gives the number of operons of the various strains.  
    # This function has application to processing of B_{noncom} matrices.
    # 29 October 2023                                          [cjw]
    
    isMat = is.matrix(C)
    n = nrow(C);   m = ncol(C)-1
    if (!isMat){ # Necessarily, C must have only 1 row
      Bmod = C[-length(C)]   # remove the last component
    } else if ((m==1) || (n==1)) { # B is either a column vector (so multiple strains but only a single ASV set to consider), or a row vector
      if (m==1){ # Column vector. Randomly choose exactly 1 strain
        Bmod = sample(1:n,1)
        strains = Bmod
      } else{# n==1
        Bmod = B
        strains = 1
      }
    } else {
      n = nrow(C);   m = ncol(C)-1
      B = matrix(as.integer(C[1:n,1:m]),nrow=n, ncol=m); Nops = C[,m+1]
      Bmod = B
      if (is.matrix(B)){
        Yp = minBset2(B)
        solved = Yp$solved
        Bmod = Yp$goodnvecs
        strains = Yp$Bindices
      }
    }
    out = Bmod
  }
  
  
###################################################################################################################

restoreLong2 = function(D0){
  # Restores full-length genus and species names of entries in the dataframes D16, D23
  # using the tables generaShortLongMap and speciesShortLongMap.
  # 13 July 2023  
  D = D0
  sapply(1:nrow(D),function(j1){
    m1 = which(generaShortLongMap[,1]==D[j1,"genus"])
    if (length(m1)>0){
      D[j1,"genus"] = generaShortLongMap[m1,2]
    }
    m2 = which(speciesShortLongMap[,1]==D[j1,"species"])
    if (length(m2)>0){
      D[j1,"species"] = speciesShortLongMap[m2,2]
    }
  })
  out = D
}
  
###################################################################################################################

remove_backslash = function(s1){
  # s11 = gsub("\\\\","x",s1)
  # spos = str_locate_all(s1,regex("\\"))
  s1c = s2c(s1)
  s1cn = sapply(1:length(s1c),function(j){charToRaw(s1c[j])})
  ib = which(sapply(1:length(s1c),function(j){identical(s1cn[j],charToRaw("\\"))}))
}

###################################################################################################################

makeShortStrainNames = function(Dg,species){
  Ns = length(Dg)
  shortNames = rep("",Ns)
  for (k in 1:Ns){
    t1 = unlist(strsplit(Dg[k],split="_"))
    m = which(sapply(1:length(t1),function(j){t1[j]==species}))
    if (m==2){shortNames[k] = "D6322"}
    if (m==3){
      nom = t1[1]
      shortNames[k] = substr(nom,start=1,stop=min(8,nchar(nom)))
    }
    if (m==4){
      nom = ifelse(nchar(t1[1])>nchar(t1[2]),t1[1],t1[2])
      shortNames[k] = substr(nom,start=1,stop=min(8,nchar(nom)))
    }
  }
  out = shortNames
}

###################################################################################################################

getSpeciesCount = function(StrData){
  # NOTE: Currently incomplete (16Sept2023)
  # StrData is a two-component list of parameters characterising each strain's cellular abundance determination.
  # From this the species cellular abundance is determined, assuming no ASV is optimally aligned
  # to 2 different strains of one species. (Note that this might need to be revised later).
  #
  # 14 September 2023                                                          [cjw]
  specSet = unique(StrData$StCounts[,"Species"]);    Ns = length(specSet);
  strainSet = StrData$StCounts[,"StrainID"];         Nstr = length(strainSet)
  for (k in 1:Ns){
    species = specSet[k]
    jj = which(sapply(1:nrow(strainsASVset.df),function(j){
      t1=unlist(strsplit(strainsASVset.df[j,"GenusSpecies"],split="_"))[2]
      t1==species} ) ) 
    # First identify the strain that has the most ASVs associating with it.
    m = which.max(sapply(1:length(jj),function(j1){length(StrData[[2]][[jj[j1]]])}))
    notm = setdiff(1:length(jj),m)
    print(StrData[[2]][[jj[m]]])
  }
}

###################################################################################################################

extract_species_from_strains = function(Csu){
  # Returns a lower and upperbound for the cellular abundance of each species identified 
  # through Csu.
  # Also returns the strains considered to be present for each species in making this 
  # abundance abundance estimation.
  # 27 September 2023                                         [cjw]
  specSet = unique(Csu$StCounts[,"Species"])
  Ns = length(specSet)
  specCounts = matrix(0,Ns,3)
  specCellAbund = matrix(0, nrow=Ns, ncol=4)
  strainsforSpecies = vector("list",Ns)   # To contain the set of strains considered to be present, and hence giving the cellular abundances)
  possStrSets = vector("list",Ns)
  for (jsp in 1:Ns){
    species = specSet[jsp]
    specCellAbund[jsp,1] = species
    cat("\n\n Species ",species,"\n")
    isp = which(Csu$StCounts[,"Species"] == species)
    nisp = length(isp)
    strainData = matrix(0,nrow=nisp+1,ncol=5)
    colnames(strainData) = c("Species","StrainID","RawCount","CellAbund","ASVu_sets")
    specCellAbundm = matrix(0,nrow=nisp+1,ncol=4)
    grp = unique(sapply(1:nisp,function(j){as.vector(Csu$strainASVs[[isp[j]]]$ASVCpts)}));  ngrp = length(grp)
    if (ngrp==1){ # Only 1 possible rawCount value
      cat("   Only 1 possible rawCount value. \n")
      rawCount = as.numeric(Csu$StCount[isp,"RawCount"])
      if (nisp==1){ # Only 1 strain for this species
        cat("    Only 1 strain for this species. \n")
        specCellAbund[jsp,2:4] = c(as.numeric(Csu$StCount[isp,"RawCount"]),
                                   as.numeric(Csu$StCount[isp,"Cellabund"]),
                                   as.numeric(Csu$StCount[isp,"Cellabund"]))
        strainData[1,1:4] = Csu$StCounts[isp,c("Species","StrainID","RawCount","Cellabund")]
        ASVustr = paste(unlist(Csu$strainASVs[[isp]]$ASVCpts),sep="_",collapse="_")
        strainData[1,5] = ASVustr
        strainsforSpecies[[jsp]] = strainData
        possStrSets[[jsp]] = Csu$StCounts[isp,"StrainID"]
      } else { # nisp>1 so there are multiple strains possible but total rawCounts is fixed by grp[[1]] elements
        # Identify strains and the number of operons each has.
        cat("  Multiple (",nisp,"}strains, but total rawCount fixed by abundance associated with grp[[1]] - that is ASVunm ",grp[[1]],"\n")
        colnames(specCellAbundm) = c("StrainID","RawCounts","CellAbund(lwr)","CellAbund(upp)")
        lownop = min(as.numeric(Csu$StCounts[isp,"Numops"]));  hinop = min(as.numeric(Csu$StCounts[isp,"Numops"]))
        rawCount = as.numeric(Csu$StCounts[isp[1],"RawCount"])
        for (jnu in 1:nisp){
          specCellAbundm[jnu,1] = Csu$StCounts[isp[jnu],"StrainID"]
          specCellAbundm[jnu,2] = rawCount
          nop.thisStrain = as.numeric(Csu$StCounts[isp[jnu],"Numops"])
          specCellAbundm[jnu,3] = round(10*rawCount/nop.thisStrain)/10
          specCellAbundm[jnu,4] = round(10*rawCount/nop.thisStrain)/10
          strainData[jnu,1:4] = Csu$StCounts[isp[jnu],c("Species","StrainID","RawCount","Cellabund")]
          ASVustr = paste(unlist(Csu$strainASVs[[isp[jnu]]]$ASVCpts),sep="_",collapse="_")
          strainData[jnu,5] = ASVustr
        }
        strainsforSpecies[[jsp]] = strainData
        # All strains identified by isp are allowed, but only 1 at a time.
        # Decision needed. Here we choose to return the lower and upper values of CellAbund 
        # An alternative is to  return the median of the cellular abundance, and the IQR in columns 3 and 4.
        # Note that nisp>1 in this computational block, so the median and IQR are computable.
        specCellAbundm[nisp+1,2] = rawCount
        specCellAbundm[nisp+1,3] = round(10*min(as.numeric(specCellAbundm[1:nisp,3])))/10
        specCellAbundm[nisp+1,4] = round(10*max(as.numeric(specCellAbundm[1:nisp,3])))/10
        pss = c(which.min(as.numeric(specCellAbundm[1:nisp,3])),which.max(as.numeric(specCellAbundm[1:nisp,3])))
        possStrSets[[jsp]] = Csu$StCounts[isp[pss],"StrainID"]
        #     print(specCellAbundm)
        specCellAbund[jsp,2] = specCellAbundm[nisp+1,2]
        specCellAbund[jsp,3] = round(as.numeric(specCellAbundm[nisp+1,3]))
        specCellAbund[jsp,4] = round(as.numeric(specCellAbundm[nisp+1,4]))
      }
    } else { # end    conditional nrgp==1   block
      # So nrgp>1.
      # Cannot have nisp<2 if ngrp>1
      # There may or may not be overlap of strains on ASV16u or ASV23u and - whether overlap or not -
      # there may be multiple strains with the same ASV16u/ASV23u sets but  different numbers of operons.
      cat("    Multiple strains and multiple rawCounts.  \n")
      m1 = matrix(0,nrow=nisp,ncol=ngrp)
      for (jj in 1:nisp){
        t1=Csu$strainASVs[[isp[jj]]]$ASVCpts
        kk = which(sapply(1:ngrp,function(k){setequal(unlist(t1),grp[[k]])}))
        m1[jj,kk]= 1
      } 
      jum1 = which(!(duplicated(m1)));  um1 = m1[jum1,]; nrum1 = nrow(um1)
      i2 = which(sapply(1:ngrp,function(j){length(grp[[j]])>1}))
      whichovl = matrix(0,nrow=nrum1,ncol=ngrp)
      whichovlList = vector("list",length(whichovl))
      if (length(i2)>0){
        # There may be overlaps.
        for (k in 1:length(i2)){
          # For each index in i2 check to see if any other groups are a subset
          i21 = setdiff(which(sapply(1:ngrp,function(j){all(grp[[j]] %in% grp[[i2[k]]])})),i2[k])
          if (length(i21)==0){cat("    No overlap of grp[[",i2[k],"]], which has ASVunm sets ",grp[[i2[k]]],  "\n")}
          else {
            for (kk in 1:length(i21)){
              m =  which(sapply(1:nrum1,function(k){um1[k,i21[kk]]==1}))
              whichovlList[[nrum1*(m-1)+i21[kk]]] = grp[[i21[kk]]]
              whichovl[m,i21[kk]] = 1  # Provides a presence/absence indicator.  The relevant ASVu23(or ASVu16) are stored in whichovlList
              # whichovl[m,i21[kk]] = grp[[i21[kk]]]  # This line only works if length(grp[[i21[kk]]])==1, which cannot be assumed to be generally true.
              cat("    whichovl[",m,",",i21[kk],"] corresponding to ASV<nn>u sets",grp[[i21[kk]]]," overlapping grp[[",i2[k],"]] (=",grp[[i2[k]]],"). \n")
            }   #    end    kk    loop
          }
        }       #    end    k     loop
      }         #    end    length(i2) conditional    block    loop
      is.overlap = sum(whichovl)>0
      if (!is.overlap){    #    
        cat("There is no overlap of the unique ASV sets for species ",species,"\n")
        # Total raw counts is the sum of the counts for each uASVsets.
        # However, to get cellular abundances we must check to see whether the strains that have 
        # the same ASVs have the same number of operons.
        for (jnu in 1:nisp){
          specCellAbundm[jnu,1] = Csu$StCounts[isp[jnu],"StrainID"]
          lownop = min(as.numeric(Csu$StCounts[isp[jnu],"Numops"]));  hinop = min(as.numeric(Csu$StCounts[isp[jnu],"Numops"]))
          rawCount = as.numeric(Csu$StCounts[isp[jnu],"RawCount"])
          specCellAbundm[jnu,2] = rawCount
          specCellAbundm[jnu,3] = round(10*rawCount/hinop)/10
          specCellAbundm[jnu,4] = round(10*rawCount/lownop)/10
          strainData[jnu,1:4] = Csu$StCounts[isp[jnu],c("Species","StrainID","RawCount","Cellabund")]
          ASVustr = paste(unlist(Csu$strainASVs[[isp[jnu]]]$ASVCpts),sep="_",collapse="_")
          strainData[jnu,5] = ASVustr
          grpMatch = which(sapply(1:length(grp),function(j){setequal(grp[[j]],unlist(Csu$strainASVs[[isp[jnu]]]$ASVCpts))}))
        }
        strainsforSpecies[[jsp]] = strainData
        # 
        inz = which(t(whichovl)>0)
        rovl = which(rowSums(whichovl)>0)
        # Form the matrix, m2, whose columns are the individual ASVu sets relevant for this species 
        # while the rows are binary vectors indicating whether a particular ASVu set is present
        # for the strain that row is derived from.
        mASVu = NULL
        for (j in 1:ngrp){ mASVu = append(mASVu,grp[[j]])}
        mASVu = unique(mASVu)
        B = matrix(0,nrow=nrow(m1), ncol=length(mASVu))
        colnames(B) = paste("ASVu",mASVu,sep="_")
        for (j in 1:nisp){
          i1 = which(m1[j,]==1)
          m2ASVu=NULL;  for (jj in 1:length(i1)){ m2ASVu = append(m2ASVu,grp[[i1[jj]]])};  m2ASu = unique(m2ASVu)
          for (kk in 1:length(m2ASVu)){B[j,which(mASVu == m2ASVu[kk])] = 1}  
        }
        print(B)
        # Now process this B matrix to determine a minimal set of strains.
        Y = minBset2(B)   # jt=16;  print(Y[[jt]]$solution[[Y[[jt]]$solution$nK+1]]);  print(Y[[jt]]$solution$task)
        if (Y$is.solved){
          possibleStrainSets = Y$solution[[Y$solution$nK + 1]]  # Each row indexes a possible strain set
          nposs = ifelse(is.vector(possibleStrainSets),1,nrow(possibleStrainSets))
          if (nposs==1){
            specCellAbundm[nisp+1,2] = round(sum(as.numeric(specCellAbundm[possibleStrainSets,2])))
            specCellAbundm[nisp+1,3] = round(10*sum(as.numeric(specCellAbundm[possibleStrainSets,3])))/10
            specCellAbundm[nisp+1,4] = round(10*sum(as.numeric(specCellAbundm[possibleStrainSets,4])))/10
            specCellAbund[jsp,2:4] = specCellAbundm[nisp+1,2:4]
            pss = possibleStrainSets
            possStrSets[[jsp]] = Csu$StCounts[isp[pss],"StrainID"]
            #       print(specCellAbundm)
          } else {
            rawCount = sum(as.numeric(Csu$StCounts[isp[possibleStrainSets[1,]],"RawCount"]))
            # Compute cellular abundances for each such set of possible strains.
            cellab = rep(0,nposs)
            for (k in 1:nposs){
              inds = possibleStrainSets[k,]
              cellab[k] = sum(as.numeric(specCellAbundm[inds,3]))
            }
            specCellAbundm[nisp+1,2] = rawCount
            specCellAbundm[nisp+1,3] = round(10*min(cellab))/10
            specCellAbundm[nisp+1,4] = round(10*max(cellab))/10
            k = which.min(cellab)
            possStrSets[[jsp]] = Csu$StCounts[isp[possibleStrainSets[k,]],"StrainID"]
            #        print(specCellAbundm)
            specCellAbund[jsp,2:4] = specCellAbundm[nisp+1,2:4]
          }    #    end    nposs     loop
        } else {
          cat("No parsimonious solution found for species ",species,"\n")
        }
      }  else {     #    end   !(is.overlap)     conditional block  - and ngrp>1
        # There is overlap to be managed.
        for (jnu in 1:nisp){
          specCellAbundm[jnu,1] = Csu$StCounts[isp[jnu],"StrainID"]
          lownop = min(as.numeric(Csu$StCounts[isp[jnu],"Numops"]));  hinop = min(as.numeric(Csu$StCounts[isp[jnu],"Numops"]))
          rawCount = as.numeric(Csu$StCounts[isp[jnu],"RawCount"])
          specCellAbundm[jnu,2] = rawCount
          specCellAbundm[jnu,3] = round(10*rawCount/hinop)/10
          specCellAbundm[jnu,4] = round(10*rawCount/lownop)/10
          strainData[jnu,1:4] = Csu$StCounts[isp[jnu],c("Species","StrainID","RawCount","Cellabund")]
          ASVustr = paste(unlist(Csu$strainASVs[[isp[jnu]]]$ASVCpts),sep="_",collapse="_")
          strainData[jnu,5] = ASVustr
          grpMatch = which(sapply(1:length(grp),function(j){setequal(grp[[j]],unlist(Csu$strainASVs[[isp[jnu]]]$ASVCpts))}))
        }
        strainsforSpecies[[jsp]] = strainData
        # inz = which(t(whichovl)>0)    - No longer used?
        #rovl = which(rowSums(whichovl)>0)    - No longer used?
        # Form the matrix, m2, whose columns are the individual ASVu sets relevant for this species 
        # while the rows are binary vectors indicating whether a particular ASVu set is present
        # for the strain that row is derived from.
        mASVu = NULL
        for (j in 1:ngrp){ mASVu = append(mASVu,grp[[j]])}
        mASVu = unique(mASVu)
        B = matrix(0,nrow=nrow(m1), ncol=length(mASVu))
        colnames(B) = paste("ASVu",mASVu,sep="_")
        for (j in 1:nisp){
          i1 = which(m1[j,]==1)
          m2ASVu=NULL;  for (jj in 1:length(i1)){ m2ASVu = append(m2ASVu,grp[[i1[jj]]])};  m2ASVu = unique(m2ASVu)
          for (kk in 1:length(m2ASVu)){B[j,which(mASVu == m2ASVu[kk])] = 1}  
        }
        #      print(B)
        # Now process this B matrix to determine a minimal set of strains.
        Y = minBset2(B)   # jt=16;  print(Y[[jt]]$solution[[Y[[jt]]$solution$nK+1]]);  print(Y[[jt]]$solution$task)
        if (Y$is.solved){
          is.mono = Y$solution$nK==1  # If TRUE each of the possibleStrainSets indexes a single strain
          if (is.mono){
            specCellAbundm[nisp+1,2] = specCellAbundm[1,2]
            specCellAbundm[nisp+1,3] = round(10*median(as.numeric(specCellAbundm[Y$solution$mono,3])))/10
            specCellAbundm[nisp+1,4] = round(10*median(as.numeric(specCellAbundm[Y$solution$mono,4])))/10
            #          print(specCellAbundm)
            specCellAbund[jsp,2:4] = specCellAbundm[nisp+1,2:4]
            possStrSets[[jsp]] = Csu$StCounts[isp[Y$solution$mono],"StrainID"]
          } else {
            possibleStrainSets = Y$solution[[Y$solution$nK + 1]]  # Each row indexes a possible strain set
            nposs = nrow(possibleStrainSets)
            sCA = vector("list",nposs)
            # Compute cellular abundances for each such set of possible strains.
            cellab = rep(0,nposs)
            for (jposs in 1:nposs){
              pSS = possibleStrainSets[jposs,]
              rawCount = sum(as.numeric(Csu$StCounts[isp[pSS],"RawCount"]))
              cellab[jposs] = sum(as.numeric(specCellAbundm[pSS,3]))
              specCellAbundm[nisp+1,2] = rawCount
              specCellAbundm[nisp+1,3] = round(10*cellab[jposs] )/10
              specCellAbundm[nisp+1,4] = round(10*cellab[jposs] )/10
              #             print(specCellAbundm)
              sCA[[jposs]] = specCellAbundm
            }
            # Decision needed.  Here we choose to return the lower and upper values of CellAbund across the nposs
            # solution sets.
            # An alternative is to  return the median of the cellular abundance, and the IQR in columns 3 and 4.
            # Note that nisp>1 in this computational block, so the median and IQR are computable.
            # Form the mean of the rawCounts.
            rawcounts = 0
            mincellAbval = rep(0,nposs)
            maxcellAbval = rep(0,nposs)
            for (jposs in 1:nposs){
              rawcounts[jposs] = as.numeric(sCA[[jposs]][nisp+1,2])
              mincellAbval[jposs] = as.numeric(sCA[[jposs]][nisp+1,3])
              maxcellAbval[jposs] = as.numeric(sCA[[jposs]][nisp+1,4])
            }
            specCellAbund[jsp,2] = median(rawcounts)
            specCellAbund[jsp,3] = round(min(mincellAbval))
            specCellAbund[jsp,4] = round(max(maxcellAbval))
            jposs = which.min(mincellAbval)
            pSS = possibleStrainSets[jposs,]
            possStrSets[[jsp]] = Csu$StCounts[isp[pSS],"StrainID"]
          }
        }else  {
          cat("No parsimonious solution found for species ",species,"\n")
          # Need to specify some strain or strains and return specCellAbund
        }    #      end     Y$is.solved   conditional  block
      }     #    end   is.overlap     conditional block
    }       #    end     ngrp>1       conditional block 
  }       #    end       jsp        loop
  iother = which(specSet=="other")
  designSpecOrder = c(setdiff(specSet,"other")[order(setdiff(specSet,"other"))],"other") 
  thisMap = sapply(1:Ns,function(j){which(specSet == designSpecOrder[j])})
  print(specCellAbund[thisMap,]) 
  out = list(Csu.df=data.frame(species=specCellAbund[thisMap,1], specCount = specCellAbund[thisMap,2],
                               CellAbundlwr=specCellAbund[thisMap,3],CellAbundupp=specCellAbund[thisMap,4]),
             strainsforSpecies=strainsforSpecies,StrainsUsed=possStrSets)
}
  
###################################################################################################################

makeBmatrixrandom = function(n,m,rate1){
  delta = 0.5 - rate1
  B = matrix(round(runif(n*m)-delta),nrow=n,ncol=m)
  j0 = which(rowSums(B)==0)
  if (length(j0)>0){
    k1 = round(0.5+m*runif(length(j0)))
    for (j in 1:length(j0)){jj = j0[j]; B[jj,k1[j]] = 1
    }
  }
  k0 = which(colSums(B)==0)
  if (length(k0)>0){
    j1 = round(0.5+n*runif(length(k0)))
    for (k in 1:length(k0)){kk = k0[k]; B[j1[k],kk] = 1}
  }
  out = B
}
  
###################################################################################################################

minrowset2 = function(B,jcom){
  # B is a binary matrix with at least one 1 in each column. jcom is a row index for 
  # one or more rows to be retained.
  # Task is to return the indices - which must include jcom if length(jcom)>0 and j
  # com > 0 -  into rows of B that give a minimal set of rows such that the colSums(Bm)
  # are all >0, where Bm is the matrix  whose rows are this minimal set.
  # 
  # 19 October 2023                                                            [cjw]
  # Check whether jcom solves the problem.
  nr = nrow(B);    nc = ncol(B)
  solved = FALSE
  ncom = length(jcom)
  if (ncom==1){
    if (sum(B[jcom,])==ncol(B)){
      soln = jcom;  solved = TRUE
    } 
  } else {
    if (length(which(colSums(B[jcom,])>0))==ncol(B)){
      soln = jcom;  solved = TRUE
    }
  }
  # Check whether there  is only a single column not accounted for by rows c(jcom,jrid) 
  # and ,if so, process this case directly.
  jrid = NULL
  krm = NULL
  for (k in 1:length(jcom)){
    t1 = setdiff(which(sapply(1:nr,function(j){identical(B[j,],B[jcom[k],])})),jcom[k])
    if (length(t1)>0){jrid = append(jrid,t1)}
    t2  = which(B[jcom[k],]==1)
    if (length(t2)>0){krm = append(krm,t2)}
  }
  krm = sort(unique(krm))
  kBp = setdiff(1:nc,krm)
  jBp = setdiff(1:nr,c(jcom,jrid))
  if (length(kBp)==1){
    jsol = which(B[jBp,kBp]==1)  # There may be multiple rows. but we take the first only.
    if (length(jsol)>1){
      # Take a row that has the fewest number of entries.
      jsolm = which.min(rowSums(B[jsol,]))
      jsol = jsol[jsolm]
    }
    solved = TRUE
    soln = order(unique(c(jcom,jrid,jsol)))
  }
  if (!solved){
    Bred = makeReducedBmatrix2(B,jcom)
    if (is.matrix(Bred$Bp)){ # Are there single rows/strains that solve the problem?
      j1 = which(rowSums(Bred$Bp)==ncol(Bred$Bp))
      if (length(j1)>0){
        solved = TRUE
        soln = c(jcom,Bred$jretained[j1[1]])
      }
    } else { # Bred$Bp is a vector
      if ( length(Bred$Bp>0)>0){
        solved = TRUE
        soln = c(jcom,Bred$jretained)
      } else {
        soln="FAIL"
      }
    } 
    if (!solved){
      Yp = minBset2(Bred$Bp) 
      # Returns  list(goodnvecs=B[Bind,],solved=solved,Bindices=Bind)
      jvec = setdiff(1:nrow(B),jcom)
      if (Yp$solved){
        solved = TRUE
        soln = c(jcom,jvec[Yp$Bindices])
      } else { # So minBset2 could not find a solution.
        cat("Failed with jcom ",jcom, " B \n");  print(B)
        soln = "FAILED with minBset2() call in  minrowset2"
      }
    }
  }
  out = soln
}

###################################################################################################################

solveOverlapB = function(B,jcom){
  # Given the binary matrix B and one or more preferred rows, jcom, create a 
  # modified, minimal rows binary matrix of the same number of columns such 
  # that 1. All columns have a single non-zero entry;
  #      2. Preference in retaining non-zero entries is given to the jcom row(s).
  # 15 October 2023                                                  [cjw]
  ncom = length(jcom)
  j1 = which(colSums(B)==1)
  if (ncom>1){
    jcom1 = which(colSums(B[jcom,])>0)
  } else {jcom1 = which(B[jcom,]==1)}
  jnotcom1 = setdiff(j1,jcom1)   
  # These  jnotcom1  rows must be kept. But any non-zero entries that overlap non-zero entries of B[jcom,]
  # must be set to zero.
}
###################################################################################################################
################################                                                   ################################
################################              CODE DEVELOPMENT AREA                ################################
################################                                                   ################################
###################################################################################################################
# Case 1: An attempt at an exhaustive search for minBset solutions
#  5 Nov 2023
n = 6;  m = 4;   rate1 = 0.2
B = makeBmatrixrandom(n,m,rate1)

minBsetSearch = function(B){
  nr = nrow(B);    nc = ncol(B)
  foundSoln = FALSE
  solnSets = NULL
  Bsub = NULL
  kch = 0
  while ((!foundSoln) && (kch<nr)){
    kch = kch+1
    if (kch==nr){ # Can directly check if using all possible rows is a valid solution or not.
      validAll = identical(colSums(B),rep(1,nc))
      if (validAll){
        foundSoln = TRUE
        out = list(foundSoln=foundSoln,  solutionSets=1:nr, task=B, Bcompromise=Bsub)
      }
    } else {
      nsets = choose(nr,kch)
      resids = vector("list",kch)
      for (j in 1:kch){
        resids[[j]] = j:(nr-kch+j)
      }
      indicesSets = recursive.strainSet3(kch,nr,resids,duplicates=FALSE)
      if (kch==1){
        allowed = which(rowSums(B)==m)
        foundSoln = length(allowed)>0
        if (foundSoln){
          nallowed = length(allowed)
          solnSets = vector("list", nallowed)
          for (ks in 1:nallowed){
            solnSets[[ks]] = allowed[ks]
          }
        }
      } else {
        check = rep(FALSE,nsets)
        for (js in 1:nsets){
          check[js] = identical(colSums(B[indicesSets$goodnvecs[js,],]),rep(1,nc))
        }
        allowed = which(check)
        foundSoln = length(allowed)>0
        if (foundSoln){
          nallowed = length(allowed)
          solnSets = vector("list", nallowed)
          for (ks in 1:nallowed){
            solnSets[[ks]] = indicesSets$goodnvecs[allowed[ks],]
          }
        }
      }
    }
  }
  if (!foundSoln){
    cat("No ideal solution found. \n Compromise solution being generated.\n")
    kch = 1
    while ((!foundSoln) && (kch<nr)){
      kch = kch+1
      if (kch==nr){ # Can directly check if using all possible rows is a valid solution or not.
        Bsub = B
      } else {
        # Form all kch-size sets of rows.  For each set form the colSums(B[,kch-set>,])
        # vector, cvec. If all elements of cvec>0 then a solution will be created for 
        # this value of kch, as follows. Identify any set that has the minimum of cvec
        # components equal to 2, and randomly choose one of these sets. 
        nsets = choose(nr,kch)
        resids = vector("list",kch)
        for (j in 1:kch){
          resids[[j]] = j:(nr-kch+j)
        }
        indicesSets = recursive.strainSet3(kch,nr,resids,duplicates=FALSE)
        check = rep(FALSE,nsets)
        allposSums = NULL;  possibleSets = matrix(0,nrow=100,ncol=kch);  nposSets = 0
        for (js in 1:nsets){
          t1 = colSums(B[indicesSets$goodnvecs[js,],])
          if (length(which(t1>0))==nc){
            allposSums = append(allposSums,js)
            nposSets = nposSets + 1
            possibleSets[nposSets,] = indicesSets$goodnvecs[js,]
          }
        }
        if (!is.null(allposSums)){foundSoln = TRUE}
      }
    }
    
    cat("Solution can be derived for",kch, " strains. \n")
    if (kch<nr){
      possibleSets = possibleSets[1:nposSets]
      # Randomly choose one of these index sets.
      jch = ifelse(nposSets==1,allposSums,sample(allposSums,1))
      Bsub = B[indicesSets$goodnvecs[jch,],]
    }
    cvec = colSums(Bsub)
    ctrim = which(cvec>1)
    BsubScoreMat = matrix(0,nrow=nrow(Bsub),ncol=ncol(Bsub))
    for (jj in 1:nrow(BsubScoreMat)){
      BsubScoreMat[jj,] = sapply(1:ncol(Bsub),function(j){Bsub[jj,j]*(ncol(Bsub)-j)})
    }
    for (jctr in ctrim){
      if ((jctr==1) && (length(which(Bsub[,1]>0)))){
        # Need to resolve ambiguity.  Choose the row of BsubScoreMat that gives the maximum cumsum.
        it0 = which(Bsub[,1]==1)
        im = which.max(sapply(it0,function(j){cumsum(BsubScoreMat[j,])[ncol(Bsub)]}))
      } else{
        it0 = which(Bsub[,jctr]==1)
        im = which.max(sapply(it0,function(j){cumsum(BsubScoreMat[j,])[jctr-1]}))
      }
      it1 = setdiff(it0,it0[im])
      Bsub[it1,jctr] = rep(0,length(it1))
      for (jj in 1:nrow(BsubScoreMat)){
        BsubScoreMat[jj,] = sapply(1:ncol(Bsub),function(j){Bsub[jj,j]*(ncol(Bsub)-j)})
      }
    }
    solnSets = indicesSets$goodnvecs[jch,]
    foundSoln = identical(colSums(Bsub),rep(1,ncol(Bsub)))
  }
  out = list(foundSoln=foundSoln,  solutionSets=solnSets, task=B, Bcompromise=Bsub)
}

recursive.strainSet3 = function(p,m,resids,duplicates=FALSE){
  # Given a list, resids, of length >=m whose components are integer vectors, constructs
  # all possible combinations of n integers where each component of resids provides one
  # integer, component j providing the j^th vector component. If duplicates=FALSE
  # any combination accepted for output has no duplicate elements.
  # The sets of integers that constitute resids must be ordered 
  # 6 November 2023                                                    [cjw]
  if (p==1){
    x = array(1:m,dim=c(m,1))
    soln = list(goodnvecs=x, is.solved=TRUE,Bindices=x)
    colnames(soln$goodnvecs) = "C1"
  } else {
    AA = recursive.strainSet3(p-1, m, resids, duplicates)
    A0 = AA$goodnvecs
    nr = dim(A0)[1]
    if (length(resids[[p]])==1){
      vec = sapply(1:nr,function(j){sample(setdiff(resids[[j]],A0[j,]),1)})
      A1 = cbind(A0,vec)
    } else {
      A1 = A0
      for (j in 2:length(resids[[p]])){A1 = rbind(A1,A0)}
      
      vec0 = NULL 
      for(jj in 1:length(resids[[p]])){
        vec1 = rep(resids[[p]][jj],nr)
        vec0=append(vec0,vec1)
      }
      A1 = cbind(A1,vec0)
    }
    allCases = unique(t(sapply(1:nrow(A1),function(j){sort(A1[j,])})))
    if (!duplicates){
      t0 = which(sapply(1:nrow(allCases),
                        function(j){length(unique(allCases[j,]))==ncol(allCases)} ) )
      t1 = allCases[t0, ]
      soln = list(goodnvecs=t1,is.solved=TRUE,Bindices=t0)
    } else {soln = A1}
    colnames(soln$goodnvecs) = paste("C",1:p,sep="")
  }
  out = soln
}



# Testing Code:-
nfail1 = 0; nfail2 = 0; nsuccess = 0;  ncompromise = 0
AllSolns = vector("list",100);  ksol = 0
for (n in c(4,6,10)){
  for (m in c(n-2,n,n+2)){
    for (rate1 in c(0.1,0.2,0.3,0.4)){
      goodResult = FALSE
      ksol = ksol+1
      B = makeBmatrixrandom(n,m,rate1)
      zz = minBsetSearch(B)
      if (is.null(zz$Bcompromise)){
        if (zz$foundSoln){
          nsoln = length(zz$solutionSets)
          for (jsol in 1:nsoln){
            if (length(zz$solutionSets[[jsol]])<2){
              check = sum(B[unlist(zz$solutionSets[[jsol]]),]) == m
            } else {
              check = sum(colSums(B[unlist(zz$solutionSets[[jsol]]),])) == m
            }
            if (!check){
              nfail1 = nfail1+1
            } else {nsuccess = nsuccess + 1; goodResult=TRUE}
          }
        } else {
          nfail2 = nfail2 + 1
          goodResult = FALSE
          cat("Failed on case ",c(m,n,rate1), " matrix \n");   print(B)
        }
        if (is.null(zz$solutionSets)){
          soln = NULL  
          goodResult = FALSE
        } else {
          soln = B[unlist(zz$solutionSets[[1]]),] 
          goodResult = TRUE
        }
        AllSolns[[ksol]] = list(Solution=zz$solutionSets,Task=B, Result = goodResult, Bcompromise=zz$Bcompromise)
      } else {
        goodResult = TRUE
        ncompromise = ncompromise + 1
        AllSolns[[ksol]] = list(Solution=zz$solutionSets,Task=B, Result = goodResult, Bcompromise=zz$Bcompromise)
      }
    }
  }
}
AllSolns = AllSolns[1:ksol]
cat("Fails: ",nfail2,"InvalidSolutions: ", nfail1,"   Compromises: ",ncompromise,"    Successes: ",nsuccess,"Case: ",c(m,n,rate1),"\n")

###################################################################################################################
# Case 2: Estimate rate constant underlying data distributed according to a Poisson process.
N = 4
distrib = data.frame(S16 = c(0.572, 0.296, 0.074, 0.040, 0.018),S23 = c(0.615, 0.220, 0.116, 0.037, 0.008))
lambda = rep(0,ncol(distrib));   error_rate = lambda
Len = c(1500, 2500)
count = 0:N
for (k in 1:ncol(distrib)){
  s1 = sum(distrib[,k])
  s2 = sum(sapply(1:length(count), function(j){out = count[j]*distrib[j,k]}))
  lambda[k] = s2/s1
  error_rate[k] = lambda[k]/Len[k]
}
print(cbind(lambda, error_rate),digits=2)




###################################################################################################################
################################                                                   ################################
################################                CODE DUMPING AREA                  ################################
################################                                                   ################################
###################################################################################################################

apportionBmat = function(B,jspSub,jcom){
  # Identify those columns contributing to the jcom strain(s). 
  # If any other strains have a non-zero entry for such a column, set that strain's vector entry to zero.
  # For any columns still summing to >1 split the non-zero values equally to add to 1.
  # 22 October 2023                                                     [cjw]
  Bmod = B
  ncom = length(jcom)
  jnotcom = setdiff(1:nrow(B),jcom)
  if (ncom==1){ # Only 1 common strain so no apportioning over common strains.  
    # Remove any non-zero entries in strains other than the common strain for those ASVs the common strain has.
    # Presume jnotcom is of length at least 1
    kBcom = which(B[jcom,]==1)
    Bmod[jnotcom,kBcom] = rep(0,length(jnotcom)*length(kBcom))
    # Check column sums over the jnotcom rows.  Apportion as necessary.
    if (length(jnotcom)>1){ # It is possible that apportioning is required - not possible if length(jnotcom)==1.
      kk1 = colSums(B[jnotcom,])
      kg1 = which(kk1>1)
      if (length(kg1)>0){
        for (k2 in kg1){B[jnotcom,k2] = B[jnotcom,k2]/kk1[k2]}
      }
    }
  } else { 
    # Remove any non-zero entries in strains other than the common strain for those ASVs the common strain has.
    if (length(jnotcom)>0){
      for (jj in 1:ncom){
        kBcom = which(B[jcom[jj],]==1) 
        if (length(jnotcom)>0){
          Bmod[jnotcom,kBcom] = rep(0,length(jnotcom)*length(kBcom))
        }
      }
    }
    # Apportion over common strains for those columns for which column sum is greater than 1.
    csum = colSums(Bmod[jspSub,])
    kb = which(csum>1)
    if (length(kb)>0){
      for (kk in 1:length(kb)){
        Bmod[jspSub,kb[kk]] = Bmod[jspSub,kb[kk]]/csum[kb[kk]]
      }
    }
  }
  out = Bmod
}

###################################################################################################################



###################################################################################################################



###################################################################################################################



###################################################################################################################

###########################################################################################
###########################################################################################
###########################################################################################
# Set of calls to get material for summary of each dataset.
# From this the material in the diary entry for 11 January 2024 was prepared.
whichDesign = 10   # Or 11, 12, 13, 14
currentrunID = "D6322" 
for (which_subunit in c("16S","23S")){
  whichDN = 20+(whichDesign-10)*2 + ifelse(which_subunit=="16S",1,2)
  out_dateString = "24012024";  in_dateString = out_dateString
  if ((whichDesign %% 5)==0){
    outRDatapath = "/stornext/Bioinf/data/lab_speed/cjw/microbiome/paper_microbiome_strain_check/output/RData" 
  } else {
    outRDatapath = "/stornext/Bioinf/data/lab_speed/cjw/microbiome/paper_microbiome_strain_check/Sub/output/RData" 
  }# Add /Sub  for subsampled datasets
  
  
  outname = paste("strain_operon_multiplicity_All16S23S_",currentrunID,"_",out_dateString,".RData",sep="")
  load(file=file.path(outRDatapath,outname))
  # This gives
  #     D16,JD16,abund16,S16,ASVtoSpeciesMap16,ASVtoStrainMap16,iab16,nlost16,
  #     D23,JD23,abund23,S23,ASVtoSpeciesMap23,ASVtoStrainMap23,iab23,nlost23
  
  #species = "faecalis"
  #print(unique(D16[which(D16[,"species"]==species),c("ASV","strain","Abund")]))
  #print(unique(D23[which(D23[,"species"]==species),c("ASV","strain","Abund")]))
  
  # Now load gss16, strainsASVset16.df, Allstrains16.df, ASV16u, specCounts16, props16 (or 23S equivalents)
  inname = paste("strainAnalysis_",which_subunit,"_",whichDN,"_",currentrunID,"_",out_dateString,".RData",sep="")
  if (which_subunit=="16S"){
    load(file=file.path(outRDatapath,inname))
  } else {
    load(file=file.path(outRDatapath,inname))
  }
  # With D16 and ASV16u and D23 ASV23u all the data required for unmerged or merged strain identification and strain and species abundance
  # is available.
}
print(strainsASVset16.df)
print(strainsASVset23.df)
  
load("/stornext/Bioinf/data/lab_speed/cjw/microbiome/paper_microbiome_strain_check/output/RData/merged_final_D6322_results.RData")
load("/stornext/Bioinf/data/lab_speed/cjw/microbiome/paper_microbiome_strain_check/output/RData/strainAnalysis_16S_21_D6322_22012024.RData")
# Returns  gss16,strainsASVset16.df,  Allstrains16.df,  ASV16u,  specCounts16,  props16
load("/stornext/Bioinf/data/lab_speed/cjw/microbiome/paper_microbiome_strain_check/output/RData/strainAnalysis_23S_22_D6322_22012024.RData")
# Returns  gss23,strainsASVset23.df,  Allstrains23.df,  ASV23u,  specCounts23,  props23


print(DesignSet[,whichDesign+1])
print(mergeStrain)
print(c(length(ASV16u),length(ASV23u)))
print(dMA16S23S[[11]],digits=2)
print(dMA16,digits=2);    print(dMA23,digits=2)
rownames(countSets[[11]]) = sort(possibleSpecies[thisSpecies])
print(countSets[[11]])
print(ASV16u)
print(ASV23u)
print(gss16)
print(gss23)

# Merging outcomes.
#  Load below returns  propSets,  countSets,  dMA16S23S,  DesignSet, mergeStrain, AbundCommon 
inname1 = paste("merged_final_",currentrunID,"_results.RData",sep="")
load(file=file.path(outRDatapath,inname1))
# Note AbundCommon contents:-
# AbundCommon[[jsp]] = data.frame(Species=specid, Strain16=Str16string, Strain23=Str23string,
#                         specCellAbund16= specCellAbund16,  specCellAbund23= specCellAbund23,
#Abund16=StrAbund16string,Abund23=StrAbund23string)
#  Load below returns  propSets,  countSets,  DesignSet,  ST,  dMA16S23S,  DisAmbig 
inname_merged = "overall_summary_D6322_mock_microbiome_16S_23S.RData"
load(file.path(outRDatapath,inname_merged))






#############################################################################################################
#############################################################################################################
#   Repeat Check of /vast/projects/rrn/microbiome/papercheck/output/(RData,plots) Analyses 9-10 March 2024  #
#############################################################################################################
#############################################################################################################
for (whichDesign in c(5:9)){
  message1 = paste("\n \n COMMENCED RUN FOR whichDesign =",whichDesign," Date: ", date(),sep="")
  print(message1)
  
  # A set of initialisations follows for the specific 16S and 23S datasets selected by whichDesign.
  # denoisein_dateString  gives the file date of the RAD output which is used as input in 
  # specifying the file names for indices and names that were output from RAD.
  if (whichDesign %in% c(0:4)){
    basepath = ifelse(whichDesign ==0,"/stornext/Bioinf/data/lab_speed/cjw/microbiome/paper_microbiome_strain_check",
                      "/stornext/Bioinf/data/lab_speed/cjw/microbiome/paper_microbiome_strain_check/Sub")
    denoisein_dateString = ifelse(whichDesign ==0,"05012024","14012024") 
  } else if (whichDesign %in% c(5:9)){
    basepath = "/vast/projects/rrn/microbiome/papercheck"
    denoisein_dateString = ifelse(whichDesign ==5,"05022024","06022024")
    # Assumes primary D6322 datasets through RAD on 05022024 and sub-samples through on 06022024  
  }
  
  # Get key information characterising the multiplicity of rRNA operons of the genomes of 
  # the bacterial strains in the reference databases. 
  # Need number of operons for any strains of interest. (IU from  unique_operons_count.R )
  inPath = file.path(basepath,"workDB1")
  inname = "strain_operon_multiplicity_All_22012024.RData"
  load(file=file.path(inPath,inname))   # Loads IU, numUniqs
  # IU[[n]] has structure  list(uniques=outline, PIDmat=M,Scoremat=S).
  
  RADFastaPath = basepath
  RADpath = basepath
  outRDatapath = file.path(basepath,"output/RData")
  plotpath = file.path(basepath,"output/plots")
  textpath = file.path(basepath,"output/text")
  intextpath = basepath
  
  SubsampleTable = data.frame(species = speciesSet, fraction=DesignSet[,whichDesign+1])
  subID = subIDSet[whichDesign+1]
  currentrunID =  subID
  which_inFasta16 = 2*whichDesign+1
  which_inFasta23 = which_inFasta16 + 1
  
  # Working with reference database at  "/vast/projects/rrn/microbiome/papercheck/workDB1"
  denoisein_dateString = ifelse(whichDesign ==5,"05022024","06022024")
  #     amplicon_D6322_16S_05022024_D6322_05022024
  
  for (which_subunit in c("16S","23S")){
    inname1 = paste("RADdenoiseBM_",paste("amplicon_D6322",which_subunit,"05022024",currentrunID,denoisein_dateString,sep="_"),".RData",sep="",collapse="")
    load(file=file.path(outRDatapath,inname1))   # Gives D,ab
    whichDN = ifelse(which_subunit=="16S",2*whichDesign+1,2*whichDesign + 2)
    inname2 = paste("strainAnalysis_",which_subunit,"_",whichDN,"_",currentrunID,"_",out_dateString,".RData",sep="")
    load(file=file.path(outRDatapath,inname2))  # Gives gss16,strainsASVset16.df,Allstrains16.df,ASV16u,specCounts16,props16)
    
    if (which_subunit=="16S"){
      cat("\n \n Dataset",currentrunID,"\n 16S Key Outputs. \n Species Cellular Counts \n")
      print(specCounts16)
      print(props16,digits=3)
      print(gss16,digits=3)
      print(Allstrains16.df, digits=3)
      print(ASV16u)
      cat("ASV Raw Counts. Total raw counts ",sum(ab)," \n")
      print(ab)
    } else {
      cat("\n \n Dataset",currentrunID,"\n 23S Key Outputs. \n Species Cellular Counts \n")
      print(specCounts23)
      print(props23,digits=3)
      print(gss23,digits=3)
      print(Allstrains23.df, digits=3)
      print(ASV23u)
      cat("ASV Raw Counts. Total raw counts ",sum(ab)," \n")
      print(ab)
    }
  }     #    end   which_subunit    loop
  
  inname3 = paste("strain_operon_multiplicity_All16S23S_",currentrunID,"_",out_dateString,".RData",sep="")
  load(file=file.path(outRDatapath,inname3))   # Gives D16,abund16,S16,iab16,nlost16,D23,abund23,S23,iab23,nlost23)

  D16 = restoreLong(D16)
  D23 = restoreLong(D23)
  
}      #       end     whichDesign    loop

# Gather data for Table 9 of paper. This data is what is returned from processing ASV alignment data, without 16S,23S merging.
# It is NOT what should be in the paper - see diary item b. of 9-12 March 2024.
D6322full_16 = c(812,1343,658,1103,320,1407,744)
D6322full_23 = c(539,  1053, 633,  713,  265,  1210, 859)
Sub21_16 = c(81,  1343, 651,  1099, 320,   700 ,149)
Sub21_23 = c(76,  1018, 620,   641, 265,   593, 203)
Sub22_16 = c(406,    13,  66,  1103, 320,  1391, 744)
Sub22_23 = c(287,    63, 121,   696, 265,  1175, 838)
Sub23_16 = c(812,    27, 652,  1101, 320,  1380,  15)
Sub23_23 = c(531,    30, 629,   647, 265,  1157,  96)
Sub24_16 = c(810,  1318,  66,  1103, 320,  1398,  15)
Sub24_23 = c(524,  1013, 119,   643, 265,  1166,  98)
Tab9 = rbind(D6322full_16,D6322full_23,Sub21_16,Sub21_23,Sub22_16,Sub22_23,Sub23_16,Sub23_23,Sub24_16,Sub24_23)
colnames(Tab9) = c("Bs","Ef","Ec","Lm","Pa","Se","Sa")

# Gather data for Table 10 of paper. See detailed working and table item c. of 9-12 March 2024.


#  Dataset                     Bs     Ef     Ec     Lm     Pa     Se     Sa        Total
#  D6322Afull_16S_rRNA gene     8     13     21      9     12     36     11         110
#  D6322Afull_23S_rRNA gene    17     13     17     13      2     33     20         115
#  Sub21_16S                    3     10     21     13      5     15      7          74
#  Sub21_23S                    2     10     22      8      5     23      5          75
#  Sub22_16S                    5      1      4     11      3     31      8          63
#  Sub22_23S                   12      2      7      6      5     35     17          84 
#  Sub23_16S                    9      1     23      7      6     35      1          82
#  Sub23_23S                   16      1     18      4      3     38      5          85
#  Sub24_16S                   11     26      4      8      5     32      1          87
#  Sub24_23S                   19      7      5      5      4     36      3          79

# Merging 16S and 23S
basepath = "/vast/projects/rrn/microbiome/papercheck"
outRDatapath = file.path(basepath,"output/RData")
plotpath = file.path(basepath,"output/plots")
textpath = file.path(basepath,"output/text")
intextpath = basepath

for (whichDesign in 5:9){
  cat("\n Design Number ",whichDesign,"\n")
  design = DesignSet[,whichDesign+1]
  design2 = design[c(5,7,3,6,2,4,1)]    
  # Gives the permuted indexing to give alphabetical order of species with "other" being moved to end.)
  
  SubsampleTable = data.frame(species = speciesSet, fraction=DesignSet[,whichDesign+1])
  subID = subIDSet[whichDesign+1]
  currentrunID =  subID
  
  # Get key information characterising the multiplicity of rRNA operons of the genomes of 
  # the bacterial strains in the reference databases. 
  # Need number of operons for any strains of interest. (IU from  unique_operons_count.R )
  inPath = file.path(basepath,"workDB1")
  inname = "strain_operon_multiplicity_All_22012024.RData"
  load(file=file.path(inPath,inname))   # Loads IU, numUniqs
  # IU[[n]] has structure  list(uniques=outline, PIDmat=M,Scoremat=S).
  
  SubsampleTable = data.frame(species = speciesSet, fraction=DesignSet[,whichDesign+1])
  subID = subIDSet[whichDesign+1]
  currentrunID =  subID
  denoisein_dateString = ifelse(whichDesign ==5,"05022024","06022024")
  
  inname = paste("strain_operon_multiplicity_All16S23S_",currentrunID,"_",in_dateString,".RData",sep="")
  load(file=file.path(outRDatapath,inname))
}
  


#############################################################################################################
#############################################################################################################
##            Ab-Initio Coding of B matrix construction and Strain-Set Determination 14 March 2024          ##
#############################################################################################################
#############################################################################################################
possibleGenera = cbind(c("Bacillus","Clostridia","Enterococcus","Escherichia","Lactobacillus",
                         "Listeria", "Pseudomonas","Salmonella","Shigella","Staphylococcus",
                         "Streptococcus","OtherGenus"),
                       c("Bacillus","Clostrid","Enteroco","Escheric","Lactobac",
                         "Listeria", "Pseudomo","Salmonel","Shigella","Staphylo",
                         "Streptoc","OtherGen"))
colnames(possibleGenera) = c("Full","Trunc8")
Ngenera = nrow(possibleGenera)/2

possibleSpecies = cbind(c("subtilis","difficile","faecalis","coli","fermentum",
                          "monocytogenes","aeruginosa","enterica","boydii","aureus",
                          "agalactiae","OtherSpecies"),
                        c("subtilis","difficil","faecalis","coli","fermentu",
                          "monocyto","aerugino","enterica","boydii","aureus",
                          "agalacti","OtherSpe"))
colnames(possibleSpecies) = c("Full","Trunc8")
Nspecies = nrow(possibleSpecies)/2
thisGenera = c(1,3,4,6,7,8,10)
thisSpecies = c(1,3,4,6,7,8,10)
speciesSet = c("subtilis","faecalis","coli","monocytogenes","aeruginosa","enterica","aureus")
specSet = speciesSet
Ns = length(speciesSet)
# Note that the preferred order of species is alphabetical except that "other" is last.
initalSpeciesOrder = possibleSpecies[thisSpecies,1]
ireqOrd = order(initalSpeciesOrder)

in_dateString = "13032024"   # This variable could be removed and replaced by out_dateString everywhere.       
out_dateString = "19032024"     

# Select a species and the fraction of the reads determined to be present in the full
# D63222 dataset that are to be used for this species.
D6322fullDesign = rep(1,7)
sub01Design = c(0.1,1,1,1,1,0.5,0.2)
sub02Design = c(0.5,0.01,0.1,1,1,1,1)   
sub03Design = c(1,0.02,1,1,1,1,0.02)    
sub04Design = c(1,1,0.1,1,1,1,0.02)     

DesignSet = data.frame(D6322fullcheck = D6322fullDesign, sub11=sub01Design,sub12=sub02Design,sub13=sub03Design,sub14=sub04Design,
                       D6322fullcheckA = D6322fullDesign, sub21=sub01Design,sub22=sub02Design,sub23=sub03Design,sub24=sub04Design)
subIDSet = c("D6322","Sub11","Sub12","Sub13","Sub14",
             "D6322A","Sub21","Sub22","Sub23","Sub24")

logging = TRUE
if (logging){
  logPath = "/vast/projects/rrn/microbiome/papercheck/output/text"
  logName = paste("make_B16_B23_Matrices_log_",out_dateString,".txt",sep="")
  sink(file = file.path(logPath,logName),type="output")
}

outRDatapath =  "/vast/projects/rrn/microbiome/papercheck/output/RData"

for (whichDesign in 5:9){
  cat("\n Processing ",subIDSet[whichDesign+1],"dataset. \n")
  inname1 = paste("merged_final_",subIDSet[whichDesign+1],"_results.RData",sep="")
  load(file=file.path(outRDatapath,inname1)) 
  # Above returns  AbundCommon, C16, C23, countSets, DesignSet, dMA16S23S, MD16, MD23, mergeStrain, propSets
  inname2 = paste("strainAnalysis_16S_",2*whichDesign+1,"_",subIDSet[whichDesign+1],"_",in_dateString,".RData",sep="") 
  load(file=file.path(outRDatapath,inname2))
  # Above returns Allstrains16.df, ASV16u, gss16, specCounts16,strainsASVset16.df
  
  inname3 = paste("strainAnalysis_23S_",2*whichDesign+2,"_",subIDSet[whichDesign+1],"_",in_dateString,".RData",sep="")
  load(file=file.path(outRDatapath,inname3))
  # Above returns Allstrains23.df, ASV23u, gss23, specCounts23,strainsASVset23.df
  print(rbind(inname1,inname2,inname3))
  
  for (jspec in 1:Ns){
    cat("\n Processing Species ",possibleGenera[thisSpecies][jspec],"\n")
    t1 = substr(possibleGenera[thisSpecies][jspec],start=1,stop=4)
    for (seqtype in c("16S","23S")){
      if (seqtype == "16S"){
        t2 = sapply(1:length(Allstrains16.df[,1]),function(j){substr(Allstrains16.df[j,1],start=1,stop=4)})
        # Trim all strain names to 11 or fewer characters (to avoid some aberrant strain names!)
        for (jj in 1:length(Allstrains16.df[,1])){
          st1 = Allstrains16.df[jj,"Strain"]
          st2 = substr(st1,start=1,stop=min(nchar(st1),11))
          Allstrains16.df[jj,"Strain"] = st2
        }
        x1 = Allstrains16.df[which(t2==t1),]
        t3 = sapply(1:length(gss16[,1]),function(j){substr(gss16[j,1],start=1,stop=4)})
        x2 = gss16[which(t3==t1),]
      } else {
        t2 = sapply(1:length(Allstrains23.df[,1]),function(j){substr(Allstrains23.df[j,1],start=1,stop=4)})
        # Trim all strain names to 11 or fewer characters (to avoid some aberrant strain names!)
        for (jj in 1:length(Allstrains23.df[,1])){
          st1 = Allstrains23.df[jj,"Strain"]
          st2 = substr(st1,start=1,stop=min(nchar(st1),11))
          Allstrains23.df[jj,"Strain"] = st2
        }
        x1 = Allstrains23.df[which(t2==t1),]
        t3 = sapply(1:length(gss23[,1]),function(j){substr(gss23[j,1],start=1,stop=4)})
        x2 = gss23[which(t3==t1),]
      }
      
      thisASVSets = unique(x1[,6]);  nASVSets = length(thisASVSets)
      thisStrains = unique(x1[,2]);    nstr = length(thisStrains)
      Bn = matrix("",nrow=nstr+1, ncol = nASVSets+2)
      rownames(Bn) = c(thisStrains,"RawCounts")
      colnames(Bn) = c("shared", paste("ASet",thisASVSets,sep="_"), "NumOps")
      # Get the ASV Set raw counts.
      rawCounts = rep(0,nASVSets)
      for (ja in 1:nASVSets){
        rawCounts[ja] = x1[which(x1[,6]==thisASVSets[ja])[1],3]
      }
      Bn[nstr+1,2:(nASVSets+1)] = rawCounts
      numops = rep(0,nstr)
      for (jstr in 1:nstr){
        ix1 = which(x1[,"Strain"]==thisStrains[jstr])
        jAS = sapply(1:length(ix1),function(j){which(thisASVSets==x1[ix1[j],6])})
        Bn[jstr,c(1+jAS)] = rep("1",length(jAS))
        Bn[jstr,"NumOps"] = x2[jstr,"No.Ops of Strain"]
      }
      if (seqtype=="16S") {B16 = Bn} else {B23 = Bn}
    }
    
    # Now identify the "shared" strains and enter information into B16, B23.
    allStrains = setdiff(unique(c(rownames(B16),rownames(B23))),"RawCounts")
    map = matrix(0,nrow=length(allStrains), ncol=2)
    for (jst in 1:length(allStrains)){
      jst1 = which (rownames(B16)==allStrains[jst])
      jst2 = which (rownames(B23)==allStrains[jst])
      if (length(jst1)>0){map[jst,1] = jst1 }
      if (length(jst2)>0){map[jst,2] = jst2 }  
    }
    rownames(map) = allStrains
    bmap = matrix(map>0,nrow=length(allStrains), ncol=2)
    shared = which(sapply(1:length(allStrains),function(j){identical((map[j,]>0),rep(TRUE,2))}))
    B16[map[shared,1],"shared"] = rep("1",length(shared))
    B23[map[shared,2],"shared"] = rep("1",length(shared))
    print(B16);    print(B23);      print(map)
  }     #    end     jspec          loop
}       #    end     whichDesign    loop

if (logging) {sink()}


make_plot_legend_label2 = function(thisSpecies){
  # Generate plot legend label given the set of species - e.g. Given "coli" generate "E.coli"
  # Output is a matrix with species in column 1 and label in column 2.
  # 18 January 2024                                                  [cjw]
  Ns = length(thisSpecies)
  labels = sapply(1:Ns,function(j){
    t1 = substr((possibleGenera[thisSpecies])[j],start=1,stop=1)
    out1 = paste(t1,(possibleSpecies[thisSpecies])[j],sep=".")})
  out = cbind(possibleSpecies[thisSpecies],labels)
}


# Form expressions of the form   sum(as.numeric(c( "165",  "174",   "175",  "267")))/7
# to get cellular abundances from B16 and B23 for a strain with 7 operons.

# Manually calculated merged 16S and 23S cellular abundance values for each dataset.
#    - see diary 16 March 2024.

#                            D6322full        Sub1          Sub2          Sub3          Sub4
#                            16S   23S     16S   23S     16S   23S     16S   23S     16S   23S
#         P.aeruginosa        80    66      80    66      80    66      80    66      80    66
#         S.aureus           124   143      25    34     124   140       3    16       3    16
#         E.coli              94    90      93    89      10    17      93    90       9    17
#         S.enterica         235   112     100    85     199   100     197   165     200   167  
#         E.faecalis         336   264     337   255      10    17       7     8     330   253 
#         L.monocytogenes    184   119     184   107     184   116     184   108     184   107
#         B.subtilis          81    54       8     8      41    29      81    53      81    52

D6322fullCellAbund = cbind(c(80,124,94,235,336,184,81),c(66,143,90,112,264,119,54)) # Ordered alphabetically by species
Sub21 = cbind(c(80,25,93,100,337,184,8),c(66,34,89,85,255,107,8))                   # Ordered alphabetically by species
Sub22 = cbind(c(80,124,9,199,3,184,41),c(66,140,17,168,16,140,29))                  # Ordered alphabetically by species
Sub23 = cbind(c(80,3,93,197,7,184,81),c(66,16,90,165,8,108,53))                     # Ordered alphabetically by species
Sub24 = cbind(c(80,3,9,200,330,184,81),c(66,16,17,167,253,107,52))                  # Ordered alphabetically by species
P16 = cbind(D6322fullCellAbund[,1],Sub21[,1],Sub22[,1],Sub23[,1],Sub24[,1])
P23 = cbind(D6322fullCellAbund[,2],Sub21[,2],Sub22[,2],Sub23[,2],Sub24[,2])
speciesAlphOrd = (possibleSpecies[thisSpecies])[ireqOrd]
Props16 = matrix(0,nrow=Ns,ncol=10)
colnames(Props16) = c(sapply(5:9,function(whichDesign){paste(c("Exp","Obs"),subIDSet[whichDesign+1],sep="")}))
rownames(Props16) = speciesAlphOrd
Props23 = matrix(0,nrow=Ns,ncol=10)
colnames(Props23) = c(sapply(5:9,function(whichDesign){paste(c("Exp","Obs"),subIDSet[whichDesign+1],sep="")}))
rownames(Props23) = speciesAlphOrd

Glengths = c(4214810,2542841,2845422,4639221,2011830,2905187,6792191,4755909,4599354,2874302,1852442,3000000)
dMA16S23S = vector("list",10)

plotname = paste("mergedB_16S23S_proportions_barplots_",out_dateString, ".pdf",sep="")
pdf(file=file.path(plotpath,plotname),paper="a4")
for (whichDesign in 5:9){
  jds = whichDesign-4
  # Now construct proportions data for the species cellular abundances as stored in AbundCommon
  # which has used the minimum for each species (11 Oct 2023).  
#  SubsampleTable = data.frame(species = speciesSet, fraction=DesignSet[,whichDesign+1])
  subID = subIDSet[whichDesign+1]
  currentrunID =  subID
  design =  DesignSet[,whichDesign+1]/sum(DesignSet[,whichDesign+1])
  design2 = design[ireqOrd] 
  
  if (whichDesign %in% c(0,5,10)){
    cellAbund = 1/Glengths[(thisSpecies)[ireqOrd]]
    theoretical16 = cellAbund/sum(cellAbund)
    theoretical23 = theoretical16
  } else {
    theoretical16 = design2*D6322fullCellAbund[,1]/sum(design2*D6322fullCellAbund[,1])
    theoretical23 = design2*D6322fullCellAbund[,2]/sum(design2*D6322fullCellAbund[,2])
  } 
  
  proportions16 = cbind(theoretical16,P16[,jds]/sum(P16[,jds]))
  Props16[,(2*jds-1):(2*jds)] = proportions16
  colnames(proportions16) = c(ifelse(whichDesign==5,"Theory","Design16"),paste(subIDSet[whichDesign+1],"16S",sep="_"))
  rownames(proportions16) = speciesAlphOrd
  proportions23 = cbind(theoretical23,P23[,jds]/sum(P23[,jds]))
  colnames(proportions23) = c(ifelse(whichDesign==5,"Theory","Design23"),paste(subIDSet[whichDesign+1],"23S",sep="_"))
  rownames(proportions23) = speciesAlphOrd
  Props23[,(2*jds-1):(2*jds)] = proportions23
  proportions = cbind(proportions16,proportions23)
  MA16 = acomp(proportions16)
  MA23 = acomp(proportions23)
  MA1623 = acomp(proportions)
  dMA16S23S[[whichDesign+1]]  = list(MA16=as.matrix(dist(t(MA16))),MA23=as.matrix(dist(t(MA23))), MA1623=as.matrix(dist(t(MA1623))))
  t1 = c(round(100*dMA16S23S[[whichDesign+1]]$MA1623[1,2])/100,round(100*dMA16S23S[[whichDesign+1]]$MA1623[3,4])/100,
         round(100*dMA16S23S[[whichDesign+1]]$MA1623[2,4])/100)
  t1s = paste(t1,sep=",  ",collapse=",  ")
  bplotRef16 = ifelse(whichDesign==0,"Theory","Des16")
  bplotRef23 = ifelse(whichDesign==0,"Theory","Des23")
  xuplim = 8               #  ifelse(nci==0,5,ifelse(nci==1,8,ifelse(nci==2,11,ifelse(nci==3,15,5+3*nci))))
  names.str = rep("",4)    #  rep("",4+2*(nci-1))
  names.str[1:2] = c(bplotRef16,"16SD")
  names.str[3:4] = c(bplotRef23,"23SD")
  colnames(proportions) = names.str
  Labels23 = make_plot_legend_label2(thisSpecies)
  barplot(as.matrix(proportions), col=c(1:7),cex.names=0.8, xlim=c(0,xuplim),
          names.arg= names.str,
          legend.text=c(Labels23[sapply(1:Ns,function(j){which(Labels23==(speciesSet[ireqOrd])[j])}),2]),
          sub=paste("(Aitchison distances of other columns to Theoretical: ",t1s,")",sep=""),
          main=paste("Merged 16S and 23S for Dataset ",currentrunID,sep="  "))
  
  # The following is a trimmed down version intended for use in a paper to be submitted (9 Dec 2023).
  
  barplot(as.matrix(proportions), col=c(1:7),cex.names=1.2, xlim=c(0,xuplim), main="",
          names.arg= names.str)
  fullNames = c("Bacillus subtilis","Enterococcus faecalis","Escherichia coli",
                "Listeria monocytogenes","Pseudomonas aeruginosa", 
                "Salmonella enterica","Stapholycoccus aureus")
  fullNames.species = c("subtilis","faecalis","coli", "monocytogenes","aeruginosa", "enterica","aureus")
  fNOrder = sapply(1:Ns,function(j){which(fullNames.species==(speciesSet[ireqOrd])[j])})
  legend(x="right",legend=fullNames[fNOrder[Ns:1]], text.font= rep(3,7),col=c(Ns:1),pch=rep(15,7))
  cat("\n Aitchison distances \n");  print(dMA16S23S[[whichDesign+1]]$MA1623,digits=2)
  
  outname = paste("mergedB_final_",currentrunID,"_results.RData",sep="")
  save(file=file.path(outRDatapath,outname),proportions16,proportions23,P16,P23,dMA16S23S,DesignSet)
}   #       end     whichDesign     loop

print(Props16,digits=2);  print(Props23,digits=2)
dev.off()

# Un-merged proportions
UProps16 = matrix(0,nrow=Ns,ncol=10)
colnames(UProps16) = c(sapply(5:9,function(whichDesign){paste(c("Exp","Obs"),subIDSet[whichDesign+1],sep="")}))
rownames(UProps16) = speciesAlphOrd
UProps23 = matrix(0,nrow=Ns,ncol=10)
colnames(UProps23) = c(sapply(5:9,function(whichDesign){paste(c("Exp","Obs"),subIDSet[whichDesign+1],sep="")}))
rownames(UProps23) = speciesAlphOrd
UProps16[,c(1,3,5,7,9)] = Props16[,c(1,3,5,7,9)]
UProps23[,c(1,3,5,7,9)] = Props23[,c(1,3,5,7,9)]
UProps16[,2] = c(0.074, 0.305,  0.085, 0.167,  0.073, 0.183, 0.113)[ireqOrd]
UProps16[,4] = c(0.010, 0.407,  0.113, 0.222,  0.097, 0.121, 0.030)[ireqOrd]
UProps16[,6] = c(0.063,   0.005, 0.015, 0.287, 0.125, 0.311, 0.194)[ireqOrd]
UProps16[,8] = c(0.126,  0.011, 0.145, 0.285, 0.124, 0.306, 0.0039)[ireqOrd]
UProps16[,10] = c(0.091, 0.372,  0.011, 0.208,  0.090, 0.225,  0.0028)[ireqOrd]
UProps23[,2] = c(0.060, 0.289,  0.099, 0.131,  0.073, 0.190, 0.157)[ireqOrd]
UProps23[,4] = c(0.012,  0.396, 0.138, 0.166,  0.103, 0.132, 0.053)[ireqOrd]
UProps23[,6] = c(0.053,   0.029, 0.029, 0.211,   0.120, 0.305, 0.254)[ireqOrd]
UProps23[,8] = c(0.105,  0.015, 0.178, 0.213,  0.131, 0.327,  0.032)[ireqOrd]
UProps23[,10] = c(0.077, 0.373,  0.025, 0.158,  0.097, 0.246,  0.024)[ireqOrd]

plotname = paste("notmergedB_16S23S_proportions_barplots_",out_dateString, ".pdf",sep="")
pdf(file=file.path(plotpath,plotname),paper="a4")
for (whichDesign in 5:9){
  jds = whichDesign-4
  # Now construct proportions data for the species cellular abundances as stored in AbundCommon
  # which has used the minimum for each species (11 Oct 2023).  
  #  SubsampleTable = data.frame(species = speciesSet, fraction=DesignSet[,whichDesign+1])
  subID = subIDSet[whichDesign+1]
  currentrunID =  subID
  design =  DesignSet[,whichDesign+1]/sum(DesignSet[,whichDesign+1])
  design2 = design[ireqOrd] 
  
  if (whichDesign %in% c(0,5,10)){
    cellAbund = 1/Glengths[(thisSpecies)[ireqOrd]]
    theoretical16 = cellAbund/sum(cellAbund)
    theoretical23 = theoretical16
  } else {
    theoretical16 = design2*D6322fullCellAbund[,1]/sum(design2*D6322fullCellAbund[,1])
    theoretical23 = design2*D6322fullCellAbund[,2]/sum(design2*D6322fullCellAbund[,2])
  } 
  MA16 = acomp(UProps16[,(2*jds-1):(2*jds)])
  MA23 = acomp(UProps23[,(2*jds-1):(2*jds)])
  MA1623 = acomp(cbind(MA16,MA23))
  dMA16S23S[[whichDesign+1]]  = list(MA16=as.matrix(dist(t(MA16))),MA23=as.matrix(dist(t(MA23))), MA1623=as.matrix(dist(t(MA1623))))
  t1 = c(round(100*dMA16S23S[[whichDesign+1]]$MA1623[1,2])/100,round(100*dMA16S23S[[whichDesign+1]]$MA1623[3,4])/100,
         round(100*dMA16S23S[[whichDesign+1]]$MA1623[2,4])/100)
  t1s = paste(t1,sep=",  ",collapse=",  ")
  bplotRef16 = ifelse(whichDesign==0,"Theory","Des16")
  bplotRef23 = ifelse(whichDesign==0,"Theory","Des23")
  xuplim = 8               #  ifelse(nci==0,5,ifelse(nci==1,8,ifelse(nci==2,11,ifelse(nci==3,15,5+3*nci))))
  names.str = rep("",4)    #  rep("",4+2*(nci-1))
  names.str[1:2] = c(bplotRef16,"16SD")
  names.str[3:4] = c(bplotRef23,"23SD")
  proportions = cbind(UProps16[,(2*jds-1):(2*jds)],UProps23[,(2*jds-1):(2*jds)])
  colnames(proportions) = names.str
  Labels23 = make_plot_legend_label2(thisSpecies)
  barplot(as.matrix(proportions), col=c(1:7),cex.names=0.8, xlim=c(0,xuplim),
          names.arg= names.str,
          legend.text=c(Labels23[sapply(1:Ns,function(j){which(Labels23==(speciesSet[ireqOrd])[j])}),2]),
          sub=paste("(Aitchison distances of other columns to Theoretical: ",t1s,")",sep=""),
          main=paste("Not-Merged 16S and 23S for Dataset ",currentrunID,sep="  "))
  
  # The following is a trimmed down version intended for use in a paper to be submitted (9 Dec 2023).
  
  barplot(as.matrix(proportions), col=c(1:7),cex.names=1.2, xlim=c(0,xuplim), main="",
          names.arg= names.str)
  fullNames = c("Bacillus subtilis","Enterococcus faecalis","Escherichia coli",
                "Listeria monocytogenes","Pseudomonas aeruginosa", 
                "Salmonella enterica","Stapholycoccus aureus")
  fullNames.species = c("subtilis","faecalis","coli", "monocytogenes","aeruginosa", "enterica","aureus")
  fNOrder = sapply(1:Ns,function(j){which(fullNames.species==(speciesSet[ireqOrd])[j])})
  legend(x="right",legend=fullNames[fNOrder[Ns:1]], text.font= rep(3,7),col=c(Ns:1),pch=rep(15,7))
  cat("\n Aitchison distances \n");  print(dMA16S23S[[whichDesign+1]]$MA1623,digits=2)
  
  outname = paste("notmergedB_final_",currentrunID,"_results.RData",sep="")
  save(file=file.path(outRDatapath,outname),UProps16,UProps23,dMA16S23S,DesignSet)
}   #       end     whichDesign     loop

print(UProps16,digits=2);  print(UProps23,digits=2)
dev.off()


ix16 = order(C16$StCounts[,"Species"])
ix23 = order(C23$StCounts[,"Species"])
for (specid in possibleSpecies[thisSpecies][ireqOrd]){
  cat("\nSpecies ",specid,"\n 16S")
  ix16x = ix16[which(C16$StCounts[ix16,"Species"]==specid)]
  miniTab16 = cbind(C16$StCounts[ix16x,],MD16[ix16x,4:ncol(MD16)])
  print(miniTab16);  cat("\n 23S")
  ix23x = ix23[which(C23$StCounts[ix23,"Species"]==specid)]
  miniTab23 = cbind(C23$StCounts[ix23x,],MD16[ix23x,4:ncol(MD23)])
  print(miniTab23)
}





