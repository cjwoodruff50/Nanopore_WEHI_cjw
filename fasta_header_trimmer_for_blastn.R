# Rscript to truncate the header of each fasta sequence in a multi-sequence fasta file to no 
# more than 50 characters in order to use that fast file in makeblastdb.  It should be possible using 
# awk but I haven't found a suitable example on the internet.  Hence I am writing R code to do it!
#
# 2 February 2024                                                          [cjw]
basepath = "/vast/projects/rrn/microbiome/blastdb/workDB1output"
seqtype = "16S"         # Or "23S"   or   "rrn"   etc.
fastaPath = file.path(basepath,seqtype)
fastaName = "workDB1_23S_19012024"      # Name of the input fasta file - note no extension, but form is .fasta .
newfastaName = paste(fastaName,"_nameAdj",sep="",collapse="")
#newfastaName = paste(unlist(strsplit(fastaName,split="[.]"))[1],"Adj.fa",sep="",collapse="")
F = readLines(con=file.path(fastaPath, fastaName))
ih = which(sapply(1:length(F),function(j){grep("^>",F[j])})>0)
for (j in 1:length(ih)){
  t1 = F[ih[j]]
  if (nchar(t1)>50){
    t2 = unlist(strsplit(t1,split="_"))
    # Shorten a reconstructed t1 by one or more of the following:-
    #    a. Remove any of the words "genome", "chromosome", "complete";
    #    b. Truncate any t2[3:length(t2)] elements to no greater than 8 characters
    #    c. Ensure that the operon number is attached at the end of the reconstructed header.
    nt2 = length(t2)
    # May have t2 elements of the form "genome16S", "genome23S". Check for these and, if present,
    # remove the "genome" part.
    igen = which(sapply(3:nt2,function(j){out=(t2[j]=="genome16S") || (t2[j]=="genome23S")}))
    if (length(igen)>0){
      t2[igen+2] = substr(t2[igen+2],start=7,stop=9)
    }
    iout = which(t2 %in% c("chromosome","complete","genome","gene","operon","strain","subspp","serotype","sequence"))
    ikeep = setdiff(1:nt2,c(iout,nt2))
    for (k in 2:nt2){if (nchar(t2[k]>8)){t2[k]=substr(t2[k],start=1,stop=8)}}
    t21 = paste("Op",t2[nt2],sep="",collapse="")
    t1n = paste(paste(t2[ikeep],sep="_",collapse="_"),t21,sep="_",collapse="_")
    if (nchar(t1n)>50){ # Implement the former brute force method
      t2n = unlist(strsplit(t1n,split="_"))
      t2n1 = t2n[length(t2n)]
      t2n2 = substr(t1n,start=1,stop=44)
      t1n = paste(t2n2,t2n1,sep="_")
    }
    F[ih[j]] = t1n
    print(c(ih[j],t1n,F[ih[j]]))
  } else {print(c("         ",j, ih[j],F[ih[j]]))}
}
writeLines(F,con=file.path(fastaPath,newfastaName))
