biomart_full<-read.csv(file="Data/Annotation biomart UniProt/mart_export_fullENS.txt",header=TRUE, stringsAsFactors = FALSE, sep='\t')
biomart_uni<-read.csv(file="Data/Annotation biomart UniProt/mart_export_withUniProt.txt",header=TRUE, stringsAsFactors = FALSE, sep='\t')
uni_annot<-read.csv(file="Data/Annotation biomart UniProt/uniprot-filtered-organism__Homo+sapiens+(Human)+[9606]_+AND+review-- (1).tab",header=TRUE, stringsAsFactors = FALSE, sep='\t')

# 1 ########################################################################   annotate genes
A1<-read.csv(file="Data/FullMergedResults.txt",header=TRUE, stringsAsFactors = FALSE, sep='\t')
A1$source=sapply(A1$TIScategory, function(x) sapply(x, function(x) as.factor(as.matrix(strsplit(x, "_")[[1]][1]))))
A1$ID=1:nrow(A1)
A1$TIS_cat_order=ifelse(grepl('UniProt_aTIS',A1$TIScategory),1,ifelse(grepl('ENST_aTIS',A1$TIScategory),2,ifelse(grepl('UniProtIso_aTIS',A1$TIScategory),3,4)))

A1$UniProt.stable.ID=ifelse(A1$source!="ENST",sapply(A1$Accession , function(x) sapply(x, function(x) as.matrix(strsplit(x, "-")[[1]][1]))),NA)
A1$ENST=ifelse(A1$source=="ENST",sapply(A1$Accession , function(x) sapply(x, function(x) as.matrix(strsplit(x, "_")[[1]][1]))),NA)

A1$ENSG=ifelse(A1$source!="ENST",biomart_uni$Gene.stable.ID[match(A1$UniProt.stable.ID , biomart_uni$UniProtKB.Swiss.Prot.ID)],biomart_full$Gene.stable.ID[match(A1$ENST , biomart_full$Transcript.stable.ID)])
A1$Gene.name=ifelse(A1$source!="ENST",biomart_uni$Gene.name[match(A1$UniProt.stable.ID , biomart_uni$UniProtKB.Swiss.Prot.ID)],biomart_full$Gene.name[match(A1$ENST , biomart_full$Transcript.stable.ID)])
A1$ENSG_biotype=biomart_full$Gene.type[match(A1$ENSG , biomart_full$Gene.stable.ID)]

A1$ENST_biotype=ifelse(A1$source=="ENST",biomart_full$Transcript.type[match(A1$ENST , biomart_full$Transcript.stable.ID)],NA)
A1$ENST_TSL=ifelse(A1$source=="ENST",biomart_full$Transcript.support.level..TSL.[match(A1$ENST , biomart_full$Transcript.stable.ID)],NA)


#some proteins without genes fix manually; add gene names from UniProt annotation

A2=A1[is.na(A1$Gene.name),]
A3=A2[,colnames(A2) %in% c( "Accession", "Start", "Description", "Isoforms")]
A3$Gene.name=uni_annot$Gene.names[match(A3$Accession,uni_annot$Entry)]
A3$SecProt=sapply(A3$Isoforms, function(x) sapply(x, function(x) strsplit(x, ",")[[1]][1]))
A3$SecStart=sapply(A3$SecProt, function(x) strsplit(x, " \\(")[[1]][2])
A3$SecStart=sapply(A3$SecStart, function(x) strsplit(x, "-")[[1]][1])
A3$SecProt=sapply(A3$SecProt, function(x) strsplit(x, " \\(")[[1]][1])
A3$SecProt=ifelse(grepl("-",A3$SecProt),sapply(A3$SecProt, function(x) strsplit(x,"-")[[1]][1]),A3$SecProt)
A3$SecProt=ifelse(grepl("ENST",A3$SecProt),sapply(A3$SecProt, function(x) strsplit(x,"_")[[1]][1]),A3$SecProt)
A3$SecGene=ifelse(substr(A3$SecProt, 1,4)=="ENST", 
                  biomart_full$Gene.name[match(A3$SecProt , biomart_full$Transcript.stable.ID)],
                  biomart_uni$Gene.name[match(A3$SecProt , biomart_uni$UniProtKB.Swiss.Prot.ID)])
A3$SecENSG=ifelse(substr(A3$SecProt, 1,4)=="ENST", 
                  biomart_full$Gene.stable.ID[match(A3$SecProt , biomart_full$Transcript.stable.ID)],
                  biomart_uni$Gene.stable.ID[match(A3$SecProt , biomart_uni$UniProtKB.Swiss.Prot.ID)])

#if 2nd protein has the same start position, different gene name, keep the 2nd protein ENSG 
#if 2nd protein has the same gene name as the 1st protein's uniprot gene name, but different start position, use 2nd prot as source of ENSG
#in other cases use the 
A3$Final.Gene.name=ifelse(is.na(A3$SecGene),sapply(A3$Gene.name, function(x) strsplit(x, " ")[[1]][1]), 
                          ifelse(A3$Start==as.numeric(A3$SecStart) | 
                                   sapply(A3$Gene.name, function(x) strsplit(x, " ")[[1]][1])==A3$SecGene ,
                                 A3$SecGene,sapply(A3$Gene.name, function(x) strsplit(x, " ")[[1]][1]))) 
A3$Final.ENSG=ifelse(is.na(A3$SecENSG ),NA, 
                     ifelse(A3$Start==as.numeric(A3$SecStart) | 
                              sapply(A3$Gene.name, function(x) strsplit(x, " ")[[1]][1])==A3$SecGene ,
                            A3$SecENSG,NA)) 

#integrate the ENSG and Gene.name to A1 table
A1$Gene.name=ifelse(is.na(A1$Gene.name),A3$Final.Gene.name[match(A1$Accession, A3$Accession)],A1$Gene.name)
A1$ENSG=ifelse(is.na(A1$ENSG),A3$Final.ENSG[match(A1$Accession, A3$Accession)],A1$ENSG)
A1$ENSG_biotype=biomart_full$Gene.type[match(A1$ENSG , biomart_full$Gene.stable.ID)]

# 1B ################### map peptide to genomic positions, make bed, find out which exon they belog to
# THIS WILL CHANGE THE biomart-selected ENSG, so update it at the end

#for ensembl gene/transcript/protein coordinate change, nice to use, but doesn't work for 98 ens
#https://bioconductor.org/packages/release/bioc/manuals/ensembldb/man/ensembldb.pdf

### therefore, my own alignment method:
library(GenomicRanges)

#for UniProt map peptide if entire UniProt protein is identical to ENST protein; 
#use the custom database for this comparison; identical use script like:
library("Biostrings")
custom.fasta="Data/proteoformer_LeeGao_uniprot_plusisoform_2019_04.fasta"
c.fa=(readAAStringSet(custom.fasta))
c.na=names(c.fa)
c.na=as.data.frame(c.na[grepl("sp\\|",c.na)])
colnames(c.na)<-"header"
c.na$header=as.character(c.na$header)
c.na$uni=sapply(c.na$header, function(x) strsplit(x, "\\|")[[1]][2])
c.na$ens=sapply(c.na$header, function(x) strsplit(x, "\\[ENS")[[1]][2])
c.na$ens[!is.na(c.na$ens)]=paste0("ENS", c.na$ens[!is.na(c.na$ens)])
c.na$sel=sapply(c.na$ens, function(x) grepl("_aTIS", strsplit(x, "#")[[1]]))
c.na$ens=sapply(c.na$ens, function(x) strsplit(x, "#")[[1]])
c.na$ens2=sapply(c.na$ens, function(x) strsplit(x, "#")[[1]][1])
for (i in 1:nrow(c.na)){c.na$ens[[i]]<-c.na$ens[[i]][c.na$sel[[i]]][1]}
c.na$ens=ifelse(is.na(c.na$ens),c.na$ens2,c.na$ens)
c.na$ens=gsub("\\]","", c.na$ens)

#use results.db export for dist.to.transcript.start
codon.1<-read.csv(file="Data/export_info_Riboseq_results/TIS_info_db1_Lee.txt",header=TRUE, stringsAsFactors = FALSE, sep='\t')
codon.2<-read.csv(file="Data/export_info_Riboseq_results/TIS_info_db2_Gaoctrl.txt",header=TRUE, stringsAsFactors = FALSE, sep='\t')
codon.3<-read.csv(file="Data/export_info_Riboseq_results/TIS_info_db3_GaoAA.txt",header=TRUE, stringsAsFactors = FALSE, sep='\t')
codon.db=rbind(codon.1,codon.2,codon.3)
rm(codon.1,codon.2,codon.3)
codon.db$name=paste(codon.db$stable_id, codon.db$chr, codon.db$start, codon.db$annotation, sep="_")
codon.db=codon.db[!duplicated(codon.db$name),]
# BiocManager::install("biomaRt")
library(biomaRt)
mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl",  host = "http://sep2019.archive.ensembl.org")
listAttributes(mart, what = c("name","description","page"))
strand_additional=getBM(attributes = c("ensembl_transcript_id","strand"),
                        filters="ensembl_transcript_id",
                        values=codon.db$stable_id,
                        mart=mart)
codon.db$strand=strand_additional$strand[match(codon.db$stable_id,strand_additional$ensembl_transcript_id)]
sum(is.na(codon.db$strand))#0
codon.db$dist_to_transcript_start_corr=ifelse(codon.db$strand==1,codon.db$dist_to_transcript_start-1, codon.db$dist_to_transcript_start)

E0U=A1$Accession[A1$source!="ENST"]
E1U=c.na$ens[match(E0U, c.na$uni)]
E1U[is.na(E1U)]="none_none_1_aTIS_111db1"
E2U=sapply(E1U, function(x) strsplit(x,"_")[[1]][1])
E3U=sapply(E1U, function(x) strsplit(x,"_")[[1]][2])
E4U=sapply(E1U, function(x) strsplit(x,"_")[[1]][3])
E4U.anot=sapply(E1U, function(x) strsplit(x,"_")[[1]][4])
E234U=paste(E2U, E3U, E4U, E4U.anot, sep="_")
E5U=A1$Sequence[A1$source!="ENST"]
E6U=A1$Start[A1$source!="ENST"]
E7=codon.db[!duplicated(codon.db$name),c("name","dist_to_transcript_start_corr")]

E8U<- GRanges(E3U, IRanges(start = as.numeric(E4U), width=rep(1,length(E3U))))

E11U=list()
for (i in 1:length(E8U)) {
  E11U[i]<-ifelse(E1U[i]=="none_none_1_aTIS_111db1",NA ,E7$dist_to_transcript_start_corr[E7$name %in% E234U[i]])
}

sum(E1U=="none_none_1_aTIS_111db1")
#109 UniProt didn't have ENST in fasta headers
su=unlist(E11U)+1+ unlist(lapply(E6U, function(x) (x-1)*3))
Wu=unlist(lapply(E5U, function(x) nchar(x)*3))
Wu[is.na(su)]<-0
su[is.na(su)]<--1
nu=E2U
E12U=IRanges(start=su, width=Wu, names=nu)

#for ENST use transcript coordinates
#convert accession start from genomic to transcript coordinate
#make a Grange from the tr_start to tr_end like tr_start+ (aa peptide * 3)-1
#change this tr_range to genome range

E1=A1$Accession[A1$source=="ENST"]
E2=sapply(E1, function(x) strsplit(x,"_")[[1]][1])
E3=sapply(E1, function(x) strsplit(x,"_")[[1]][2])
E4=sapply(E1, function(x) strsplit(x,"_")[[1]][3])
E4.anot=sapply(E1, function(x) strsplit(x,"_")[[1]][4])
E234=paste(E2, E3, E4, E4.anot, sep="_")
E5=A1$Sequence[A1$source=="ENST"]
E6=A1$Start[A1$source=="ENST"]
E8<- GRanges(E3, IRanges(start = as.numeric(E4),width=rep(1,length(E3))))

E11=list()
for (i in 1:length(E8)) {
  E11[i]<-E7$dist_to_transcript_start_corr[E7$name %in% E234[i]]
}
s=unlist(E11)+1+ unlist(lapply(E6, function(x) (x-1)*3))
W=unlist(lapply(E5, function(x) nchar(x)*3))
n=E2
E12E=IRanges(start=s, width=W, names=n)
E12=c(E12U,E12E)

#http://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesHOWTOs.pdf
# 
library(GenomicFeatures)
#BiocManager::install("RMariaDB")
library(RMariaDB)
library(rtracklayer)

txdb <- makeTxDbFromEnsembl(organism="Homo sapiens",  release=98)
ex=exonsBy(txdb, by=c("tx"), use.names=TRUE)
#tr=transcriptsBy(txdb, by=c("cds"), use.names=TRUE)
ex1=ex[names(ex) %in% names(E12)]

ex.tr=IRangesList()

for (i in 1:length(ex1)){
  w.ex=list()
  s.ex=list()
  
  for (n in 1:length(ex1[[i]])){
    w.ex[n]=width(ex1[[i]])[n]
    s.ex[n]=ifelse(n==1,1,s.ex[[n-1]]+w.ex[[n-1]])
  }
  
  ex.tr[[i]]=IRanges(start=unlist(s.ex), width=unlist(w.ex), names=rep(names(ex1[i]), length(ex1[[i]])),
                     exon_name=ex1[[i]]@elementMetadata$exon_name, exon_rank=ex1[[i]]@elementMetadata$exon_rank)
  
  
}
names(ex.tr)<-names(ex1)

ov=IRangesList()
for (i in 1:length(E12)){
  if (start(E12[i])==-1){ov[[i]]=IRanges()}
  else {
    ov[[i]]=overlapsRanges(ex.tr[[names(E12)[i]]], E12[i])
    ov[[i]]=IRanges(start=start(ov[[i]]), width=width(ov[[i]]), names=rep(names(E12)[i], length(ov[[i]])),
                    exon_rank=ex.tr[[names(E12)[i]]]@elementMetadata$exon_rank[poverlaps(ex.tr[[names(E12)[i]]], E12[i])],
                    exon_name=ex.tr[[names(E12)[i]]]@elementMetadata$exon_name[poverlaps(ex.tr[[names(E12)[i]]], E12[i])])
  }
}
#table(unlist(lapply(ov, function(x) length(x)==0))) #108 didn't have overlap

gr=GRangesList()
for (i in 1:length(ov)){
  if (length(ov[[i]])!=0){
    ex2=ex1[names(ov[[i]])][[1]][elementMetadata(ov[[i]])$exon_rank]
    ex3=ex2
    tr2=ex.tr[names(ov[[i]])][[1]][elementMetadata(ov[[i]])$exon_rank]
    if(decode(strand(ex2)=="+")){
      for (n in 1:length(ex2)){
        sel=start(tr2)[n]:end(tr2)[n] %in% start(ov[[i]])[n]:end(ov[[i]])[n]
        start(ex3)[n]<-(start(ex2)[n]:end(ex2)[n])[sel][1]
        end(ex3)[n]<-tail((start(ex2)[n]:end(ex2)[n])[sel], 1)
        width(ex3)[n]<-length(start(ex3)[n]:end(ex3)[n])
      }
    }    
    
    if(decode(strand(ex2)=="-")){
      for (n in 1:length(ex2)){
        sel=end(tr2)[n]:start(tr2)[n] %in% start(ov[[i]])[n]:end(ov[[i]])[n]
        start(ex3)[n]<-(start(ex2)[n]:end(ex2)[n])[sel][1]
        end(ex3)[n]<-tail((start(ex2)[n]:end(ex2)[n])[sel], 1)
        width(ex3)[n]<-length(start(ex3)[n]:end(ex3)[n])
      }    
    }
    
    gr[[i]]=ex3
  }
  else {gr[[i]]=GRanges()}
}



#export to bam file
#names(gr)<-c(E5U,E5)
E0=c(E0U, E1)
E13U=A1$Sequence.information..longest.[A1$source!="ENST"]
E13E=A1$Sequence.information..longest.[A1$source=="ENST"]
E13=c(E13U, E13E)
E14=paste0(E13, ";Protein_accession=",E0)
E14=sapply(E14, function(x) gsub(" ","",unlist(x)))
names(gr)<-E14
bed.pep=gr[unlist(lapply(gr, length))!=0]
#problems with 0-based coordinates + strand (should be first start-1, last end -1)
#to correct it 100%, I would need to get strand to codon.db e.g. from biomaRt
#for every + strand gene, the distance to transcript start should be -1
#the rest of the code the same
export(asBED(bed.pep),"peptide.genomicPos.bed",format="bed")


#on which exon the proteoform starts?

exon_rank_TIS=lapply(ov, function(x) (elementMetadata(x)$exon_rank))
E0.db=data.frame(E0)
E0.db$ID=c(A1$ID[A1$source!="ENST"],A1$ID[A1$source=="ENST"])

E0.db$exon_rank_TIS=exon_rank_TIS


A1$exon_rank_TIS=E0.db$exon_rank_TIS[match(A1$ID,E0.db$ID)]

#on which exon the aTIS starts?
#for UniProt peptides: are all acorresponding ENST_x_x_aTIS or some are CDS UTR?
table(grepl("_aTIS", E1U)) #13 not aTIS
E1U[!grepl("_aTIS", E1U)]
E0U[!grepl("_aTIS", E1U)]

#for ENST peptides: get the aTIS position in that transcript (if there is aTIS e.g. ntr don't have one)
E9=codon.db[codon.db$annotation=="aTIS" & (codon.db$stable_id %in% E2), c("stable_id", "chr", "start", "annotation")]
E9=unique(E9)
E9=merge(as.data.frame(E2),E9, by.x="E2", by.y="stable_id", all.x=TRUE)
E9=E9[sapply(E2,function(x) grep(x,E9$E2)[1]),]
E9$chr[is.na(E9$chr)]<-"none"
E9$start[is.na(E9$start)]<-1
E10=GRanges(E9$chr, IRanges(start = E9$start,width=rep(1,nrow(E9))))


E8.aTIS=c(E8U, E10)

ov.aTIS=GRangesList()
for (i in 1:length(E8.aTIS)){
  if (as.character(seqnames(E8.aTIS[i]))=="none"){ov.aTIS[[i]]=GRanges()}
  else {
    ov.aTIS[[i]]=intersect(ex1[[names(E12)[i]]], E8.aTIS[i], ignore.strand=TRUE)
    ov.aTIS[[i]]=GRanges(seqnames=rep(names(E12)[i], length(ov.aTIS[[i]])),IRanges(start=start(ov.aTIS[[i]]), width=width(ov.aTIS[[i]]),
                                                                                   exon_rank=ex1[[names(E12)[i]]]@elementMetadata$exon_rank[pintersect(ex1[[names(E12)[i]]], E8.aTIS[i])@elementMetadata$hit],
                                                                                   exon_name=ex1[[names(E12)[i]]]@elementMetadata$exon_name[pintersect(ex1[[names(E12)[i]]], E8.aTIS[i])@elementMetadata$hit]))
  }
}

exon_rank_of_corresponding_aTIS=lapply(ov.aTIS, function(x) (elementMetadata(x)$exon_rank))
E0.db$exon_rank_of_corresponding_aTIS=exon_rank_of_corresponding_aTIS

A1$exon_rank_of_corresponding_aTIS=E0.db$exon_rank_of_corresponding_aTIS[match(A1$ID,E0.db$ID)]

#add ENST to UniProt IDs?
E0.db$ENST=c(E2U,E2)
E0.db$ENSG=biomart_full$Gene.stable.ID[match(E0.db$ENST , biomart_full$Transcript.stable.ID)]
E0.db$ENST_biotype=biomart_full$Transcript.type[match(E0.db$ENST , biomart_full$Transcript.stable.ID)]
E0.db$ENST_TSL=biomart_full$Transcript.support.level..TSL.[match(E0.db$ENST , biomart_full$Transcript.stable.ID)]
E0.db$A1ENSG=A1$ENSG[match(E0.db$ID, A1$ID)]
E0.db$test=(E0.db$ENSG==E0.db$A1ENSG)
E0.db$Gene.name=biomart_full$Gene.name[match(E0.db$ENST , biomart_full$Transcript.stable.ID)]
E0.db$A1name=A1$Gene.name[match(E0.db$ID, A1$ID)]
E0.db$test2=(E0.db$Gene.name==E0.db$A1name)
#actually 62 times a different gene comes out from Proteoofrmer header - based analysis (for peptide-genome map)
# and different from  biomart; often biomart gives gene (copies) present on the wierd chromosomes, like   Chromosome CHR_HSCHR5_2_CTG1_1: 69,673,667-69,702,549 reverse strand.                                        
#sometimes proteoformer gives strange things, like pseudogenes etc.
#12 times gene names are not the same
#

#are we having empty ENSG for 109 non-mapped? for sure! these have "none" as E0.db$ENST
E.db=E0.db[E0.db$ENST!="none",]
A1$Gene.name=ifelse(A1$ID %in% E.db$ID, E.db$Gene.name[match(A1$ID, E.db$ID)],A1$Gene.name)
A1$ENSG=ifelse(A1$ID %in% E.db$ID, E.db$ENSG[match(A1$ID, E.db$ID)],A1$ENSG)
A1$ENST=ifelse(A1$ID %in% E.db$ID ,E.db$ENST[match(A1$ID, E.db$ID)],A1$ENST)

A1$ENSG_biotype=biomart_full$Gene.type[match(A1$ENSG , biomart_full$Gene.stable.ID)]
A1$ENST_biotype=biomart_full$Transcript.type[match(A1$ENST , biomart_full$Transcript.stable.ID)]
A1$ENST_TSL=biomart_full$Transcript.support.level..TSL.[match(A1$ENST , biomart_full$Transcript.stable.ID)]

############################################################################### END peptide-genome-map

A4=aggregate(TIS_cat_order ~ Gene.name , data = A1, FUN = function(x) unique(sort(x)))

require(dplyr)
A1 %>% count(Gene.name) ->A5

A1$Multiple_TIS=paste0(A4$TIS_cat_order[match(A1$Gene.name , A4$Gene.name)])
A1$Gene.TIS.count=A5$n[match(A1$Gene.name , A5$Gene.name)]

# 2 ################################################################   find the genes with multiple TIS
#categories:
#1.annot+alt_TIS
#2.multiple_alt-TIS
#3.multiple_annot-TIS
#4.one_TIS

table(A1$Multiple_TIS)
#table(A1$Multiple_TIS~A1$Gene.TIS.count)
#aggregate(Multiple_TIS~Gene.TIS.count, data=A1, length)

# A1$Multiple_TIS=ifelse(A1$Multiple_TIS %in% c("c(1, 2, 3, 4)", "c(1, 2, 4)", "c(1, 3, 4)","c(1, 4)","c(2, 4)", "c(3, 4)"),"1.annot+alt_TIS", A1$Multiple_TIS)
# A1$Multiple_TIS=ifelse((A1$Multiple_TIS %in% c("c(1, 2)","c(1, 3)", "c(2, 3)")) & A1$Gene.TIS.count>1,"3.multiple_annot-TIS", A1$Multiple_TIS)
# A1$Multiple_TIS=ifelse((A1$Multiple_TIS %in% c("1","2", "3")) & A1$Gene.TIS.count>1,"3.multiple_annot-TIS", A1$Multiple_TIS)
# A1$Multiple_TIS=ifelse((A1$Multiple_TIS %in% c("1","2", "3", "4")) & A1$Gene.TIS.count==1,"4.one_TIS", A1$Multiple_TIS)
# A1$Multiple_TIS=ifelse((A1$Multiple_TIS %in% c( "4")) & A1$Gene.TIS.count>1,"2.multiple_alt-TIS", A1$Multiple_TIS)

#change Multiple_TIS column according to Annelies request UniProtIso_aTIS is alt TIS
A1$Multiple_TIS=ifelse(A1$Multiple_TIS %in% c("c(1, 2, 3, 4)", "c(1, 2, 4)", "c(1, 3, 4)","c(1, 4)","c(2, 4)", "c(1, 3)", "c(2, 3)"),"1.annot+alt_TIS", A1$Multiple_TIS)
A1$Multiple_TIS=ifelse((A1$Multiple_TIS %in% c("c(1, 2)")) & A1$Gene.TIS.count>1,"3.multiple_annot-TIS", A1$Multiple_TIS)
A1$Multiple_TIS=ifelse((A1$Multiple_TIS %in% c("1","2")) & A1$Gene.TIS.count>1,"3.multiple_annot-TIS", A1$Multiple_TIS)
A1$Multiple_TIS=ifelse((A1$Multiple_TIS %in% c("1","2", "3", "4")) & A1$Gene.TIS.count==1,"4.one_TIS", A1$Multiple_TIS)
A1$Multiple_TIS=ifelse((A1$Multiple_TIS %in% c( "4", "3")) & A1$Gene.TIS.count>1,"2.multiple_alt-TIS", A1$Multiple_TIS)
A1$Multiple_TIS=ifelse(A1$Multiple_TIS %in% c( "c(3, 4)"),"2.multiple_alt-TIS", A1$Multiple_TIS)


# 3  ############################################# parse UniProt squence feature annotation (for canonical proteins)

B1=uni_annot[,colnames(uni_annot) %in% c("Entry","Signal.peptide","Transit.peptide","Propeptide", "Region", "Motif", "Coiled.coil", "Compositional.bias", "Repeat", "Zinc.finger")]
library(reshape2)
B2=melt(B1, id=1)
B2=B2[B2$value!="",]

#B3=sapply(B2$value, function(x) strsplit(x, "; (?=[A-Z]{4,})", perl=TRUE))
B3=sapply(B2$value, function(x) strsplit(x, "; (?=COILED|COMPBIAS|MOTIF|PROPEP|REGION|REPEAT|SIGNAL|TRANSIT|ZN_FING)", perl=TRUE))
table(sapply(B3, length))
names(B3)<-B2$Entry
B3[sapply(B3, length)>3]
B4=melt(B3)
B4$value=as.character(B4$value)
str(B4)
B5=B4
B5$Sequence.feature.type=sapply(B5$value, function(x) strsplit(x, " ")[[1]][1])
B5$Sequence.feature.name=sapply(B5$value, function(x) strsplit(x, "/note=")[[1]][2])
B5$Sequence.feature.name=sapply(B5$Sequence.feature.name, function(x) strsplit(x, ";")[[1]][1])
B5$Sequence.feature.start=sapply(B5$value, function(x) strsplit(x, "\\..")[[1]][1])
B5$Sequence.feature.start=sapply(B5$Sequence.feature.start, function(x) strsplit(x, " ")[[1]][2])
B5$Sequence.feature.stop=sapply(B5$value, function(x) strsplit(x, "\\..")[[1]][2])
B5$Sequence.feature.stop=sapply(B5$Sequence.feature.stop, function(x) strsplit(x, ";")[[1]][1])
B5$strand="+"
table(B5$Sequence.feature.start)
table(B5$Sequence.feature.stop)
B5[B5$Sequence.feature.stop=="?",]

B5$Sequence.feature.start=ifelse(B5$Sequence.feature.start=="?",1,B5$Sequence.feature.start)
B5$Sequence.feature.start=ifelse(B5$Sequence.feature.start=="?994",994,B5$Sequence.feature.start)
B5$Sequence.feature.stop=ifelse(grepl(";",B5$Sequence.feature.start),"start",B5$Sequence.feature.stop)
B5$Sequence.feature.start=ifelse(grepl(";",B5$Sequence.feature.start),gsub(";","",B5$Sequence.feature.start),B5$Sequence.feature.start)
B5$Sequence.feature.stop=ifelse(grepl("start",B5$Sequence.feature.stop),B5$Sequence.feature.start,B5$Sequence.feature.stop)
B5$Sequence.feature.stop=ifelse(B5$Sequence.feature.stop=="?"&B5$Sequence.feature.type%in%c("SIGNAL","TRANSIT"),10,B5$Sequence.feature.stop)
B5$Sequence.feature.stop=ifelse(B5$Sequence.feature.stop=="?",B5$Sequence.feature.start,B5$Sequence.feature.stop)
B5$Sequence.feature.stop=ifelse(grepl("\\?",B5$Sequence.feature.stop),gsub("\\?","",B5$Sequence.feature.stop),B5$Sequence.feature.stop)
B5$Sequence.feature.stop=ifelse(grepl(">",B5$Sequence.feature.stop),gsub(">","",B5$Sequence.feature.stop),B5$Sequence.feature.stop)


B5$Sequence.feature.start=as.numeric(B5$Sequence.feature.start)
B5$Sequence.feature.stop=as.numeric(B5$Sequence.feature.stop)

B5$combo=paste0(B5$Sequence.feature.type, " ",B5$Sequence.feature.start,"..",B5$Sequence.feature.stop, " note=", B5$Sequence.feature.name)
B5=B5[order(B5$L1,B5$Sequence.feature.start),]
table(B5$Sequence.feature.type)

# 4 ############################# add codon and frame info to ENST 5UTR, 3UTR, CDS entries
#ad unique IDs
A10=A1[order(A1$Multiple_TIS,A1$Gene.name,A1$TIS_cat_order,A1$Start),]
A10$ID=1:nrow(A10)

A10$codon=NA
A10$frame=NA
A10$dist.aa.to.aTIS=NA
A10$dist.nt.to.aTIS=NA

#for UniProt and UniProt Isoforms:
A10$codon=ifelse(A10$source!="ENST"&A10$Start==1,"AUG",A10$codon)
A10$codon=ifelse(A10$source!="ENST"&A10$Start>1 & (A10$First.AA=="M" | A10$Preceding.AA=="M"),"AUG",A10$codon)
A10$codon=ifelse(A10$source!="ENST"&A10$Start>1 
                 & A10$First.AA!="M" & A10$Preceding.AA!="M" & A10$secondPreceding.AA!="M","non-AUG",A10$codon)
A10$codon=ifelse(A10$source!="ENST"&A10$Start>2 
                 & A10$First.AA!="M" & A10$Preceding.AA!="M" & A10$secondPreceding.AA=="M","AUG",A10$codon)

A10$dipeptidase=ifelse(A10$Start>2 & A10$First.AA!="M" & A10$Preceding.AA!="M" & A10$secondPreceding.AA=="M",paste0(A10$secondPreceding.AA,A10$Preceding.AA),NA)

A10$frame=ifelse(A10$source!="ENST","in",A10$frame)
A10$dist.aa.to.aTIS=ifelse(A10$source!="ENST",A10$Start-1,A10$dist.aa.to.aTIS)
A10$dist.nt.to.aTIS=A10$dist.aa.to.aTIS*3

#for ENST use results.db export

A11=as.data.frame(cbind(A10$Accession[A10$source=="ENST"], A10$ID[A10$source=="ENST"]))
library(stringr)
A11$name=str_sub(A11$V1, end=-8)
A11$codon=codon.db$start_codon[match(A11$name,codon.db$name)]
A11$codon2=ifelse(A11$codon=="ATG", "AUG", "non-AUG")
A11$dist=codon.db$dist_to_aTIS[match(A11$name,codon.db$name)]
A11$frame=A11$dist%%3
A11$frame=ifelse(A11$frame==0,"in",ifelse(is.na(A11$frame),NA,"out"))

#Add to A10
A10$codon=ifelse(A10$source=="ENST"&A10$Start==1,A11$codon2[ match(A10$ID,A11$V2)],A10$codon)
A10$codon=ifelse(A10$source=="ENST"&A10$Start==2 & A10$First.AA!="M", A11$codon2[ match(A10$ID,A11$V2)],A10$codon)
A10$codon=ifelse(A10$source=="ENST"&A10$Start==3 & A10$First.AA!="M" & A10$Preceding.AA!="M", A11$codon2[ match(A10$ID,A11$V2)],A10$codon)

A10$codon=ifelse(A10$source=="ENST"&A10$Start>1  & A10$First.AA=="M","AUG",A10$codon)
A10$codon=ifelse(A10$source=="ENST"&A10$Start>2  & A10$First.AA!="M" & A10$Preceding.AA=="M","AUG",A10$codon)
A10$codon=ifelse(A10$source=="ENST"&A10$Start>2  & A10$First.AA!="M" & A10$Preceding.AA!="M","non-AUG",A10$codon)


A10$frame=ifelse(A10$source=="ENST",A11$frame[ match(A10$ID,A11$V2)],A10$frame)
A10$dist.aa.to.aTIS=ifelse(A10$source=="ENST" & A10$frame=="in",(A11$dist[ match(A10$ID,A11$V2)]/3)  +  (A10$Start-1), A10$dist.aa.to.aTIS)
A10$dist.nt.to.aTIS=ifelse(A10$source=="ENST" & A10$frame=="in",A10$dist.aa.to.aTIS*3, A10$dist.nt.to.aTIS)
A10$dist.nt.to.aTIS=ifelse(A10$source=="ENST" & A10$frame=="out",A11$dist[ match(A10$ID,A11$V2)]  +  (A10$Start-1)*3, A10$dist.nt.to.aTIS)

# 5 ####################################### map via dbtoolkit
map_uni<-read.csv(file="Data/Blast & uniprot exact map/C13merge_map_to_Unicanon.txt",header=TRUE, stringsAsFactors = FALSE, sep='\t')
map_uni=map_uni[map_uni$Sequence!="",]
map_uni$Isoforms=sapply(map_uni$Isoforms, function(x) gsub("\\^A",";", x))
map_uni$all=paste0(map_uni$Accession, " (", map_uni$Start,"-",map_uni$Stop,");",map_uni$Isoforms)
map_uni$all=ifelse(map_uni$Isoforms=="",sapply(map_uni$all, function(x) gsub(";","", x)),map_uni$all)
#seems to be selected on start and alphabetical

A10$map_uni=NA
A10$map_uni=map_uni$all[match(A10$Sequence,map_uni$Sequence)]

# 6 ########################################### proteoforms loose signal peptides etc.? match linear motifs
## 6A ### divide into several tables

#For ENST truncations use mapping.UniProt.canon (via dbtoolkit)
A21=A10[!is.na(A10$map_uni)& A10$source=="ENST",]
A21$trunc.start=1
A21$trunc.end=as.numeric(map_uni$Start[match(A21$Sequence,map_uni$Sequence)])-1
A21$trunc.ref=map_uni$Accession[match(A21$Sequence,map_uni$Sequence)]
A21$strand="+"

#for 5'extensions, check if ENST has UniProt and make a range 1 to 2 (to see if some features were N-terminal but no longer)
A24=A10[A10$dist.aa.to.aTIS<0,]
A24=A24[!is.na(A24$dist.aa.to.aTIS),]
A24=A24[is.na(A24$map_uni),]
A24$trunc.start=1
A24$trunc.end=2
A24$trunc.ref=biomart_uni$UniProtKB.Swiss.Prot.ID[match(A24$ENST , biomart_uni$Transcript.stable.ID)]
A24$strand="+"
A24=A24[A24$trunc.ref!="",]


#UniProtIsoforms use mapping.UniProt.canon (via dbtoolkit)
A22=A10[!is.na(A10$map_uni)& A10$source=="UniProtIso",]
A22$trunc.start=1
A22$trunc.end=as.numeric(map_uni$Start[match(A22$Sequence,map_uni$Sequence)])-1
A22$trunc.ref=map_uni$Accession[match(A22$Sequence,map_uni$Sequence)]
A22$strand="+"
table(A22$UniProt.stable.ID==A22$trunc.ref) #check if the accessions are the same
table(A22$trunc.end<2) #make sure no UniProt Isoform peptide get actually mapped at UniProt canonical aTIS

#for "UniProt_pos>2" make range of truncation.start =1 truncation.end =A23$start-1 truncation.ref=accession
A23=A10[A10$source=="UniProt" & A10$Start>2,]
A23$trunc.start=1
A23$trunc.end=A23$dist.aa.to.aTIS
A23$trunc.ref=A23$Accession
A23$strand="+"

## 6A cd ######## join these tables to make one genomic range
trunc.start=c(A21$trunc.start,A22$trunc.start,A23$trunc.start,A24$trunc.start)
length(trunc.start)
trunc.end=c(A21$trunc.end,A22$trunc.end,A23$trunc.end,A24$trunc.end)
length(trunc.end)
strand=c(A21$strand,A22$strand,A23$strand,A24$strand)
length(strand)
trunc.ref=c(A21$trunc.ref,A22$trunc.ref,A23$trunc.ref,A24$trunc.ref)
length(trunc.ref)
peptide=c(A21$Sequence ,A22$Sequence ,A23$Sequence ,A24$Sequence )
length(peptide)
ID=c(A21$ID ,A22$ID ,A23$ID ,A24$ID )
length(ID)

A25=as.data.frame(cbind(trunc.start,trunc.end, strand, trunc.ref, peptide, ID))
A25$trunc.start=as.numeric(as.character(A25$trunc.start))
A25$trunc.end=as.numeric(as.character(A25$trunc.end))

##### 6B ######################################################################################################
#for proteoforms not mapping to UniProt via dbtoolkit, pairwise alignment to selected UniProt
#in absence of any reliable UniProt match, use only prediction ELM
#if UniProt is aligned, use ELM predicted in lost and gained sequences
#ALSO: find UniProt annotated features overlapping gained sequences (i.e. the same analysis as 6A)

#for UniProtIso not mapping to UniProt use uniprot stable id
A10$proteoform.sequence=substring(A10$protein.sequence, A10$Start)

#for ENST not mapping to UniProt use uniprot for ESNG. If none UniProt found (e.g. ntr), report all ELM results
biomart_uni_simple=subset(biomart_uni, UniProtKB.Swiss.Prot.ID!="", select=c(Gene.stable.ID,UniProtKB.Swiss.Prot.ID) )
get.elm=data.frame("proteoform.sequence"= c(A10[is.na(A10$map_uni) & A10$source=="UniProtIso","proteoform.sequence"],
                                            A10[is.na(A10$map_uni) & A10$source=="ENST","proteoform.sequence"]),
                   "ID" = c(A10[is.na(A10$map_uni) & A10$source=="UniProtIso","ID"],
                            A10[is.na(A10$map_uni) & A10$source=="ENST","ID"]),
                   "sequence.feature.reference.ID" = c(A10[is.na(A10$map_uni) & A10$source=="UniProtIso","UniProt.stable.ID"],
                                                       biomart_uni_simple$UniProtKB.Swiss.Prot.ID[match(A10[is.na(A10$map_uni) & A10$source=="ENST","ENSG"],biomart_uni_simple$Gene.stable.ID )])
)
get.elm$proteoform.sequence=as.character(get.elm$proteoform.sequence)
get.elm$ID=as.character(get.elm$ID)
#load UniProt fasta to get  sequence.feature.reference  
library(Biostrings)
uniprot.fasta=readAAStringSet("Data/Uniprot_human_2019_04_canonical.fasta")
names(uniprot.fasta)<-sapply(strsplit(names(uniprot.fasta), split="\\|"), function(x) x[2])
get.elm$sequence.feature.reference=match(get.elm$sequence.feature.reference.ID, names(uniprot.fasta))
get.elm$sequence.feature.reference=ifelse(is.na(get.elm$sequence.feature.reference.ID),1,get.elm$sequence.feature.reference)
get.elm$sequence.feature.reference=as.character(uniprot.fasta[get.elm$sequence.feature.reference])
get.elm$sequence.feature.reference=ifelse(is.na(get.elm$sequence.feature.reference.ID),NA,get.elm$sequence.feature.reference)

## 6C ### ELM prediction: write Linux command to get prediction from ELM API
fileConn<-file("wget.proteoforms.sh")
writeLines(con=fileConn ,paste0("wget -O ",get.elm$ID,".elm.txt http://elm.eu.org/start_search/", get.elm$proteoform.sequence, " && sleep 2m &&"))
close(fileConn)
fileConn1<-file("wget.ref.sh")
writeLines(con=fileConn1 ,paste0("wget -O ",unique(get.elm$sequence.feature.reference.ID),".elm.txt http://elm.eu.org/start_search/", unique(get.elm$sequence.feature.reference), " && sleep 2m &&"))
close(fileConn1)

##### 6B cd ###  align proteoform to selected reference (use gapOpening -9.5, gapExtension -0.5 to prefer longer but fewer gaps)
#find gained in proteoform = insertions to reference (subject)
#find lost in proteoform = deletions in protoform (pattern)

library(Biostrings)
data(BLOSUM62)

#https://stat.ethz.ch/pipermail/bioconductor/2013-September/055194.html
# >It's a questionable design choice but the ranges describing the
#  >     deletions are reported with respect to the "original" (aka unaligned)
#  >     pattern, not to the aligned pattern. For example the 1st range you
#  >     see in 'deletion(pa)' means there is a deletion of 3 letters (when
#  >     going from subject to pattern) and that this deletion starts at position
#  >     4 in the original pattern. 
# >     As you noticed, those ranges can be shifted to make them refer to the
# >     aligned pattern:
# >
pa.aa=list()
offset.del=list()
DEL=list()
offset.ins=list()
INS=list()
lost.ranges=list()
gain.ranges=list()

for (i in 1:nrow(get.elm)){
  if (is.na(get.elm$sequence.feature.reference[i]))
  {
    pa.aa[i]=NA
    offset.del[i]=NA
    DEL[i]=NA
    offset.ins[i]=NA
    INS[i]=NA
    lost.ranges[i]=NA
    gain.ranges[i]=NA
  }
  else
  {
    pa.aa[[i]]<-pairwiseAlignment(pattern = get.elm$proteoform.sequence[i], 
                                  subject = get.elm$sequence.feature.reference[i],
                                  substitutionMatrix = BLOSUM62, gapOpening = -9.5, gapExtension =-0.5, 
                                  scoreOnly = FALSE, type="overlap")
    offset.del[[i]] <- c(0L, cumsum(head(width(unlist(deletion(pa.aa[[i]]))), n=-1)))
    DEL[[i]]=shift(unlist(deletion(pa.aa[[i]])), offset.del[[i]])
    offset.ins[[i]] <- c(0L, cumsum(head(width(unlist(insertion(pa.aa[[i]]))), n=-1)))
    INS[[i]]=shift(unlist(insertion(pa.aa[[i]])), offset.ins[[i]])
    lost.ranges[[i]]=setdiff(DEL[[i]], INS[[i]])
    gain.ranges[[i]]=setdiff(INS[[i]] ,DEL[[i]])
    
  }
}
#Explanation: DEL= deletion(range vs. alignment); INS=insertion(range vs. alignment)
#Explanation: lost.ranges=DEL-INS; gain.ranges=INS-DEL; most likely not needed 
#(as these never overlap i.e. region of alignment cannot be an ins AND del at the same time; added just in case)

for (i in 1:length(pa.aa)){
  if (!is.na(pa.aa[[i]])){
    writePairwiseAlignments(pa.aa[[i]],
                            file=paste0("prot_align/",
                                        paste0("proteoformID",get.elm$ID, "to", get.elm$sequence.feature.reference.ID,".aln"))[i], 
                            block.width = 100)
  }
}

#### 6B cd ### #get ranges of lost sequence in UniProt reference and gained sequence in proteoform
library(GenomicRanges)
lost.ref.seq=GRangesList()
gain.pform.seq=GRangesList()
for (i in 1:nrow(get.elm)){
  if (is.na(get.elm$sequence.feature.reference[i]))
  {
    lost.ref.seq[[i]]=GRanges(IRanges(start=0, end=0), seqnames="None")
    gain.pform.seq[[i]]=GRanges(ranges=IRanges(start=1, width=nchar(get.elm$proteoform.sequence[i])),
                                seqnames=get.elm$ID[i])
  }
  else
  {
    #PATTERN = PROTEOFORM part included in alignment
    pattern=as.character(pattern(pa.aa[[i]]))
    pat.split=strsplit(pattern, "")[[1]]
    count.pat=0
    pat.pos=list()
    for (n in 1:nchar(pattern)){
      if (pat.split[n]=="-"){
        pat.pos[n]="-"
      }
      else{
        pat.pos[n]=count.pat+1
        count.pat=count.pat+1
      }
    }
    if (length(gain.ranges[[i]])!=0){
      
      gain.pform.seq[[i]]=GRanges(IRanges(start=unlist(pat.pos[start(gain.ranges[[i]])]), width=width(gain.ranges[[i]])),
                                  seqnames=get.elm$ID[i])
    }
    else {gain.pform.seq[[i]]=GRanges(IRanges(start=0, end=0), seqnames=get.elm$ID[i])}
    
    #SUBJECT = REF UNIPROT part included in alignment
    subject=as.character(subject(pa.aa[[i]]))
    sub.split=strsplit(subject, "")[[1]]
    count.sub=0
    sub.pos=list()
    for (n in 1:nchar(subject)){
      if (sub.split[n]=="-"){
        sub.pos[n]="-"
      }
      else{
        sub.pos[n]=count.sub+1
        count.sub=count.sub+1
      }
    }
    if (length(lost.ranges[[i]])!=0){
      lost.ref.seq[[i]]=GRanges(IRanges(start=unlist(sub.pos[start(lost.ranges[[i]])]), width=width(lost.ranges[[i]])),
                                seqnames=get.elm$sequence.feature.reference.ID[i])
    }
    else {lost.ref.seq[[i]]=GRanges(IRanges(start=0, end=0), seqnames=get.elm$sequence.feature.reference.ID[i])}
    
  }
}

#we know that N-terminus (not matching in proteoform vs. ref) is not included in alignment
#test if C-terminus in matching
library(stringr)
table(str_sub(get.elm$proteoform.sequence[!is.na(get.elm$sequence.feature.reference)], start= -5)
      ==str_sub(get.elm$sequence.feature.reference[!is.na(get.elm$sequence.feature.reference)], start= -5))
#C-terminus is also sometimes excluded from alignment

## 6B cd ### add N-terminal ranges  and C-terminal ranges excluded from alignment to gained  and lost
lost.ref.seq.backup=lost.ref.seq
gain.pform.seq.backup=gain.pform.seq
for (i in 1:nrow(get.elm)){
  if (!is.na(pa.aa[[i]])){
    #N- and C- terminal parts of gained proteoform seq
    st.pat=start(pa.aa[[i]]@pattern@range)
    end.pat=end(pa.aa[[i]]@pattern@range)
    
    if (st.pat!=1){
      gain.pform.seq[[i]]=c(GRanges(seqnames=seqnames(gain.pform.seq[[i]])[1],
                                    IRanges(start=1, end=st.pat-1)), shift(gain.pform.seq[[i]], st.pat-1 ))
    }
    
    if (end.pat!=nchar(get.elm$proteoform.sequence[i])){
      gain.pform.seq[[i]]=c(gain.pform.seq[[i]], GRanges(seqnames=seqnames(gain.pform.seq[[i]])[1], 
                                                         IRanges(start=end.pat+1, end=nchar(get.elm$proteoform.sequence[i]))))
    }
    #N- and C- terminal parts of lost ref seq
    st.sub=start(pa.aa[[i]]@subject@range)
    end.sub=end(pa.aa[[i]]@subject@range)
    
    if (st.sub!=1){
      lost.ref.seq[[i]]=c(GRanges(seqnames=seqnames(lost.ref.seq[[i]])[1],
                                  IRanges(start=1, end=st.sub-1)), shift(lost.ref.seq[[i]], st.sub-1 ))
    }
    
    if (end.sub!=nchar(get.elm$sequence.feature.reference[i])){
      lost.ref.seq[[i]]=c(lost.ref.seq[[i]], GRanges(seqnames=seqnames(lost.ref.seq[[i]])[1],
                                                     IRanges(start=end.sub+1, end=nchar(get.elm$sequence.feature.reference[i]))))
    }
    
  }
}

names(lost.ref.seq)<-as.character(get.elm$ID)
names(gain.pform.seq)<-as.character(get.elm$ID)

######################################################################################################
## 6C ### ELM prediction: gather results, transform to GRanges
path <- "Data/ELM/with C13/WITH-C13-ready/"

#load files
lst=list()
library("stringr")
files=list.files(path,pattern="\\.elm.txt$",full.names=F) #Make a list of files
filen=str_extract(files, '.*(?=\\.elm.txt)') #Pretty the file names for object names
for (i in 1:length(files)){
  lst[[i]]=read.table(paste(path,files[i], sep=''), header=TRUE, sep="\t")#Load the files
  names(lst)[i]<-filen[i] #Name the entries
  lst[[i]]=lst[[i]][,order(names(lst[[i]]))]
}
#add human readable names
elm.names=read.csv(file="Data/ELM/elm_classes.tsv", sep="\t", header = TRUE, comment.char = "#")
for (i in 1:length(lst)){
  lst[[i]]$elm_name=elm.names$FunctionalSiteName[match(lst[[i]]$elm_identifier, elm.names$ELMIdentifier)]
}
#make GRanges
elm.gr=GRangesList()
for (i in 1:length(lst)){
  if(length(lst[[i]]$start)==0){elm.gr[[i]]=GRanges()}
  else
    elm.gr[[i]]=GRanges(seqnames=names(lst[i]), IRanges(start=lst[[i]]$start, end=lst[[i]]$stop, 
                                                        mcols=data.frame(elm_name=lst[[i]]$elm_name, elm_ID=lst[[i]]$elm_identifier)))
}
elm.ranges=unlist(elm.gr)

## overlap lost with elm
lost.fo=list()
lost.df=list()
for (i in 1:length(lost.ref.seq)){
  lost.fo[[i]]=findOverlaps(lost.ref.seq[[i]], elm.ranges)
  lost.ov <- pintersect(lost.ref.seq[[i]][queryHits(lost.fo[[i]])], elm.ranges[subjectHits(lost.fo[[i]])])
  lost.perc.ov <- width(lost.ov) / width(elm.ranges[subjectHits(lost.fo[[i]])])
  lost.fo[[i]] <- lost.fo[[i]][lost.perc.ov ==1]
  if (length(lost.fo[[i]])==0) {get.elm$Lost_ELM[i]=NA }
  else {
    lost.df[[i]]=as.data.frame(lost.fo[[i]])
    lost.df[[i]]$elm_name=elm.ranges$mcols.elm_name[lost.df[[i]]$subjectHits]
    lost.df[[i]]$elm_ID=elm.ranges$mcols.elm_ID[lost.df[[i]]$subjectHits]
    lost.df[[i]]$elm_from_protein=as.character(seqnames(elm.ranges[lost.df[[i]]$subjectHits]))
    lost.df[[i]]$start=start(elm.ranges[lost.df[[i]]$subjectHits])
    lost.df[[i]]$end=end(elm.ranges[lost.df[[i]]$subjectHits])
    lost.df[[i]]$combo=paste0(lost.df[[i]]$start,"..",lost.df[[i]]$end)
    lost.agg=aggregate(combo~elm_name+elm_from_protein,FUN= function(x) paste(x,collapse=";"),data=lost.df[[i]])
    get.elm$Lost_ELM[i]=paste0(get.elm$sequence.feature.reference.ID[i],": ",paste(paste0(lost.agg$elm_name," " , lost.agg$combo), collapse="#"))
  }
}

## overlap gained with elm

gain.fo=list()
gain.df=list()
for (i in 1:length(gain.pform.seq)){
  gain.fo[[i]]=findOverlaps(gain.pform.seq[[i]], elm.ranges)
  gain.ov <- pintersect(gain.pform.seq[[i]][queryHits(gain.fo[[i]])], elm.ranges[subjectHits(gain.fo[[i]])])
  gain.perc.ov <- width(gain.ov) / width(elm.ranges[subjectHits(gain.fo[[i]])])
  gain.fo[[i]] <- gain.fo[[i]][gain.perc.ov ==1]
  if (length(gain.fo[[i]])==0) {get.elm$Gained_ELM[i]=NA }
  else {
    gain.df[[i]]=as.data.frame(gain.fo[[i]])
    gain.df[[i]]$elm_name=elm.ranges$mcols.elm_name[gain.df[[i]]$subjectHits]
    gain.df[[i]]$elm_ID=elm.ranges$mcols.elm_ID[gain.df[[i]]$subjectHits]
    gain.df[[i]]$elm_from_protein=as.character(seqnames(elm.ranges[gain.df[[i]]$subjectHits]))
    gain.df[[i]]$start=start(elm.ranges[gain.df[[i]]$subjectHits])
    gain.df[[i]]$end=end(elm.ranges[gain.df[[i]]$subjectHits])
    gain.df[[i]]$combo=paste0(gain.df[[i]]$start,"..",gain.df[[i]]$end)
    gain.agg=aggregate(combo~elm_name+elm_from_protein,FUN= function(x) paste(x,collapse=";"),data=gain.df[[i]])
    get.elm$Gained_ELM[i]=paste0("proteoformID ", get.elm$ID[i],": ",paste(paste0(gain.agg$elm_name," " , gain.agg$combo), collapse="#"))
  }
}

#######################################################################################################
## 6A + 3 ### overlap truncated regions with protein sequence features parsed from UniProt database
library(GenomicRanges)
motif=makeGRangesFromDataFrame(B5,seqnames.field = "L1", start.field="Sequence.feature.start", end.field="Sequence.feature.stop",strand.field="strand")
#add metacolumns
mcols(motif)<-B5[,c("value", "Sequence.feature.type", "Sequence.feature.name", "combo","Sequence.feature.stop" )]

A25GR=makeGRangesFromDataFrame(A25,seqnames.field = "trunc.ref", start.field="trunc.start", end.field="trunc.end",strand.field="strand")
mcols(A25GR)<-A25[,c("ID", "peptide", "trunc.end", "trunc.ref")]


fo=findOverlaps(A25GR, motif)
fo.df=as.data.frame(fo)

overlaps <- pintersect(A25GR[queryHits(fo)], motif[subjectHits(fo)])
percentOverlap <- width(overlaps) / width(motif[subjectHits(fo)])
hits <- fo[percentOverlap > 0.5]
hits.df=as.data.frame(hits)
hits.df$combo=motif$combo[hits.df$subjectHits]
hits.df$ID=A25GR$ID[hits.df$queryHits]
hits.df$sequence.feature.reference.ID=A25GR$trunc.ref[hits.df$queryHits]
lost.features=aggregate(combo ~ ID+sequence.feature.reference.ID , data = hits.df, paste, collapse="#")

A20=A10
A20$lost.Nterm.sequence.feature=lost.features$combo[match(A20$ID,lost.features$ID)]
#change to report all reference IDs, even if no feature found
#A20$sequence.feature.reference.ID=lost.features$sequence.feature.reference.ID[match(A20$ID,lost.features$ID)]
A20$sequence.feature.reference.ID=as.character(A25$trunc.ref[match(A20$ID,A25$ID)])

#find if peptide can be direct effect of N-terminal processing (removed signal, transit or pro- peptide)
motif.spt=motif[motif$Sequence.feature.type %in% c("SIGNAL", "PROPEP", "TRANSIT")]
spt.end=findOverlaps(A25GR, motif.spt, type="end")
length(spt.end) #1 
#find if peptide can be approximately an effect of Nt processing (+/-1aa)
spt=findOverlaps(A25GR, motif.spt)
spt.df=as.data.frame(spt)
spt.df$ID=A25GR$ID[spt.df$queryHits]
spt.df$feature.end=motif.spt$Sequence.feature.stop[spt.df$subjectHits]
spt.df$trunc.end=A25GR$trunc.end[spt.df$queryHits]
spt.df$feature.type=motif.spt$Sequence.feature.type[spt.df$subjectHits]
spt.df$dist.signal.transit.propeptide=spt.df$feature.end-spt.df$trunc.end
spt.df$potential.signal.transit.propeptide.processing=ifelse(spt.df$dist.signal.transit.propeptide %in% c(-1,0,1), paste0(spt.df$feature.type, "_",spt.df$dist.signal.transit.propeptide),NA)


A20$signal.transit.propeptide.processing_distance=spt.df$potential.signal.transit.propeptide.processing[match(A20$ID,spt.df$ID)]
#motif.spt[seqnames(motif.spt)=="Q8WUF8"]
#motif.spt[3761]
#A25GR[seqnames(A25GR)=="Q8WUF8"]
#A25GR[628]

## 6D #####
## for proteoforms not mapping via dbtoolkit, take alignment-defined lost UniProt regions and overlap with UniProt features
unilost.fo=list()
for (i in 1:length(lost.ref.seq)){
  unilost.fo[[i]]=findOverlaps(lost.ref.seq[[i]], motif)
  unilost.ov <- pintersect(lost.ref.seq[[i]][queryHits(unilost.fo[[i]])], motif[subjectHits(unilost.fo[[i]])])
  unilost.perc.ov <- width(unilost.ov) / width(motif[subjectHits(unilost.fo[[i]])])
  unilost.fo[[i]] <- unilost.fo[[i]][unilost.perc.ov ==1]
  if (length(unilost.fo[[i]])==0) {get.elm$Lost_UniProt_feature[i]=NA }
  else {
    get.elm$Lost_UniProt_feature[i]=paste(motif$combo[subjectHits(unilost.fo[[i]])], collapse="#")
  }
}

A20$lost.Nterm.sequence.feature=ifelse(is.na(A20$lost.Nterm.sequence.feature),
                                       get.elm$Lost_UniProt_feature[match(A20$ID,get.elm$ID)],
                                       A20$lost.Nterm.sequence.feature)
A20$sequence.feature.reference.ID=ifelse(is.na(A20$sequence.feature.reference.ID),
                                         as.character(get.elm$sequence.feature.reference.ID[match(A20$ID,get.elm$ID)]),
                                         A20$sequence.feature.reference.ID)

## 6E ###
# add information from ELMto A20 table
A20$lost_ELM=get.elm$Lost_ELM[match(A20$ID, get.elm$ID)]
A20$gained_ELM=get.elm$Gained_ELM[match(A20$ID, get.elm$ID)]
#

# 7 ############################# TopFIND - look for evidence of degradation
A26=data.frame(sapply(A20$map_uni, function(x) strsplit(x," ")[[1]][1]), A20$Sequence)
A26=A26[!is.na(A26[,1]),]
write.table(A26, file="C13_TopFind_input.txt", row.names = FALSE,  quote=FALSE, sep="\t", col.names = FALSE)

top<-read.csv(file="Data/TopFIND/2020_07_07_C13_dbtoolkitIDs0aaprecision/Full_Table.txt",header=TRUE, stringsAsFactors = FALSE, sep='\t')
top$Cleaving.proteases[top$Cleaving.proteases==""] <- NA
top$Cleaving.proteases[top$Cleaving.proteases=="MAP2"] <- NA
top$Cleaving.proteases[top$Cleaving.proteases=="MAP12"] <- NA

A20$Cleaving.proteases=top$Cleaving.proteases[match(A20$Sequence,top$Input.Sequence)]
A20$Distance.to.last.transmembrane.domain.shed=top$Distance.to.last.transmembrane.domain..shed.[match(A20$Sequence,top$Input.Sequence)]

# 8 ############################ Blast 



library(tidyverse)
library(seqinr)
seqs = as.list(dplyr::pull(A20[is.na(A20$map_uni),], Sequence))
names = dplyr::pull(A20[is.na(A20$map_uni),], Sequence)
write.fasta(seqs, names, "blast-C13.fasta",
            open = "w", as.string = FALSE)

#use cutoff like openprot:
#over 80% of protein sequence identity over 50% of the length (Basic Local Alignment Search Tool (BLAST)
# Program	blastp
# Word size	6
# Expect value	10
# Hitlist size	100
# Gapcosts	11,1
# Matrix	BLOSUM62
# Filter string	F
# Genetic Code	1
# Window Size	40
# Threshold	21
# Composition-based stats	2 
# 
# Database either human UniProt or human non-redundant
# UniProt:
# Posted date	Jul 5, 2020 4:28 AM
# Number of letters	11,364,533
# Number of sequences	20,386
# Entrez query	
# Includes: Homo sapiens (taxid:9606)
# ntr:
# Database
# Posted date	Jul 5, 2020 1:22 AM
# Number of letters	105,805,817,297
# Number of sequences	293,844,859
# Entrez query	
# Includes: human (taxid:9606)

#download Hits (csv) table from pblast website and load it to R

blast.names=c( "query acc.ver",	 "subject acc.ver",	 "% identity",	 "alignment length",	 "mismatches",	 "gap opens",	 "q. start",	 "q. end",	 "s. start",	 "s. end",	 "evalue",	 "bit score",	 "% positives")

blast.uni<-read.csv(file="Data/Blast & uniprot exact map/blast_vs_uniprot_human/G8RK09SM014-Alignment-HitTable.csv",header=FALSE, stringsAsFactors = FALSE, sep=',')

colnames(blast.uni)<-blast.names

blast.uni$perc.cover=round(100*(blast.uni$`q. end`-blast.uni$`q. start`+1) /sapply(blast.uni$`query acc.ver`, nchar),0)
blast.uni$perc.ide=round(blast.uni$`% identity`,0)
blast.uni$summary=paste0(blast.uni$`subject acc.ver`,"_", blast.uni$perc.ide,"%ide_",  blast.uni$perc.cover, "%cov_",  "Eval=",blast.uni$evalue)

blast.uni.agg=aggregate(.~`query acc.ver`, blast.uni, FUN=head,1)
blast.uni.agg=blast.uni.agg[as.numeric(as.character(blast.uni.agg$perc.cover))>=50 & as.numeric(as.character(blast.uni.agg$perc.ide))>=80,]


blast.nr<-read.csv(file="Data/Blast & uniprot exact map/blast_vs_nonredund_human/G8S3NUDH016-Alignment-HitTable.csv",header=FALSE, stringsAsFactors = FALSE, sep=',')

colnames(blast.nr)<-blast.names

blast.nr$perc.cover=round(100*(blast.nr$`q. end`-blast.nr$`q. start`+1) /sapply(blast.nr$`query acc.ver`, nchar),0)
blast.nr$perc.ide=round(blast.nr$`% identity`,0)
blast.nr$summary=paste0(blast.nr$`subject acc.ver`,"_", blast.nr$perc.ide,"%ide_",  blast.nr$perc.cover, "%cov_",  "Eval=",blast.nr$evalue)

blast.nr.agg=aggregate(.~`query acc.ver`, blast.nr, FUN=head,1)
blast.nr.agg=blast.nr.agg[as.numeric(as.character(blast.nr.agg$perc.cover))>=50 & as.numeric(as.character(blast.nr.agg$perc.ide))>=80,]

A20$blast.vs.uniprot=blast.uni.agg$summary[match(A20$Sequence,blast.uni.agg$`query acc.ver`)]
A20$blast.vs.nonredundant=blast.nr.agg$summary[match(A20$Sequence,blast.nr.agg$`query acc.ver`)]

# 10 ################## BioGrid
biogrid<-read.csv(file="Data/biogrid/BIOGRID-ORGANISM-Homo_sapiens-3.5.182.tab2.txt",header=TRUE, stringsAsFactors = FALSE, sep='\t')

#unique interactions per row
#get tabla for aggregation (each interaction once for A and once for B)

biogrid2=cbind.data.frame(c(biogrid$Official.Symbol.Interactor.A, biogrid$Official.Symbol.Interactor.B),
                          c(biogrid$Entrez.Gene.Interactor.A,biogrid$Entrez.Gene.Interactor.B),
                          c(biogrid$Official.Symbol.Interactor.B, biogrid$Official.Symbol.Interactor.A),
                          c(biogrid$Entrez.Gene.Interactor.B,biogrid$Entrez.Gene.Interactor.A))
colnames(biogrid2)<-c("Official.Symbol.Interactor.A","Entrez.Gene.Interactor.A","Official.Symbol.Interactor.B","Entrez.Gene.Interactor.B")
biogrid2$Official.Symbol.Interactor.A=as.character((biogrid2$Official.Symbol.Interactor.A))
biogrid2$Official.Symbol.Interactor.B=as.character((biogrid2$Official.Symbol.Interactor.B))

biogrid3=aggregate(Official.Symbol.Interactor.B ~ Official.Symbol.Interactor.A+Entrez.Gene.Interactor.A , data = biogrid2, FUN = function(x) unique(sort(x)))
biogrid3$result=paste0("bioGRID-interactions: ",lapply(biogrid3$Official.Symbol.Interactor.B, length),"; partners: ",biogrid3$Official.Symbol.Interactor.B)
biogrid3$result=gsub("c\\(","",biogrid3$result)
biogrid3$result=gsub("\\)","",biogrid3$result)
biogrid3$result=gsub("\"","",biogrid3$result)


A20$biogrid=biogrid3$result[match(A20$Gene.name, biogrid3$Official.Symbol.Interactor.A)]

# # BiocManager::install("biomaRt")
  library(biomaRt)
#  mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl",  host = "http://sep2019.archive.ensembl.org")
listAttributes(mart, what = c("name","description","page"))
biomart_additional=getBM(attributes = c("ensembl_gene_id","mim_gene_description","mim_morbid_description", "entrezgene_id"),
                         filters="ensembl_gene_id",
                         values=A20$ENSG,
                         mart=mart)
# 72                              mim_gene_description
# 73                                mim_gene_accession
# 74                            mim_morbid_description
# 75                              mim_morbid_accession
# 81                                     entrezgene_id

biogrid3$ENSG=biomart_additional$ensembl_gene_id[match(biogrid3$Entrez.Gene.Interactor.A, biomart_additional$entrezgene_id)]
A20$biogrid[is.na(A20$biogrid)]=biogrid3$result[match(A20$ENSG[is.na(A20$biogrid)], biogrid3$ENSG)]

# 10 ################# OMIM

biomart_additional2=aggregate(mim_morbid_description ~ ensembl_gene_id+entrezgene_id , data = biomart_additional, FUN = function(x) unique(x))
biomart_additional2$omim=paste0(biomart_additional2$mim_morbid_description)
biomart_additional2$omim=gsub("c\\(","",biomart_additional2$omim)
biomart_additional2$omim=gsub("\","," #",biomart_additional2$omim)
biomart_additional2$omim=gsub("\\)","",biomart_additional2$omim)
biomart_additional2$omim=gsub("\"","",biomart_additional2$omim)

A20$omim=biomart_additional2$omim[match(A20$ENSG, biomart_additional2$ensembl_gene_id)]
sum(is.na(A20$omim)) #25
sum(A20$omim=="", na.rm=TRUE) #2504
A20$omim[is.na(A20$omim)] <-""

#memory issues
rm(c.fa)

#12 ################### data Petra PMID: 24623590 
mcp<-read.csv(file="Data/data_Petra/peptides_Petra.txt",header=FALSE, stringsAsFactors = FALSE, sep='\t')
#A20$VanDamme.2014=A20$Sequence %in% mcp$V1
mcp.add1=substring(mcp$V1[substring(mcp$V1, 1,1)=="M"],2)
mcp.add2=paste0("M",mcp$V1[substring(mcp$V1, 1,1)!="M"])
mcp.join=c(mcp$V1,mcp.add1,mcp.add2)
#find the position of match(every peptide in MCP vs. every proteoform.sequenca Annelies)
mcp.match=sapply(mcp.join, regexpr, A20$proteoform.sequence, ignore.case=TRUE)
#match has to occur at position 1 to say that MCP found the same TIS
mcp.12=sapply(as.data.frame(mcp.match), function(x) x %in% c(1))
mcp.12=rowSums(mcp.12)>0
sum(mcp.12)#1584
sum(A20$VanDamme.2014)#1001
A20$VanDamme.2014=mcp.12


# 13 ################### molecular mass
#install.packages("Peptides")
library("Peptides")
A20$MolWeigth.kDa=round(mw(substring(A20$protein.sequence, A20$Start))/1000, digits = 1)
table(strsplit(A20$protein.sequence[725],'')[[1]])
#U=selenocysteine
A20$MW.difference.aTIS.kDa=ifelse(A20$source!="ENST",round(mw(A20$protein.sequence)/1000, digits = 1) - A20$MolWeigth.kDa,NA)
A20$MW.difference.aTIS.kDa=ifelse(grepl("ENST_aTIS",A20$TIScategory),round(mw(A20$protein.sequence)/1000, digits = 1) - A20$MolWeigth.kDa, A20$MW.difference.aTIS.kDa)
A20$MW.difference.aTIS.kDa=ifelse(grepl("ENST_5UTR",A20$TIScategory),round(mw(substring(A20$protein.sequence, -(A20$dist.aa.to.aTIS)+1))/1000, digits = 1) - A20$MolWeigth.kDa, A20$MW.difference.aTIS.kDa)
A20$MW.difference.aTIS.kDa=ifelse(grepl("ENST_CDS",A20$TIScategory),round(A20$dist.aa.to.aTIS*0.11, digits=1), A20$MW.difference.aTIS.kDa)
A20$MW.difference.aTIS.kDa=ifelse(grepl("ENST_5UTR",A20$TIScategory) & is.na(A20$dist.aa.to.aTIS), NA, A20$MW.difference.aTIS.kDa)

# 14 ############### cytosolic
cyt<-read.csv(file="Data/cytosol_shotgun/Cytosolic Map-GOCC-With Unique.txt",header=TRUE, stringsAsFactors = FALSE, sep='\t')
cyt$sel=grepl("cytosol", cyt$C..GOCC.slim.name) | grepl("cytoplasm", cyt$C..GOCC.slim.name)
cyt$sel.unique=ifelse(grepl("organelle",cyt$C..GOCC.slim.name),FALSE,cyt$sel )
cyt$prot=sapply(cyt$T..T..Protein.IDs, function(x) strsplit(x,";")[[1]][1])
cyt$ENSG=A20$ENSG[match(cyt$prot, A20$UniProt.stable.ID)]

cyt.sel=cyt$ENSG[cyt$sel]
cyt.sel<- cyt.sel[!is.na(cyt.sel)]
cyt.sel.unique=cyt$ENSG[cyt$sel.unique]
cyt.sel.unique<- cyt.sel.unique[!is.na(cyt.sel.unique)]

#human protein atlas
hpa<-read.csv(file="Data/HumanProteinAtlas/subcellular_location.tsv",header=TRUE, stringsAsFactors = FALSE, sep='\t')
hpa.sel=hpa$Gene[grepl("Cytosol",hpa$GO.id)]
cyt.fin.sel=intersect( cyt.sel, hpa.sel)
A20$cytosolic.in.HEK.proteome=ifelse(A20$ENSG %in% cyt.sel.unique, TRUE, ifelse(A20$ENSG %in% cyt.fin.sel, TRUE, FALSE))

# 15 ############################ make score for more interesting genes
#fix confidence
A20$TISconfidence[grepl("extEvid",A20$TIScategory)]<-"high-conf-TIS"


A20$interestig.gene.score=0

A20$interestig.gene.score=ifelse(A20$TISconfidence=="high-conf-TIS", A20$interestig.gene.score+1, A20$interestig.gene.score)
A20$interestig.gene.score=ifelse(A20$perc.Ace.TIS>=50, A20$interestig.gene.score+1, A20$interestig.gene.score)
A20$interestig.gene.score=ifelse(A20$sum.spectral.count.TIS >1, A20$interestig.gene.score+1, A20$interestig.gene.score)
A20$interestig.gene.score=ifelse(A20$Shorter.Sequences.found!="", A20$interestig.gene.score+1, A20$interestig.gene.score)
A20$interestig.gene.score=ifelse(A20$TIS_in_anyDB!="", A20$interestig.gene.score+1, A20$interestig.gene.score)
A20$interestig.gene.score=ifelse(grepl(";",A20$Enzyme), A20$interestig.gene.score+1, A20$interestig.gene.score)
A20$interestig.gene.score=ifelse(A20$codon=="non-AUG" & grepl("ENST_5UTR",A20$TIScategory), A20$interestig.gene.score+1, A20$interestig.gene.score)
A20$interestig.gene.score=ifelse(A20$VanDamme.2014, A20$interestig.gene.score+1, A20$interestig.gene.score)

#truncation =<50% protein length and domain lost = +3 points
A20$interestig.gene.score=ifelse(!is.na(A20$lost.Nterm.sequence.feature) & !is.na(A20$dist.aa.to.aTIS) & (A20$dist.aa.to.aTIS/nchar(A20$protein.sequence))<=0.5 , A20$interestig.gene.score+3, A20$interestig.gene.score)
# ELM lost or gained +1
A20$interestig.gene.score=ifelse(!is.na(A20$lost_ELM) | !is.na(A20$gained_ELM), A20$interestig.gene.score+1, A20$interestig.gene.score)
#score drops to 0 if:
A20$interestig.gene.score=ifelse(!is.na(A20$dipeptidase), 0, A20$interestig.gene.score)
A20$interestig.gene.score=ifelse(!is.na(A20$signal.transit.propeptide.processing_distance), 0, A20$interestig.gene.score)
A20$interestig.gene.score=ifelse(!is.na(A20$Cleaving.proteases), 0, A20$interestig.gene.score)
#keep scores for alternative TIS only
A20$interestig.gene.score=ifelse(A20$TIS_cat_order %in% c(3,4),  A20$interestig.gene.score,0)
#keep scores for genes 1.annot+alt_TIS or 4_one_TIS, 2_multiple_altTIS if cytosolic
A20$interestig.gene.score=ifelse(A20$Multiple_TIS=="1.annot+alt_TIS" |
                                   (A20$Multiple_TIS %in% c("2.multiple_alt-TIS", "4.one_TIS")
                                    & A20$cytosolic.in.HEK.proteome), A20$interestig.gene.score, 0)
A20$interestig.TIS.score=A20$interestig.gene.score


A27=aggregate(interestig.gene.score~Gene.name, A20, FUN=sum)

A20$interestig.gene.score=A27$interestig.gene.score[match(A20$Gene.name,A27$Gene.name)]

#add to gene score if there is OMIM disease or if there are no interactions known in bioGRID
A20$interestig.gene.score=ifelse(A20$interestig.gene.score>0 & A20$omim!="",A20$interestig.gene.score+1, A20$interestig.gene.score )
A20$interestig.gene.score=ifelse(A20$interestig.gene.score>0 & is.na(A20$biogrid),A20$interestig.gene.score+1, A20$interestig.gene.score )



A20.copy=A20
colnames(A20)[colnames(A20)=="lost.Nterm.sequence.feature"]<- "lost.UniProtDB.sequence.feature"
A20$exon_rank_TIS=sapply(A20$exon_rank_TIS, function(x) paste(x[1]))
A20$exon_rank_of_corresponding_aTIS=sapply(A20$exon_rank_of_corresponding_aTIS, function(x) paste(x[1]))
table(paste0(sapply(A1$exon_rank_TIS, function(x) paste(x[1])), A1$Sequence) %in% paste0(A20$exon_rank_TIS, A20$Sequence))
#TRUE 
#3305 

A30=A20[A20$interestig.gene.score>=1,]
A30=A30[order(A30$interestig.gene.score, decreasing = TRUE),]
hist(log2(A20$interestig.gene.score))
hist((A20$interestig.gene.score), breaks=0:max(A20$interestig.gene.score))

write.table(A20, file="07072020-C13-FullMergedResults-R_output.txt", row.names = FALSE,  quote=TRUE, sep="\t", col.names = TRUE)
write.table(A30, file="07072020-MOST-INTERESTING-C13-FullMergedResults-R_output.txt", row.names = FALSE,  quote=TRUE, sep="\t", col.names = TRUE)








