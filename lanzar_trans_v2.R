library(pureseqtmr)
library(Gviz)
library(GenomicRanges)
library(data.table)
library(stringr)
library("msa")
library("ape")
library("seqinr")
library(ggtree)
library(ggplot2)


options(ucscChromosomeNames=FALSE)

isEmpty <- function(x) { #This function checks if a data frame is empty or not
  return(length(x)==0)
}

direct="/storage/parras/CMV/DataSet"

for(pedazos in 1:length(list.files(direct)))
{
# Use an example proteome

#pedazos=1

MiCMV=read.delim(paste(direct,list.files(direct)[pedazos],sep="/"),header = F)

fasta_filename=paste(direct,list.files(direct)[pedazos],sep="/")
# Predict the topology
topology <- predict_topology(fasta_filename)

topology2=topology
o=0
res_location_trans=NULL


for(i in 1:length(topology$topology))
{
	if(!isEmpty(grep("1",topology[i,2])))
	{
		#print(topology[i,2])
		o=o+1
		topology2[o,]=topology[i,]
		
	
		dist_aux=strsplit(strsplit(as.character(topology[i,1]), "\\location=")[[1]][2],"\\]")[[1]][1]
		dist_aux=str_replace_all(dist_aux,"<","")
		dist_aux=str_replace_all(dist_aux,">","")

		gene_name=NULL
		seq_1=NULL

		if(!isEmpty(grep("gene=",as.character(topology[i,1]))))
		{
			gene_name=strsplit(strsplit(as.character(topology[i,1]),"gene=")[[1]][2],"\\]")[[1]][1]
		}else
		{
			gene_name=paste(strsplit(strsplit(as.character(topology[i,1]),"protein=")[[1]][2],"\\]")[[1]][1],strsplit(strsplit(as.character(topology[i,1]),"protein_id=")[[1]][2],"\\]")[[1]][1],sep=";")

		}

  		if(!isEmpty(grep("complement",dist_aux)))
  		{
    			if(!isEmpty(grep("join",dist_aux)))
    			{ 
     				dist_aux=strsplit(strsplit(as.character(dist_aux), "\\.")[[1]][1],"\\(")[[1]][3] 
    			}else
    			{
      				dist_aux=strsplit(strsplit(as.character(dist_aux), "\\.")[[1]][1],"\\(")[[1]][2]
    			}
  		}else
  		{
    			if(!isEmpty(grep("join",dist_aux)))
    			{ 
      				dist_aux=strsplit(strsplit(as.character(dist_aux), "\\.")[[1]][1],"\\(")[[1]][2]
    			}else
    			{
      				dist_aux=strsplit(as.character(dist_aux), "\\.")[[1]][1]
    			}
  		}

		dist_aux=as.numeric(dist_aux)

		seq_1=sort(c(((grep("1",strsplit(as.character(topology[i,2]),"")[[1]]))*3),((grep("1",strsplit(as.character(topology[i,2]),"")[[1]]))*3)+1,((grep("1",strsplit(as.character(topology[i,2]),"")[[1]]))*3)+2))

		band_dist=0

		save_trans=data.frame(matrix(ncol=3))
		e=1

		for(u in 1:(length(seq_1)-1))
		{
			if(band_dist==0)
			{
				total_inicial=dist_aux+seq_1[u]
				band_dist=1
			}

			if(!(dist_aux+seq_1[u+1]==(dist_aux+seq_1[u])+1))
			{
				total_final=dist_aux+seq_1[u]
				print(paste(i,gene_name,total_inicial,total_final))
				save_trans[e,]=c(gene_name,total_inicial,total_final)
				e=e+1
				band_dist=0	
			}	
		}

		total_final=dist_aux+seq_1[length(seq_1)]
		print(paste(gene_name,total_inicial,total_final))
		save_trans[e,]=c(gene_name,total_inicial,total_final)
		e=e+1

		res_location_trans=rbind(res_location_trans,save_trans)
	}
}

write.table(res_location_trans,paste("/storage/parras/CMV/Salida/",list.files(direct)[pedazos],"_res_location_trans.txt",sep=""),sep="\t",quote = F,row.names = F,col.names = F)


topology2=topology2[1:o,]
topology2=topology2[order(topology2$name),]


sub_MiCMV=data.frame(matrix(ncol=1))
e=1

for(i in 1:length(MiCMV[,1]))
{
	if(sub('.', '', as.character(MiCMV[i,1])) %in% as.character(topology2$name))
	{
		#print(as.character(MiCMV[i,1]))
		sub_MiCMV[e,1]=as.character(MiCMV[i,1])
		sub_MiCMV[e+1,1]=as.character(MiCMV[i+1,1])
		e=e+2
	}
}

#sub_MiCMV=sub_MiCMV[1:e,]

write.table(sub_MiCMV,paste("/storage/parras/CMV/Salida/",list.files(direct)[pedazos],"_seqs_transm.txt",sep=""),quote = F,row.names = F,col.names = F)


#topology2[,1]=as.character(topology2[,1])

for(top in 1:length(topology2$name))
{
	if(!isEmpty(grep("gene=",as.character(topology2[top,1]))))
	{
		topology2[top,1]=strsplit(strsplit(as.character(topology2[top,1]),"gene=")[[1]][2],"\\]")[[1]][1]
	}else
	{
		topology2[top,1]=paste(strsplit(strsplit(as.character(topology2[top,1]),"protein=")[[1]][2],"\\]")[[1]][1],strsplit(strsplit(as.character(topology2[top,1]),"protein_id=")[[1]][2],"\\]")[[1]][1],sep=";")
	}
}




#topology2$name=sapply(strsplit(as.character(topology2$name), " "), "[[", 2)



pdf(paste("/storage/parras/CMV/Salida/",list.files(direct)[pedazos],"_transm2.pdf",sep=""),width=15,height=50)

print(plot_topology(topology2))

dev.off()


####Diamond
#system(paste('/storage/parras/CMV/EggNOG/A_ver/eggnog-mapper-master/emapper.py -i ',paste("/storage/parras/CMV/Salida/",list.files(direct)[pedazos],"_seqs_transm.txt",sep=""),' -o ',paste("/storage/parras/CMV/Salida/",list.files(direct)[pedazos],"_resultadosEggNOG",sep=""),' -m diamond --data_dir /storage/parras/CMV/EggNOG/A_ver/eggnog-mapper-master/data/',sep=""))

####HMMER
system(paste('/storage/parras/CMV/EggNOG/A_ver/eggnog-mapper-master/emapper.py -i ',paste("/storage/parras/CMV/Salida/",list.files(direct)[pedazos],"_seqs_transm.txt",sep=""),' -o ',paste("/storage/parras/CMV/Salida/",list.files(direct)[pedazos],"_resultadosEggNOG",sep=""),' --cpu 30 -m hmmer --data_dir /storage/parras/CMV/EggNOG/A_ver/eggnog-mapper-master/data/ -d Virus',sep=""))


################################################################# Plot

multifa_original=read.delim(paste(direct,list.files(direct)[pedazos],sep="/"),header = F)

multifa_original=as.data.frame(multifa_original[ seq(1, nrow(multifa_original), by = 2),1])

multifa_original$info=sapply(strsplit(sapply(strsplit(as.character(multifa_original[,1]), "\\location="), "[[", 2), "\\]"), "[[", 1)

multifa_original$info=str_replace_all(multifa_original$info,"<","")
multifa_original$info=str_replace_all(multifa_original$info,">","")

for(i in 1:length(multifa_original[,1]))
{
  if(!isEmpty(grep("complement",multifa_original[i,1])))
  {
    if(!isEmpty(grep("join",multifa_original[i,1])))
    { 
      multifa_original[i,3]=strsplit(strsplit(as.character(multifa_original[i,2]), "\\.")[[1]][1],"\\(")[[1]][3]
      multifa_original[i,4]=strsplit(strsplit(as.character(multifa_original[i,2]), "\\.")[[1]][length(strsplit(as.character(multifa_original[i,2]), "\\.")[[1]])],"\\)")[[1]][1]
      multifa_original[i,5]="-"
    }else
    {
      multifa_original[i,3]=strsplit(strsplit(as.character(multifa_original[i,2]), "\\.")[[1]][1],"\\(")[[1]][2]
      multifa_original[i,4]=strsplit(strsplit(as.character(multifa_original[i,2]), "\\.")[[1]][3],"\\)")[[1]][1]
      multifa_original[i,5]="-"
    }
  }else
  {
    if(!isEmpty(grep("join",multifa_original[i,1])))
    { 
      multifa_original[i,3]=strsplit(strsplit(as.character(multifa_original[i,2]), "\\.")[[1]][1],"\\(")[[1]][2]
      multifa_original[i,4]=strsplit(strsplit(as.character(multifa_original[i,2]), "\\.")[[1]][length(strsplit(as.character(multifa_original[i,2]), "\\.")[[1]])],"\\)")[[1]][1]
      multifa_original[i,5]="+"
    }else
    {
      multifa_original[i,3]=strsplit(as.character(multifa_original[i,2]), "\\.")[[1]][1]
      multifa_original[i,4]=strsplit(as.character(multifa_original[i,2]), "\\.")[[1]][3]
      multifa_original[i,5]="+"
    }
  }
}



		if(length(grep("gene=",as.character(multifa_original[,1])))==length(grep("protein=",as.character(multifa_original[,1]))))
		{
			multifa_original$info=sapply(strsplit(sapply(strsplit(as.character(multifa_original[,1]), "\\gene="), "[[", 2), "\\]"), "[[", 1)
		}else
		{
			multifa_original$info=paste(sapply(strsplit(sapply(strsplit(as.character(multifa_original[,1]), "\\protein="), "[[", 2), "\\]"), "[[", 1),sapply(strsplit(sapply(strsplit(as.character(multifa_original[,1]), "\\protein_id="), "[[", 2), "\\]"), "[[", 1),sep=";")

		}




multifa=read.delim(paste("/storage/parras/CMV/Salida/",list.files(direct)[pedazos],"_seqs_transm.txt",sep=""),header = F)


genome_name=strsplit(strsplit(as.character(multifa[1,1]),"\\|")[[1]][2],"\\_")[[1]][1]

multifa=as.data.frame(multifa[ seq(1, nrow(multifa), by = 2),1])

multifa$info=sapply(strsplit(sapply(strsplit(as.character(multifa[,1]), "\\location="), "[[", 2), "\\]"), "[[", 1)

multifa$info=str_replace_all(multifa$info,"<","")
multifa$info=str_replace_all(multifa$info,">","")

for(i in 1:length(multifa[,1]))
{
  if(!isEmpty(grep("complement",multifa[i,1])))
  {
    if(!isEmpty(grep("join",multifa[i,1])))
    { 
      multifa[i,3]=strsplit(strsplit(as.character(multifa[i,2]), "\\.")[[1]][1],"\\(")[[1]][3]
      multifa[i,4]=strsplit(strsplit(as.character(multifa[i,2]), "\\.")[[1]][length(strsplit(as.character(multifa[i,2]), "\\.")[[1]])],"\\)")[[1]][1]
      multifa[i,5]="-"
    }else
    {
      multifa[i,3]=strsplit(strsplit(as.character(multifa[i,2]), "\\.")[[1]][1],"\\(")[[1]][2]
      multifa[i,4]=strsplit(strsplit(as.character(multifa[i,2]), "\\.")[[1]][3],"\\)")[[1]][1]
      multifa[i,5]="-"
    }
  }else
  {
    if(!isEmpty(grep("join",multifa[i,1])))
    { 
      multifa[i,3]=strsplit(strsplit(as.character(multifa[i,2]), "\\.")[[1]][1],"\\(")[[1]][2]
      multifa[i,4]=strsplit(strsplit(as.character(multifa[i,2]), "\\.")[[1]][length(strsplit(as.character(multifa[i,2]), "\\.")[[1]])],"\\)")[[1]][1]
      multifa[i,5]="+"
    }else
    {
      multifa[i,3]=strsplit(as.character(multifa[i,2]), "\\.")[[1]][1]
      multifa[i,4]=strsplit(as.character(multifa[i,2]), "\\.")[[1]][3]
      multifa[i,5]="+"
    }
  }
}



		if(length(grep("gene=",as.character(multifa[,1])))==length(grep("protein=",as.character(multifa[,1]))))
		{
			multifa$info=sapply(strsplit(sapply(strsplit(as.character(multifa[,1]), "\\gene="), "[[", 2), "\\]"), "[[", 1)
		}else
		{
			multifa$info=paste(sapply(strsplit(sapply(strsplit(as.character(multifa[,1]), "\\protein="), "[[", 2), "\\]"), "[[", 1),sapply(strsplit(sapply(strsplit(as.character(multifa[,1]), "\\protein_id="), "[[", 2), "\\]"), "[[", 1),sep=";")
		}




multifasta_backup=multifa

anno_mcv=fread(paste("/storage/parras/CMV/Salida/",list.files(direct)[pedazos],"_resultadosEggNOG.emapper.annotations",sep=""))
anno_mcv=as.data.frame(anno_mcv)
anno_mcv[,1]=sapply(strsplit(anno_mcv[,1], "\\_"), "[[", 4)
multifasta_backup[,1]=sapply(strsplit(sapply(strsplit(as.character(multifasta_backup[,1]), "\\_"), "[[", 4), " "), "[[", 1)


for(oju in 1:length(multifasta_backup[,1]))
{
  if(!isEmpty(anno_mcv[multifasta_backup[oju,1]==anno_mcv[,1],8]))
  {
    multifasta_backup[oju,6]=anno_mcv[multifasta_backup[oju,1]==anno_mcv[,1],8]
  }else
  {
    multifasta_backup[oju,6]="-"
  }
}

#anno_mcv$gene_name=multifasta_backup[multifasta_backup[,1] %in% anno_mcv[,1],2]
#anno_mcv$start=multifasta_backup[multifasta_backup[,1] %in% anno_mcv[,1],3]
#anno_mcv$end=multifasta_backup[multifasta_backup[,1] %in% anno_mcv[,1],4]
#anno_mcv$stand=multifasta_backup[multifasta_backup[,1] %in% anno_mcv[,1],5]

#tab_model_CMV=as.data.frame(cbind(genome_name,anno_mcv$start,anno_mcv$end,(as.numeric(anno_mcv$end)-as.numeric(anno_mcv$start)),anno_mcv$stand,paste(anno_mcv$gene_name," (",anno_mcv$Description,")",sep = "")))
tab_model_CMV=as.data.frame(cbind(genome_name,multifasta_backup$V3,multifasta_backup$V4,(as.numeric(multifasta_backup$V4)-as.numeric(multifasta_backup$V3)),multifasta_backup$V5,paste(multifasta_backup$info," (",multifasta_backup$V6,")",sep = "")))

colnames(tab_model_CMV)=c("chromosome","start","end","width","strand","symbol")
tab_model_CMV$start=as.numeric(as.character(tab_model_CMV$start))
tab_model_CMV$end=as.numeric(as.character(tab_model_CMV$end))
tab_model_CMV$width=as.numeric(as.character(tab_model_CMV$width))
tab_model_CMV[,1]=as.character(tab_model_CMV[,1])

multifa_original[,1]=genome_name

colnames(multifa_original)=c("chr","gene","start","end","strand")

anno_location=as.data.frame(fread(paste("/storage/parras/CMV/Salida/",list.files(direct)[pedazos],"_res_location_trans.txt",sep=""),header=F,sep = "\t"))

tab_model_CMV2=as.data.frame(cbind(genome_name,anno_location$V2,anno_location$V3,(as.numeric(anno_location$V3)-as.numeric(anno_location$V2)),
                                   anno_location$V1))
colnames(tab_model_CMV2)=c("chromosome","start","end","width","symbol")
tab_model_CMV2$start=as.numeric(as.character(tab_model_CMV2$start))
tab_model_CMV2$end=as.numeric(as.character(tab_model_CMV2$end))
tab_model_CMV2$width=as.numeric(as.character(tab_model_CMV2$width))
tab_model_CMV2[,1]=as.character(tab_model_CMV2[,1])

genome_GRange=makeGRangesFromDataFrame(multifa_original)
atrack <- AnnotationTrack(genome_GRange, name = genome_name)
gtrack <- GenomeAxisTrack()
grtrack <- GeneRegionTrack(tab_model_CMV, name = "Transm. ORF",transcriptAnnotation = "symbol",background.panel = "azure",
                           background.title = "darkblue")
funtrack<- GeneRegionTrack(tab_model_CMV2, name = "Transm. region",transcriptAnnotation = "symbol",background.panel = "#fbf3f4",
                           background.title = "#9884b9")
#plotTracks(list(gtrack, atrack, grtrack,funtrack))


#library(BSgenome.Vcholerae.NCBI.N16961)
#options(ucscChromosomeNames=FALSE)
#strack <- SequenceTrack(Vcholerae,chromosome = "FJ527563.1",name = "FJ527563.1")
#plotTracks(list(gtrack, atrack, grtrack,strack))


pdf(paste("/storage/parras/CMV/Salida/",list.files(direct)[pedazos],"_genome_plot.pdf",sep=""),width = 100)
plotTracks(list(gtrack, atrack, grtrack,funtrack))
dev.off()



databaCMV=read.delim("ViralRefSeq_DB_OUT.txt",header = F)
databaCMV[,1]=as.character(databaCMV[,1])
databaCMV[ seq(1, nrow(databaCMV), by = 2),2]=sapply(strsplit(sapply(strsplit(databaCMV[ seq(1, nrow(databaCMV), by = 2),], " "), "[[", 1), ">"), "[[", 2)

#MiCMV=read.delim("beta5_prots_NCBI.txt_OUT.txt",header = F)

system("blastp -query ",paste(direct,list.files(direct)[pedazos],sep="/")," -db ViralRefDB -out Salida_last -outfmt 6")


BlastSalida=read.delim("Salida_last",header = F)
sub_CMV=BlastSalida[BlastSalida[,11]<0.001,]


pdf(paste("/storage/parras/CMV/Salida/",list.files(direct)[pedazos],"_salida_arbol.pdf",sep=""))
for(i in 1:length(MiCMV[,1]))
{
  if(!isEmpty(grep(">",MiCMV[i,1])))
  {
    print(strsplit(strsplit(as.character(MiCMV[i,1])," ")[[1]][1],">")[[1]][2])
    
    acc_interest=as.character(sub_CMV[sub_CMV[,1]==strsplit(strsplit(as.character(MiCMV[i,1])," ")[[1]][1],">")[[1]][2],][,2])
  
    tab_res_blast=databaCMV[sort(c(which(databaCMV[,2] %in% acc_interest),which(databaCMV[,2] %in% acc_interest)+1)),]
    
    seq_mias=as.data.frame(as.character(MiCMV[c(i,i+1),1]))
    seq_mias[,1]=as.character(seq_mias[,1])
    colnames(seq_mias)="V1"
    colnames(tab_res_blast)="V1"
    
    
    lol=c(seq_mias[,1],tab_res_blast[,1])
    
    write.table(as.data.frame(lol),paste("/storage/parras/CMV/Salida/todas_prots_arboles/CMV_",strsplit(strsplit(strsplit(as.character(MiCMV[i,1])," ")[[1]][1],">")[[1]][2],"\\|")[[1]][2],".txt",sep = ""),quote = F,row.names = F,col.names = F)  

    if(dim(data.frame(lol))[1]>4)
    {
      mySequences <- readAAStringSet(paste("/storage/parras/CMV/Salida/todas_prots_arboles/CMV_",strsplit(strsplit(strsplit(as.character(MiCMV[i,1])," ")[[1]][1],">")[[1]][2],"\\|")[[1]][2],".txt",sep = ""))
      myFirstAlignment <- msa(mySequences)
      hemoAln2 <- msaConvert(myFirstAlignment, type="seqinr::alignment")
      d <- dist.alignment(hemoAln2, "identity")
      hemoTree <- njs(d)
      plot(hemoTree, main=paste("Phylogenetic Tree of",strsplit(strsplit(strsplit(as.character(MiCMV[i,1])," ")[[1]][1],">")[[1]][2],"\\|")[[1]][2]))
      print(ggtree(hemoTree,layout = 'circular') + geom_tiplab(size=2.5) + ggtitle(paste("Phylogenetic Tree of",strsplit(strsplit(strsplit(as.character(MiCMV[i,1])," ")[[1]][1],">")[[1]][2],"\\|")[[1]][2])))
    }
  }
  
}
dev.off()




}

#################################################################################################


####HMMER TODO
system('cat /storage/parras/CMV/Salida/todas_prots_arboles/* > /storage/parras/CMV/Salida/todas_prots_arboles/todas_prots_analizar.fasta')
system(paste('/storage/parras/CMV/EggNOG/A_ver/eggnog-mapper-master/emapper.py -i /storage/parras/CMV/Salida/todas_prots_arboles/todas_prots_analizar.fasta -o /storage/parras/CMV/Salida/todas_prots_arboles/todas_prots_analizar_SALIDA --cpu 30 -m hmmer --data_dir /storage/parras/CMV/EggNOG/A_ver/eggnog-mapper-master/data/ -d Virus',sep=""))




