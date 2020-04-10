gbf3.raw <- read.table(file="../../../../BRC1 network paper/from Fran Romero/network/GBF3_ABA_targetGenes_5column.txt",as.is=T,header=F)
head(gbf3.raw)

gbf3.targets <- unique(gbf3.raw$V5)
length(gbf3.targets)

network.data <- read.table(file="brc1_network.tsv",header = TRUE,as.is=TRUE,sep="\t",quote = "",comment.char = "%")
rownames(network.data) <- network.data$names
gbf3.network.targets <- intersect(network.data$names,gbf3.targets)
length(gbf3.network.targets)

network.data[gbf3.network.targets,"GBF3"] <- 1

write.table(x = network.data,file = "brc1_network.tsv",quote = F,sep = "\t",row.names = F)
