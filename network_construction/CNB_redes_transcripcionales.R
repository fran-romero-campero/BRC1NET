############################################
## Construcci?n y An?lisis de             ##
## Redes Transcripcionales.               ## 
##                                        ##
## Autor: Francisco J. Romero-Campero     ##
## Email: fran@us.es                      ##
## Fecha: Febrero 2018                    ##
############################################

## Las redes transcripcionales son redes dirigidas donde los nodos representan genes 
## y se traza una arista de un nodo g1 al node g2 cuando g1 codifica por un factor de 
## transcripci?n se une al promoter del g2. T?picamente, est?s redes se construyen a 
## partir de datos de ChIP-seq y suponen un refinameinto de la redes de co-expresi?n 
## g?nica. 

## Este script recibe como entrada un fichero donde a su vez se guardan los 
## identificadores de los genes para los cuales a partir de datos de ChIP-seq se han 
## determinado sus dianas potenciales.

tfs.data <- read.table(file = "targets.txt",as.is=TRUE)
tfs.data
tfs <- tfs.data[,1]
tfs
length(tfs)

## Para obtener el conjunto total de genes de la red leemos todos los ficheros 
## correspondientes a la dianas potenciales de los factores de transcripci?n y los
## vamos almacenando en un vector. 
genes.in.network <- c()

for(i in 1:length(tfs))
{
  current.tf.file <- paste(tfs[i],".txt",sep="")
  current.tf.data <- read.table(file=current.tf.file,as.is=TRUE)
  current.tf.targets <- current.tf.data[,1]
  print(length(current.tf.targets))
  genes.in.network <- c(genes.in.network,current.tf.targets)
}


tfs[3] # AT4G34000 ABF3 11562
tfs[7] # AT2G46270 GBF3 17542
tfs[8] # AT4G36740 HB40 14554
genes.in.network <- unique(genes.in.network)
genes.in.network
length(genes.in.network)

## Generamos incialmente una matriz de adyacencia con todos los valores 0 y vamos
## rellen?ndola con los valores 1 son las dianas de los correspondientes factores
## de transcripci?n. 
adj.matrix <- matrix(0,nrow=length(genes.in.network),ncol=length(genes.in.network))

rownames(adj.matrix) <- genes.in.network
colnames(adj.matrix) <- genes.in.network

adj.matrix[1:4,1:4]
adj.matrix[1:4,1:4]

for(i in 1:length(tfs))
{
  current.tf <- tfs[i]
  current.tf.file <- paste(current.tf, ".txt",sep="")
  current.tf.data <- read.table(file=current.tf.file,as.is=TRUE)
  current.tf.targets <- current.tf.data[,1]
  adj.matrix[current.tf,current.tf.targets] <- 1
}

## Finalmente, construimos la red transcripcional a partir de su matriz de adyacencia
## y la guardamos en formato gml. 
library(igraph)
install.packages("igraph")
gene.transcriptional.network <- graph.adjacency(adj.matrix, mode="directed")
write.graph(gene.transcriptional.network,file="transcriptional_gene_network.gml",format="gml")

gene.transcriptional.network

## Esta red es excesivamente grande y por motivos de recursos computacionales y tiempo
## vamos a trabajar s?lo con la red inducida que contiene los factores de transcripci?n
## y las aristas entre ellos.
tfs.network <- induced.subgraph(gene.transcriptional.network,tfs)
plot.igraph(x = tfs.network,vertex.size=8,edge.arrow.size=0.5,vertex.label="",vertex.color="blue")
write.graph(tfs.network,file="red_transcripcional_tf.gml",format="gml")

## Uno de los primeros an?lisis que se realizan sobre redes transcripcionales es la 
## b?queda de subgrafos o patrones no aleatorios. Estos subgrafos se denominan 
## motivos de red.
## Formalmente, un motivo de red es un subgrafo que aparece un n?mero significativamente
## de veces mayor en la red de inter?s que en redes aleatorias que cumplen las 
## mismas propiedades.

## Primero, realizamos un an?lisis r?pido para comprobar que esta red inducida
## NO es libre de escala. 

# Calculo del grado de los nodos
network.degrees <- degree(tfs.network)
hist(network.degrees,col="blue",xlab="Node degree", ylab="Probability",main="Degree distribution")

# Calculo de la frecuencia absoluta del grado de los nodos
degree.frequencies <- table(network.degrees)
# Eliminamos nodos con grado 0 para poder aplicar log10

# Transformaci?n logar?tmica
log10.degrees.frequencies <- log10(degree.frequencies)
log10.node.degrees <- log10(as.numeric(names(degree.frequencies)))

# Regresi?n lineal
lm.r <- lm(log10.degrees.frequencies ~ log10.node.degrees)
summary(lm.r)

## Por lo tanto, las redes aleatorias que cumplen las mismas propiedas que nuestra
## red de inter?s son las generadas seg?n el modelo de Erdos-Renyi. La funci?n
## erdos.renyi.game genera redes aleatorias que siguen una distribuci?n de Poisson
## con un n?mero de nodos y aristas dado. 
tfs.network

random.graph <- erdos.renyi.game(n=21, p.or.m=115, type="gnm", directed=TRUE, loops=TRUE)
plot.igraph(x = random.graph,vertex.size=8,edge.arrow.size=0.5,vertex.label="",vertex.color="blue")

## Empezamos viendo si es un motivo de red la autorregulaci?n.
## Para ver el n?mero de genes autorregulados en una red basta con sumar los
## elementos de la diagonal principal. 

## N?mero de genes autorregulados en la red de factores de transcripci?n. 
tfs.adjacency <- as.matrix(get.adjacency(tfs.network))
tfs.adjacency
diag(tfs.adjacency)
autorregulation.in.tfs <- sum(diag(tfs.adjacency))
autorregulation.in.tfs

## N?mero de genes autorregulados en la red 
autorregulation.in.random <- sum(diag(as.matrix(get.adjacency(random.graph))))
autorregulation.in.random

## Para comprobar la signficancia de la autorregulaci?n generamos 1000 redes
## aleatorias y calculamos el n?mero de autorregulaciones en cada una almacen?ndola
## en un vector.
autorregulation.random.graphs <- vector(length=1000, mode="numeric")

for (i in 1:1000)
{
  random.graph <- erdos.renyi.game(n=21, p.or.m=115, type="gnm", directed=TRUE, loops=TRUE)
  autorregulation.random.graphs[i] <- sum(diag(as.matrix(get.adjacency(random.graph))))
}

mean(autorregulation.random.graphs)
sd(autorregulation.random.graphs)

sum(autorregulation.random.graphs > autorregulation.in.tfs)/10000


## Motivos de red de tres nodos

## La funci?n graph.motifs recibe como entrada una red y un tama?o de subgrafo k
## (en la actualidad s?lo puede recibir tama?os 3 o 4) y devuelve el n?mero de
## veces que se encuentra cada subgrafo con k nodos en la red. 

occurrency.subgraph.three.nodes <- graph.motifs(tfs.network, size=3)
occurrency.subgraph.three.nodes
length(occurrency.subgraph.three.nodes)

## La funci?n graph.isocreate genera todos los grafos posibles de un tama?o dado
## empezando a enumerarlos por 0. 
plot.igraph(graph.isocreate(size=3, number=0))
plot.igraph(graph.isocreate(size=3, number=1))
plot.igraph(graph.isocreate(size=3, number=2))
plot.igraph(graph.isocreate(size=3, number=3))
plot.igraph(graph.isocreate(size=3, number=4))
plot.igraph(graph.isocreate(size=3, number=5))
plot.igraph(graph.isocreate(size=3, number=6))
plot.igraph(graph.isocreate(size=3, number=7))
plot.igraph(graph.isocreate(size=3, number=8))
plot.igraph(graph.isocreate(size=3, number=9))
plot.igraph(graph.isocreate(size=3, number=10))
plot.igraph(graph.isocreate(size=3, number=11))
plot.igraph(graph.isocreate(size=3, number=12))
plot.igraph(graph.isocreate(size=3, number=13))

## Para ver si cada uno de los 16 posibles subgrafos es o no motivo de red
## realizamos le procedimiento anterior basado en la generaci?n de redes aleatorias.
## En este caso el acumulador es matricial. 
motifs.3.random.graph <- matrix(0,nrow=1000, ncol=16)
motifs.3.random.graph[1:3,]

for (i in 1:1000)
{
  random.graph <- erdos.renyi.game(n=21, p.or.m=115, type="gnm", directed=TRUE, loops=TRUE)
  motifs.3.random.graph[i,] <- graph.motifs(random.graph, size=3)
}

motifs.3.random.graph[1:3,]

## subgrafo 3
plot(graph.isocreate(size=3, number=2))
occurrency.subgraph.three.nodes[3]

mean(motifs.3.random.graph[,3])
sd(motifs.3.random.graph[,3])
sum(motifs.3.random.graph[,3] >71)/1000
##No significativo

## subgrafo 14
plot(graph.isocreate(size=3, number=13))
occurrency.subgraph.three.nodes[14]

mean(motifs.3.random.graph[,14])
sd(motifs.3.random.graph[,14])
sum(motifs.3.random.graph[,14] > 30)/1000
##No significativo

feedback.loop.with.output <- graph.isocreate(size=3, number=13)

## Para encontrar todas las instancias de este subgrafo en nuestra red usamos 
## la siguiente instrucci?n.
maps.feedback.loop.with.output <- graph.subisomorphic.lad(feedback.loop.with.output, tfs.network, all.maps=TRUE)[["maps"]]
maps.feedback.loop.with.output

