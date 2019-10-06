## carrega os pacotes
library(pd.hg.u133.plus.2)
library(HDF5Array)
library(affyio)
library(limma)
library(R.utils)
library(dplyr)
library(Rcpp)
library(doParallel)
library(foreach)
library(future)

#http://dirk.eddelbuettel.com/code/rcpp/html/Rmath_8h_source.html

## cria um objeto que armazena os nomes dos dados que serão carregados
files.name <- as.character(list.files(all.files = TRUE, pattern = "[gz]", recursive = FALSE, include.dirs = FALSE, full.names = FALSE, no.. = TRUE))

## testa se todos os CELs são do mesmo pacote
input.kind <- lapply(files.name, read.celfile.header)
same.kind <- length(unique(sapply(input.kind, "[[", "cdfName"))) == 1
if (!same.kind)
  stop("All files need to be from the same package.")

## extrai as dimensões do banco de dados
nrows <- prod(input.kind[[1]][["CEL dimensions"]])
dim <- as.integer(c(nrows, length(files.name)))

## cria conexão com banco de dados do pacote e extrai as informações dos perfect match
pkganno <- "pd.hg.u133.plus.2"  ### Why is not input.kind[[1]][["cdfName"]] ?
library(pkganno, character.only = TRUE)
conn <- db(get(pkganno))
pminfo <- dbGetQuery(conn, "SELECT * FROM pmfeature")

## cria um arquivo HDF5 com um grupo e um espaço para armazenar os dados
h5createFile("myhdf5file.h5")
h5createGroup("myhdf5file.h5", "analysis")
h5createDataset(file = "myhdf5file.h5", dataset = "analysis/dataset", dims = dim, maxdims = dim, storage.mode = "double", chunk = dim, level = 5, showWarnings = FALSE)

## cria objetos da classe arrayRealizationSink para processamentos em blocos, assim como o "bloco" que será utilizado
real.sink <- RealizationSink(dim = dim, dimnames = NULL, type = "double")
block <- colGrid(real.sink, ncol = 100)

## salva os dados em HDF5
begin <- 1
end <- ncol(block[[1L]])
for (i in seq_along(block)) {
  aux <- read.celfiles(files.name[begin:end])
  aux <- exprs(aux)
  h5write(aux, "myhdf5file.h5", "analysis/dataset", index = list(NULL, begin:end))
  begin <- begin + ncol(block[[i]])
  ifelse(test = i == length(block)-1, yes = end <- end + ncol(block[[length(block)]]), no = end <- end + ncol(block[[i]]))
}

# Fecha o arquivo HDF5 e abre novamente
h5closeAll()
h5f <- H5Fopen("myhdf5file.h5")
h5g <- H5Gopen(h5f, "analysis")

## Reserva um espaço de memória para salvar os PM
h5createDataset(file = h5g, dataset = "PM", dims = c(length(pminfo$fid), length(files.name)), maxdims = c(length(pminfo$fid), length(files.name)), storage.mode = "double", chunk = c(length(pminfo$fid), length(files.name)), level = 5, showWarnings = FALSE)

## Salva os PM em HDF5
begin <- 1
end <- ncol(block[[1L]])
for (i in seq_along(block)) {
  aux <- h5read(h5g, "dataset", index = list(NULL, begin:end))
  aux <- apply(aux, 2, "[", pminfo$fid)
  h5write(aux, h5g, "PM", index = list(NULL, begin:end))
  begin <- begin + ncol(block[[i]])
  ifelse(test = i == length(block)-1, yes = end <- end + ncol(block[[length(block)]]), no = end <- end + ncol(block[[i]]))
}

## Reserva um espaço de memória para salvar os PM com background corrigido
h5createDataset(file = h5g, dataset = "background(PM)", dims = c(length(pminfo$fid), length(files.name)), maxdims = c(length(pminfo$fid), length(files.name)), storage.mode = "double", chunk = c(length(pminfo$fid), length(files.name)), level = 5, showWarnings = FALSE)

## cria objetos da classe arrayRealizationSink para processamentos em blocos, assim como o "bloco" que será utilizado
real.sink <- RealizationSink(dim = c(length(pminfo$fid), length(files.name)), dimnames = NULL, type = "double")
block <- colGrid(real.sink, ncol = 100) 

## Salva os PM com background corrigido em HDF5
begin <- 1
end <- ncol(block[[1L]])
for (i in seq_along(block)) {
  aux <- h5read(h5g, "PM", index = list(NULL, begin:end))
  aux <- apply(aux, 2, backgroundCorrect.matrix, normexp.method = "rma", verbose = TRUE)
  h5write(aux, h5g, "background(PM)", index = list(NULL, begin:end))
  begin <- begin + ncol(block[[i]])
  ifelse(test = i == length(block)-1, yes = end <- end + ncol(block[[length(block)]]), no = end <- end + ncol(block[[i]]))
} 

# Define que os dados do disco serão reconhecidos como HDF5
setRealizationBackend("HDF5Array")

# Reconhece os dados dos PM com background corrigido como HDF5 no environment
#data <- realize(h5g$`background(PM)`)

## Obtém os índices dos PM com background corrigido, assim como os armazena no objeto index
begin <- 1
end <- ncol(block[[1L]])
index <- as(matrix(nrow = nrow(data), ncol = ncol(data)), "HDF5Array")
for (i in seq_along(block)) {
  aux <- h5read(h5g, "dataset", index = list(NULL, begin:end))
  #aux <- read_block(data, block[[i]])
  aux <- apply(aux, 2, order)
  index <- as(write_block(index, block[[i]], aux), "HDF5Array")   
  begin <- begin + ncol(block[[i]])
  ifelse(test = i == length(block)-1, yes = end <- end + ncol(block[[length(block)]]), no = end <- end + ncol(block[[i]]))
}

## Reserva um espaço de memória para salvar os PM com background corrigido e normalizados
h5createDataset(file = h5g, dataset = "normalize(PM)", dims = c(length(pminfo$fid), length(files.name)), maxdims = c(length(pminfo$fid), length(files.name)), storage.mode = "double", chunk = c(length(pminfo$fid), length(files.name)), level = 5, showWarnings = FALSE)  

## cria objetos da classe arrayRealizationSink para processamentos em blocos, assim como os "blocos" que serão utilizados
real.sink <- RealizationSink(dim = c(length(pminfo$fsetid), length(files.name)), dimnames = NULL, type = "double")
block.col <- colGrid(real.sink, ncol = 100)
block.row <- rowGrid(real.sink, nrow = 100000)

## Ordena do menor para o maior os PM com background corrigido
begin <- 1
end <- ncol(block.col[[1L]])
for (i in seq_along(block.col)) {
  aux <- h5read(h5g, "background(PM)", index = list(NULL, begin:end))
  aux <- apply(aux, 2, sort)
  h5write(aux, h5g, "normalize(PM)", index = list(NULL, begin:end))
  begin <- begin + ncol(block.col[[i]])
  ifelse(test = i == length(block.col)-1, yes = end <- end + ncol(block.col[[length(block.col)]]), no = end <- end + ncol(block.col[[i]]))
}

## Obtém as médias das linhas da matrix dos PM com background corrigido e ordenados e as insere em todas as colunas, assim como armazena no espaço de memória selecionado acima
begin <- 1
end <- nrow(block.row[[1L]])
for (i in seq_along(block.row)) {
  aux <- h5read(h5g, "normalize(PM)", index = list(begin:end, NULL))
  aux <- matrix(data = rep(rowMeans(aux), each = length(files.name)), nrow = nrow(block.row[[i]]), ncol = length(files.name), byrow = TRUE )
  h5write(aux, h5g, "normalize(PM)", index = list(begin:end, NULL))
  begin <- begin + nrow(block.row[[i]])
  ifelse(test = i == length(block.row)-1, yes = end <- end + nrow(block.row[[length(block.row)]]), no = end <- end + nrow(block.row[[i]]))
}

## Reordena os PM em sua ordem original e armazena em .h5
begin <- 1
end <- ncol(block.col[[1L]])
for (i in seq_along(block.col)) {
  aux <- h5read(h5g, "normalize(PM)", index = list(NULL, begin:end))
  aux1 <- index[,begin:end]
  for(j in seq(ncol(aux))){
    aux[,j] <- aux[,j][aux1[,j]]
  }
  h5write(aux, h5g, "normalize(PM)", index = list(NULL, begin:end))
  begin <- begin + ncol(block.col[[i]])
  ifelse(test = i == length(block.col)-1, yes = end <- end + ncol(block.col[[length(block.col)]]), no = end <- end + ncol(block.col[[i]]))
}

# Reconhece os dados como HDF5 no environment
#data <- realize(h5g$dataset)
#normalized.pm <- realize(h5g$`normalize(PM)`)

## cria objetos da classe arrayRealizationSink para processamentos em blocos, assim como o "bloco" que será utilizado
real.sink <- RealizationSink(dim = dim, dimnames = NULL, type = "double")
real.sink.pm <- RealizationSink(dim = c(length(pminfo$fsetid), length(files.name)), dimnames = NULL, type = "double")
block <- colGrid(real.sink, ncol = 100)
block.pm <- colGrid(real.sink.pm, ncol = 100)   

## cria o objeto no qual serão salvos os raw data com os PM normalizados
raw.with.pm <- as(matrix(nrow = nrow(real.sink), ncol = ncol(real.sink)), "HDF5Array")

## cria uma função em C++ que insere valores em uma matriz dado um vetor de índices
cppFunction('
            NumericMatrix insert_normalized_pm(NumericMatrix data, NumericMatrix normalized_pm, NumericVector index){
              
              int ncols = normalized_pm.ncol(), nrows = normalized_pm.nrow();

              for(int i = 0; i < ncols; i++){
                for(int j = 0; j < nrows; j++){
                  
                  data(index[j],i) = normalized_pm(j,i);

                }
              }
            return data;
            }
            ')

## armazena os dados com PM com background corrigido e normalizado no banco de dados original
begin <- 1
end <- ncol(block.pm[[1L]])
for (i in seq_along(block.pm)) {
  aux1 <- h5read(h5g, "normalize(PM)", index = list(NULL, begin:end))
  aux <- h5read(h5g, "test", index = list(NULL, begin:end))
  #aux1 <- read_block(normalized.pm, block.pm[[i]])
  #aux <- read_block(data, block[[i]])
  aux <- insert_normalized_pm(data = as.matrix(aux), normalized_pm = as.matrix(aux1), index = c(pminfo$fid - 1))
  #for(j in seq(ncol(aux))){
  #  aux[,j] <- insert(aux[,j], pminfo$fid, aux1[,j])[-((pminfo$fid)+c(1:length(pminfo$fid)))]
  #}
  raw.with.pm <- as(write_block(raw.with.pm, block[[i]], aux), "HDF5Array")
  begin <- begin + ncol(block.pm[[i]])
  ifelse(test = i == length(block.pm)-1, yes = end <- end + ncol(block.pm[[length(block.pm)]]), no = end <- end + ncol(block.pm[[i]]))
}

rm(aux,aux1)

## ordena os valores dos indices dos PM em relação aos identificadores dos genes (subtrai 1 para ser usado em um loop em C++)
by.gene.id <- pminfo[order(pminfo$fsetid),-c(2:4)] - 1 

## cria uma função em C++ que retorna o agrupamento os dados de PM por gene
cppFunction('
            NumericMatrix genes_matrix_pm(NumericMatrix data, NumericMatrix genes_pm, NumericVector by_gene_id){

              int nrows = genes_pm.nrow(), ncols = genes_pm.ncol();

              for(int i=0; i<nrows; i++){
                for(int j=0; j<ncols; j++){

                  genes_pm(i,j) = data(by_gene_id[i],j);

                }
              }
              return genes_pm;
            }
            ') 


## armazena a quantidade de valores de expressão medidos para cada indentificador
counting.genes.id <- as.matrix(table(pminfo$fsetid))

## cria o objeto no qual serão salvos os PM dos genes ordenados pelo identificador
ord.genes.pm <- as(matrix(nrow = nrow(real.sink.pm), ncol = ncol(real.sink.pm)), "HDF5Array")

## armazena os dados dos PM dos genes ordenados pelo identificador
begin <- begin.pm <- 1
end <- end.pm <- ncol(block.pm[[1L]])
for (i in seq_along(block.pm)) {
  
  aux <- read_block(raw.with.pm, block[[i]])
  
  ## cria o objeto em que são armazenadas as matrizes com os PM dos genes ordenadas pelo identificador
  genes.pm <- genes_matrix_pm(as.matrix(aux), matrix(nrow = nrow(pminfo), ncol = ncol(aux)), as.matrix(by.gene.id))
  
  ord.genes.pm <- as(write_block(ord.genes.pm, block.pm[[i]], genes.pm), "HDF5Array")
  
  begin.pm <- begin <- begin + ncol(block.pm[[i]])
  ifelse(test = i == length(block.pm)-1, yes = end.pm <- end <- end + ncol(block.pm[[length(block.pm)]]), no = end.pm <- end <- end + ncol(block.pm[[i]]))
}

rm(raw.with.pm, by.gene.id, genes.pm, aux, begin, begin.pm, end, end.pm)

by.gene.id <- pminfo[order(pminfo$fsetid),-c(3:4)]
splited.index <- split(by.gene.id[,1], by.gene.id[,2])

## Reserva um espaço de memória para salvar os PM com background corrigido, normalizados e sumarizados
h5createDataset(file = h5g, dataset = "summarize(PM)", dims = c(length(unique(pminfo$fsetid)), length(files.name)), maxdims = c(length(unique(pminfo$fsetid)), length(files.name)), storage.mode = "double", chunk = c(length(unique(pminfo$fsetid)), length(files.name)), level = 5, showWarnings = FALSE)    

## cria um objeto no qual serão salvos os dados de expressão
exprs.values <- c()
begin <- 1
end <- 0

## aplica a sumarização e salva os valores de expressão no objeto exprs.value
for (i in seq_along(counting.genes.id)) {
  end <- end + counting.genes.id[i]  
  gene.specific.matrix <- genes.pm[begin:end,]
  gene.specific.matrix <- log(colMeans(sweep(gene.specific.matrix, 2, medpolish(gene.specific.matrix)$col)),2)
  exprs.values <- rbind(exprs.values, gene.specific.matrix)
  begin <- begin + counting.genes.id[i]
}


## realiza sumarização em paralelo

# foreach package

no_cores <- detectCores(logical = F)

cl <- makeCluster(no_cores-1)

registerDoParallel(cl)

exprs.values <- foreach(i = seq_along(counting.genes.id), .combine = rbind, .verbose = T, .packages = c("DelayedArray", "stats")) %dopar% {
                  end <- end + counting.genes.id[i]  
                  gene.specific.matrix <- ord.genes.pm[begin:end,]
                  begin <- begin + counting.genes.id[i] 
                  gene.specific.matrix <- log(colMeans(sweep(gene.specific.matrix, 2, medpolish(gene.specific.matrix)$col)),2)
                }
gene.expression <- function(index){
  log(abs(colMeans(sweep(raw.with.pm[index,],2,medpolish(raw.with.pm[index,])$col))),2)
}

clusterExport(cl, "raw.with.pm")

exprs.values <- parSapplyLB(cl, splited.index, gene.expression)

stopCluster(cl)
# doFuture package

registerDoFuture()
plan(multiprocess)


## cria objetos da classe arrayRealizationSink para processamentos em blocos, assim como o "bloco" que será utilizado
real.sink <- RealizationSink(dim = dim(raw.with.pm), dimnames = NULL, type = "double")
real.sink.pm <- RealizationSink(dim = dim(genes.pm), dimnames = NULL, type = "double")
block <- colGrid(real.sink, ncol = 100)
block.pm <- colGrid(real.sink.pm, ncol = 100)

# Armazena a matriz dos PM para um gene específico em .h5
for (i in seq_along(block)) {
  aux1 <- c()
  aux <- read_block(raw.with.pm, block[[i]])
  for (j in seq(ncol(aux))) {
    aux1 <- cbind(aux1, aux[,j][pminfo$fid])
  }
  genes.pm <- write_block(genes.pm, block.pm[[i]], aux1)
}









# Fecha o arquivo .h5
h5closeAll()