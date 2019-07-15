## carrega os pacotes
library(pd.hg.u133.plus.2)
library(HDF5Array)
library(affyio)
library(limma)
library(R.utils)
library(dplyr)


###### FUNÇÃO QUE IMPORTA OS DADOS COMO HDF5

## cria um objeto que armazena os nomes dos dados que serão carregados
setwd("~/Documents/GSE25507")
files.name <- as.character(list.files(all.files = TRUE, pattern = "[gz]", recursive = FALSE, include.dirs = FALSE, full.names = FALSE, no.. = TRUE))

## testa se todos os CELs sao do mesmo tipo
input.kind <- lapply(files.name, read.celfile.header)
same.kind <- length(unique(sapply(input.kind, "[[", "cdfName"))) == 1
if (!same.kind)
  stop("All files need to be from the same package.")

## extrai as dimensões do banco de dados
nrows <- prod(input.kind[[1]][["CEL dimensions"]])
dim <- c(nrows, length(files.name))

pkganno <- "pd.hg.u133.plus.2"  # Why is not input.kind[[1]][["cdfName"]] ?
library(pkganno, character.only = TRUE)
conn <- db(get(pkganno))
## dbListTables(conn) retorna as tabelas possíveis de serem acessadas através da conexão gerada acima
pminfo <- dbGetQuery(conn, "SELECT * FROM pmfeature")
featinfo <- dbGetQuery(conn, 'SELECT * FROM feature')

## cria um arquivo HDF5 com um grupo e um espaço para armazenar os dados
h5createFile("myhdf5file.h5")
h5createGroup("myhdf5file.h5", "analysis")
h5createDataset(file = "myhdf5file.h5", dataset = "analysis/dataset", dims = dim, maxdims = dim, storage.mode = "double", chunk = dim, level = 5, showWarnings = FALSE)

## cria objeto da classe arrayRealizationSink para processamento em bloco,
## assim como o "bloco" que será utilizado
real.sink <- RealizationSink(dim = dim, dimnames = NULL, type = "double")
block <- colGrid(real.sink, ncol = 100)

## determina um limite para a quantidade de dados que seram processadas,
## bem como será o tipo de leitura dos dados (HDF5)
makeCappedVolumeBox(maxvol = 50000000, maxdim = c(500000,100), shape = "first-dim-grows-first")
setRealizationBackend("HDF5Array")

## armazena os dados no arquivo criado acima já com o background corrigido pelo método rma

### ARMAZENAR A MATRIZ COM OS DADOS E PASSAR AS DEMAIS INFORMAÇÕES COMO ATRIBUTO PARA O ARQUIVO HDF5

#importação por blocos de 100 colunas
#num.blocks <- length(files.name) %/% 100
#length.last.block <-length(files.name) %% 100
#i <- 0
#dim(files.name) <- c(1, length(files.name))
#apply(files.name, 2, function(name.celfile){
#  i <- i+1
#  h5write(read.celfiles(name.celfile), "myhdf5file.h5", "analysis/fulldata", index = list(NULL,i))
#}) 

#importação por amostra com background corrigido
for (i in seq_along(files.name)) {
  aux <- read.celfile(files.name[i], intensity.means.only = TRUE)[["INTENSITY"]][["MEAN"]] ## armazena uma amostra
  newaux <- aux[pminfo$fid] ## seleciona apenas PM
  newaux <- backgroundCorrect.matrix(newaux, normexp.method = "rma", verbose = TRUE) ## realiza correção de ruído de fundo
  dim(newaux) <- c(NULL) 
  aux <- insert(aux, pminfo$fid, newaux) ## devolve PM com ruído de fundo corrigido
  aux <- aux[-((pminfo$fid)+c(1:length(pminfo$fid)))] ## extrai PM crus
  h5write(aux, "myhdf5file.h5", "analysis/dataset", index = list(NULL,i))
}

## seleciona outro espaço de memória para armazenar os dados dos índices 
#h5createDataset(file = "myhdf5file.h5", dataset = "analysis/index0", dims = c(length(pminfo$fid), length(files.name)), maxdims = dim, storage.mode = "double", chunk = dim, level = 5, showWarnings = FALSE)


h5closeAll()
h5f <- H5Fopen("myhdf5file.h5")
h5g <- H5Gopen(h5f, "analysis")
h5d <- H5Dopen(h5g, "dataset")
#h5id <- H5Dopen(h5g, "index")

## localiza os dados de expressão na memória
data <- realize(h5g$dataset)

## seleciona apenas PM
pm <- matrix(nrow = length(pminfo$fid), ncol = length(files.name))
for (i in 1:length(files.name)) {
  aux <- data[,i]
  pm[,i] <- aux[pminfo$fid]
}

## obtem os índices assim como os armazena no espaço selecionado acima
for (i in 1:length(files.name)) {
#  aux <- h5read("myhdf5file.h5", "analysis/dataset", index = list(NULL, i))
  aux <- order(pm[,i])
  h5write(aux, "myhdf5file.h5", "analysis/index0", index = list(NULL,i))
}

pm <- apply(pm, 2, sort) ## ordena-os PM do menor para o maior
mu <- apply(pm, 1, mean) ## calcula a média das linhas dos PM

## atribui a cada linha o valor de sua respectiva média
for (i in 1:length(mu)) {
  pm[i, ] <- rep(mu[i], length(files.name))
}

## localiza os dados dos índices e reordena os dados de expressão
id <- realize(h5g$index0)
for(i in 1:length(files.name)){
    aux <- id[ ,i]
    newaux <- pm[,i]
    newaux <- newaux[order(aux)]
    pm[ ,i] <- newaux
}

## armazena os dados com PM com background corrigido e normalizado 
for (i in 1:length(files.name)) {
  aux <- data[,i]
  newaux <- pm[,i]
  aux <- insert(aux, pminfo$fid, newaux) ## devolve PM com background corrigido e normalizado
  aux <- aux[-((pminfo$fid)+c(1:length(pminfo$fid)))] ## extrai PM apenas com background corrigido
  h5write(aux, "myhdf5file.h5", "analysis/test", index = list(NULL,i))  
}

## cria objeto onde serão armazenadas as matrizes de PM para cada gene 
pminfo <- as_tibble(pminfo) 
genes.id <- unique(pminfo$fsetid) ## armazena todos os indentificadores dos genes 
genes.pm <- c()
exprs.values <- c()

## armazena as matrizes com os PM de cada gene para todas as amostras
for (j in 1:length(genes.id)) {
  newaux <- filter(pminfo, pminfo$fsetid == genes.id[j])
  for (i in 1:length(files.name)) {
  aux <- data[,i]
  genes.pm <- cbind(genes.pm, aux[newaux$fid])
  }
  column.effect <- medpolish(genes.pm)$col
  genes.pm <- sweep(genes.pm, 2, column.effect)
  genes.pm <- apply(genes.pm, 2, mean)
  exprs.values <- rbind(exprs.values, genes.pm)
  }

## cria "grades" para processamento em blocos  
original.grid <- blockGrid(data, block.length = 5000000, block.shape = "first-dim-grows-first")
medpolish.grid <- blockGrid(data, block.length = 5000000, block.shape = "first-dim-grows-first")

## aplica o polimento de mediana em blocos 
for(i in seq_along(original.grid)){
  original.block <- read_block(data, original.grid[[i]])
  medpolish.block <- medpolish(original.block)$residuals
  return.block <- original.block - medpolish.block
  data <- write_block(data, original.grid[[i]], return.block)
}

## salva os dados processados 
for (i in 1:length(files.name)) {
  h5write(data[, i], "myhdf5file.h5", "analysis/dataset", index = list(NULL,i))
}

h5closeAll()



#run <- by(pminfo$fid, pminfo$fsetid, function(pm.index){
#  for (i in 1:length(files.name)) {
#    genes.pm <- cbind(genes.pm, data[,i]pm.index    
#  }
#})

gene.pm.id <- as(as(by(pminfo$fid, pminfo$fsetid, rbind, simplify = T), "DataFrame"), "HDF5Array")


#  k <- k + 1    
# apply(data, 2, function(amostra){
#        genes.pm[k] <- cbind(genes.pm[k], amostra[pm.index])
#      })
#    })
