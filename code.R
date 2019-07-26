## carrega os pacotes
library(pd.hg.u133.plus.2)
library(HDF5Array)
library(affyio)
library(limma)
library(R.utils)
library(dplyr)
library(sqldf)


###### FUNÇÃO QUE IMPORTA OS DADOS COMO HDF5

## Vai até o diretório no qual estão os arquivos .CEL
setwd("~/Documents/GSE25507")

## cria um objeto que armazena os nomes dos dados que serão carregados
files.name <- as.character(list.files(all.files = TRUE, pattern = "[gz]", recursive = FALSE, include.dirs = FALSE, full.names = FALSE, no.. = TRUE))

## testa se todos os CELs são do mesmo tipo
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
mminfo <- dbGetQuery(conn, "SELECT * FROM mmfeature")
featinfo <- dbGetQuery(conn, 'SELECT * FROM featureSet')
tableinfo <- dbGetQuery(conn, "SELECT * FROM table_info")


## cria um arquivo HDF5 com um grupo e um espaço para armazenar os dados
h5createFile("myhdf5file.h5")
h5createGroup("myhdf5file.h5", "analysis")
h5createDataset(file = "myhdf5file.h5", dataset = "analysis/dataset", dims = dim, maxdims = dim, storage.mode = "double", chunk = dim, level = 5, showWarnings = FALSE)


## determina um limite para a quantidade de dados que seram processadas,
## bem como será o tipo de leitura dos dados (HDF5)
#makeCappedVolumeBox(maxvol = 500000000, maxdim = c(5000000,100), shape = "first-dim-grows-first")
setRealizationBackend("HDF5Array")

## armazena os dados no arquivo criado acima já com o background corrigido pelo método rma

### ARMAZENAR A MATRIZ COM OS DADOS E DEMAIS INFORMAÇÕES PARA O ARQUIVO HDF5 EM ESPAÇOS DE MEMÓRIA SEPARADOS

#importação por blocos de 100 colunas


#importação por amostra sem background corrigido
for (i in seq_along(files.name)) {
  aux <- read.celfile(files.name[i], intensity.means.only = T)[["INTENSITY"]][["MEAN"]]
  h5write(aux, "myhdf5file.h5", "analysis/index1", index = list(NULL,i))
}


#importação por amostra com background corrigido
#for (i in seq_along(files.name)) {
#  aux <- read.celfile(files.name[i], intensity.means.only = TRUE)[["INTENSITY"]][["MEAN"]] ## armazena uma amostra
#  newaux <- aux[pminfo$fid] ## seleciona apenas PM
#  newaux <- backgroundCorrect.matrix(newaux, normexp.method = "rma", verbose = TRUE) ## realiza correção de ruído de fundo
#  dim(newaux) <- c(NULL) 
#  aux <- insert(aux, pminfo$fid, newaux) ## devolve PM com ruído de fundo corrigido
#  aux <- aux[-((pminfo$fid)+c(1:length(pminfo$fid)))] ## extrai PM crus
#  h5write(aux, "myhdf5file.h5", "analysis/dataset", index = list(NULL,i))
#}


h5closeAll()
h5f <- H5Fopen("myhdf5file.h5")
h5g <- H5Gopen(h5f, "analysis")
#h5d <- H5Dopen(h5g, "dataset")
#h5id <- H5Dopen(h5g, "index1")

## localiza os dados de expressão na memória
#data <- realize(h5g$test)

## Reserva um espaço de memória para salvar os PM
h5createDataset(file = "myhdf5file.h5", dataset = "analysis/PM", dims = c(length(pminfo$fid), length(files.name)), maxdims = c(length(pminfo$fid), length(files.name)), storage.mode = "double", chunk = c(length(pminfo$fid), length(files.name)), level = 5, showWarnings = FALSE)

## cria objetos da classe arrayRealizationSink para processamentos em bloco,
## assim como os "blocos" que serão utilizados
dim <- as.integer(dim)
real.sink <- RealizationSink(dim = dim, dimnames = NULL, type = "double")
real.sink.pm <- RealizationSink(dim = c(nrow(pminfo), length(files.name)), dimnames = NULL, type = "double")
block.pm <- colGrid(real.sink.pm, ncol = 100)
block.pm.row <- rowGrid(real.sink.pm, nrow = 30000)
block <- colGrid(real.sink, ncol = 30)

## Salva os PM em HDF5
begin <- 1
end <- ncol(block[[1L]])
for (i in seq_along(block)) {
  aux <- h5read(h5g, "test", index = list(NULL, begin:end))
  aux <- apply(aux, 2, "[", pminfo$fid)
  h5write(aux, "myhdf5file.h5", "analysis/PM", index = list(NULL, begin:end))
  begin <- begin + ncol(block[[i]])
  ifelse(test = i == length(block)-1, yes = end <- end + ncol(block[[length(block)]]), no = end <- end + ncol(block[[i]]))
  }


## Reserva um espaço de memória para salvar os PM com background corrigido
h5createDataset(file = "myhdf5file.h5", dataset = "analysis/background(PM)", dims = c(length(pminfo$fid), length(files.name)), maxdims = c(length(pminfo$fid), length(files.name)), storage.mode = "double", chunk = c(length(pminfo$fid), length(files.name)), level = 5, showWarnings = FALSE)


## Salva os PM com background corrigido em HDF5
begin <- 1
end <- ncol(block.pm[[1L]])
for (i in seq_along(block.pm)) {
  aux <- h5read(h5g, "PM", index = list(NULL, begin:end))
  aux <- apply(aux, 2, backgroundCorrect.matrix, normexp.method = "rma", verbose = TRUE)
  h5write(aux, "myhdf5file.h5", "analysis/background(PM)", index = list(NULL, begin:end))
  begin <- begin + ncol(block.pm[[i]])
  ifelse(test = i == length(block.pm)-1, yes = end <- end + ncol(block.pm[[length(block.pm)]]), no = end <- end + ncol(block.pm[[i]]))
}


## Reserva um espaço de memória para salvar os índices dos PM com background corrigido
h5createDataset(file = "myhdf5file.h5", dataset = "analysis/indexbg(PM)", dims = c(length(pminfo$fid), length(files.name)), maxdims = c(length(pminfo$fid), length(files.name)), storage.mode = "double", chunk = c(length(pminfo$fid), length(files.name)), level = 5, showWarnings = FALSE)

## Obtém os índices dos PM com background corrigido, assim como os armazena no espaço de memória selecionado acima
begin <- 1
end <- ncol(block.pm[[1L]])
for (i in seq_along(block.pm)) {
  aux <- h5read(h5g, "background(PM)", index = list(NULL, begin:end))
  aux <- apply(aux, 2, order)
  h5write(aux, "myhdf5file.h5", "analysis/indexbg(PM)", index = list(NULL, begin:end))
  begin <- begin + ncol(block.pm[[i]])
  ifelse(test = i == length(block.pm)-1, yes = end <- end + ncol(block.pm[[length(block.pm)]]), no = end <- end + ncol(block.pm[[i]]))
}


## Reserva um espaço de memória para salvar os PM com background corrigido e normalizados
h5createDataset(file = "myhdf5file.h5", dataset = "analysis/normalize(PM)", dims = c(length(pminfo$fid), length(files.name)), maxdims = c(length(pminfo$fid), length(files.name)), storage.mode = "double", chunk = c(length(pminfo$fid), length(files.name)), level = 5, showWarnings = FALSE)

## Ordena do menor para o maior os PM com background corrigido
begin <- 1
end <- ncol(block.pm[[1L]])
for (i in seq_along(block.pm)) {
  aux <- h5read(h5g, "background(PM)", index = list(NULL, begin:end))
  aux <- apply(aux, 2, sort)
  h5write(aux, "myhdf5file.h5", "analysis/normalize(PM)", index = list(NULL, begin:end))
  begin <- begin + ncol(block.pm[[i]])
  ifelse(test = i == length(block.pm)-1, yes = end <- end + ncol(block.pm[[length(block.pm)]]), no = end <- end + ncol(block.pm[[i]]))
}


## Obtém as médias das linhas da matrix dos PM com background corrigido e ordenados e as insere em todas as colunas, assim como armazena no espaço de memória selecionado acima
begin <- 1
end <- nrow(block.pm.row[[1L]])
for (i in seq_along(block.pm.row)) {
  aux <- h5read(h5g, "normalize(PM)", index = list(begin:end, NULL))
  aux <- apply(aux, 1, mean)
  aux <- rep(aux, each = length(files.name))
  aux <- matrix(data = aux, nrow = nrow(block.pm.row[[i]]), ncol = length(files.name), byrow = TRUE)
  h5write(aux, "myhdf5file.h5", "analysis/normalize(PM)", index = list(begin:end, NULL))
  begin <- begin + nrow(block.pm.row[[i]])
  ifelse(test = i == length(block.pm.row)-1, yes = end <- end + nrow(block.pm.row[[length(block.pm.row)]]), no = end <- end + nrow(block.pm.row[[i]]))
}


## Reordena os PM em sua ordem original
begin <- 1
end <- ncol(block.pm[[1L]])
for (i in seq_along(block.pm)) {
  aux <- h5read(h5g, "normalize(PM)", index = list(NULL, begin:end))
  aux1 <- h5read(h5g, "indexbg(PM)", index = list(NULL, begin:end))
  for(j in seq(ncol(aux))){
    aux[,j] <- aux[,j][aux1[,j]]
  }
  h5write(aux, "myhdf5file.h5", "analysis/normalize(PM)", index = list(NULL, begin:end))
  begin <- begin + ncol(block.pm[[i]])
  ifelse(test = i == length(block.pm)-1, yes = end <- end + ncol(block.pm[[length(block.pm)]]), no = end <- end + ncol(block.pm[[i]]))
}



#pm <- apply(pm, 2, sort) ## ordena-os PM do menor para o maior
#mu <- apply(pm, 1, mean) ## calcula a média das linhas dos PM

## atribui a cada linha o valor de sua respectiva média
#for (i in 1:length(mu)) {
#  pm[i, ] <- rep(mu[i], length(files.name))
#}

## localiza os dados dos índices e reordena os dados de expressão
#id <- realize(h5g$index0)
#for(i in 1:length(files.name)){
#    aux <- id[ ,i]
#    newaux <- pm[,i]
#    newaux <- newaux[order(aux)]
#    pm[ ,i] <- newaux
#}

## armazena os dados com PM com background corrigido e normalizado no banco de dados original
begin <- 1
end <- ncol(block.pm[[1L]])
for (i in seq_along(block.pm)) {
  aux1 <- h5read(h5g, "normalize(PM)", index = list(NULL, begin:end))
  aux <- h5read(h5g, "test", index = list(NULL, begin:end))
  for(j in seq(ncol(aux))){
    aux[,j] <- insert(aux[,j], pminfo$fid, aux1[,j])[-((pminfo$fid)+c(1:length(pminfo$fid)))]
  }
  h5write(aux, "myhdf5file.h5", "analysis/test", index = list(NULL, begin:end))
  begin <- begin + ncol(block.pm[[i]])
  ifelse(test = i == length(block.pm)-1, yes = end <- end + ncol(block.pm[[length(block.pm)]]), no = end <- end + ncol(block.pm[[i]]))
}


#for (i in 1:length(files.name)) {
#  aux <- data[,i]
# newaux <- pm[,i]
#  aux <- insert(aux, pminfo$fid, newaux) ## devolve PM com background corrigido e normalizado
#  aux <- aux[-((pminfo$fid)+c(1:length(pminfo$fid)))] ## extrai PM apenas com background corrigido
#  h5write(aux, "myhdf5file.h5", "analysis/test", index = list(NULL,i))  
#}

## Reserva um espaço de memória para salvar os PM com background corrigido, normalizados e sumarizados
h5createDataset(file = "myhdf5file.h5", dataset = "analysis/summarize(PM)", dims = c(length(unique(pminfo$fsetid)), length(files.name)), maxdims = c(length(unique(pminfo$fsetid)), length(files.name)), storage.mode = "double", chunk = c(length(unique(pminfo$fsetid)), length(files.name)), level = 5, showWarnings = FALSE)


## cria objeto onde serão armazenadas as matrizes de PM para cada gene 
pminfo <- as_tibble(pminfo) 
genes.id <- unique(pminfo$fsetid) ## armazena todos os indentificadores dos genes 
genes.pm <- c()
exprs.values <- c()
exprs.values.matrix <- c()

begin.row <- 1
end.row <- nrow(block.pm.row[[1L]])
for (h in seq_along(block.pm.row)) {
  sapply(genes.id, invisible(function(id){
    exprs.values <- c()
    begin <- 1
    end <- ncol(block.pm[[1L]])
    pminfo.specific.gene <- filter(pminfo, pminfo$fsetid == id)
    for (i in seq_along(block.pm)) {
      aux <- h5read(h5g, "test", index = list(NULL, begin:end))
      genes.pm <- sapply(X = seq(ncol(aux)), FUN = invisible(function(j){
                    genes.pm <- cbind(genes.pm, aux[,j][pminfo.specific.gene$fid])
                  }))
      print(genes.pm)
      column.effect <- medpolish(genes.pm)$col
      genes.pm <- sweep(genes.pm, 2, column.effect)
      genes.pm <- apply(genes.pm, 2, mean)
      exprs.values <- c(exprs.values, log(genes.pm,2))
      begin <- begin + ncol(block.pm[[i]])
      ifelse(test = i == length(block.pm)-1, yes = end <- end + ncol(block.pm[[length(block.pm)]]), no = end <- end + ncol(block.pm[[i]]))
      genes.pm <- c()
    }
    exprs.values.matrix <- rbind(exprs.values.matrix, exprs.values)
  }))
  h5write(exprs.values.matrix, "myhdf5file.h5", "analysis/summarize(PM)", index = list(begin.row:end.row, NULL))
  begin.row <- begin.row + ncol(block.pm.row[[h]])
  ifelse(test = h == length(block.pm.row)-1, yes = end.row <- end.row + ncol(block.pm.row[[length(block.pm.row)]]), no = end.row <- end.row + ncol(block.pm.row[[h]]))
}


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
#original.grid <- blockGrid(data, block.length = 5000000, block.shape = "first-dim-grows-first")
#medpolish.grid <- blockGrid(data, block.length = 5000000, block.shape = "first-dim-grows-first")

## aplica o polimento de mediana em blocos 
#for(i in seq_along(original.grid)){
#  original.block <- read_block(data, original.grid[[i]])
#  medpolish.block <- medpolish(original.block)$residuals
#  return.block <- original.block - medpolish.block
#  data <- write_block(data, original.grid[[i]], return.block)
#}

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
