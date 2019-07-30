## carrega os pacotes
library(pd.hg.u133.plus.2)
library(HDF5Array)
library(affyio)
library(limma)
library(R.utils)
library(dplyr)


# IMPORTAÇÃO DOS DADOS EM HDF5

import.h5 <- function(pathdir){
  
  ## Vai até o diretório no qual estão os arquivos .CEL
  setwd(pathdir)
  
  ## cria um objeto que armazena os nomes dos dados que serão carregados
  files.name <- as.character(list.files(all.files = TRUE, pattern = "[gz]", recursive = FALSE, include.dirs = FALSE, full.names = FALSE, no.. = TRUE))

  ## testa se todos os CELs são do mesmo pacote
  input.kind <- lapply(files.name, read.celfile.header)
  same.kind <- length(unique(sapply(input.kind, "[[", "cdfName"))) == 1
  if (!same.kind)
    stop("All files need to be from the same package.")
  
  ## extrai as dimensões do banco de dados
  nrows <- prod(input.kind[[1]][["CEL dimensions"]])
  dim <- c(nrows, length(files.name))
  
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
  dim <- as.integer(dim)
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
  
  # Define que os dados do disco serão reconhecidos como HDF5
  setRealizationBackend("HDF5Array")
  
  # Reconhece os dados como HDF5 no environment
  data <- realize(h5g$dataset)
  
  # Fecha o arquivo .h5
  h5closeAll()
  
  # Retorna o banco de dados em formato .h5
  return(data)
}

####################################################################################################################

# SALVA OS PERFECT MATCH EM HDF5

importPM.h5 <- function(){
  
  ## cria conexão com banco de dados do pacote e extrai as informações dos perfect match
  pkganno <- "pd.hg.u133.plus.2"  ### Why is not input.kind[[1]][["cdfName"]] ?
  library(pkganno, character.only = TRUE)
  conn <- db(get(pkganno))
  pminfo <- dbGetQuery(conn, "SELECT * FROM pmfeature")
  
  ## cria um objeto que armazena os nomes dos dados que serão carregados
  files.name <- as.character(list.files(all.files = TRUE, pattern = "[gz]", recursive = FALSE, include.dirs = FALSE, full.names = FALSE, no.. = TRUE))
  
  ## extrai as dimensões do banco de dados
  input.kind <- read.celfile.header(files.name[1])
  nrows <- prod(input.kind[["CEL dimensions"]])
  dim <- c(nrows, length(files.name))
  
  # Abre o arquivo HDF5 
  h5f <- H5Fopen("myhdf5file.h5")
  h5g <- H5Gopen(h5f, "analysis")
  
  ## Reserva um espaço de memória para salvar os PM
  h5createDataset(file = h5g, dataset = "PM", dims = c(length(pminfo$fid), length(files.name)), maxdims = c(length(pminfo$fid), length(files.name)), storage.mode = "double", chunk = c(length(pminfo$fid), length(files.name)), level = 5, showWarnings = FALSE)

  ## cria objetos da classe arrayRealizationSink para processamentos em blocos, assim como o "bloco" que será utilizado
  dim <- as.integer(dim)
  real.sink <- RealizationSink(dim = dim, dimnames = NULL, type = "double")
  block <- colGrid(real.sink, ncol = 100)  
    
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
  
  # Define que os dados do disco serão reconhecidos como HDF5
  setRealizationBackend("HDF5Array")
  
  # Reconhece os PM como HDF5 no environment
  data <- realize(h5g$PM)
  
  # Fecha o arquivo .h5
  h5closeAll()
  
  # Retorna os PM em formato .h5
  return(data)
}


######################################################################################################################################


# SALVA OS PERFECT MATCH COM BACKGROUND CORRIGIDO EM HDF5

background.h5  <- function(){
  
  ## cria conexão com banco de dados do pacote e extrai as informações dos perfect match
  pkganno <- "pd.hg.u133.plus.2"  ### Why is not input.kind[[1]][["cdfName"]] ?
  library(pkganno, character.only = TRUE)
  conn <- db(get(pkganno))
  pminfo <- dbGetQuery(conn, "SELECT * FROM pmfeature")
  
  ## cria um objeto que armazena os nomes dos dados que serão carregados
  files.name <- as.character(list.files(all.files = TRUE, pattern = "[gz]", recursive = FALSE, include.dirs = FALSE, full.names = FALSE, no.. = TRUE))
  
  # Abre o arquivo HDF5 
  h5f <- H5Fopen("myhdf5file.h5")
  h5g <- H5Gopen(h5f, "analysis")
  
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
  
  # Reconhece os dados como HDF5 no environment
  data <- realize(h5g$`background(PM)`)
  
  # Fecha o arquivo .h5
  h5closeAll()
  
  # Retorna o banco de dados em formato .h5
  return(data)  
}


#######################################################################################################################################


# FUNÇÃO QUE RETORNA UM OBJECTO .h5 COM OS ÍNDICES DOS PM COM BACKGROUND CORRIGIDO COMO SE ESTIVESSEM EM ORDEM CRESCENTE

ascend.indexPM <- function(){
  
  # Abre o arquivo HDF5 
  h5f <- H5Fopen("myhdf5file.h5")
  h5g <- H5Gopen(h5f, "analysis")
  
  # Define que os dados do disco serão reconhecidos como HDF5
  setRealizationBackend("HDF5Array")
  
  # Reconhece os dados como HDF5 no environment
  data <- realize(h5g$`background(PM)`)
  
  ## cria objetos da classe arrayRealizationSink para processamentos em blocos, assim como o "bloco" que será utilizado
  real.sink <- RealizationSink(dim = dim(data), dimnames = NULL, type = "double")
  block <- colGrid(real.sink, ncol = 100)   
  
  ## Obtém os índices dos PM com background corrigido, assim como os armazena no espaço de memória selecionado acima
  begin <- 1
  end <- ncol(block[[1L]])
  index.h5 <- as(matrix(nrow = nrow(data), ncol = ncol(data)), "HDF5Array")
  for (i in seq_along(block)) {
    aux <- read_block(data, block[[i]])
    aux <- apply(aux, 2, order)
    index.h5 <- as(write_block(index.h5, block[[i]], aux), "HDF5Array")   
    begin <- begin + ncol(block[[i]])
    ifelse(test = i == length(block)-1, yes = end <- end + ncol(block[[length(block)]]), no = end <- end + ncol(block[[i]]))
  }
  
  # Fecha o arquivo .h5
  h5closeAll()
  
  # Retorna o banco de dados em formato .h5
  return(index.h5)  
  
}


###################################################################################################################################################  
  

# SALVA OS PERFECT MATCH COM BACKGROUND CORRIGIDO E NORMALIZADOS EM HDF5

normalize.h5 <- function(){
  
  ## cria conexão com banco de dados do pacote e extrai as informações dos perfect match
  pkganno <- "pd.hg.u133.plus.2"  ### Why is not input.kind[[1]][["cdfName"]] ?
  library(pkganno, character.only = TRUE)
  conn <- db(get(pkganno))
  pminfo <- dbGetQuery(conn, "SELECT * FROM pmfeature")
  
  ## cria um objeto que armazena os nomes dos dados que serão carregados
  files.name <- as.character(list.files(all.files = TRUE, pattern = "[gz]", recursive = FALSE, include.dirs = FALSE, full.names = FALSE, no.. = TRUE))
  
  # Abre o arquivo HDF5 
  h5f <- H5Fopen("myhdf5file.h5")
  h5g <- H5Gopen(h5f, "analysis")
  
  ## Reserva um espaço de memória para salvar os PM com background corrigido e normalizados
  h5createDataset(file = h5g, dataset = "normalize(PM)", dims = c(length(pminfo$fid), length(files.name)), maxdims = c(length(pminfo$fid), length(files.name)), storage.mode = "double", chunk = c(length(pminfo$fid), length(files.name)), level = 5, showWarnings = FALSE)  

  ## cria objetos da classe arrayRealizationSink para processamentos em blocos, assim como os "blocos" que serão utilizados
  real.sink <- RealizationSink(dim = c(length(unique(pminfo$fsetid)), length(files.name)), dimnames = NULL, type = "double")
  block.col <- colGrid(real.sink, ncol = 146)
  block.row <- rowGrid(real.sink, nrow = 54675)
  
    
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
  index <- ascend.indexPM()
  h5closeAll()
  h5f <- H5Fopen("myhdf5file.h5")
  h5g <- H5Gopen(h5f, "analysis")
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
  
  # Define que os dados do disco serão reconhecidos como HDF5
  setRealizationBackend("HDF5Array")
  
  # Reconhece os dados como HDF5 no environment
  data <- realize(h5g$`normalize(PM)`)
  
  # Fecha o arquivo .h5
  h5closeAll()
  
  # Retorna o banco de dados em formato .h5
  return(data)    
  
}


#######################################################################################################################################################################


# CRIA UM OBJETO EM HDF5 COMO O RAW DATA, MAS SUBSTITUI OS PM RAW POR PM NORMALIZADOS  

raw.with.pm.h5 <- function(){

  ## cria conexão com banco de dados do pacote e extrai as informações dos perfect match
  pkganno <- "pd.hg.u133.plus.2"  ### Why is not input.kind[[1]][["cdfName"]] ?
  library(pkganno, character.only = TRUE)
  conn <- db(get(pkganno))
  pminfo <- dbGetQuery(conn, "SELECT * FROM pmfeature")  
    
  # Abre o arquivo HDF5 
  h5f <- H5Fopen("myhdf5file.h5")
  h5g <- H5Gopen(h5f, "analysis")
  
  # Define que os dados do disco serão reconhecidos como HDF5
  setRealizationBackend("HDF5Array")
  
  # Reconhece os dados como HDF5 no environment
  data <- realize(h5g$dataset)
  normalized.pm <- realize(h5g$`normalize(PM)`)
  
  ## cria objetos da classe arrayRealizationSink para processamentos em blocos, assim como o "bloco" que será utilizado
  real.sink <- RealizationSink(dim = dim(data), dimnames = NULL, type = "double")
  real.sink.pm <- RealizationSink(dim = dim(normalized.pm), dimnames = NULL, type = "double")
  block <- colGrid(real.sink, ncol = 100)
  block.pm <- colGrid(real.sink.pm, ncol = 100)   
  
  ## cria o objeto no qual serão salvos os raw data com os PM normalizados
  raw.with.pm <- as(matrix(nrow = nrow(real.sink), ncol = ncol(real.sink)), "HDF5Array")
  
  ## armazena os dados com PM com background corrigido e normalizado no banco de dados original
  begin <- 1
  end <- ncol(block.pm[[1L]])
  for (i in seq_along(block.pm)) {
    aux1 <- read_block(normalized.pm, block.pm[[i]])
    aux <- read_block(data, block[[i]])
    for(j in seq(ncol(aux))){
      aux[,j] <- insert(aux[,j], pminfo$fid, aux1[,j])[-((pminfo$fid)+c(1:length(pminfo$fid)))]
    }
    raw.with.pm <- as(write_block(raw.with.pm, block[[i]], aux), "HDF5Array")
    begin <- begin + ncol(block.pm[[i]])
    ifelse(test = i == length(block.pm)-1, yes = end <- end + ncol(block.pm[[length(block.pm)]]), no = end <- end + ncol(block.pm[[i]]))
  }
  
  # Fecha o arquivo .h5
  h5closeAll()
  
  # Retorna os dados em formato .h5
  return(raw.with.pm)      
}


###########################################################################################################################################################


# CRIA UMA MATRIZ COM APENAS OS PM DE UM GENE ESPECÍFICO

gene.specific.matrix <- function(raw.with.pm.h5.result, id){
  
  ## cria conexão com banco de dados do pacote e extrai as informações dos perfect match
  pkganno <- "pd.hg.u133.plus.2"  ### Why is not input.kind[[1]][["cdfName"]] ?
  library(pkganno, character.only = TRUE)
  conn <- db(get(pkganno))
  pminfo <- dbGetQuery(conn, "SELECT * FROM pmfeature")
  pminfo <- as_tibble(pminfo) 
  pminfo <- filter(pminfo, pminfo$fsetid == id)
  
  ## cria o objeto no qual será armazenada a matriz com os pm de um gene especifico
  genes.pm <- as(matrix(nrow = nrow(pminfo), ncol = ncol(raw.with.pm.h5.result)), "HDF5Array")
  
  ## cria objetos da classe arrayRealizationSink para processamentos em blocos, assim como o "bloco" que será utilizado
  real.sink <- RealizationSink(dim = dim(raw.with.pm.h5.result), dimnames = NULL, type = "double")
  real.sink.pm <- RealizationSink(dim = dim(genes.pm), dimnames = NULL, type = "double")
  block <- colGrid(real.sink, ncol = 50)
  block.pm <- colGrid(real.sink.pm, ncol = 50)
  
  # Armazena a matriz dos PM para um gene específico em .h5
  for (i in seq_along(block)) {
    aux1 <- c()
    aux <- read_block(data, block[[i]])
    for (j in seq(ncol(aux))) {
      aux1 <- cbind(aux1, aux[,j][pminfo$fid])
    }
    genes.pm <- write_block(genes.pm, block.pm[[i]], aux1)
  }
  
  # Retorna os dados em formato .h5
  return(genes.pm)
}
  

################################################################################################################################################################

# REALIZA A SUMARIZAÇÃO DOS DADOS

summarize.by.block <- function(gene.specific.matrix.result){

  ## cria conexão com banco de dados do pacote e extrai as informações dos perfect match
  pkganno <- "pd.hg.u133.plus.2"  ### Why is not input.kind[[1]][["cdfName"]] ?
  library(pkganno, character.only = TRUE)
  conn <- db(get(pkganno))
  pminfo <- as_tibble(dbGetQuery(conn, "SELECT * FROM pmfeature"))
  featinfo <- as_tibble(dbGetQuery(conn, "SELECT * FROM featureSet"))
  
  ## Abre o arquivo HDF5 
  h5f <- H5Fopen("myhdf5file.h5")
  h5g <- H5Gopen(h5f, "analysis")
  
  ## Reserva um espaço de memória para salvar os PM com background corrigido, normalizados e sumarizados
  h5createDataset(file = h5g, dataset = "summarize(PM)", dims = c(length(unique(pminfo$fsetid)), length(files.name)), maxdims = c(length(unique(pminfo$fsetid)), length(files.name)), storage.mode = "double", chunk = c(length(unique(pminfo$fsetid)), length(files.name)), level = 5, showWarnings = FALSE)    

  ## armazena todos os indentificadores dos genes
  genes.id <- unique(pminfo$fsetid)
  
  ## salva o resultado de raw.with.pm.h5 em um objeto
  raw.with.pm.result <- raw.with.pm.h5()
  
  ## cria um objeto no qual serão salvos os dados de expressão
  exprs.values <- c()
    #as(matrix(nrow = length(genes.id), ncol = ncol(raw.with.pm.result)), "HDF5Array")  

  ## aplica a sumarização e salva os valores de expressão no objeto exprs.value
  for (i in seq_along(genes.id)) {
    genes.pm <- gene.specific.matrix(raw.with.pm.result, genes.id[i])
    genes.pm <- log(colMeans(sweep(genes.pm, 2, medpolish(genes.pm)$col)),2)
    exprs.values <- rbind(exprs.values, genes.pm)
    #exprs.values[i,] <- genes.pm
  }
  
  ## Salva os valores de expressão em HDF5
  h5writeDataset(exprs.values, h5g, "summarize(PM)")
  
  ## Define que os dados do disco serão reconhecidos como HDF5
  setRealizationBackend("HDF5Array")
  
  ## Reconhece os valores de expressão no disco como HDF5 e salva no objeto normalized
  summarized <- realize(h5g$`summarize(PM)`)
  
  ## Nomeia as linhas e colunas do conjunto de dados
  genes.id <- order(genes.id)
  genes.id <- featinfo$man_fsetid[genes.id]
  colnames(summarized) <- files.name
  row.names(summarized) <- genes.id
  
  # Fecha o arquivo .h5
  h5closeAll()
  
  # Retorna os dados em formato .h5
  return(summarized)    
}
