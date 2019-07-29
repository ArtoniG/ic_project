## carrega os pacotes
library(pd.hg.u133.plus.2)
library(HDF5Array)
library(affyio)
library(limma)
library(R.utils)
library(dplyr)


# IMPORTAÇÃO DOS DADOS EM HDF5

import.h5 <- function(pathdir, pkganno){
  
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
  real.sink <- RealizationSink(dim = c(length(pminfo$fid), length(files.name)), dimnames = NULL, type = "double")
  block <- colGrid(real.sink, ncol = 100)   
  
  ## Obtém os índices dos PM com background corrigido, assim como os armazena no espaço de memória selecionado acima
  begin <- 1
  end <- ncol(block[[1L]])
  index.h5 <- as(matrix(nrow = length(pminfo$fid), ncol = length(files.name)), "HDF5Array")
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
  block.col <- colGrid(real.sink.pm, ncol = 146)
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
    aux <- apply(aux, 1, mean)
    aux <- rep(aux, each = length(files.name))
    aux <- matrix(data = aux, nrow = nrow(block.pm.row[[i]]), ncol = length(files.name), byrow = TRUE)
    h5write(aux, h5g, "normalize(PM)", index = list(begin:end, NULL))
    begin <- begin + nrow(block.row[[i]])
    ifelse(test = i == length(block.row)-1, yes = end <- end + nrow(block.row[[length(block.row)]]), no = end <- end + nrow(block.row[[i]]))
  }
  
  ## Reordena os PM em sua ordem original e armazena em .h5
  begin <- 1
  end <- ncol(block.col[[1L]])
  index <- ascend.indexPM()
  for (i in seq_along(block.col)) {
    aux <- h5read(h5g, "normalize(PM)", index = list(NULL, begin:end))
    #aux1 <- h5read(h5g, "indexbg(PM)", index = list(NULL, begin:end))
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


raw.with.pm <- function(){
  
  # Abre o arquivo HDF5 
  h5f <- H5Fopen("myhdf5file.h5")
  h5g <- H5Gopen(h5f, "analysis")
  
  # Define que os dados do disco serão reconhecidos como HDF5
  setRealizationBackend("HDF5Array")
  
  # Reconhece os dados como HDF5 no environment
  data <- realize(h5g$`normalize(PM)`)
  
  ## cria objetos da classe arrayRealizationSink para processamentos em blocos, assim como o "bloco" que será utilizado
  real.sink <- RealizationSink(dim = c(length(pminfo$fid), length(files.name)), dimnames = NULL, type = "double")
  block <- colGrid(real.sink, ncol = 100)   
  
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
  
}


















