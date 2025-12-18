# LAINFS - Controle de Qualidade e Trimagem
![Badge em Desenvolvimento](http://img.shields.io/static/v1?label=STATUS&message=EmDesenvolvimento&color=YELLOW&style=for-the-badge)

## Sistema Operacional Linux-Ubuntu 

Este projeto utiliza o pacote **ShortRead** (Bioconductor). \
Em sistemas Linux (especialmente Ubuntu), é comum ocorrer erro durante a instalação por falta de dependências do sistema operacional necessárias para compilar pacotes a partir do código-fonte. \
O ShortRead depende de vários pacotes do Bioconductor, que utilizam código em C/C++ e precisam de bibliotecas externas do sistema, se essas bibliotecas não estiverem instaladas, o R falha ao compilar os pacotes.

Caso tenha encontrado problemas para rodar o projeto e o ShortRead ou tidyverse não estejam instalados, abra um terminal e execute:

```bash
sudo apt update

sudo apt install -y \
  build-essential \
  gfortran \
  libblas-dev \
  liblapack-dev \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev \
  libpng-dev \
  libjpeg-dev \
  zlib1g-dev \
  libbz2-dev \
  liblzma-dev
  
  libfribidi-dev \
  libfontconfig1-dev \
  libharfbuzz-dev \
  libfreetype6-dev \
  libpng-dev \
  libtiff5-dev \
  libjpeg-dev

```

Esse comando instala todas as bibliotecas necessárias para compilar os pacotes usados pelo ShortRead. 

Após instalar as dependências do sistema, reinicie o RStudio para garantir que o R reconheça as bibliotecas recém-instaladas. 

Então é possível rodar o script "dependencies.R" para instalar todas as bibliotecas necessárias!
