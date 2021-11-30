# Container of R, dplyr, ggrepel, BiocManager, pheatmap, RColorBrewer, tidyverse, EnhancedVolcano, factoextra, DEGreport, apeglm, FactoMineR

# https://cran.r-project.org/

#Image base
FROM r-base:3.6.3

#downloading the required packages

RUN apt-get update
RUN apt-get install -y procps libssl-dev libcurl4-gnutls-dev curl 
RUN apt-get install -y git libopenmpi-dev openmpi-bin openmpi-doc 
RUN apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)){install.packages('BiocManager')}" \
   && R -e "BiocManager::install('DESeq2')" \
   && R -e "BiocManager::install('EnhancedVolcano')" \
   && R -e "url='https://cran.r-project.org/src/contrib/ggplot2_3.3.5.tar.gz';install.packages(url,repos=NULL, type='source')" \ 
   && R -e "url='https://cran.r-project.org/src/contrib/tidyverse_1.3.1.tar.gz';install.packages(url, repos=NULL, type='source')" \ 
   && R -e "url='https://cran.r-project.org/src/contrib/pheatmap_1.0.12.tar.gz';install.packages(url, repos=NULL, type='source')" \ 
   && R -e "url='https://cran.r-project.org/src/contrib/RColorBrewer_1.1-2.tar.gz';install.packages(url,repos=NULL, type='source')" \
   && R -e "url='https://cran.r-project.org/src/contrib/matrixStats_0.61.0.tar.gz';install.packages(url,repos=NULL, type='source')" \  
   && R -e "BiocManager::install('factoextra')" \
   && R -e "BiocManager::install('FactoMineR')" 


ENTRYPOINT ["/usr/bin/R"] 

