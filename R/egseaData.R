# Gene set collections for the Ensemble of Gene Set Enrichment Analyses (EGSEA) package 
# 
# Author: Monther Alhamdoosh, E:m.hamdoosh@gmail.com

#' @name EGSEAdata-package
#' @aliases EGSEAdata
#' @docType package 
#' @title Gene Set Collections for the EGSEA package
#' @author Monther Alhamdoosh, Yifang Hu and Gordon K. Smyth 
#' @description This package provides R data objects for the  gene set collections that 
#' are used by the EGSEA package.
NULL

#' EGSEA analysis results on the human IL-13 dataset  
#' 
#' EGSEA analysis was performed on the il13.data from the EGSEAdata package using
#' the KEGG pathways, c2 and c5 gene set collections.    
#' @name il13.gsa
#' @docType data
#' @format An object of class EGSEAResults 
#' @source The dataset of this analysis is available in the EGSEAdata package.
NULL

#' Human IL-13 dataset  
#' 
#' The voom object calculated from the TMM normalized count matrix of RNA-seq 
#' performed on samples of human normal PBMCs, IL-13 stimulated PBMCs and IL-13R 
#' antagnonist PBBMCs. It also contains the contrast matrix of this experiment.  
#' @name il13.data
#' @docType data
#' @format A List object with two components: voom and contra. 
#' @source The count matrix of this experiment is vailable from the GEO database
#' \url{www.ncbi.nlm.nih.gov/geo/} as series GSE79027.
NULL

#' Human IL-13 dataset - Raw Counts
#' 
#' It contains the raw count matrix of RNA-seq performed on samples of human normal
#' PBMCs, IL-13 stimulated PBMCs and IL-13R antagnonist PBBMCs. It also contains the 
#' contrast and design matrices of this experiment. The gene symbols mapping is also included.  
#' @name il13.data.cnt
#' @docType data
#' @format A List object with five components: counts, group, design, contra and genes. 
#' @source The FASTQ files of this experiment are vailable from the GEO database
#' \url{www.ncbi.nlm.nih.gov/geo/} as series GSE79027.
NULL

#' Mouse mammary cell dataset
#' 
#' The voom object calculated from TMM normalized count matrix of RNA-seq
#' performed on samples of the epithelial cells of the mouse mammary glands 
#' from three populations: basal, luminal progenitor and mature luminal. It also
#' contains the contrast matrix of this experiment.  
#' @name mam.data
#' @docType data
#' @format A List object with two components: voom and contra.
#' @source The count matrix of this experiment is vailable from the GEO database
#' \url{www.ncbi.nlm.nih.gov/geo/} as series GSE63310.
NULL 

#' @title GeneSetDB Human Collections
#' 
#' @description  Human gene set collections from the GeneSetDB
#' @name gsetdb.human
#' @docType data
#' @format list 
#' @source
#' Araki Hiromitsu,Knapp Christoph,Tsai Peter and Print Cristin(2012), 
#' GeneSetDB: A comprehensive meta-database, statistical and visualisation framework 
#' for gene set analysis, FEBS Open Bio, 2, doi: 10.1016/j.fob.2012.04.003 
#' 
NULL

#' @title GeneSetDB Mouse Collections
#' 
#' @description  Mouse gene set collections from the GeneSetDB
#' @name gsetdb.mouse
#' @docType data
#' @format list 
#' @source
#' Araki Hiromitsu,Knapp Christoph,Tsai Peter and Print Cristin(2012), 
#' GeneSetDB: A comprehensive meta-database, statistical and visualisation framework 
#' for gene set analysis, FEBS Open Bio, 2, doi: 10.1016/j.fob.2012.04.003 
#' 
NULL

#' @title GeneSetDB Rat Collections
#' 
#' @description  Rat gene set collections from the GeneSetDB
#' @name gsetdb.rat
#' @docType data
#' @format list 
#' @source
#' Araki Hiromitsu,Knapp Christoph,Tsai Peter and Print Cristin(2012), 
#' GeneSetDB: A comprehensive meta-database, statistical and visualisation framework 
#' for gene set analysis, FEBS Open Bio, 2, doi: 10.1016/j.fob.2012.04.003 
#' 
NULL


#' @title KEGG Pathways Collections
#' 
#' @description  Human, Mouse and Rat gene set collections from the KEGG database
#' @name kegg.pathways 
#' @docType data
#' @format list 
#' @source
#' Luo, W., Friedman, M., Shedden K., Hankenson, K. and Woolf, P GAGE: Generally Applicable
#' Gene Set Enrichment for Pathways Analysis. BMC Bioinformatics 2009, 10:161
#' 
NULL

#' @title MSigDB Gene Set Collections 
#' 
#' @description  Gene set collections from the MSigDB database Version 5
#' @name msigdb
#' @docType data
#' @format list 
#' @source
#' Aravind Subramanian, Pablo Tamayo, Vamsi K. Mootha, Sayan Mukherjee, Benjamin L. Ebert, Michael A. Gillette, Amanda Paulovich, Scott L. Pomeroy, Todd R. Golub, Eric S. Lander, and Jill P. Mesirov
#' Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles
#' PNAS 2005 102 (43) 15545-15550
#' 
NULL


#' @title Mouse H MSigDB Gene Set Collections
#' 
#' @description  Mouse orthologs gene set collections from the MSigDB database Version 5
#' @name Mm.H
#' @docType data
#' @format list 
#' @source http://bioinf.wehi.edu.au/software/MSigDB/
#' 
NULL


#' @title Mouse C2 MSigDB Gene Set Collections
#' 
#' @description  Mouse orthologs gene set collections from the MSigDB database Version 5
#' @name Mm.c2
#' @docType data
#' @format list 
#' @source http://bioinf.wehi.edu.au/software/MSigDB/
#' 
NULL

#' @title Mouse C3 MSigDB Gene Set Collections
#' 
#' @description  Mouse orthologs gene set collections from the MSigDB database Version 5
#' @name Mm.c3
#' @docType data
#' @format list 
#' @source http://bioinf.wehi.edu.au/software/MSigDB/
#' 
NULL

#' @title Mouse C4 MSigDB Gene Set Collections
#' 
#' @description  Mouse orthologs gene set collections from the MSigDB database Version 5
#' @name Mm.c4
#' @docType data
#' @format list 
#' @source http://bioinf.wehi.edu.au/software/MSigDB/
#' 
NULL


#' @title Mouse C5 MSigDB Gene Set Collections
#' 
#' @description  Mouse orthologs gene set collections from the MSigDB database Version 4
#' @name Mm.c5
#' @docType data
#' @format list 
#' @source http://bioinf.wehi.edu.au/software/MSigDB/
#' 
NULL

#' @title Mouse C6 MSigDB Gene Set Collections
#' 
#' @description  Mouse orthologs gene set collections from the MSigDB database Version 5
#' @name Mm.c6
#' @docType data
#' @format list 
#' @source http://bioinf.wehi.edu.au/software/MSigDB/
#' 
NULL

#' @title Mouse C7 MSigDB Gene Set Collections
#' 
#' @description  Mouse orthologs gene set collections from the MSigDB database Version 5
#' @name Mm.c7
#' @docType data
#' @format list 
#' @source http://bioinf.wehi.edu.au/software/MSigDB/
#' 
NULL


