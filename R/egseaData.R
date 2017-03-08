# Gene set collections for the Ensemble of Gene Set Enrichment Analyses (EGSEA) package 
# 
# Author: Monther Alhamdoosh, E:m.hamdoosh@gmail.com

#' @name EGSEAdata-package
#' @aliases EGSEAdata
#' @docType package 
#' @title Gene Set Collections for the EGSEA package
#' @author Monther Alhamdoosh, Yifang Hu and Gordon K. Smyth 
#' @description This package includes gene set collections that are used for 
#' the Ensemble of Gene Set Enrichment Analyses (EGSEA) method for gene set testing. 
#' It includes Human and Mouse versions of the MSidDB (Subramanian, et al. (2005) 
#' PNAS, 102(43):15545-15550) and GeneSetDB (Araki, et al. (2012) FEBS Open Bio, 2:76-82) 
#' collections.
#' @details While the gene set collections in MSigDB and GeneSetDB have different names and 
#' purposes, some of these collections overlap. For example, both databases contain a Gene Ontology
#' collection but MSigDB's collection aimed for a higher level of abstraction for the GO terms.  
#' @references 
#' Monther Alhamdoosh, Milica Ng, Nicholas J. Wilson, Julie M. Sheridan, Huy Huynh, 
#' Michael J. Wilson, Matthew E. Ritchie; Combining multiple tools outperforms 
#' individual methods in gene set enrichment analyses. Bioinformatics 2017; 
#' 33 (3): 414-424. doi: 10.1093/bioinformatics/btw623
NULL

#' EGSEA analysis results on the human IL-13 dataset  
#' 
#' EGSEA analysis was performed on the il13.data from the EGSEAdata package using
#' the KEGG pathways, c2 and c5 gene set collections. Type show(il13.gsa) to see
#' the version of datasets/packages that were used.  
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
#' @details Procedure \cr 
#' 1. The Human GMT file was downloaded from the website. \cr 
#' 2. The gene set sources and categories were manually compiled from the Help page. \cr 
#' 3. An R list was created for the gene set categories.  \cr 
#' 4. An annotation data frame was created for the gene sets. \cr 
#' 5. An R data object was written using save(). 
#' @name gsetdb.human
#' @docType data
#' @format list 
#' @seealso Invoke egsea.data() to find out the current version and latest download/update date.
#' @source
#' Araki Hiromitsu,Knapp Christoph,Tsai Peter and Print Cristin(2012), 
#' GeneSetDB: A comprehensive meta-database, statistical and visualisation framework 
#' for gene set analysis, FEBS Open Bio, 2, doi: 10.1016/j.fob.2012.04.003 
#' Downloaded from http://www.genesetdb.auckland.ac.nz/
#' 
NULL

#' @title GeneSetDB Mouse Collections
#' 
#' @description  Mouse gene set collections from the GeneSetDB
#' @details Procedure \cr 
#' 1. The Mouse GMT file was downloaded from the website. \cr 
#' 2. The gene set sources and categories were manually compiled from the Help page. \cr 
#' 3. An R list was created for the gene set categories.  \cr 
#' 4. An annotation data frame was created for the gene sets. \cr 
#' 5. An R data object was written using save(). 
#' @name gsetdb.mouse
#' @docType data
#' @format list 
#' @seealso Invoke egsea.data() to find out the current version and latest download/update date.
#' @source
#' Araki Hiromitsu,Knapp Christoph,Tsai Peter and Print Cristin(2012), 
#' GeneSetDB: A comprehensive meta-database, statistical and visualisation framework 
#' for gene set analysis, FEBS Open Bio, 2, doi: 10.1016/j.fob.2012.04.003 
#' Downloaded from http://www.genesetdb.auckland.ac.nz/
#' 
NULL

#' @title GeneSetDB Rat Collections
#' 
#' @description  Rat gene set collections from the GeneSetDB
#' @details Procedure \cr 
#' 1. The Rat GMT file was downloaded from the website. \cr 
#' 2. The gene set sources and categories were manually compiled from the Help page. \cr 
#' 3. An R list was created for the gene set categories.  \cr 
#' 4. An annotation data frame was created for the gene sets. \cr 
#' 5. An R data object was written using save(). 
#' @name gsetdb.rat
#' @docType data
#' @format list 
#' @seealso Invoke egsea.data() to find out the current version and latest download/update date.
#' @source
#' Araki Hiromitsu,Knapp Christoph,Tsai Peter and Print Cristin(2012), 
#' GeneSetDB: A comprehensive meta-database, statistical and visualisation framework 
#' for gene set analysis, FEBS Open Bio, 2, doi: 10.1016/j.fob.2012.04.003 
#' Downloaded from http://www.genesetdb.auckland.ac.nz/
#' 
NULL


#' @title KEGG Pathways Collections
#' 
#' @description  Human, Mouse and Rat gene set collections from the KEGG database
#' @details The collections were generated using the following R script \cr 
#'  library(gage) \cr 
#' 	kegg.pathways = list() \cr 
#' 	species.all = c("human", "mouse", "rat") \cr 
#' 	for (species in species.all){ \cr 
#' 		kegg = kegg.gsets(species = species, id.type = "kegg") \cr 
#' 		kegg.pathways[[species]] = kegg \cr 
#' 	} \cr 
#' 	save(kegg.pathways, file='kegg.pathways.rda')
#' @name kegg.pathways 
#' @docType data
#' @format list 
#' @seealso Invoke egsea.data() to find out the current version and latest download/update date.
#' @source
#' Luo, W., Friedman, M., Shedden K., Hankenson, K. and Woolf, P GAGE: Generally Applicable
#' Gene Set Enrichment for Pathways Analysis. BMC Bioinformatics 2009, 10:161
#' Obtained from \pkg{gage} using the function \code{kegg.gsets()}
#' 
NULL

#' @title MSigDB Gene Set Collections 
#' 
#' @description  Gene set collections from the MSigDB database
#' @details Procedure \cr 
#' 1. The current msigdb_vx.xml file was downloaded. \cr 
#' 2. It was parsed using xmlParse() and then converted to list using xmlToList()
#' 3. The resulting list was written to an RData file using save(). \cr 
#' This dataset is mainly used to extract MSigDB gene set annotation and the human gene set 
#' collections.
#' @name msigdb
#' @docType data
#' @format list 
#' @seealso Invoke egsea.data() to find out the current version and latest download/update date.
#' @source
#' Aravind Subramanian, Pablo Tamayo, Vamsi K. Mootha, Sayan Mukherjee, Benjamin L. Ebert, Michael A. Gillette, Amanda Paulovich, Scott L. Pomeroy, Todd R. Golub, Eric S. Lander, and Jill P. Mesirov
#' Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles
#' PNAS 2005 102 (43) 15545-15550  
#' 
NULL


#' @title Mouse H MSigDB Gene Set Collections
#' 
#' @description  Mouse orthologs gene set collections from the MSigDB database
#' @details Procedure \cr 
#' 1. The current msigdb_vx.xml file was downloaded. \cr 
#' 2. Human Entrez Gene IDs were mapped to Mouse Entrez Gene IDs, using the HGNC 
#' Comparison of Orthology Predictions (HCOP). \cr 
#' 3. The collection was converted to a list in R, and written to a RData file using save().
#' @name Mm.H
#' @docType data
#' @format list 
#' @seealso Invoke egsea.data() to find out the current version and latest download/update date.
#' @source Downloaded from http://bioinf.wehi.edu.au/software/MSigDB/. 
#' Extracted from
#' Aravind Subramanian, Pablo Tamayo, Vamsi K. Mootha, Sayan Mukherjee, Benjamin L. Ebert, Michael A. Gillette, Amanda Paulovich, Scott L. Pomeroy, Todd R. Golub, Eric S. Lander, and Jill P. Mesirov
#' Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles
#' PNAS 2005 102 (43) 15545-15550
NULL


#' @title Mouse C2 MSigDB Gene Set Collections
#' 
#' @description  Mouse orthologs gene set collections from the MSigDB database
#' @details Procedure \cr 
#' 1. The current msigdb_vx.xml file was downloaded. \cr 
#' 2. Human Entrez Gene IDs were mapped to Mouse Entrez Gene IDs, using the HGNC 
#' Comparison of Orthology Predictions (HCOP). \cr 
#' 3. The collection was converted to a list in R, and written to a RData file using save().
#' @name Mm.c2
#' @docType data
#' @format list 
#' @seealso Invoke egsea.data() to find out the current version and latest download/update date.
#' @source Downloaded from http://bioinf.wehi.edu.au/software/MSigDB/. 
#' Extracted from
#' Aravind Subramanian, Pablo Tamayo, Vamsi K. Mootha, Sayan Mukherjee, Benjamin L. Ebert, Michael A. Gillette, Amanda Paulovich, Scott L. Pomeroy, Todd R. Golub, Eric S. Lander, and Jill P. Mesirov
#' Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles
#' PNAS 2005 102 (43) 15545-15550
#' 
NULL

#' @title Mouse C3 MSigDB Gene Set Collections
#' 
#' @description  Mouse orthologs gene set collections from the MSigDB database
#' @details Procedure \cr 
#' 1. The current msigdb_vx.xml file was downloaded. \cr 
#' 2. Human Entrez Gene IDs were mapped to Mouse Entrez Gene IDs, using the HGNC 
#' Comparison of Orthology Predictions (HCOP). \cr 
#' 3. The collection was converted to a list in R, and written to a RData file using save().
#' @name Mm.c3
#' @docType data
#' @format list 
#' @seealso Invoke egsea.data() to find out the current version and latest download/update date.
#' @source Downloaded from http://bioinf.wehi.edu.au/software/MSigDB/. 
#' Extracted from 
#' Aravind Subramanian, Pablo Tamayo, Vamsi K. Mootha, Sayan Mukherjee, Benjamin L. Ebert, Michael A. Gillette, Amanda Paulovich, Scott L. Pomeroy, Todd R. Golub, Eric S. Lander, and Jill P. Mesirov
#' Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles
#' PNAS 2005 102 (43) 15545-15550
#' 
NULL

#' @title Mouse C4 MSigDB Gene Set Collections
#' 
#' @description  Mouse orthologs gene set collections from the MSigDB database
#' @details Procedure \cr 
#' 1. The current msigdb_vx.xml file was downloaded. \cr 
#' 2. Human Entrez Gene IDs were mapped to Mouse Entrez Gene IDs, using the HGNC 
#' Comparison of Orthology Predictions (HCOP). \cr 
#' 3. The collection was converted to a list in R, and written to a RData file using save().
#' @name Mm.c4
#' @docType data
#' @format list 
#' @seealso Invoke egsea.data() to find out the current version and latest download/update date.
#' @source Downloaded from http://bioinf.wehi.edu.au/software/MSigDB/
#' Extracted from
#' Aravind Subramanian, Pablo Tamayo, Vamsi K. Mootha, Sayan Mukherjee, Benjamin L. Ebert, Michael A. Gillette, Amanda Paulovich, Scott L. Pomeroy, Todd R. Golub, Eric S. Lander, and Jill P. Mesirov
#' Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles
#' PNAS 2005 102 (43) 15545-15550
#' 
NULL


#' @title Mouse C5 MSigDB Gene Set Collections
#' 
#' @description  Mouse orthologs gene set collections from the MSigDB database
#' @details Procedure \cr 
#' 1. The current msigdb_vx.xml file was downloaded. \cr 
#' 2. Human Entrez Gene IDs were mapped to Mouse Entrez Gene IDs, using the HGNC 
#' Comparison of Orthology Predictions (HCOP). \cr 
#' 3. The collection was converted to a list in R, and written to a RData file using save().
#' @name Mm.c5
#' @docType data
#' @format list 
#' @seealso Invoke egsea.data() to find out the current version and latest download/update date.
#' @source Downloaded from http://bioinf.wehi.edu.au/software/MSigDB/. 
#' Extracted from
#' Aravind Subramanian, Pablo Tamayo, Vamsi K. Mootha, Sayan Mukherjee, Benjamin L. Ebert, Michael A. Gillette, Amanda Paulovich, Scott L. Pomeroy, Todd R. Golub, Eric S. Lander, and Jill P. Mesirov
#' Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles
#' PNAS 2005 102 (43) 15545-15550
#' 
NULL

#' @title Mouse C6 MSigDB Gene Set Collections
#' 
#' @description  Mouse orthologs gene set collections from the MSigDB database
#' @details Procedure \cr 
#' 1. The current msigdb_vx.xml file was downloaded.\cr 
#' 2. Human Entrez Gene IDs were mapped to Mouse Entrez Gene IDs, using the HGNC 
#' Comparison of Orthology Predictions (HCOP).\cr 
#' 3. The collection was converted to a list in R, and written to a RData file using save().
#' @name Mm.c6
#' @docType data
#' @format list 
#' @seealso Invoke egsea.data() to find out the current version and latest download/update date.
#' @source Downloaded from http://bioinf.wehi.edu.au/software/MSigDB/. 
#' Extracted from
#' Aravind Subramanian, Pablo Tamayo, Vamsi K. Mootha, Sayan Mukherjee, Benjamin L. Ebert, Michael A. Gillette, Amanda Paulovich, Scott L. Pomeroy, Todd R. Golub, Eric S. Lander, and Jill P. Mesirov
#' Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles
#' PNAS 2005 102 (43) 15545-15550
#' 
NULL

#' @title Mouse C7 MSigDB Gene Set Collections
#' 
#' @description  Mouse orthologs gene set collections from the MSigDB database
#' @details Procedure \cr 
#' 1. The current msigdb_vx.xml file was downloaded. \cr 
#' 2. Human Entrez Gene IDs were mapped to Mouse Entrez Gene IDs, using the HGNC 
#' Comparison of Orthology Predictions (HCOP).\cr 
#' 3. The collection was converted to a list in R, and written to a RData file using save().
#' @name Mm.c7
#' @docType data
#' @format list 
#' @seealso Invoke egsea.data() to find out the current version and latest download/update date.
#' @source Downloaded from http://bioinf.wehi.edu.au/software/MSigDB/. 
#' Extracted from
#' Aravind Subramanian, Pablo Tamayo, Vamsi K. Mootha, Sayan Mukherjee, Benjamin L. Ebert, Michael A. Gillette, Amanda Paulovich, Scott L. Pomeroy, Todd R. Golub, Eric S. Lander, and Jill P. Mesirov
#' Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles
#' PNAS 2005 102 (43) 15545-15550
#' 
NULL

#' @title EGSEAdata databases information
#' @description It displays information about the available gene set collections in 
#' EGSEAdata for a species.
#' @details It prints out for each database: the database name, version, update/download
#' date, data source, supported species, gene set collections, the names of the related
#'  R data objects and the number of gene sets in each collection. 
#' @param species character, a species name and used to retreive the number of gene sets
#' for that particular species. Default is "human". Accepted values are  
#' "human", "homo sapiens", "hs", "mouse", "mus musculus", "mm",
#' "rat", "rattus norvegicus" or "rn".
#' @param simple logical, whether to display the number of gene sets in each 
#' collection or not.
#' @param returnInfo, logical, whether to print out the databases information or 
#' return it as a list.
#' 
#' @return nothing.
#' @export
#' @examples
#' # Example of egsea.data 
#' egsea.data()
#' 

egsea.data <- function(species = "human", simple=FALSE, returnInfo=FALSE){
    db.info = egseaData.info(species, simple)
    if (!returnInfo){
        cat(paste0("The following databases are available in EGSEAdata for ", 
                        normalizeSpecies(species), ":\n\n"))
        for (db in names(db.info)){
            cat(paste0("Database name: ", db.info[[db]]$info$name, "\n"))
            cat(paste0("Version: ", db.info[[db]]$info$version, "\n"))
            cat(paste0("Download/update date: ", db.info[[db]]$info$date, "\n"))
            cat(paste0("Data source: ", db.info[[db]]$info$source, "\n"))
            cat(paste0("Supported species: ", db.info[[db]]$info$species, "\n"))
            cat(paste0("Gene set collections: ", paste(db.info[[db]]$info$collections, 
                                    collapse=", "), "\n"))
            cat(paste0("Related data objects: ", paste(db.info[[db]]$info$data, 
                                    collapse=", "), "\n"))   
            if (!simple){
                cat("Number of gene sets in each collection for ", 
                        normalizeSpecies(species),": \n")
                for (col in names(db.info[[db]]$size)){
                    cat(paste0("\t", col, ": ", db.info[[db]]$size[[col]], "\n"))
                }
            }
            cat("\n")
        }
        cat("Type ?<data object name> to get a specific information 
				about it, e.g., ?kegg.pathways.")
    }
    else    
        return(db.info)
}

#TODO: If you update any database, make sure you update the fields of egseaData.info()
egseaData.info <- function(species = "human", simple=TRUE){
    db.info = list()
    # load preliminary data names/descriptions 
    fullToShort = list()
    fullToShort[["homo sapiens"]] = "human"
    fullToShort[["mus musculus"]] = "mouse"
    fullToShort[["rattus norvegicus"]] = "rat"
    
    msigdb.names = list(c1="c1 Positional Gene Sets", 
            c2="c2 Curated Gene Sets", 
            c3="c3 Motif Gene Sets", 
            c4="c4 Computational Gene Sets", c5="c5 GO Gene Sets", 
            c6="c6 Oncogenic Signatures", c7="c7 Immunologic Signatures",
            h="h Hallmark Signatures")
    
    gsetdb.names = list("gsdbdrug" = "GeneSetDB Drug/Chemical",   
            "gsdbdis" = "GeneSetDB Disease/Phenotype",
            "gsdbgo" = "GeneSetDB Gene Ontology",
            "gsdbpath" = "GeneSetDB Pathway",
            "gsdbreg" = "GeneSetDB Gene Regulation")
    # normalize species name
    normSpecies = normalizeSpecies(species)
    shortSpecies = fullToShort[[tolower(normSpecies)]]    
    # strore databases information and collection sizes
    # KEGG database
    keggdb.info = list(name = "KEGG Pathways", version = "NA", 
            date = "07 March 2017", source = "gage::kegg.gsets()",
            species= "human, mouse, rat", 
            collections= c("Signaling","Metabolism", "Disease"),
            data = "kegg.pathways")
    if (!simple){
        data(kegg.pathways)
        keggdb.sizes = list(Signaling=length(kegg.pathways[[shortSpecies]]$sig.idx),
                Metabolism=length(kegg.pathways[[shortSpecies]]$met.idx),
                Disease=length(kegg.pathways[[shortSpecies]]$dise.idx))
    }
    else
        keggdb.sizes = NULL
    db.info[["kegg"]] = list(info=keggdb.info, size=keggdb.sizes)
    # MSIGDB database     
    if (shortSpecies != "rat"){
        msigdb.info = list(name = "Molecular Signatures Database (MSigDB)", version = "5.2", 
                date = "07 March 2017", source = "http://software.broadinstitute.org/gsea",
                species= "human, mouse", 
                collections= c("h", "c1", "c2", "c3", "c4", "c5", 
                        "c6","c7"),
                data = c("msigdb", "Mm.H", "Mm.c2", "Mm.c3", "Mm.c4", "Mm.c5", 
                        "Mm.c6", "Mm.c7"))
        if (!simple){
            msigdb.sizes = list()
            if (shortSpecies == "mouse"){
                for (geneSet in c("h", "c2", "c3", "c4", "c5", 
                        "c6","c7")){
                    geneSet1 = ifelse(geneSet == "h", "H", geneSet)               
                    data(list=paste0("Mm.", geneSet1))                
                    gsets =get(paste0("Mm.", geneSet1))
                    msigdb.sizes[[msigdb.names[[geneSet]]]] = length(gsets)
                }
            }else{
                data(msigdb)
                gsc.all = msigdb  
                gsc.all = gsc.all[names(gsc.all) == "GENESET"]  
                organisms = sapply(gsc.all, function(x) x["ORGANISM"])
                gsc.all = gsc.all[organisms == "Homo sapiens"]
                types = tolower(sapply(gsc.all, function(x) x["CATEGORY_CODE"]))         
                for (geneSet in c("h", "c1", "c2", "c3", "c4", "c5", 
                        "c6","c7")){
                    gsets = gsc.all[types == geneSet]
                    msigdb.sizes[[msigdb.names[[geneSet]]]] = length(gsets)
                }
            }
        }
        else
            msigdb.sizes = NULL
        db.info[["msigdb"]] = list(info=msigdb.info, size=msigdb.sizes)
    }else{
        cat("WARNING: MSigDB gene sets do not support the species Rattus Norvegicus\n")
    }
    # GeneSetDB database
    gsetdb.info = list(name = "GeneSetDB Database", version = "NA", 
            date = "15 January 2016", source = "http://www.genesetdb.auckland.ac.nz/",
            species= "human, mouse, rat", 
            collections= c("gsdbdis", "gsdbgo", "gsdbdrug", "gsdbpath", "gsdbreg"),
            data = c("gsetdb.human", "gsetdb.mouse", "gsetdb.rat"))    
    if (!simple){
        data(list=paste0('gsetdb.', shortSpecies))
        gsetdb.all = get(paste0('gsetdb.', shortSpecies))
        categories = unique(as.character(gsetdb.all$anno$Category))
        genesetdb.gs.labels = loadGeneSetDBCategoryLabels()
        gsetdb.sizes = list()
        for (cat in categories){
            label = genesetdb.gs.labels[[cat]]
            gsets = gsetdb.all$anno[
                    gsetdb.all$anno[,"Category"] == cat,"GeneSet"]
            gsetdb.sizes[[gsetdb.names[[label]]]] = length(gsets)        
        }
    }else
        gsetdb.sizes = NULL
    db.info[["gsetdb"]] = list(info=gsetdb.info, size=gsetdb.sizes)
    
    return(db.info) 
}


loadGeneSetDBCategoryLabels <- function(){
    genesetdb.gs.labels = list("GeneSetDB Drug/Chemical"="gsdbdrug",   
            "GeneSetDB Disease/Phenotype" = "gsdbdis",
            "GeneSetDB Gene Ontology" = "gsdbgo",
            "GeneSetDB Pathway" = "gsdbpath",
            "GeneSetDB Gene Regulation" = "gsdbreg")
    return(genesetdb.gs.labels)
}

normalizeSpecies <- function(species){
    human.names = c("human", "homo sapiens", "hs")
    mouse.names = c("mouse", "mus musculus", "mm")
    rat.names = c("rat", "rattus norvegicus" , "rn")
    species = tolower(species)
    if (species %in% human.names){
        species = "Homo sapiens"
    }else if (species %in% mouse.names){
        species = "Mus musculus"
    }
    else if (species %in% rat.names){
        species = "Rattus norvegicus"
    }
    else{
        stop("Unrecognized species.")
    }
    return(species)
}


