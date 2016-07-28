#
## Author: monther
################################################################################
#
### Reads GeneSetDB genesets ###############
##data.dir = "~/Documents/Workspace/Pathway_analysis/EGSEA/EGSEA/data/"
##root.path = "~/Documents/Workspace/Pathway_analysis/EGSEA/"
##setwd(root.path)
##gsdb.human = updateGeneSetDB(gmt.file.name = "download-gmt_h", species='human', data.dir= data.dir)
##gsdb.mouse = updateGeneSetDB(gmt.file.name = "download-gmt_m", species='mouse', data.dir= data.dir)
##gsdb.rat = updateGeneSetDB(gmt.file.name = "download-gmt_r", species='rat', data.dir= data.dir)
#
### KEGG
##updateKEGGpathways(data.dir = data.dir)
#
#updateGeneSetDB <- function(gmt.file.name, species="human", data.dir = "./"){
#	fc <- file(gmt.file.name)
#	gsets.raw <- strsplit(readLines(fc), "\t")
#	close(fc)
#	gset.category = list()
#	gset.category[["Biocarta"]] = c("GeneSetDB Pathway", "Biocarta", "http://www.biocarta.com/")
#	gset.category[["EHMN"]] = c("GeneSetDB Pathway", "Edinburgh Human Metabolic Network", "http://www.ehmn.bioinformatics.ed.ac.uk/")
#	gset.category[["HumanCyc"]] = c("GeneSetDB Pathway", "Encyclopedia of Human Genes and Metabolism", "http://humancyc.org/")
#	gset.category[["INOH"]] = c("GeneSetDB Pathway", "Integrating Network Objects with Hierarchies", "http://www.inoh.org/")
#	gset.category[["NetPath"]] = c("GeneSetDB Pathway", "NetPath", "http://www.netpath.org/")
#	gset.category[["PID"]] = c("GeneSetDB Pathway", "Pathway Interaction Database", "http://pid.nci.nih.gov/")
#	gset.category[["Reactome"]] = c("GeneSetDB Pathway", "Reactome", "http://www.reactome.org/")
#	gset.category[["WikiPathways"]] = c("GeneSetDB Pathway", "WikiPathways", "http://www.wikipathways.org/")
#
#	gset.category[["CancerGenes"]] = c("GeneSetDB Disease/Phenotype", "Cancer Genes", "http://cbio.mskcc.org/CancerGenes/")
#	gset.category[["KEGG(Disease)"]] = c("GeneSetDB Disease/Phenotype", "KEGG Disease", "http://www.genome.jp/kegg/disease/")
#	gset.category[["HPO"]] = c("GeneSetDB Disease/Phenotype", "Human Phenotype Ontology", "http://www.human-phenotype-ontology.org/")
#	gset.category[["MethCancerDB"]] = c("GeneSetDB Disease/Phenotype", "Cancer type and Methylation status", "http://www.methcancerdb.net/")
#	gset.category[["MethyCancer"]] = c("GeneSetDB Disease/Phenotype", "Cancer type", "http://methycancer.psych.ac.cn/")
#	gset.category[["MPO"]] = c("GeneSetDB Disease/Phenotype", "Mammalian Phenotype Ontology", "http://www.informatics.jax.org/searches/MP_form.shtml")
#	gset.category[["SIDER"]] = c("GeneSetDB Disease/Phenotype", "SIDe Effect Resource", "http://sideeffects.embl.de/")
#	
#	gset.category[["CTD"]] = c("GeneSetDB Drug/Chemical", "Comparative Toxicogenomics Database", "http://ctd.mdibl.org/")
#	gset.category[["DrugBank"]] = c("GeneSetDB Drug/Chemical", "DrugBank", "http://www.drugbank.ca/")
#	gset.category[["MATADOR"]] = c("GeneSetDB Drug/Chemical", "Manually Annotated Targets and Drugs Online Resource", "http://matador.embl.de/")
#	gset.category[["SMPDB"]] = c("GeneSetDB Drug/Chemical", "Small Molecular Pathway DataBase", "http://www.smpdb.ca/")
#	gset.category[["STITCH"]] = c("GeneSetDB Drug/Chemical", "Search Tool for Interactions of Chemicals", "http://stitch.embl.de/")
#	gset.category[["T3DB"]] = c("GeneSetDB Drug/Chemical", "Toxin and Toxin Target Database", "http://www.t3db.org/")
#	
#	gset.category[["MicroCosm Targets"]] = c("GeneSetDB Gene Regulation", "MicroCosm Targets", "http://www.ebi.ac.uk/enright-srv/microcosm/")
#	gset.category[["miRTarBase"]] = c("GeneSetDB Gene Regulation", "miRTarBase", "http://mirtarbase.mbc.nctu.edu.tw/")
#	gset.category[["TFactS"]] = c("GeneSetDB Gene Regulation", "TFactS", "http://www.tfacts.org/")
#	gset.category[["Rel/NF-kappaB target genes"]] = c("GeneSetDB Gene Regulation", "Rel/NF-kappaB target genes", "http://bioinfo.lifl.fr/NF-KB/")
#	
#	gset.category[["GO_BP"]] = c("GeneSetDB Gene Ontology", "GO Biological Processes", "http://www.geneontology.org/")
#	gset.category[["GO_MF"]] = c("GeneSetDB Gene Ontology", "GO Molecular Functions", "http://www.geneontology.org/")
#	gset.category[["GO_CC"]] = c("GeneSetDB Gene Ontology", "GO Cellular Components", "http://www.geneontology.org/")
#
#
#
#	gsetdb.all = list()
#	gsetdb.all$species = species	
#	gsetdb.all$category = gset.category
#	gsetdb.all$original = lapply(gsets.raw, function(x) x[-c(1,2)])
#	names(gsetdb.all$original) = sapply(gsets.raw, function(x) x[1])
#	sourceDB = sapply(gsets.raw, function(x) x[2])
#	gsetdb.all$anno = data.frame(ID=paste0("gsdb", seq(1, length(gsets.raw))),
#			GeneSet=names(gsetdb.all$original),
#			NumGenes= sapply(gsetdb.all$original, length),
#			Category=as.character(sapply(sourceDB, function(x) gset.category[[x]][1])),
#			SourceDB=sourceDB,
#			SourceURL=as.character(sapply(sourceDB, function(x) gset.category[[x]][3]))
#	)
#	if (!is.null(data.dir)){
#		save(gsetdb.all, file=paste0(data.dir, 'genesetdb.',species,'.rda'))
#	}
#	print(paste0('A file has been written to ',data.dir, 'genesetdb.',species,'.rda'))
#	return(gsetdb.all)
#}
#
#
#updateKEGGpathways <-function(data.dir="./"){
#	
#	kegg.pathways = list()
#	species.all = c("human", "mouse", "rat")
#	for (species in species.all){
#		kegg = kegg.gsets(species = species, id.type = "kegg")
#		kegg.pathways[[species]] = kegg
#	}
#	save(kegg.pathways, file=paste0(data.dir, 'kegg.pathways.rda'))
#	print(paste0("A file has been written to ", data.dir, 'kegg.pathways.rda' ))
#}
#
#updateMSigDBAnno <- function(xml.file, out.dir="./"){
#
#
#	msigdb.xml.file="../msigdb_v5.0.xml"
#	msigdb = xmlToList(xmlParse(msigdb.xml.file))
#	save(msigdb, file="msigdb.rda")
#}
#
#updateMsigDBMouseC5 <- function(out.dir="./"){	
#	
#	data(msigdb.mouse_c5_v4)
#	go.terms = sapply(names(Mm.c5), function(x) tolower(gsub("_", " ", x)))
#	xx = GOTERM	
#	t = sapply(1:length(xx), function(x) tolower(Term(xx[x])))
#	t = sapply(t, function(x) gsub("-", " ", x))
#	t.splits = sapply(t, function(x) strsplit(x, " ", fixed=TRUE)[1])
#	t.splits = sapply(t.splits, function(x) as.vector(sapply(x, function(y) sub(",$", "", y))))
#	
#	ontology = character(0)
#	goID = character(0)
#	for (i in 1:length(go.terms)){
#		temp = xx[go.terms[i] == t]
#		if (length(temp) == 1){
#			ontology = c(ontology, Ontology(temp))
#			goID = c(goID, GOID(temp))
#			next
#		}
#		temp = xx[grep(go.terms[i], t)]
#		if (length(temp) == 0){
#			tokens = strsplit(go.terms[i], " ", fixed=TRUE)[[1]]
#			temp = xx[sapply(t.splits, function(x) length(intersect(tokens, x)) / length(tokens) == 1)]
##			print(go.terms[i])
##			print(length(temp))			
#		}
#		if (length(temp) == 1){
#			ontology = c(ontology, Ontology(temp))
#			goID = c(goID, GOID(temp))
#			next
#		}
#		if (length(temp) > 1){				
#			max.overlapping = 0
#			best = NULL
#			genes.original = Mm.c5[[i]]
#			for (j in 1:length(temp)){
#				tmp = temp[j]
#				tryCatch({genes = get(GOID(tmp), org.Mm.egGO2ALLEGS)}, 
#						error=function(e){genes=c()})
#				overlapping = length(intersect(genes.original, genes))
#				if ( overlapping > max.overlapping){
#					max.overlapping = overlapping
#					best = tmp
#				}
#			}
#			if (!is.null(best)){				
#				ontology = c(ontology, Ontology(best))
#				goID = c(goID, GOID(best))	
##				print(Term(best))
#				next	
#			}
#		}		
#		ontology = c(ontology, "N/A")
#		goID = c(goID, "N/A")
#	}
#	c5.go.anno = cbind(Ontology=ontology, GOID=goID)
#	rownames(c5.go.anno) = names(Mm.c5)
#	file.path = paste0(out.dir, "msigdb.mouse_c5_v4_anno.rda")
#	save(c5.go.anno, file=file.path)
#}