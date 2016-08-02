# Weighted correlation network analysis

> MBL Neurobiology 2016
> 
> Joon An, Sanders Lab, UCSF (<JoonYong.An@ucsf.edu>)


![Willsey 2013](https://cbs.asu.edu/sites/default/files/images/Octopus_MBLTeal.jpg)*Image credit from MBL*

Co-exression network is a systems biology method to describe pattern(s) among genes across samples with multi-dimension information. This tutorial will use a R package, called Weighted correlation network analysis (WGCNA), for co-expression network analysis.

## Step0. Set-up the current run 

To do WGCNA, you first need to install WGCNA on your machine. In the R Studio console, please type the following command. The installation of WGCNA requires several R packages and dependency. Also, you need to install some packages for this script. This step will take some time (and often ERRORS).

~~~R
source("http://bioconductor.org/biocLite.R") 

biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 

install.packages("WGCNA")
~~~

Now load packages.

~~~R
library(WGCNA)
library(igraph) ## Plot networks
#library(GEOquery) ## Automatically get GEO data
#library(biomaRt) ## Annotations from biomaRt (host is set to use Gencode v10)
~~~

Once you finished the installation, then you can set up the current run. First, specifiy your working directory, where **input and output data** will be located. Also, set some working enrivorment parameters.

~~~R
options(stringsAsFactors=FALSE)
# Number of processes for this script. It depends on your computing power.
allowWGCNAThreads(nThreads = 2) 
~~~

### Download tutorial files
The tutorial data is located in `Bioinformatics/WGCNA` on the MBL Dropbox. Please download following files in your working directory:


* `20160608_kallisto_tpm.txt`
* `sleuth_sample_example.txt`
* `20160605_kallisto_tpm.geneData.txt`
* `data_expression.txt`



## Step1. Load your data 

WGCNA requires a matrix format (aka table). In this tutorial, you will use three matrix data. 

### Load expression profiles
First, load an expression profile. This data contains an expression profile, obtained from your RNAseq. After loading this matrix, check the dimension of your data. You can also see the first few lines of the data using the `head` function.

~~~R
datExpr <- read.csv('20160608_kallisto_tpm.txt',header=T,sep='\t') 
dim(datExpr)
head(datExpr)
~~~

### Load a meta data
Then, load a meta data. This data contains external information regarding samples or experimental conditions - mutant or wile-type, developmental stages, tissue types, etc. Again, this is a matrix format.

~~~R
datMeta <- read.csv('sampleInfo.txt',header=T,sep='\t') 
head(datMeta)
~~~

### Load a gene data
Lastly, load a gene data. This data contains information about genes. Since the expression dataset is based on transcripts, we will collapse the expression data into the gene level.

~~~R
# Check the number of transcripts and genes

length(geneData$target_id); length(unique(geneData$gene_name))

# Collapse expression values among transcripts from the same gene.

rownames(datExpr) = geneData$target_id

collapseDat <- collapseRows(datET = datExpr,method="MaxMean",
                            rowGroup = geneDat$gene_name,rowID = rownames(datExpr))
~~~

To do faster tutorial, we will use the pre-made dataset. Please load the data.

~~~R
geneData <- read.csv('20160605_kallisto_tpm.geneData.txt',header=T,sep='\t')

head(geneData)

~~~

## Step2. Explore Topology of the Network

A biological system is a kind of scale-free network, where interactions are not randomly distributed across genes, proteins or any form of interactors, and thus exibits hierachical and structural characteristics. 

Gene-expression can be used to describe this using the network concept, in which nodes are genes and edges between nodes exist when gene expression profiles are significantly correlated (coexpressed) across all samples. 

Okay. I will stop writing the book. Perhaps you might want to read some interesting literatures: 

[Network biology: understanding the cell's functional organization, Barabási and Oltval, 2004](http://www.ncbi.nlm.nih.gov/pubmed/14735121)

[Evolution of synapse complexity and diversity. Emes and Grant, 2012](http://www.ncbi.nlm.nih.gov/pubmed/22715880)

[Systems biology and gene networks in neurodevelopmental and neurodegenerative disorders, Parikshak et al., 2015](http://www.ncbi.nlm.nih.gov/pubmed/26149713)


### Find a soft-thresholding power for network construction

This step will find scale free topology of the network given a range of multiple soft thresholding powers. The aim of this step is to help you choose an appropriate power for network construction.

~~~R
# set the powers to examine network topoloy
powers <- seq(1,10,by=1)

sft <- pickSoftThreshold(t(datExpr),
                         powerVector=powers,
                         corFnc="bicor",networkType="signed")


# This will return you a table like below (numbers may not match):
# Power SFT.R.sq slope truncated.R.sq mean.k. median.k. max.k.
# 1    0.635 12.400          0.869 12000.0  1.28e+04  15600
# 2    0.731 13.000          0.651  3630.0  3.74e+03   6140
# 3    0.461  7.000          0.862  1270.0  1.21e+03   2790
# 4    0.130 -0.221          0.691   526.0  4.33e+02   1910
# 5    0.821 -1.580          0.773   266.0  1.70e+02   1580
# 6    0.828 -1.510          0.839   161.0  7.41e+01   1410
# 7    0.816 -1.380          0.876   112.0  3.51e+01   1310
# 8    0.804 -1.280          0.903    85.9  1.82e+01   1240
# 9    0.801 -1.200          0.918    70.3  1.00e+01   1180
# 10    0.792 -1.110          0.946    60.1  5.90e+00   1140
~~~

Now you plot your result and pick up the power. 

~~~R
par(mfrow = c(1,2)) # Set the plot screen 
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
	xlab="Soft Threshold (power)", 
	ylab="Scale Free Topology Model Fit, signed R^2",type="n",
	main = paste("Scale independence"))

text(sft$fitIndices[,1], 
	-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	labels=powers,cex=cex1,col="red");

# Mean connectivity as a function of the soft-thresholding power

plot(sft$fitIndices[,1], sft$fitIndices[,5],
	xlab="Soft Threshold (power)",
	ylab="Mean Connectivity", type="n", 
	main = paste("Mean connectivity"))

text(sft$fitIndices[,1], sft$fitIndices[,5], 
	labels=powers, cex=cex1,col="red")
~~~

Choosing the optimal power is not easy but the general suggestion is to choose a number you see the saturation of the graph. You can read the developer's suggestions. Please have the 6th question: [WGCNA FAQ](https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/faq.html)  

Now you can set the power for your network construction

~~~R
power <- 9
~~~

## Step3. Construct a co-expression network

#### Run WGCNA using the semi-automated function

~~~R
net <- blockwiseModules(t(datExpr),
		power=power,
		deepSplit=2,
		minModuleSize=200,
		minKMEtoStay=0,
		mergeCutHeight=0.50,
		detectCutHeight=0.99995,
		impute=TRUE,
		corType="bicor",
		networkType="signed",
		pamStage=FALSE,
		verbose=3,saveTOMs=FALSE,
		maxBlockSize=1000) 
~~~

You can see the details about parameters. Type this: 

~~~R
?blockwiseModules
~~~

Now you can check the number of modules in your network and the number of module genes. 

~~~R
table(net$colors)

# The result will be like:
# black      blue     brown     green      grey       red turquoise    yellow 
# 415      1914       887       749      9519       661      2494       782 
~~~

Module colours do not have any meaning, except a visual/labelling purpose. 

Note that **GREY IS RESERVED to color genes that are not part of any module**. Do you have the large number of Grey module genes? What would be a reason for this?

## Step4. Visualize the relationship of modules 

### Get a adjacency matrix for the hierarchical clustering function

This step creates an adjacency matrix from an expression profile. This matrix will be used to compute network topology (Topological Overlap Matrix) ([Zhang et al., 2005](https://labs.genetics.ucla.edu/horvath/GeneralFramework/WeightedNetwork2005.pdf))

~~~R
adjacency = adjacency(t(datExpr), power = power)
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average")

# Plot
plotDendroAndColors(geneTree, 
						net$colors, 
						"Dynamic Tree Cut",
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05,
                    main = "Gene Dendrogram and Module 						Colors")
~~~

You can also plot a multidimensional scaling plot can show similar information about module relationships in two dimensions. This plots the first two principal components of the distance matrix.

~~~R
par(mfrow = c(1,1),mar=c(5,4,2,2))

MEs = net$MEs #the module eigengenes

# Choose module colors
moduleColors_idx = unique(net$colors) 

colnames(MEs) = moduleColors_idx 

MDS = cmdscale(as.dist(distance),2) 

plot(MDS, 
	col=moduleColors_idx, 
	xlab='1st Principle Component', 
	ylab='2nd Principle Component',
	xlim=c(-0.5,0.5),ylim=c(-0.4,0.4))
	
text(MDS[,1], MDS[,2], 
	moduleColors_idx, 
	pos= 3, cex=0.8)
~~~

## Step 5. Characterize modules by external information 

Now, you will examine which external information is related to network modules.

~~~R
# List your external information
colnames(datMeta)

# Choose the genotype for your analysis
# Create a categorical value and add to a new column
datMeta$mut <- ifelse(datMeta$genotype!="ctr_wt",'red', 'darkblue')

# Plot
barplot(MEs$green,col=datMeta$mut,
        width=0.1, xlab='Samples', ylab='Eigen Value', 
        names.arg=datMeta$sample,
        las=2,beside = TRUE,
        main='Green Module') 

legend("topright",inset=c(0.10,0), 
       fill=c("darkblue","red"), 
       legend=c("Controls","Mutants"))
~~~


Here you will test association of network modules with the expression data of head tissues. Create a new column containing categorical values for your external information.

~~~R
# Again, set your comparison set from external information

datMeta$exp <- ifelse(datMeta$genotype!="ctr_wt",1, 0)

# Calculate correlation between module eigen values and external information

# net$MEs is the column contains module eigen values across samples

head(net$MEs)

moduleTraitCor = cor(net$MEs, datMeta$exp, use = "p")

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, 96)

textMatrix = paste(signif(moduleTraitCor, 2),
				"(", signif(moduleTraitPvalue, 1), 
				")", sep = " ")

dim(textMatrix) = dim(moduleTraitCor)

par(mar = c(3, 5, 3, 3))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = " ", 
               yLabels = names(MEs),
               ySymbols = names(MEs), 
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.9,
               zlim = c(-1,1),
               main = paste("Module Association with Motor-neuron sample "))
~~~

## Step6. Relationship of top module genes

Here you will do gene-level visualization. One example is to look at the relationship of the top genes with the highest eigen-values of the module. 

~~~R
# Retrieve module colours
table(net$colors)

kMEdat <- signedKME(t(datExpr), MEs, corFnc="cor")

colnames(kMEdat) <- colnames(MEs)

colnames(kMEdat) # module colours

mod_col = 'turquoise'

mod_kMEdat <- kMEdat[net$colors==mod_col,]

# Select the top 100 genes
topGenes=rank(mod_kMEdat[,mod_col],
				ties.method="first")<=100

gene.names= rownames(mod_kMEdat)[topGenes]

# Subset the adjacency matrix
adj <- adjacency[net$colors==mod_col,
					net$colors==mod_col]

adj <- adj[topGenes,topGenes]

quantile(adj)

# This is optional. I discard few edges with a low correlation value by replacing with 0 value, This will highlight few nodes with high connectivity. If you want to use this option, please remove the hashtag(#) from the sentence below.
# adj[adj<0.001] <- 0

# Define a graph object for a network plot
g1 <- graph.adjacency(as.matrix(adj),
						mode="undirected",
						weighted=TRUE,
						diag=FALSE) 

# Create an object for plotting graph
par(mar = c(1, 1, 1, 1))
plot.igraph(g1,
            vertex.size=10, 
            vertex.label.dist=0.5,
            vertex.label.font=0.2,
            vertex.color=mod_col,
            edge.color="grey",
            edge.width=E(g1)$weight*100)
~~~

## Step7. Annotate modules genes with Drosophila Gene Ontology (Web)

Here you will annotate your module genes by pathway database. Please go to the [Gene Ontology Website](http://geneontology.org/page/go-enrichment-analysis) and copy-and-paste your genes into the box (left) for enrichment text. Alternatively, you can use a R package for this (GO.db). 

~~~R
# Save your gene list into "topGenes.txt"
writeLines(gene.names, 'topGenes.txt',sep='\n')

# Save your module genes and do Gene ontology annotation to see which genes are represented in your module.

# Turquoise
geneturquoise <- rownames(datExpr)[net$colors=="turquoise"]

writeLines(geneturquoise,'modGenes_turquoise.txt',sep='\n')

# Red Module
geneRed <- rownames(datExpr)[net$colors=="red"]

writeLines(geneRed, 'modGenes_Red.txt',sep='\n')
~~~

## What you might want to do next

* Gene set enrichment with the gene list of your interest (e.g. synatic genes)