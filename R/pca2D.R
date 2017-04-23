#' PCA2D function
#'
#' This function allows you to draw 2D PCA plot.
#' @param eset 
#' @param logTransform If the data should be log transformed? Default is FALSE.
#' @param property Property or the variable of interest
#' @export
#' @examples
#' PCA2D(eset,logTransform=FALSE, property="ER")

######### PCA functions #########

	PCA2D<- function (eset,logTransform=FALSE, property) {

		exprs(eset)=exprs(eset)[complete.cases(exprs(eset)),]
		if (logTransform) {
			exps=log2(exprs(eset)+1)
		} 
		else {
			exps= as.matrix(exprs(eset))
		}
		
		pseudoCount= as.matrix(exps)
			  
		b.PCA=prcomp(t(pseudoCount),retx=TRUE,center=TRUE)
		pca.mat= cbind(b.PCA$x[,1],b.PCA$x[,2],b.PCA$x[,3],b.PCA$x[,4],b.PCA$x[,5])
		groups <- pData(eset)[,which(colnames(pData(eset))==property)] #input$property=='subtype'
		colcols= colorRampPalette(brewer.pal(11,"Spectral"))(length(unique(groups)))

		pairsD3(pca.mat, group = groups, labels =c("PC1","PC2","PC3","PC4","PC5"), cex = 3,
		width = 1200, col = colcols[factor(as.character(groups),exclude=NULL)], theme = "colour",opacity = 1)

	}
	
