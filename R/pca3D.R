#' PCA3D function
#'
#' This function allows you to draw 3D PCA plot.
#' @param eset 
#' @param logTransform If the data should be log transformed? Default is FALSE.
#' @param property Property or the variable of interest
#' @export
#' @examples
#' PCA3D(eset,logTransform=FALSE, property="ER")

	PCA3D<- function (eset,logTransform=FALSE, property) {

		exprs(eset)=exprs(eset)[complete.cases(exprs(eset)),]
		if (logTransform) {
			exps=log2(exprs(eset)+1)
		} 
		else {
			exps= as.matrix(exprs(eset))
		}

		pseudoCount= as.matrix(exps)      
		b.PCA=prcomp(t(pseudoCount),retx=TRUE,center=TRUE)
		groups <- pData(eset)[,which(colnames(pData(eset))==property)] #input$property=='subtype'
		colcols= colorRampPalette(brewer.pal(11,"Spectral"))(length(unique(groups)))

		par3d(windowRect = c(10, 10, 800, 800),cex=1.2)
		plot3d(b.PCA$x[,1],b.PCA$x[,2],b.PCA$x[,3],type="s",col = colcols[factor(as.character(groups),exclude=NULL)],  xlab ="PC1", ylab = "PC2 ", size=1.5,
		zlab = "PC3 ", asp=.5, main=paste("3D PCA plot of the ", property, sep=""))
		rglwidget()

	}

