#' Heatmap scroll function
#'
#' This function allows you to draw heatmap that let you scroll down with the top most highly variable genes using ComplexHeatmap package.
#' @param property Property or the variable of interest
#' @param eset 
#' @param precomp_dir 
#' @param datasetname 
#' @param clustering_methods Which method to use for clustering.
#' @export
#' @examples
#' heatmap_fit(property="ER",eset, precomp_dir,datasetname="EDS-1013", clustering_methods="Clustering by groups")

######### Heatmap scroll function #########

	heatmap_scroll =function (property,eset, precomp_dir,datasetname, clustering_methods) {
	
		if(!is.null(property)) {

			if (dir.exists(paste(precomp_dir,datasetname,"/","gimmOutCol/",sep=""))) {

				load(paste(precomp_dir,datasetname,"/","gimmOutCol/",datasetname,"_geneClustTopGenes1000_gimmOutCol.rda",sep=""))
				load(paste(precomp_dir,datasetname,"/","gimmOutRow/",datasetname,"_geneClustTopGenes1000_gimmOutRow.rda",sep=""))
				toCluster_R <-gimmOutRow$clustData

				forHeatmap<-(data.matrix(toCluster_R[,-(1:2)]))
				forHeatmap[forHeatmap>2]<-2
				forHeatmap[forHeatmap<(-2)]<-(-2)
				df <- pData(eset)[,property,drop=FALSE]
				df[is.na(df)] <- "NA"
				df[] <- lapply( df, factor)
				colnames(df) <- property
				
				if(abs(12-(ncol(forHeatmap)/30))<2) {
					col_fontsize <- abs(12-(ncol(forHeatmap)/30)) +2
				} else {
					col_fontsize <- abs(12-(ncol(forHeatmap)/30))
				}
				
				load("allColors.rda")
				colr= vector(mode="list", length=ncol(df))
				names(colr)= colnames(df)
				for(i in seq_len(ncol(df))) {
					if (length(unique(as.character(df[,i])))>105) {
						colr[[i]]=setNames(colorRampPalette(brewer.pal(11,"Spectral"))(length(unique(as.character(df[,i])))), unique(as.character(df[,i])))
					} else {
						colr[[i]]=setNames(allColors[seq_along(unique(as.character(df[,i])))], unique(as.character(df[,i])))
					}
				}
				
				if (clustering_methods=="Clustering by groups") {
					df <- df[do.call(order, c(data.frame(df[,1:ncol(df)]))),,drop=FALSE]
					forHeatmap <- forHeatmap[,rownames(df)]
					
					ha <- HeatmapAnnotation(df, which="column", width = unit(1,"mm"), col=colr,
					annotation_legend_param=list(title_gp = gpar(fontsize = 14),ncol=1))
					
					ht <- Heatmap(forHeatmap,  name = "Expression", col = colorRamp2(c(-2,0, 2), c("blue", "black","yellow")),
					cluster_columns = FALSE, cluster_rows = as.dendrogram(gimmOutRow$hGClustData), 
					show_row_names = TRUE, row_names_side = "right", row_names_gp = gpar(fontsize = 10),
					row_names_max_width = unit(8, "cm"), show_column_names = if(length(colnames(forHeatmap))<100) {TRUE} else {FALSE},column_names_max_height=unit((4/10)*max(nchar(colnames(forHeatmap))), "cm"),
					column_names_gp= gpar(fontsize = col_fontsize), row_dend_reorder=FALSE, top_annotation = ha,
					top_annotation_height= unit(0.7*ncol(df), "cm"),
					heatmap_legend_param = list(color_bar = "continuous",title_gp = gpar(fontsize = 13), legend_direction = "horizontal",nrow=1,
					legend_width = unit(5, "cm"), title_position = "topcenter"))
					draw(ht, heatmap_legend_side = "top", annotation_legend_side = "right")
					
					for(an in colnames(df)) {
						decorate_annotation(an, {
							# annotation names on the right
							grid.text(an, unit(1, "npc") + unit(.25, "cm"), 0.5, default.units = "npc", just = "left", gp=gpar(fontsize = 14))
						})
					}
				} else {
					ha <- HeatmapAnnotation(df, which="column", width = unit(1,"mm"), col=colr,
					annotation_legend_param=list(title_gp = gpar(fontsize = 14),ncol=1))
					
					ht <- Heatmap(forHeatmap,  name = "Expression", col = colorRamp2(c(-2,0, 2), c("blue", "black","yellow")),
					cluster_columns = as.dendrogram(gimmOutCol$hGClustData), cluster_rows = as.dendrogram(gimmOutRow$hGClustData), 
					show_row_names = TRUE, row_names_side = "right", row_names_gp = gpar(fontsize = 10),
					row_names_max_width = unit(8, "cm"), show_column_names = if(length(colnames(forHeatmap))<100) {TRUE} else {FALSE},column_names_max_height=unit((4/10)*max(nchar(colnames(forHeatmap))), "cm"),
					column_names_gp= gpar(fontsize = col_fontsize), row_dend_reorder=FALSE, top_annotation = ha,
					top_annotation_height= unit(0.7*ncol(df), "cm"),
					heatmap_legend_param = list(color_bar = "continuous",title_gp = gpar(fontsize = 13), legend_direction = "horizontal",nrow=1,
					legend_width = unit(5, "cm"), title_position = "topcenter"))
					draw(ht, heatmap_legend_side = "top", annotation_legend_side = "right")
					
					for(an in colnames(df)) {
						decorate_annotation(an, {
							# annotation names on the right
							grid.text(an, unit(1, "npc") + unit(.25, "cm"), 0.5, default.units = "npc", just = "left", gp=gpar(fontsize = 14))
						})
					}
				}
				
			} else {
			
				exps= as.matrix(exprs(eset)) - rowMeans(as.matrix(exprs(eset)))
				medAbsDev<-apply(exps,1,function(x) median(abs(x)))

				topGenes= function(exps, medAbsDev) {
					if (dim(exps)[1]>= 1000 & dim(exps)[1]<=1500) {
						topGenes<-order(medAbsDev,decreasing=T)[1:1500]
					} else if (dim(exps)[1]> 1500 ) {
						topGenes<-order(medAbsDev,decreasing=T)[1:1000]
					} else if (dim(exps)[1]< 1000) {
						topGenes<-order(medAbsDev,decreasing=T)
					} else {
						topGenes=medAbsDev
					}
					return(topGenes)
				}
				topGenes=topGenes(exps, medAbsDev)
				topGenes= topGenes[!is.na(topGenes)]
				
				forHeatmap<- data.matrix(exps[topGenes,])
				forHeatmap[forHeatmap>2]<-2
				forHeatmap[forHeatmap<(-2)]<-(-2)
				df <- pData(eset)[,property,drop=FALSE]
				df[is.na(df)] <- "NA"
				df[] <- lapply( df, factor)
				colnames(df) <- property

				if(abs(12-(ncol(forHeatmap)/30))<2) {
					col_fontsize <- abs(12-(ncol(forHeatmap)/30)) +2
				} else {
					col_fontsize <- abs(12-(ncol(forHeatmap)/30))
				}

				load("allColors.rda")
				colr= vector(mode="list", length=ncol(df))
				names(colr)= colnames(df)
				for(i in seq_len(ncol(df))) {
					if (length(unique(as.character(df[,i])))>105) {
						colr[[i]]=setNames(colorRampPalette(brewer.pal(11,"Spectral"))(length(unique(as.character(df[,i])))), unique(as.character(df[,i])))
					} else {
						colr[[i]]=setNames(allColors[seq_along(unique(as.character(df[,i])))], unique(as.character(df[,i])))
					}
				}
				
				if (clustering_methods=="Clustering by groups") {
					df <- df[do.call(order, c(data.frame(df[,1:ncol(df)]))),,drop=FALSE]
					forHeatmap <- forHeatmap[,rownames(df)]

					ha <- HeatmapAnnotation(df, which="column", width = unit(1,"mm"), col=colr,
					annotation_legend_param=list(title_gp = gpar(fontsize = 14),ncol=1))

					ht <- Heatmap(forHeatmap,  name = "Expression", col = colorRamp2(c(-2,0, 2), c("blue", "black","yellow")),
					cluster_columns = FALSE, cluster_rows = TRUE, clustering_distance_rows = "pearson", 
					show_row_names = TRUE, row_names_side = "right", row_names_gp = gpar(fontsize = 10),
					row_names_max_width = unit(8, "cm"), show_column_names = if(length(colnames(forHeatmap))<100) {TRUE} else {FALSE},column_names_max_height=unit((4/10)*max(nchar(colnames(forHeatmap))), "cm"),
					column_names_gp= gpar(fontsize = col_fontsize), row_dend_reorder=FALSE, top_annotation = ha,
					top_annotation_height= unit(0.7*ncol(df), "cm"),
					heatmap_legend_param = list(color_bar = "continuous",title_gp = gpar(fontsize = 13), legend_direction = "horizontal",nrow=1,
					legend_width = unit(5, "cm"), title_position = "topcenter"))
					draw(ht, heatmap_legend_side = "top", annotation_legend_side = "right")
					
					for(an in colnames(df)) {
						decorate_annotation(an, {
							# annotation names on the right
							grid.text(an, unit(1, "npc") + unit(.25, "cm"), 0.5, default.units = "npc", just = "left", gp=gpar(fontsize = 14))
						})
					}
				} else {
					ha <- HeatmapAnnotation(df, which="column", width = unit(1,"mm"), col=colr,
					annotation_legend_param=list(title_gp = gpar(fontsize = 14),ncol=1))

					ht <- Heatmap(forHeatmap,  name = "Expression", col = colorRamp2(c(-2,0, 2), c("blue", "black","yellow")),
					cluster_columns = TRUE, clustering_distance_columns = "pearson", cluster_rows = TRUE, clustering_distance_rows = "pearson",
					show_row_names = TRUE, row_names_side = "right", row_names_gp = gpar(fontsize = 10),
					row_names_max_width = unit(8, "cm"), show_column_names = if(length(colnames(forHeatmap))<100) {TRUE} else {FALSE},column_names_max_height=unit((4/10)*max(nchar(colnames(forHeatmap))), "cm"),
					column_names_gp= gpar(fontsize = col_fontsize), row_dend_reorder=FALSE, top_annotation = ha,
					top_annotation_height= unit(0.7*ncol(df), "cm"),
					heatmap_legend_param = list(color_bar = "continuous",title_gp = gpar(fontsize = 13), legend_direction = "horizontal",nrow=1,
					legend_width = unit(5, "cm"), title_position = "topcenter"))
					draw(ht, heatmap_legend_side = "top", annotation_legend_side = "right")
					
					for(an in colnames(df)) {
						decorate_annotation(an, {
							# annotation names on the right
							grid.text(an, unit(1, "npc") + unit(.25, "cm"), 0.5, default.units = "npc", just = "left", gp=gpar(fontsize = 14))
						})
					}
				}
			}
		}
		else {
			warning("select a property first")
		}			

	}
	
