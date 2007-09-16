
####################
# HarshExt
####################

HarshExt <- function(affy.object, my.ErrorImage = NULL, extended.radius = 10){

Harshlight(affy.object, my.ErrorImage, extended.radius, compact.quant.bright = 0, compact.quant.dark = 0, compact.size.limit = 15, compact.connect = 8, compact.pval = 0.01, diffuse.bright = 0, diffuse.dark = 0, diffuse.pval = 0.001, diffuse.connect = 8, diffuse.radius = 10, diffuse.size.limit = (3*314), percent.contiguity = 50, report.name = '', na.sub = FALSE, interpolate = TRUE, diffuse.close = TRUE)

}

####################
# HarshComp
####################

HarshComp <- function(affy.object, my.ErrorImage = NULL, extended.radius = 10, compact.quant.bright = 0.025, compact.quant.dark = 0.025, compact.size.limit = 15, compact.connect = 8, compact.pval = 0.01, percent.contiguity = 50, report.name = 'R.report.ps', na.sub = FALSE, interpolate = TRUE){

Harshlight(affy.object, my.ErrorImage, extended.radius, compact.quant.bright, compact.quant.dark, compact.size.limit, compact.connect, compact.pval, diffuse.bright = 0, diffuse.dark = 0, diffuse.pval = 0.001, diffuse.connect = 8, diffuse.radius = 10, diffuse.size.limit = (3*314), percent.contiguity, report.name, na.sub, interpolate, diffuse.close = TRUE)

}

####################
# Harshlight
####################

Harshlight <- function(affy.object, my.ErrorImage = NULL, extended.radius = 10, compact.quant.bright = 0.025, compact.quant.dark = 0.025, compact.size.limit = 15, compact.connect = 8, compact.pval = 0.01, diffuse.bright = 40, diffuse.dark = 35, diffuse.pval = 0.001, diffuse.connect = 8, diffuse.radius = 10, diffuse.size.limit = (3*3.14*(diffuse.radius**2)), percent.contiguity = 50, report.name = 'R.report.ps', na.sub = FALSE, interpolate = TRUE, diffuse.close = TRUE){

	####### Check: affy.object must be a valid affy object

	if(class(affy.object) != 'AffyBatch'){
		print("Error: affy.object must be a valid AffyBatch object")
		flush.console()
		return()
	}

	if(length(affy.object) < 2){
		print("Error: affy.object must have at least 2 chips")
		flush.console()
		return()
	}
	NROW<- nrow(affy.object)
	NCOL <- ncol(affy.object)
	NROW_new <- .GetBound(affy.object)$y2
	#print (NROW_new)
	#print (NROW)		

	###########################################################################

	compact.connect <- (compact.connect == 8) #1 = 8-neighbours connectivity, 0 = 4-neighbours connectivity
	diffuse.connect <- (diffuse.connect == 8) #1 = 8-neighbours connectivity, 0 = 4-neighbours connectivity

	if(.CheckParam(diffuse.dark = diffuse.dark, diffuse.bright = diffuse.bright, compact.quant.bright = compact.quant.bright, compact.quant.dark = compact.quant.dark, compact.size.limit = compact.size.limit, diffuse.size.limit = diffuse.size.limit, compact.connect = compact.connect, diffuse.connect = diffuse.connect, diffuse.radius = diffuse.radius, diffuse.pval = diffuse.pval, extended.radius = extended.radius, compact.pval = compact.pval, NROW = NROW, NCOL = NCOL, percent.contiguity = percent.contiguity, na.sub = na.sub, interpolate = interpolate, diffuse.close = diffuse.close))
		return()

	#######################################################

	int.chips <- intensity(affy.object)
	int.subs <- int.chips



	
	if(!is.null(my.ErrorImage)){
		if(class(my.ErrorImage) != 'matrix'){
			print("Error: my.ErrorImage must be a valid matrix object")
			flush.console()
			return()
		}
		ErrorInt <- my.ErrorImage
	} else {
		print("Generating Error Images")
		ErrorInt <- .ErrorIntensity(int.chips, transfo = log2)
		if(is.character(ErrorInt[1])){
			print("Couldn't generate Error Images: bailing out")
			flush.console()
			return()
		}
	}

	############################## SET CHIP SIZE ###############################

	#Set the number of rows and columns of the chip for the analysis
	print("Initializing Harshlight")
	trash <- .C("init_harshlight",
		as.integer(NROW),
		as.integer(NCOL),
		as.integer(length(affy.object)),
		as.integer(0),
		PACKAGE = "Harshlight")
	status <- trash[[4]]

	if(status){
		return()
	}

	############################################

	#Start creating the ps file
	#Header of the report
	if(report.name == ''){
		report.name <- 0
	}

	if(report.name != 0){
		chip.name <- sampleNames(affy.object)

		chip.header <- paste(chip.name, collapse = ", ")
		if(nchar(chip.header) > 40){
			chip.header <- paste(substr(chip.header, 0, 40), "...", sep = ", ")
		}

		trash <- .C("report_overall_header",
     			as.character(report.name),
			as.integer(extended.radius),
			as.double(compact.quant.bright),
			as.double(compact.quant.dark),
			as.integer(compact.size.limit),
			as.integer(compact.connect),
			as.double(compact.pval),
			as.double(diffuse.bright),
			as.double(diffuse.dark),
			as.double(diffuse.pval),
			as.integer(diffuse.connect),
			as.integer(diffuse.size.limit),
			as.integer(diffuse.radius),
			as.double(percent.contiguity),
			as.integer(na.sub),
			as.integer(interpolate),
			as.integer(diffuse.close),
			as.character(chip.header),
			as.integer(0),
			PACKAGE = "Harshlight")
		status <- trash[[19]]

		if(status){
			print("An error occurre while trying to write to file")
			print("The report file will not be created")
			flush.console()
			report.name <- 0
		}
	}

	########################### Choose the best simulations for the chip ##########################

	allowed.dim <- c(534, 640, 712)	#dimensions of the chips used for the simulations
	side <- (NROW + NCOL)/2
	min.size <- abs(side - allowed.dim)

	allowed.pval <- c(0.01, 0.02, 0.05, 0.10, 0.20, 0.25, 0.30, 0.40)

	for(i in 1:2){
		if(i == 1){
			min.pval <- abs(compact.quant.bright - allowed.pval)
		}else{
			min.pval <- abs(compact.quant.dark - allowed.pval)
		}

		best.dim.score <- c(0, 1, 2)
		best.dim <- best.dim.score[which(min.size == min(min.size))[1]]	#simulation closest in size to the chip in analysis

		if(!min(min.pval)){
	
			if(compact.connect){
				nbhood <- 4 + best.dim
			} else {
				nbhood <- 1 + best.dim
			}

			data(sim, package = "Harshlight")
			if(!exists("sim")){
				print("Couldn't find simulations")
				flush.console()
				trash <- .C("closepsfile", PACKAGE = "Harshlight")
				return()
			}

			sim.pval <- sim[,nbhood]
			#Chooses the best p.values for the given quantile
			if(i == 1){
				p.values.bright <- .RetrievePvalues_exact(compact.quant.bright, min.pval, sim.pval, length(ErrorInt[,1])/2)
			}else{
				p.values.dark <- .RetrievePvalues_exact(compact.quant.dark, min.pval, sim.pval, length(ErrorInt[,1]/2))
			}
			rm("sim", envir = environment(Harshlight))
	
		}else{
			if(interpolate == FALSE){
				if(i == 1){
					simulation.bright <- seq(length = 409600/2, 0, 0)
					trash <- .C("simulations", simulation.bright <- as.integer(simulation.bright), as.double(compact.quant.bright), as.integer(compact.connect), status <- as.integer(0), DUP = FALSE, PACKAGE = "Harshlight")
					if(status){
						print("	An error occurrred while running the simulations. The null distribution will be recovered through interpolation")
						flush.console()
						interpolate <- TRUE
					}
					p.values.bright <- simulation.bright/100000
				}else{
					if(compact.quant.dark != compact.quant.bright){
						simulation.dark <- seq(length = 409600/2, 0, 0)
						trash <- .C("simulations", simulation.dark <- as.integer(simulation.dark), as.double(compact.quant.dark), as.integer(compact.connect), status <- as.integer(0),  DUP = FALSE, PACKAGE = "Harshlight")
						if(status){
							print("	An error occurrred while running the simulations. The null distribution will be recovered through interpolation")
							flush.console()
							interpolate <- TRUE
						}
						p.values.dark <- simulation.dark/100000
					}else{
						p.values.dark <- p.values.bright
					}
				}
			}
			if(interpolate == TRUE){

				if(compact.connect){
					nbhood <- 7 + best.dim*2
				} else {
					nbhood <- 1 + best.dim*2
				}

				data(sim.int, package = "Harshlight")
				if(!exists("sim.int")){
					print("Couldn't find simulations")
					flush.console()
					trash <- .C("closepsfile", PACKAGE = "Harshlight")
					return()
				}

				#Chooses the best p.values for the given quantile
				if(i == 1){
					p.values.bright <- .RetrievePvalues(compact.quant.bright, sim.int[,nbhood], sim.int[,nbhood + 1], length(ErrorInt[,1])/2)
				}else{
					p.values.dark <- .RetrievePvalues(compact.quant.dark, sim.int[,nbhood], sim.int[,nbhood + 1], length(ErrorInt[,1]/2))
				}
				rm("sim.int", envir = environment(Harshlight))
			}

		}
	}

	for(i in 1:dim(ErrorInt)[2]){
		print(paste("Analyzing chip number", i, sep = " "))
		flush.console()

	if(report.name != 0){
		# Chip header
		if(nchar(chip.name[i]) > 10){
			chip.name[i] <- paste(substr(chip.name[i], 0, 10), "...", sep = "")
		}
		trash <- .C("chip_overall_header",
			as.integer(i),
			as.character(chip.name[i]),
			PACKAGE = "Harshlight")
	
		# Original image
		trash <- .C("chip_image",
			as.integer(50),
			as.integer(475),
			as.integer(.ScaleInt(log2(int.chips[,i]))),
			as.character("Original image"),
			as.integer(0),
			as.integer(0),
			as.integer(0),
			as.integer(0),
			PACKAGE = "Harshlight")
	
		trash <- .C("chip_image",
			as.integer(300),
			as.integer(475),
			as.integer(.ScaleInt(ErrorInt[,i])),
			as.character("Error image"),
			as.integer(0),
			as.integer(0),
			as.integer(0),
			as.integer(0),
			PACKAGE = "Harshlight")
		}


		result <- .AnalyzeChip(ErrorInt[,i], extended.radius = extended.radius, compact.quant.bright = compact.quant.bright, compact.quant.dark = compact.quant.dark, compact.size.limit = compact.size.limit, compact.connect = compact.connect, compact.pval = compact.pval, diffuse.bright = diffuse.bright, diffuse.dark = diffuse.dark, diffuse.pval = diffuse.pval, diffuse.connect = diffuse.connect, diffuse.radius = diffuse.radius, diffuse.size.limit = diffuse.size.limit, percent.contiguity = percent.contiguity, NROW = NROW, NCOL = NCOL, diffuse.close = diffuse.close, p.values.bright = p.values.bright, p.values.dark = p.values.dark, report.name = report.name)

		if(!na.sub){
			int.chips[,i][is.na(result)] <- apply(int.subs[is.na(result),-i], 1, median)
		}
		else {
			print("substitue na here")
			int.chips[,i][is.na(result)] <- NA
		}

	}

	###################### CLEAR THE MEMORY #################

	trash <- .C("closepsfile", PACKAGE = "Harshlight")
	trash <- .C("free_memory", PACKAGE = "Harshlight")

	print("Substituting values")
	flush.console()

	rm(ErrorInt)
	
	#from yupu: change to assign the matrix all together to solve the memory problem Jan 10,2007
	#for(i in 1:length(affy.object)){
	#	intensity(affy.object)[,i] <- int.chips[,i]
	#}
	intensity(affy.object) = int.chips
	rm(int.chips)

	return(affy.object)
}




.ErrorIntensity <- function(object,transfo=log2){
	# object	intensity of an affy object
	#d0 = date()
	if (is.function(transfo)) { object <- transfo(object) }

	ncol = ncol(object)
	nrow = nrow(object)



	for(i in 1: nrow){
		my_row <- object[i,]
		bool_NA <- is.na(my_row)
		trash <- .C('ErrorInt_row',
		    as.double(my_row),
		    as.integer(ncol),
		    as.integer(0),
		    NAOK = TRUE,
		    PACKAGE = "Harshlight")
		my_row <- trash[[1]]
		status <- trash[[3]]

		object[i,] <- my_row
		object[i,][bool_NA] <- NA
		
		if(status){
			return(object <- "a")
		}

	}

	for(i in 1:ncol){
		my_col <- object[,i]
		bool_NA <- is.na(my_col)
		trash <- .C('norm',
			  as.double(my_col),
			  as.integer(nrow),
			  as.integer(0),
			  NAOK = TRUE,
			  PACKAGE = "Harshlight")
		my_col <- trash[[1]]
		status <- trash[[3]]

		if(status){
			return(object <- "a")
		}

		object[,i] <- my_col
                    object[,i][bool_NA] <- NA
	}
	#d1 = date()
	return(object)
	#res <- list(E=E,d0 = d0,d1 =d1) 
}

#this is the function that get the non-empty bound of chip--yupu April 19

.GetBound <-function(data){
	cdf = getCdfEnvAffy(data) 
	index = indexProbes.CdfEnvAffy(cdf,which=c("pm","mm"))
	min_x = ncol(data)
	min_y = nrow(data)
	max_x = 0
	max_y = 0

	len = length(index)

	for(i in 1:len){
		xy = index2xy(cdf,index[[i]])
		mi_x = min(xy[,1])
		ma_x = max(xy[,1])

		if(mi_x < min_x){min_x = mi_x}
		if(ma_x > max_x){max_x = ma_x}

		mi_y = min(xy[,2])
		ma_y = max(xy[,2])

		if(mi_y < min_y){min_y = mi_y}
		if(ma_y > max_y){max_y = ma_y}
	}
	rm(cdf)
	rm(index)
	list(x1=min_x,x2=max_x,y1=min_y,y2=max_y)
}




###############
# AnalyzeChip --- purpose is to define the defected areas
###############

.AnalyzeChip <- function(chip, extended.radius = 10, compact.quant.bright = 0.025, compact.quant.dark = 0.025, compact.size.limit = 15, compact.connect = 8, compact.pval = 0.01, diffuse.bright = 40, diffuse.dark = 35, diffuse.pval = 0.001, diffuse.connect = 8, diffuse.radius = 10, diffuse.size.limit = (3*3.14*(diffuse.radius**2)), percent.contiguity = 50, NROW, NCOL, diffuse.close = TRUE, p.values.bright, p.values.dark, report.name)

{
	#sink("bug.txt")

	image1 <- chip

	################### BIG DEFECTS ########################

	#Calculates the difference between the standard deviation of the error image and the median image.
	#perc.var > 35% seems a good cutoff for big defected chips.

	dbg <- 0

	if(dbg){
		print("enter EXT")
		flush.console()
	}

	extended.median <- .FindExtended(chip, extended.radius) #return the image convoluted with a median kernel

	if(dbg){
		print("exit EXT")
		flush.console()
	}
	
        sd.Im <- sd(chip)
        sd.med <- sd(extended.median)
        perc.var <- 100*(sd.med^2/sd.Im^2)
	#perc.var <- 20

	if(perc.var >= 30){
		print("The chip has an extended defect.")
		print("It is recommended that the chip be eliminated from the batch and the analysis repeated.")
		flush.console()
		if(report.name != 0){
			trash <- .C("extended_stop",
				as.double(perc.var),
				PACKAGE = "Harshlight")
		}
		return(image1)
	}

	#################### SMALL DEFECTS #####################
	#Eliminates the small defects from X

	#Returns the original matrix without the small defects

	if(dbg){
		print("enter COM")
		flush.console()
	}

	chip <- .FindCompact(chip, compact.quant.bright, compact.quant.dark, compact.size.limit, compact.connect, p.values.bright, p.values.dark, compact.pval, radius = diffuse.radius, type = "small", percent.contiguity)

	if(dbg){
		print("exit COM")
		flush.console()
	}

	image2 <- chip$comb.mat

	if(report.name != 0){
		stat.c <- .ClusterNumber(image2)

		cluster.stat.c <- chip$size
		exist.cluster.c <- which(cluster.stat.c != 0)

		trash <- .C("chip_image",
			as.integer(50),
			as.integer(225),
			as.integer(.ScaleErr(image2)),
			as.character("Compact defects"),
			as.integer(exist.cluster.c),
			as.integer(cluster.stat.c[exist.cluster.c]),
			as.integer(length(exist.cluster.c)),
			as.integer(1),
			NAOK = TRUE,
			PACKAGE = "Harshlight")
	}

	chip <- chip$original

	#################### DIFFUSE DEFECTS ###################

	#Cutoff value to identify outliers in diffuse defects

	diffuse.bright <- log2(1 + diffuse.bright/100) #transform the percent into a log2 value
	diffuse.dark <- log2(1 + diffuse.dark/100) #transform the percent into a log2 value

	if(!diffuse.bright){
		diffuse.bright <- max(chip)
	}
	if(!diffuse.dark){
		diffuse.dark <- -min(chip)
	}

	q0 <- qnorm(diffuse.pval, lower.tail = FALSE)
	thres.dark <- mean((chip < (-diffuse.dark)), na.rm = TRUE)
	thres.bright <- mean((chip >  diffuse.bright ), na.rm = TRUE)

	#Returns two binary matrices, one for bright and one for dark outliers, based on the fraction of outliers present in a certain area

	if(dbg){
		print("enter DIF")
		flush.console()
	}

	diff.flag <- .FindDiffuse(img = chip, diffuse.bright = diffuse.bright, diffuse.dark = diffuse.dark, radius = diffuse.radius, quant = q0, thres.dark = thres.dark, thres.bright = thres.bright)

	if(dbg){
		print("exit DIF")
		flush.console()
	}

	#Clusters the outliers. Remember: ClusterDefects returns the outliers having negative numbers

	if(dbg){
		print("en cluster DIF")
		flush.console()
	}

	diffuse.bright <- .ClusterDefects(diff.flag$out.bright, diffuse.size.limit, diffuse.connect, p.values.bright, compact.pval, type = "diffuse")
	diffuse.dark <- .ClusterDefects(diff.flag$out.dark, diffuse.size.limit, diffuse.connect, p.values.dark, compact.pval, type = "diffuse")

	if(dbg){
		print("ex cluster DIF")
		flush.console()
	}

##########################################

	cluster.stat.d1 <- diffuse.bright$sz
	cluster.stat.d2 <- diffuse.dark$sz
	cluster.stat.d <- cluster.stat.d1 + cluster.stat.d2
	exist.cluster.d <- which(cluster.stat.d != 0)

	diffuse.bright <- diffuse.bright$img
	diffuse.dark <- diffuse.dark$img

##########################################

	stat.d1 <- .ClusterNumber(diffuse.bright)
	diffuse.bright[diffuse.bright < 0] <- -1
	stat.d2 <- .ClusterNumber(diffuse.dark)
	diffuse.dark[diffuse.dark < 0] <- -1
	
	cluster.stat.d1.p <- (sum(which(cluster.stat.d1 != 0)*cluster.stat.d1[which(cluster.stat.d1 != 0)])/length(cluster.stat.d1))*100
	cluster.stat.d2.p <- (sum(which(cluster.stat.d2 != 0)*cluster.stat.d2[which(cluster.stat.d2 != 0)])/length(cluster.stat.d2))*100

	#Puts together the diffuse.dark and diffuse.bright arrays (one image for both dark and bright values)

	if(!diffuse.close){
		diff.overlap <- diffuse.dark + diffuse.bright
		diff.overlap <- as.integer(diff.overlap == -2)

		diffuse.bright <- (-diffuse.bright)

	} else {
		#Closes the holes to better define the areas
	
		if(dbg){
			print("enter IMG CLOSE")
			flush.console()
		}


		diffuse.bright <- .ImageClose(diffuse.bright * (-1), radius = diffuse.radius)
		diffuse.dark <- .ImageClose(diffuse.dark * (-1), radius = diffuse.radius)

		if(dbg){
			print("exit IMG CLOSE")
			flush.console()
		}
	
		#Combines the two arrays again in one single image
	
		diff.overlap <- diffuse.dark + diffuse.bright
		diff.overlap <- as.integer(diff.overlap == 2)
		diffuse.dark <- (-diffuse.dark)

		stat.d1$perc <- (sum(diffuse.bright != 0)/length(diffuse.bright))*100
		stat.d2$perc <- (sum(diffuse.dark != 0)/length(diffuse.dark))*100
		cluster.stat.d1.p <- (sum(diffuse.bright != 0)/length(diffuse.bright))*100
		cluster.stat.d2.p <- (sum(diffuse.dark != 0)/length(diffuse.dark))*100

	}

	cluster.stat.overlap <- (sum(diff.overlap == 1)/length(diff.overlap))*100
	diffuse.comb <- diffuse.bright + diffuse.dark + diff.overlap/2

	if(report.name != 0){

		stat.d <- as.list(c())
		stat.d$num <- stat.d1$num + stat.d2$num
		stat.d$perc <- stat.d1$perc + stat.d2$perc - cluster.stat.overlap
		cluster.stat.d.perc <- cluster.stat.d1.p + cluster.stat.d2.p - cluster.stat.overlap

          	trash <- .C( "chip_image",
			as.integer(300),
			as.integer(225),
			as.integer(.ScaleErr(diffuse.comb)),
			as.character("Diffuse defects"),
			as.integer(exist.cluster.d),
			as.integer(cluster.stat.d[exist.cluster.d]),
			as.integer(length(exist.cluster.d)),
			as.integer(2),
			NAOK = TRUE,
			PACKAGE = "Harshlight")

		trash <- .C("chip_summary",
			as.double(perc.var),
			as.integer(stat.c$num),
			as.integer(stat.d$num),
			as.double(stat.c$perc),
			as.double(stat.d$perc),
			PACKAGE = "Harshlight")
	
		#trash <- .C("chip_summary",
		#	as.double(perc.var),
		#	as.integer(sum(cluster.stat.c)),
		#	as.integer(sum(cluster.stat.d)),
		#	as.double((sum(exist.cluster.c*cluster.stat.c[exist.cluster.c])/length(cluster.stat.c))*100),
		#	as.double(cluster.stat.d.perc),
		#	NAOK = TRUE,
		#	PACKAGE = "Harshlight")
	}

	image1[image2 != 0] <- NaN
	image1[diffuse.comb != 0] <- NaN

	return(image1)

	#sink()
}

######################
# CheckParam --- purpose is to verify that the input data are valid
######################

.CheckParam <- function(diffuse.dark, diffuse.bright, compact.quant.bright, compact.quant.dark, compact.size.limit, diffuse.size.limit, compact.connect, diffuse.connect, diffuse.radius, diffuse.pval, extended.radius, compact.pval, NROW, NCOL, percent.contiguity, na.sub, interpolate, diffuse.close){

	if(!(NROW == round(NROW)) || !(NCOL == round(NCOL)) || NROW <= 0 || NCOL <= 0){
		print("Error: the chip has no valid row or column length")
		return(1)
	}

	if(diffuse.radius <= 0 || diffuse.radius > NROW || diffuse.radius > NCOL){
		######### less than 20 for 640X640 ###########
		print("Error: diffuse.radius must be a positive number")
		print("and less than the dimensions of the chip")
		return(1)
	}

	if(extended.radius <= 0 || extended.radius > NROW || extended.radius > NCOL){
		########### less than 1/10 ###############
		print("Error: extended.radius must be a positive number")
		print("and less than the dimensions of the chip")
		return(1)
	}

	if(compact.quant.bright > 0.40 || compact.quant.bright < 0 || compact.quant.dark > 0.40 || compact.quant.dark < 0){
		print("Error: compact.quant.bright/compact.quant.dark must be between 0 and 0.40")
		return(1)
	}

	if(compact.size.limit < 0 ||compact.size.limit >= NROW*NCOL || diffuse.size.limit < 0 || diffuse.size.limit >= NROW*NCOL){
		print("Error: compact.size.limit/diffuse.size.limit out of bound")
		return(1)
	}

	if(!is.numeric(diffuse.bright) || !is.numeric(diffuse.dark)){
		print("Error: diffuse.bright/diffuse.dark cannot be bigger than 100%")
		return(1)
	}

	if((compact.connect != 0 & compact.connect != 1) || (diffuse.connect != 0 & diffuse.connect != 1)){
		print("Error: compact.connect/diffuse.connect out of bound")
		return(1)
	}

	if(diffuse.pval <= 0 || diffuse.pval > 1){
		print("Error: diffuse.pval out of bounds")
		return(1)
	}

	if(compact.pval <= 0 || compact.pval > 1){
		print("Error: compact.pval out of bounds")
		return(1)
	}

	if(percent.contiguity < 0 || percent.contiguity > 100){
		print("Error: percent.contiguity out of bounds")
		return(1)
	}
	if(!is.logical(na.sub)){
		print("Error: na.sub must be either TRUE or FALSE")
		return(1)
	}

	if(!is.logical(interpolate)){
		print("Error: interpolate must be either TRUE or FALSE")
		return(1)
	}
	if(!is.logical(diffuse.close)){
		print("Error: diffuse.close must be either TRUE or FALSE")
		return(1)
	}

	return(0)
}

#######################
# RetrievePvalues_exact --- purpose is to assigne the best simulation according to the quantile used
#######################

.RetrievePvalues_exact <- function(quantile.small, min.pval, sim.pval, size){

	#Choose the best p.values for the quantile used
	#allowed.pval <- c(0.01, 0.02, 0.05, 0.10, 0.20, 0.25, 0.30, 0.40)
	#min.pval <- abs(quantile.small - allowed.pval)

	p.values <- unlist(sim.pval[which(min.pval == min(min.pval))[1]])
	p.values <- c(p.values, seq(length = (size - length(p.values) + 1), 0, 0))/100000
	
	return(p.values)
}

#######################
# RetrievePvalues --- purpose is to assigne the best simulation according to the quantile used
#######################

.RetrievePvalues <- function(quantile.small, a, b, size){

	#Choose the best p.values for the quantile used

	new.a <- a[1]*(quantile.small - a[2])^2 # k*(new.d-c)^2
	new.b <- b[1]*(quantile.small - b[2])^2 # k*(new.d-c)^2
	new.p <- c(0, 1/(1+exp((c(1:size) - new.b)/new.a)))

	return(new.p)
}

####################
# FindExtended --- purpose is to convolve an image with a median filter applying a circular kernel
####################

.FindExtended <- function(img, radius){

	#x error intensity matrix
	#radius of the kernel to use

	med.obs <- seq(length = length(img), 0, 0) #matrix in which to store the convolved matrix
	
	f<-.C("extended_defects",
		as.double(img),
		med.obs <- as.double(med.obs),
		as.integer(radius),
		status <- as.integer(0),
		NAOK = TRUE,
		DUP = FALSE,
		PACKAGE = "Harshlight")

	if(status){
		print("	Memory problem: cannot analyze extended defects")
		flush.console()
	}

	return(median.obs = med.obs)
}

########################
# FindCompact --- purpose is to identify the small defects in a chip. Returns the original matrix without the small defects - try splitting dark and bright again
#		        but remembering the cluster.id from the previous call
########################

.FindCompact <- function(img, quant.bright = 0.025, quant.dark = 0.025, size.limit = 15, connect = 1, simulation.bright, simulation.dark, compact.pval = 0.01, radius = 10, type = "small", percent.contiguity = 50){

	#img		error intensity matrix
	#quant.bright	quantile to define the bright outliers
	#quant.dark	quantile to define the dark outliers
	#size.limit	limit of the cluster size to be considered of an acceptable size
	#connect		0 = 4-point neighborhood, 1 = 8-point nrighborhood
	#simulation	cluster size from random matrices
	#compact.pval	significance of cluster size

	original <- img

	img_dr <- img

	#######################################################
	#Calculates the quantiles to define the outliers
	
	q.bright <- quantile(img, 1 - quant.bright, na.rm = TRUE)
	q.dark <- quantile(img, quant.dark, na.rm = TRUE)

	#######################################################
	#Substitutes the values in the matrix: 1 if >= q.bright or <= q.dark, 0 otherwise

	img <- as.integer(img > q.bright)
	img_dr <- as.integer(img_dr < q.dark)

	return(.Contiguity(original, img, img_dr, size.limit, connect, simulation.bright, simulation.dark, compact.pval, radius, type = type, percent.contiguity))
}

########################
# Called by FindCompact
# Contiguity --- purpose is to find the small brigth and dark defects in the error intensity matrix
########################

.Contiguity <- function(original, mat.bright, mat.dark, size.limit = 15, connect = 1, simulation.bright, simulation.dark, compact.pval = 0.01, radius = 10, type, percent.contiguity = 50){

	#Cluster the small defects - returns a matrix in which flagged cells are clustered; each cluster is assigned to a negative integer
	mat.bright <- .ClusterDefects(mat.bright, size.limit, connect, simulation.bright, compact.pval, type = type)
	mat.dark <- .ClusterDefects(mat.dark, size.limit, connect, simulation.dark, compact.pval, type = type)

##############################
	
	size.brdr <- mat.bright$sz + mat.dark$sz
	mat.bright <- mat.bright$img
	mat.dark <- mat.dark$img

##############################

	#Combined image from the closed images
	comb.mat <- mat.bright*(-1) + mat.dark

	if(percent.contiguity){

		#Close the areas covered by small defects - returns a matrix in which the clusters = 1 and the background = 0
		mat.br.closed <- .ImageClose(mat.bright * (-1), radius = radius)
		mat.dr.closed <- .ImageClose(mat.dark * (-1), radius = radius)

		#Since the closing procedure might slightly change the boundaries of the clusters, we need to add all the points
		#that were flagged to the closed image, so that the cleaning procedure for the small clusters that belong to
		#diffuse defects includes all the flagged points inside the closed image - returns a matrix in which closed clusters are assigned a negative integer
		mat.br.closed[mat.bright != 0] <- 1
		mat.dr.closed[mat.dark != 0] <- 1

		#Assign to each closed region a cluster negative integer
		mat.br.closed <- .ClusterDefects(mat.br.closed, size.limit, connect, simulation.bright, compact.pval, type = type)
		mat.dr.closed <- .ClusterDefects(mat.dr.closed, size.limit, connect, simulation.dark, compact.pval, type = type)

		mat.br.closed <- mat.br.closed$img
		mat.dr.closed <- mat.dr.closed$img

		#Includes overlapping among closed regions of dark and bright defects
		n <- min(mat.dr.closed) - 1
		lng <- length(mat.dr.closed)
		for(i in 1: lng){
			if((mat.br.closed[i] != 0 && mat.dr.closed[i] != 0) && (mat.br.closed[i] != mat.dr.closed[i])){
				mat.br.closed[mat.br.closed == mat.br.closed[i]] <- n
				mat.dr.closed[mat.dr.closed == mat.dr.closed[i]] <- n
	
				mat.br.closed[mat.dr.closed == n] <- n
	
				n <- n - 1
			}
		}
	
		mat.dr.closed[mat.br.closed != 0] <- 0
		comb.mat.closed <- mat.br.closed + mat.dr.closed
	
		f <- function(x){
			n <- max(-x)
			s <- seq(length = n, 0, 0)
			for(i in 1:n){
				s[i] <- sum(x == -i)
			}
			return(s)
		}
		s <- f(comb.mat.closed)
	
		comb.mat.closed[comb.mat.closed%in%c(-(which(s < 100 & s > 0)))] <- 0
		image3 <- comb.mat.closed * (comb.mat == 0)
	
		if(sum(image3 != 0)){
			for(i in 1:max(-image3)){
				n <- (100*(1 - sum(image3 == -i)/sum(comb.mat.closed == -i)))
				if(!is.nan(n) && n < percent.contiguity){
					comb.mat[comb.mat.closed == -i] <- 0
				}
			}
		}
	}

	#Returns the original matrix in which the small defects are erased

	original[comb.mat != 0] <- 0
	return(list(original = original, comb.mat = comb.mat, size = size.brdr))

}

###################
# FindDiffuse --- purpose is to define the diffuse defects in the error matrix x,
# calculating the fraction of outliers present in a sliding window, considering the significance of this fraction
###################

.FindDiffuse <- function(img, diffuse.bright, diffuse.dark, radius = 10, quant, thres.dark, thres.bright){

	#img			error intensity matrix
	#diffuse.bright, diffuse.dark	threshold value to determine outliers in the error matrix
	#radius			radius of the sliding window in which to count the fraction of outliers
	#quant		
	#thres.dark, thres.bright	values used to calculate the significance of the fraction of outliers found in the sliding window
				#in the c function diffuse_defects
	#diff.bright and diff.dark	store the matrices returned from the c function, in which the pixels in defected areas have value 1

	diff.bright <- seq(length = length(img), 0, 0)
	diff.dark <- seq(length = length(img), 0, 0)

	trash <- .C("diffuse_defects",
		img <- as.double(img),
		as.double(diffuse.bright),
		as.double(-diffuse.dark),
		as.integer(radius),
		diff.bright <- as.double(diff.bright),
		diff.dark <- as.double(diff.dark),
		as.double(quant),
		as.double(thres.dark),
		as.double(thres.bright),
		status <- as.integer(0),
		DUP = FALSE,
		NAOK = TRUE,
		PACKAGE = "Harshlight")

	if(status){
		print("	An error occurred while analyzing diffuse defects")
		flush.console()
	}

	return(list(out.bright = diff.bright, out.dark = diff.dark));
}

###################
# ClusterDefects --- purpose is to cluster the outliers - up, left, down, right
###################

.ClusterDefects <- function(img, size.limit, connect, simul.pval, compact.pval, type)
{

	#img		matrix (0,1) in which to look for clusters
	#limit		limit of cluster size below which the cluster is not considered as such
	#connect		0 = 4-point neighborhood, 1 = 8-point neighborhood
	#simul.pval	cluster size from random arrays
	#compact.pval	significance of cluster size
	#type		'small' or 'diffuse'; defines the function that called ClusterDefects

	type <- 1 - (type == "small")

	array.size <- seq(length = length(img) + 1, 0, 0) #stores the information about the different cluster size
	trash <- .C("cluster_defects", 
      		img <- as.integer(img), 
		array.size <- as.integer(array.size), 
		as.integer(size.limit), 
		as.integer(connect),
		as.double(simul.pval),
		as.double(compact.pval), 
		as.integer(type),
		status <- as.integer(0),
		DUP = FALSE,
		PACKAGE = "Harshlight")

	array.size <- array.size[2:length(array.size)]
	if(status){
		print("A problem occurred while analyzing compact or diffuse defects")
		flush.console()
		img <- (abs(img))*(-1)
	}

	return(list(img = img, sz = array.size))
}

###################
# ScaleInt --- scales the original and error images
###################

.ScaleInt <- function(img){
	
	#q <- max(-quantile(img, .01, na.rm=T), quantile(img, .99, na.rm=T))
	q <- max(-quantile(img, .001, na.rm = TRUE), quantile(img, .999, na.rm = TRUE))
	#q <- max(-quantile(img, .0001, na.rm=T), quantile(img, .9999, na.rm=T))
	img[img < (-q)]<- (-q)
	img[img > ( q)]<- ( q)
	image.min <- min(img)
	image.max <- max(img)
	image.range <- image.max - image.min
	img <- round(((img - image.min)/image.range)*255)

	return(img)

}

###################
# ScaleErr --- scales the different defect images
###################

.ScaleErr <- function(img){

	img[img > 0.5] <- 255
	img[img == 0.5] <- 185
	img[img == 0] <- 128
	img[img < 0] <- 0

	return(img)

}

###################
# ClusterNumber --- statistics for clustered image (total number of clusters and percentage of area covered by the defects)
###################

.ClusterNumber <- function(img){
	n <- 0
	for(i in min(img):max(img)){
		if(sum(img[img == i]) & i){
			n <- n + 1
		}
	}
	return(list(num = n, perc = (sum(img != 0)/length(img))*100))
}

###################
# ImageClose --- purpose is to fill the blank spaces between features of an image (defines areas better)
###################

.ImageClose <- function(img, radius = 10){

	#img		image to which to apply the dilatation and erosion
	#radius		radius of the kernel to use

     	result <- seq(length = length(img), 0, 0)

	trash <- .C("image_dilation", as.double(img), result <- as.double(result), as.integer(radius), status <- as.integer(0), DUP = FALSE, PACKAGE = "Harshlight")
	if(status){
		print("	A problem occurred while dilating the image")
		flush.console()
		return(abs(img))
	}

	trash <- .C("image_erosion", as.double(result), img <- as.double(img), as.integer(radius), status <- as.integer(0), DUP = FALSE, PACKAGE = "Harshlight")
	if(status){
		print("	A problem occurred while eroding the image")
		flush.console()
		return(abs(img))
	}

	return(img)

}

