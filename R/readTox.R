readTox <- function(File, light = TRUE){
	# # read concentration-response data from txt file
	getExtension <- function(File){ 
		ex <- strsplit(basename(File), split="\\.")[[1]]
		return(ex[-1])
		}
	
	extension = getExtension(File)

	if(extension == 'txt') sept <- "" else sept <- ","
	if (light == TRUE){
		# Read dose-response data in a light format
		ds <- read.table(File, fill = T, sep = sept, row.names = NULL, header = F, comment.char = "", quote = "")
		conc <- as.vector(ds[-c(1), 1])
		mode(conc) <- "numeric"
		tier <- as.vector(unlist(ds[1, 2 : length(ds[1, ])]))
		expr <- as.matrix(ds[-c(1), 2 : ncol(ds)])
		mode(expr) <- "numeric"
		colnames(expr) <- tier
		list(x = conc, y = expr)

	} else{
		# Read dose-response data in a detailed format
		ds <- read.table(File, fill = T, sep = sept, row.names = NULL, header = F, comment.char = "", quote = "")
		cmpd <- as.vector(ds[1, 1])
		concNum <- as.numeric(as.character(ds[1, 2]))
		tierNum <- as.numeric(as.character(ds[1, 3]))
		type <- as.vector(ds[1, 4])
		conc <- as.vector(ds[-c(1, 2), 2])
		mode(conc) <- "numeric"
		tier <- as.vector(unlist(ds[2, 3 : length(ds[2, ])]))
		expr <- as.matrix(ds[-c(1, 2), 3 : ncol(ds[3, ]) ])
		mode(expr) <- "numeric"
		colnames(expr) <- tier
		list(x = conc, y = expr, name = cmpd, concNum = concNum, tierNum = tierNum, type = type)
	}
}
