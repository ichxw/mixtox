readModel <- function(File){
	# # read concentration-response data from txt file
	getExtension <- function(File){ 
		ex <- strsplit(basename(File), split="\\.")[[1]]
		return(ex[-1])
		}

	extension = getExtension(File)

	if(extension == 'txt') sept <- "" else sept <- ","

    ds <- read.table(File, fill = T, sep = sept, row.names = NULL, header = TRUE, comment.char = "", quote = "")

    model <- as.vector(ds[, 2])
    param <- as.matrix(ds[, c(3 : 5)])
    rownames(param) <- ds[, 1]

    list(model = model, param = param)
}