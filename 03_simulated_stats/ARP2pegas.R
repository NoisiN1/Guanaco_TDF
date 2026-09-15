#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#						Claudio S. Quilodrán | 24.06.2025
#							claudio.quilodran@unige.ch
# 						   Importing ARP files to pegas
# 	   This function works for diploid DNA sequences or SNP data from SPLATCHE simulations
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
require(pegas)

#' @title ARP files to pegas
#'
#' @description This function generates an object of class \code{"loci"} from .arp files generated from SPLATCHE simulations. 
#' @param filename arlequin file generated from SPLATCHE
#' @param otherpackage generates an object for use with either the 'adegenet' or 'hierfstat' package. By default (NULL), it is imported into 'pegas'.
#' @details This function works for diploid DNA sequences or SNP data from SPLATCHE simulations.
#' @return An object of class \code{"loci"}, \code{"genind"} or a \code{"data.frame"}.
#' @seealso \code{\link{loci2genind}}
ARP2pegas <- function(filename, otherpackage=NULL){

if (!is.null(otherpackage) && !(otherpackage %in% c("pegas", "adegenet", "hierfstat"))) { stop("This function only imports .arp files to packages 'pegas', 'adegenet' or 'hierfstat'.", call. = FALSE) }

arlequin_file<-readLines(filename)

###############
# Verify data #
###############	
dataType <- grep("DataType=", arlequin_file, value = T)
diploid<-as.numeric(gsub("[^0-9]", "", grep("GenotypicData=", arlequin_file, value = T)))

if (!grepl("RFLP|DNA", dataType)) stop("This function only works for diploid sequences or SNP data from SPLATCHE.", call.=F)
if (diploid!=1) stop("This function only works for diploid data from SPLATCHE.", call.=F)
###############

dataType <- ifelse(grepl("RFLP", dataType), "SNP", "DNA")

start <- grep("SampleData=", arlequin_file)+1
end <- grep("}", arlequin_file)-2

popnames <- (gsub('.*SampleName="([^"]+)".*', '\\1', grep("SampleName", arlequin_file, value = TRUE)))


	switch(dataType, DNA ={

sample<-sapply(1:length(start), function(i){
	
	samp<-read.table(textConnection(arlequin_file[start[i]:end[i]]), fill=TRUE)
	samp2<-sapply(1:nrow(samp), function(j){ grep("^[^N]*[TACG][^N]*$", samp[j,], value = TRUE) })
		
	even<-which(1:length(samp2) %% 2 == 0)
	odd<-which(1:length(samp2) %% 2 == 1)

	samp3<-t(sapply(1:length(even), function(k){ 
		
	evenseq<-unlist(strsplit(samp2[even[k]], ""))
	oddseq<-unlist(strsplit(samp2[odd[k]], ""))

	paste(oddseq, evenseq, sep = "-")

	}) )
	
	cbind(popnames[i], samp3)
}, simplify=F)


samples<-as.data.frame(do.call(rbind, lapply(seq_along(sample), function(i) {
  df <- sample[[i]]
  df
})), stringsAsFactors = F)
colnames(samples)<-c("population", paste("p", 1:(ncol(samples)-1), sep=""))		
			
		
	}, 
	SNP = {
	
sample<-sapply(1:length(start), function(i){
	
	samp<-read.table(textConnection(arlequin_file[start[i]:end[i]]), fill=TRUE)
	
	evenseq<-samp[which(1:nrow(samp) %% 2 == 0),1:(ncol(samp)-2)]
	oddseq<-samp[which(1:nrow(samp) %% 2 == 1),-c(1:2)]

	evenseq<-apply(evenseq, c(1, 2), as.numeric)	+1
	oddseq<-apply(oddseq, c(1, 2), as.numeric)+1
	
	evenseq<-apply(evenseq, c(1, 2), as.character)	
	oddseq<-apply(oddseq, c(1, 2), as.character)

	
	samp3<-t(sapply(1:nrow(evenseq), function(k){ 		
	paste(oddseq[k,], evenseq[k,], sep = "-")
	}) )
	
	cbind(popnames[i], samp3)
}, simplify=F)


samples<-as.data.frame(do.call(rbind, lapply(seq_along(sample), function(i) {
  df <- sample[[i]]
  df
})), stringsAsFactors = F)
colnames(samples)<-c("population", paste("p", 1:(ncol(samples)-1), sep=""))
		
	})	

temp <- tempfile()
  write.table(samples, sep = " ", file = temp, col.names = TRUE, quote = FALSE, row.names = FALSE)

data<-pegas::read.loci(temp, col.pop=1, allele.sep="-")

if (!is.null(otherpackage)) {
if (otherpackage == "adegenet") { 
  #require(adegenet)
  data <- loci2genind(data)

} else if (otherpackage == "hierfstat") {
  #require(adegenet)
  require(hierfstat)
  data <- loci2genind(data)
  data <- genind2hierfstat(data)
}	
}	
	return(data)
	 
}