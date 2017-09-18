\name{submit_to_david}
\alias{submit_to_david}
\title{
Doing DAVID analysis
}
\description{
Doing DAVID analysis
}
\usage{
submit_to_david(genes, email, catalog = c("GOTERM_CC_FAT", "GOTERM_BP_FAT", "GOTERM_MF_FAT", "KEGG_PATHWAY"),
    idtype = "ENSEMBL_GENE_ID", species = "Homo sapiens")
}
\arguments{

  \item{genes}{a vector of genes}
  \item{email}{the email user registered on David web service}
  \item{catalog}{a vector of functional catelogs}
  \item{idtype}{id types}
  \item{species}{species if the id type is not unique to a species}

}
\examples{
# There is no example
NULL

}
