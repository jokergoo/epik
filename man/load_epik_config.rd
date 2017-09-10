\name{load_epik_config}
\alias{load_epik_config}
\title{
Load and validate configuration file
}
\description{
Load and validate configuration file
}
\usage{
load_epik_config(config_file, export_env = parent.frame(), validate = TRUE)
}
\arguments{

  \item{config_file}{path of configuration file}
  \item{export_env}{environment where to export variables}
  \item{validate}{whether do validation}

}
\details{
To run functions in epic package smoothly, users can validate their data by \code{\link{load_config}} function.
All necessary variables are initialized in a configuration file.
The configuration file should provide following variables:

\code{SAMPLE}: a data frame, row names must be sample id and there must be a
  'class' column which corresponds to classes of samples. There can also
  be other annotation columns.

\code{COLOR}: a list of color settings corresponding to annotation column in 
  'SAMPLE'. The value of each list element must be either a named vector
  or a color mapping function. 'COLOR$class' must be defined or random color
  will be assigned. Names of other color settings should be same as
  corresponding columns in 'SAMPLE'.

\code{TXDB} (optional): a \code{\link[GenomicFeatures]{TxDb}} object.

\code{GTF_FILE} (optional): GTF file which is used to built \code{TXDB}. If it is not specified, the function tries to extract from \code{TXDB}.

\code{EXPR} (optional): a matrix which contains expression values. Column names 
  should be sample id and row names should be gene ids. Note gene ids in the 
  expression matrix should be same type as genes in \code{\link[GenomicFeatures]{TxDb}}.

\code{CHROMOSOME}: a vector of chromosome names.

\code{GENOME}: abbreviation of GENOME.

\code{GENOMIC_FEATURE_LIST}: a list of genomic features as GRanges objects. There
  must be a element named 'cgi'.

\code{MARKS} (optional): a vector of ChIP-Seq markers.

\code{methylation_hooks()} must be defined.

\code{chipseq_hooks()} is optional unless you want to do integrative analysis.

\code{CGI_SHORE_EXTEND}: extension of cgi, by default it is 2kb both upstream and downstream.

\code{\link{register_global_var}} registers optional global variables.
}
\value{
No value is returned.
}
\seealso{
\code{\link[epik.cmd]{epik}}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
