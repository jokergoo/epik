\name{systemdf}
\alias{systemdf}
\title{
Wrapper of system calls in which input and output are all table-like files
}
\description{
Wrapper of system calls in which input and output are all table-like files
}
\usage{
systemdf(cmd, envir = parent.frame(), verbose = FALSE)
}
\arguments{

  \item{cmd}{shell command}
  \item{envir}{environment where to look for variables encoded in \code{cmd}}
  \item{verbose}{whether print messages}

}
\details{
This function (system + data frame) provides a convinient way to invoke
system calls in R. Since most of system calls expect tables as inputs and outputs, 
\code{systemdf()} does following step by step:

\itemize{
  \item use backtick to mark variables which are data frames or other variables which can be converted to data frames by \code{\link[base]{as.data.frame}}
  \item extract the variable names and look for data frames in \code{envir}
  \item write data frames into temporary files 
  \item replace variables names with paths that correspond to temporary files
  \item make the system call 
  \item finally send back the output by piping back to R
}

A simple example is as follows:

  \preformatted{
  bed1 = circlize::generateRandomBed(nr = 1000)
  bed2 = circlize::generateRandomBed(nr = 1000)
  df = systemdf("bedtools closest -a `bed1` -b `bed2` | awk '$1==\"chr1\"'")  }
}
\value{
A data frame. Column names may be lost.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
if(Sys.info()["sysname"] \%in\% c("Linux", "Darwin")) {
    df = data.frame(x = sample(1:10, 10), y = sample(11:20, 10))
    systemdf("sort -k1,1n `df`")
}
}
