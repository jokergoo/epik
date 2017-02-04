\name{run_pipeline}
\alias{run_pipeline}
\title{
run pipeline in HPC through job scheduling system
}
\description{
run pipeline in HPC through job scheduling system
}
\usage{
run_pipeline(config_file, prefix = "", email = NULL, enforce = FALSE, Rscript_binary = "Rscript", submit_by = "qsub")
}
\arguments{

  \item{config_file}{path of configuration script, check \code{\link{load_config}} for all configurations.}
  \item{prefix}{prefix of the job name}
  \item{email}{email address if you want to be notified of your jobs}
  \item{enforce}{enforce to run all the steps even for those successfully finished jobs.}
  \item{Rscript_binary}{path of Rscript binary}
  \item{submit_by}{which job scheduling you are using}

}
\details{
A workflow will be submitted to HPC.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
