\name{register_global_var}
\alias{register_global_var}
\title{
Register global variables
}
\description{
Register global variables
}
\usage{
register_global_var(var_name)
}
\arguments{

  \item{var_name}{a vector of variable names}

}
\details{
Besides mandatory global variables checked by \code{\link{load_epik_config}}, the optional global variables
can be set by \code{\link{register_global_var}}. These optional variables will be exported to the working environment.
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
