\name{epik}
\alias{epik}
\title{
Run pre-defiend scripts
}
\description{
Run pre-defiend scripts
}
\usage{
epik()
}
\details{
There are several R scripts which can be run directly. The path of all scripts can be obtained by

  \preformatted{
   dir(system.file("pipeline", "script", package = "epik"), pattern = "\\.R$")  }

You can either directly run these R scripts by:

  \preformatted{
   Rscript full-path-of-cmd.R [options]  }

or use the short command:

  \preformatted{
   Rscript -e "epik::epik()" cmd [options]  }

to get the list of commands:

  \preformatted{
   Rscript -e "epik::epik()"  }

For each cmd, use \code{Rscript -e "epik::epik()" cmd --help} to get help.

Basically all scripts need \code{--config} option which corresponds to a configuration R file
that defines how to import data. A template for the configuration file is at:

  \preformatted{
  system.file("pipeline", "pipeline_configure_template.R", package = "epik")  }

Available commands (\code{cmd}) are:

\describe{
  \item{chromatin_states_transitions}{visualize chromatin states transitions by chord diagram}
  \item{correlated_enriched}{visualize enrichment of correlated regions on tss/cgi/tfbs/enhancers}
  \item{correlated_regions}{find regions in which methylation is correlated to expression of the associated gene}
  \item{correlated_regions_downstream}{visualize statistics of correlated regions, visualize genome-wide distribution of correlated regions by Hilbert curve}
  \item{correlated_regions_filter}{only keep correlation regions with significant correlations. Subgroup specificity for each class is calculated if needed. }
  \item{correlated_regions_gviz}{visualize correlated regions and other associated information by Gviz package}
  \item{correlated_regions_reduce}{reduce correlated regions}
  \item{differential_methylation_in_cgi_and_shore}{visualize differentially methylated regions in cgi and shores}
  \item{differential_methylation_in_genomic_features}{visualize differentially methylated regions in a set of genomic features.}
  \item{general_methylation_distribution}{use heatmap to visualize methylation distribution}
  \item{methylation_subtype_classification_in_cgi_and_shore}{classify subgroups by methylation in cgi and shores}
}
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
