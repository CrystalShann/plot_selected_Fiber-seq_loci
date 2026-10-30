\name{cluster_colors}
\alias{cluster_colors}
\title{Original Fiber-seq cluster colors}
\description{Return the cluster palette used by locus heatmaps, accessibility
profiles and read annotations, matching the original Fiber-seq plotting helpers.}
\usage{cluster_colors(clusters)}
\arguments{
\item{clusters}{Character vector of cluster labels, or a factor whose levels
define the complete set of clusters. Labels must be nonempty and nonmissing.}
}
\value{Named character vector of colors, one per distinct label or factor level.}
\details{Labels of the form cluster1, cluster2, and so on use their numeric
suffix as a palette index. Thus cluster3 is green4 even when cluster1 and
cluster2 are absent. Other labels use their order of appearance. The original
25-color palette is interpolated when more colors are needed; pass the same
complete factor levels to keep interpolated colors consistent across subsets.
Empty input returns an empty named character vector.}
\examples{
cluster_colors(c("cluster1", "cluster2", "cluster3"))
cluster_colors(c("cluster3", "cluster1"))
}
