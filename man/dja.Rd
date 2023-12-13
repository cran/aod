\name{dja}
\alias{dja}
\docType{data}

\title{ Mortality of Djallonke Lambs in Senegal }

\description{
  Field trial to assess the effect of ewes deworming (prevention of gastro-intestinal parasitism) on the mortality
  of their offspring (age < 1 year). This data set is extracted from a large database on small ruminants production
  and health in Senegal (Lancelot et al., 1998). Data were collected in a sample of herds in Kolda (Upper Casamance,
  Senegal) during a multi-site survey (Faugère et al., 1992). See also the references below for a presentation of the
  follow-up survey (Faugère and Faugère, 1986) and a description of the farming systems (Faugère et al., 1990).
  }

\usage{data(dja)}

\format{
  A data frame with 21 observations on the following 4 variables.
  \describe{
    \item{group}{a factor with 2 levels: \code{CTRL} and \code{TREAT}, indicating the treatment.}
    \item{village}{a factor indicating the village of the herd.}
    \item{herd}{a factor indicating the herd.}
    \item{n}{a numeric vector: the number of animals exposed to mortality.}
    \item{trisk}{a numeric vector: the exposition time to mortality (in year).}
    \item{y}{a numeric vector: the number of deaths.}
  }
}

\references{
  Faugère, O., Faugère, B., 1986. \emph{Suivi de troupeaux et contrôle des performances individuelles des petits
  ruminants en milieu traditionnel africain. Aspects méthodologiques.} Rev. Elev. Méd. vét. Pays trop., 39 (1): 29-40.\cr
  Faugère, O., Dockès, A.-C., Perrot, C., Faugère, B., 1990. \emph{L'élevage traditionnel des petits ruminants
  au Sénégal. I. Pratiques de conduite et d'exploitation des animaux chez les éleveurs de la région de Kolda.} Revue
  Elev. Méd. vét. Pays trop. 43: 249-259.\cr
  Faugère, O., Tillard, E., Faugère, B., 1992. \emph{Prophylaxie chez les petits ruminants au Sénégal : régionalisation
  d'une politique nationale de protection sanitaire}. In: B. Rey, S. H. B. Lebbie, L. Reynolds (Ed.), First biennial
  conference of the African Small Ruminant Research Network, ILCA, 1990, ILRAD, Nairobi, pp. 307-314.\cr
  Lancelot, R., Faye, B., Juanès, X., Ndiaye, M., Pérochon, L., Tillard, E., 1998. \emph{La base de données BAOBAB:
  un outil pour modéliser la production et la santé des petits ruminants dans les systèmes d'élevage traditionnels
  au Sénégal.} Revue Elev. Méd. vét. Pays trop., 51 (2): 135-146.\cr
  }

\keyword{datasets}
