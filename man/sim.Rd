\name{sim}
\docType{data}
\alias{sim}
\alias{sim.int}
\alias{a}
\alias{b}
\alias{simulations}
\title{blemish simulations from 100.000 random chips}
\description{This data set contains the probability distribution of the cluster size under the assumption of spatially randomly distributed outliers. The distribution depends on the chip size, the density of outliers and the definition of connectivity (see help(Harshlight)).
The sets sim\_n4\_/sim\_n8\_ contain the results from 100.000 simulations (sim.pval), while the sets sim\_4\_/sim\_8\_ contain the parameters that are used to interpolate the probability for values of density of outliers and chip sizes not simulated (a, b). The data set is used by the package Harshlight and is not intended for direct use by the user.
}
\usage{
 data(sim_n4_534)
 data(sim_n4_640)
 data(sim_n4_712)
 data(sim_n8_534)
 data(sim_n8_640)
 data(sim_n8_712)
 data(sim_4_534)
 data(sim_4_640)
 data(sim_4_712)
 data(sim_8_534)
 data(sim_8_640)
 data(sim_8_712)
}
\format{The number of occurrences of at least one cluster of a certain size in 100.000 chips with randomly distributed outliers.}
\keyword{datasets}
