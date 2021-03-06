% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mass_action_prediction.R
\name{set_colours}
\alias{set_colours}
\title{Sets colours based on an attribute}
\usage{
set_colours(attribute, color_list, ncolors, minimum, maximum)
}
\arguments{
\item{attribute}{the attribute to be mapped, a vector of values}

\item{color_list}{a list of colours to be made into a colour ramp}

\item{ncolors}{number of colours in the ramp}

\item{minimum}{minimum value in the colour ramp}

\item{maximum}{maximum value in the colour ramp}
}
\value{
a list of colours, linearly mapped to the colour ramp on the interval between minimum and maximum
}
\description{
Sets colours based on an attribute
}

