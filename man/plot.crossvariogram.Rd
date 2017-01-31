\name{plot.crossvariogram}
\alias{plot.crossvariogram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Crossvariogram and lmc fitting plot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
It displays both crossvariogram and lmc fitting.
}
\usage{
plot.crossvariogram(CV, RES = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CV}{
%%     ~~Describe \code{CV} here~~
Crossvariogram resulting from \code{crossvariogram}
}
  \item{RES}{
%%     ~~Describe \code{RES} here~~
lmc fitting resulting from \code{lmc}. If NULL, it only shows the crossvariogram plot.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Victor Vicente Palacios
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{crossvariogram},\link{lmc}}
}
\examples{
data("iris")
Versicolor <- iris[which(iris$Species=='versicolor'),-5]
##Data Standarization
means_vers <- apply(Versicolor,2,mean)
sd_vers <- apply(Versicolor,2,sd)
Versicolor_st <- Versicolor

for (i in 1:length(Versicolor[1,]))
{Versicolor_st[,i] <- (Versicolor[,i]-means_vers[i])/sd_vers[i]}

##PrComp
PC_train <- princomp(Versicolor_st)
## CrossVariogram Calculation
CV_vers <- crossvariogram(as.data.frame(PC_train$scores[,1:2]),as.data.frame(Versicolor_st),11)
## Plot
plot.crossvariogram(CV_vers)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dplot }% use one of  RShowDoc("KEYWORDS")

