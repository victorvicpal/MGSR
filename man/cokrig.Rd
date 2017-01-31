\name{cokrig}
\alias{cokrig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Cokriging
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Subspatial Cokriging estimation.
}
\usage{
cokrig(RES, xygrid)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RES}{
%%     ~~Describe \code{RES} here~~
Resulting lmc fitting
}
  \item{xygrid}{
%%     ~~Describe \code{xygrid} here~~
Grid
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
Dataframe containing coordinates and their estimated values.
}
\references{
%% ~put references to the literature/web site here ~
Myers D.E. Matrix Formulation of Co-Kriging. Math Geol. 1982;14(3):249-257.
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
\code{\link{lmc},\link{GRID_MGSR}}
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
## lmc 
RES_vers <- lmc(CV_vers,'Pow',1.6)

## Cokriging (standarized values)
Z_vers_st <- cokrig(RES_vers,xygrid)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
