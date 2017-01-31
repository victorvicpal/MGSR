\name{CrossValidation}
\alias{CrossValidation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Cross-Validation function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The cross-validation method is based on a common resampling. Each sample value Z(uα) is removed in turn from the data set and a value Z(u[α]) is estimated at that location using the N − 1 other samples.
}
\usage{
CrossValidation(RES)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RES}{
%%     ~~Describe \code{RES} here~~
Resulting lmc fitting.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Large datasets can last a very long time.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Object containing RMSE, pseudo-R^2, residuals and variance of residuals
}
\references{
%% ~put references to the literature/web site here ~
Wackernagel H. Multivariate Geostatistics: An Introduction with Applications. Springer-Verlag 2003.
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
\code{\link{lmc}}
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
#CrossValidation
Cr_val <- CrossValidation(RES_vers)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")

