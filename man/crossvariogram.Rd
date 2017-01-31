\name{crossvariogram}
\alias{crossvariogram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Subspatial CrossVariogram Calculation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function calculates all possible crossvariograms within the values you referred. 
}
\usage{
CrossVariogram(coord,values,n=15)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coord}{
%%     ~~Describe \code{x} here~~
Dataframe containing the coordinates from Factorial Techniques.
}
\item{values}{Dataframe containing values from data.}
\item{n}{Number of pairs.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Only 2 dimension allowed.
Usually standarized values keep better performance.
Use \code{n_pairs_opt} to evaluate it.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Object containing all possible CrossVariograms
}
\references{
%% ~put references to the literature/web site here ~
Lark R.M. Robust estimation of the pseudo cross-variogram for cokriging soil properties. Eur J Soil Sci. 2002
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
\code{\link{n_pairs_opt},\link{plot.crossvariogram}}
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
