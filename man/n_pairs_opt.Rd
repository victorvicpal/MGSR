\name{n_pairs_opt}
\alias{n_pairs_opt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Help Function.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Help Function to choose the proper number of pairs for the crossvariogram.
}
\usage{
n_pairs_opt(coord, values, n_min = 12, n_max = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coord}{
%%     ~~Describe \code{coord} here~~
Dataframe containing the coordinates from Factorial Techniques.
}
  \item{values}{
%%     ~~Describe \code{values} here~~
Dataframe containing values from data.
}
  \item{n_min}{
%%     ~~Describe \code{n_min} here~~
minimal number of pairs to be checked
}
  \item{n_max}{
%%     ~~Describe \code{n_max} here~~
maximal number of pairs to be checked
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Plotting 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Object containing all (max - min) values of pairs and the suggested one.
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
This function just helps, n values can be freely choosen.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{crossvariogram}}
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
n_pairs <- n_pairs_opt(as.data.frame(PC_train$scores[,1:2]),as.data.frame(Versicolor_st),6,15)
plot(n_pairs$dif_pairs)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{}% use one of  RShowDoc("KEYWORDS")
¡
