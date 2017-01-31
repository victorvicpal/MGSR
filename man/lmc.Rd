\name{lmc}
\alias{lmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Linear Model of Corregionalization
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Linear Model of Corregionalization fitting
}
\usage{
lmc(CV, fun, a, tol = 0.001, mode = "aut")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CV}{
%%     ~~Describe \code{CV} here~~
Crossvariogram resulting from \code{crossvariogram}
}
  \item{fun}{
%%     ~~Describe \code{fun} here~~
Vector of functions. Use 'Sph' for spherical, 'Gau' for Gaussian, 'Exp' for Exponential, 'Pow' for Power and 'Lin' for lineal.
}
  \item{a}{
%%     ~~Describe \code{a} here~~
Vector of ranges. Same length as 'fun'. 
If you use the Power function which is more suitable for subspatial fitting (0<a<2).
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
Tolerance
}
  \item{mode}{
%%     ~~Describe \code{mode} here~~
Automatic 'aut' or Manual 'man'. Manual mode allows for chosing sill values individually.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Power distribution is suitable only for small distances.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Object containing sill matrices and WSS results.
}
\references{
%% ~put references to the literature/web site here ~
Pelletier B, Dutilleul P, Larocque G, Fyles J.W. Fitting the Linear Model of Coregionaliza- tion by Generalized Least Squares. Math Geol. 2004;36(3):323-343.
}
\author{
%%  ~~who you are~~
Victor Vicente Palacios
}
\note{
%%  ~~further notes~~
Combination of multiple ranges and functions are also allowed.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{crossvariogram},\link{plot.crossvariogram}}
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
## plot
plot.crossvariogram(CV_vers,RES_vers)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
