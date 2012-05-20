library(plyr);

countries <-
structure(1:51, .Label = c("Albania", "Andorra", "Armenia", "Austria", 
"Azerbaijan", "Belarus", "Belgium", "Bosnia-Herzegovina", "Bulgaria", 
"Croatia", "Cyprus", "Czech Republic", "Denmark", "Estonia", 
"Finland", "France", "Georgia", "Germany", "Greece", "Hungary", 
"Iceland", "Ireland", "Israel", "Italy", "Latvia", "Lithuania", 
"Luxembourg", "Macedonia", "Malta", "Moldova", "Monaco", "Montenegro", 
"Morocco", "Netherlands", "Norway", "Poland", "Portugal", "Romania", 
"Russia", "San Marino", "Serbia", "Serbia and Montenegro", "Slovakia", 
"Slovenia", "Spain", "Sweden", "Switzerland", "Turkey", "Ukraine", 
"United Kingdom", "Yugoslavia"), class = "factor");

n.countries <- length(countries);

cs <- function (x) c(0, cumsum(x));

gibbs <- function (prefs, clist, iters=100) {
  rank <- sample(clist);
  for (i in 1:iters) {
    n <- sample.int(length(rank), 1);
    item <- rank[n];
    rank <- rank[-n];
    pr <- prefs[rank] - prefs[item];
    pos <- sample.int(length(rank) + 1, 1, prob=exp(cs(pr) - rev(cs(rev(pr)))));
    rank <- append(rank, item, pos-1);
  }
  return(rank);
}

build.prefs <- function (data, alpha=1) {
  l <- list()
  for (giver in countries) {
    l[[giver]] <- f.giver(subset(data, Giver==giver), alpha);
  }
  return(l);
}

simulate.vote <- function(entrants, voter, prefs) {
  scores <- numeric(n.countries);
  names(scores) <- countries;
  pr <- prefs[,voter];
  rank <- gibbs(pr, entrants[entrants != voter]);
  scores[rank[1:10]] <- c(12, 10, 8, 7, 6, 5, 4, 3, 2, 1);
  return(scores[entrants]);
}

simulate.contest <- function (entrants, voters, prefs) {
  scores <- matrix(0, n.countries, n.countries, dimnames=list(countries, countries));
  for (v in voters) {
    scores[entrants, v] <- simulate.vote(entrants, v, prefs);
  }
  return(scores[entrants, voters]);
}

simulate.twostage <- function (semi1, semi1.v, semi2, semi2.v, prefs) {
  extras <- countries[c(semi1.v[!(semi1.v %in% semi1)], semi2.v[!(semi2.v %in% semi2)])];
  quals.1 <- semi1[rank(-rowSums(simulate.contest(semi1, semi1.v, prefs)), ties.method='random') <= 10];
  quals.2 <- semi2[rank(-rowSums(simulate.contest(semi2, semi2.v, prefs)), ties.method='random') <= 10];
  entrants <- countries[c(quals.1, quals.2, extras)];
  voters <- countries[c(semi1.v, semi2.v)];
  return(simulate.contest(entrants, voters, prefs));
}

comparisons <- function (data) {
  ddply(data, .(Year, Giver), function (d) {
      gts <- which(outer(d$Score, d$Score, '>'));
      Country.A <- rep(d$Country, nrow(d))[gts];
      Country.B <- rep(d$Country, each=nrow(d))[gts];
      return(data.frame(Country.A=Country.A, Country.B=Country.B));
    });
}

fit.jags <- function (data, ...) {
  library(rjags);
  
  comps <- comparisons(data);
  n <- nrow(comps);
  comps$Year <- substr(comps$Year, 1, 4);
  comps$Song.A <- paste(comps$Country.A, comps$Year, sep='/');
  comps$Song.B <- paste(comps$Country.B, comps$Year, sep='/');
  songs <- sort(unique(c(comps$Song.A, comps$Song.B)));
  comps$Song.A <- factor(comps$Song.A, songs);
  comps$Song.B <- factor(comps$Song.B, songs);
  cat("Fitting model:\n");
  cat(n, "comparisons,", length(songs), "songs", length(songs) + n.countries^2 + 2, "parameters\n");
  j <- jags.model('ev.bugs', list(n=n, n.countries=n.countries, n.songs=length(songs), y=rep(1, n), country.a=comps$Country.A, country.b=comps$Country.B, giver=comps$Giver, song.a=comps$Song.A, song.b=comps$Song.B), ...);
  return(list(songs=songs, model=j))
}

write.graph <- function(data, filename) {
  f <- file(filename, 'w');
  cat('digraph G {\n', file=f);
  cat('graph [splines=true, overlap=false, pack=true];\n', file=f);
  nodes <- unique(c(data$Giver, data$Country));
  for (n in nodes) {
    cat(n, '[label="', as.character(countries[n]), '"];\n', sep='', file=f);
  }
  for (i in 1:nrow(data)) {
    cat(data$Giver[i], '-> ', data$Country[i], ';\n', sep='', file=f);
  }
  cat('}\n', file=f);
  close(f);
}