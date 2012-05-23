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

gibbs <- function (prefs, clist, iters=50) {
  rank <- sample(clist);
  len <- length(rank);
  for (i in 1:iters) {
    n <- sample.int(length(rank), 1);
    item <- rank[n];
    rank <- rank[-n];
    pr <- prefs[rank] - prefs[item];
    pos <- sample.int(len, 1, prob=exp(cs(pr) - rev(cs(rev(pr)))));
    rank <- append(rank, item, pos-1);
  }
  return(rank);
}

simulate.vote <- function(entrants, voter, prefs) {
  scores <- numeric(n.countries);
  names(scores) <- countries;
  pr <- prefs[,voter];
  rank <- gibbs(pr, entrants[entrants != voter]);
  scores[rank[1:10]] <- c(12, 10, 8, 7, 6, 5, 4, 3, 2, 1);
  return(scores[entrants]);
}

simulate.contest <- function (entrants, voters, friends, quality) {
  names(quality) <- countries;
  colnames(friends) <- countries;
  rownames(friends) <- countries;
  prefs <- friends + quality;
  scores <- matrix(0, n.countries, n.countries, dimnames=list(countries, countries));
  for (v in voters) {
    scores[entrants, v] <- simulate.vote(entrants, v, prefs);
  }
  results <- list(
      scores=scores[entrants, voters],
      totals=rowSums(scores[entrants, voters]),
      entrants=entrants,
      voters=voters,
      friends=friends,
      quality=quality
    );
  return(results);
}

simulate.twostage <- function (semi1.e, semi1.v, semi2.e, semi2.v, friends, quality) {
  extras <- countries[c(semi1.v[!(semi1.v %in% semi1.e)], semi2.v[!(semi2.v %in% semi2.e)])];
  semi1 <- simulate.contest(semi1.e, semi1.v, friends, quality);
  quals.1 <- qualifiers(semi1);
  semi2 <- simulate.contest(semi2.e, semi2.v, friends, quality);
  quals.2 <- qualifiers(semi2);
  
  entrants <- countries[c(quals.1, quals.2, extras)];
  voters <- countries[c(semi1.v, semi2.v)];
  
  results <- simulate.contest(entrants, voters, friends, quality)
  output <- list(
      final=results,
      semis=list(semi1, semi2)
    );
  return(output);
}

qualifiers <- function (contest, n=10) {
  return(contest$entrants[rank(-contest$totals, ties.method='random') <= n]);
}

winner <- function (contest) {
  return(contest$entrants[which.max(contest$totals)]);
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

sample.model <- function (data, n=100, thin=10) {
  m <- fit.jags(data, inits=list(sigma.friend=2, sigma.song=2), n.adapt=100);
  samps <- jags.samples(m$model, c('b.friend', 'b.song', 'sigma.friend', 'sigma.song'), n*thin, thin=thin);
  id <- function(x) x;
  friends <- alply(samps$b.friend, c(3,4), id);
  songs.quality <- alply(samps$b.song, c(2,3), id);
  sigma.friend <- alply(samps$sigma.friend, c(2,3), id);
  sigma.song <- alply(samps$sigma.song, c(2,3), id);
  return(list(friends=friends, songs=m$songs, songs.quality=songs.quality, sigma.friend=sigma.friend, sigma.song=sigma.song));
}

random <- function (x) sample(x, 1)[[1]];

winning.probs <- function (results) {
  table(unlist(llply(results, function (x) winner(x$final)))) / length(results);
}

qualify.probs <- function (results) {
  v <- rowMeans(sapply(results, function (x) countries %in% x$final$entrants));
  names(v) <- countries;
  return(v);
}
simulate.final.12 <- function (samps) {
  i <- sample(length(samps$friends), 1);
  friends <- samps$friends[[i]]
  quality <- rnorm(n.countries, 0, samps$sigma.song[[i]]);
  simulate.twostage(semi12.1.e, semi12.1.v, semi12.2.e, semi12.2.v, friends, quality);
}

simulate.given.s1 <- function (results, s1, n) {
  extras <- factor(c('Azerbaijan', 'France', 'Germany', 'Italy', 'Spain', 'United Kingdom'), countries);
  s1 <- factor(s1, countries);
  qualities <- llply(s1, function (x) quality.if.qualify(results, x));
  names(qualities) <- s1;
  res <- list();
  for (i in 1:n) {
    s2 <- results[[i]]$semis[[2]];
    friends <- s2$friends;
    quality <- s2$quality;
    for (j in 1:10) {
      quality[s1[j]] <- sample(qualities[[j]], 1);
    }
    res[[i]] <- simulate.contest(countries[c(s1, qualifiers(s2), extras)], countries[c(semi12.1.v, semi12.2.v)], friends, quality);
  }
  return(res);
}
if.qualify <- function (results, country) {
  Filter(function (x) country %in% x$final$entrants, results);
}

quality.if.qualify <- function (results, country) {
  sapply(if.qualify(results, country), function (x) x$final$quality[country])
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