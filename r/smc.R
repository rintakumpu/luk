# SMC-algoritmi
# Argumentit
# y: Havaintotaulukko (data.frame)
# N: haluttu partikkelien määrä
# t: aika-askelten määrä
# resampling: uudelleenotantastrategia
# proposal: ehdotusjakauma
# locator_n: kuinka montaa paikanninta käyttää. Vakiona käytetään kaikkia.
# plotting: Luodaanko partikkelikartat?
# max_plots: Jos asetettu, luodaan kartat ainoastaan ensimmäisellä max_plots aika-askeleella
# treshold: Jos asetettu, pysäytetään algoritmi, jos paikannusvirheen muutos alle tämän arvon.
# verbose: Tulostetaanko algoritmin tila?
# grid: Käytetäänkö alustuksessa satunnaista vai determinististä gridiä.

# Palauttaa list(e,x,w,p_maps,v)
# e: Sijaintivirhevektori (t-pituinen)
# x: x (partikkeli) matriisi
# w: w (paino) matriisi
# p_maps: partikkelikartat listana
# v: keskivirhevektori (t-pituinen)

smc <- function(y, N, t = NA, resampling=c("adaptive", "none", "each-step"), proposal=c("likelihood", "prior"), locator_n=NA, plotting=T, max_plots=NA, treshold=NA, verbose=F, grid=c("random", "deterministic")) {
  
  # Sovitetaan vektorimuotoiset arvot
  grid <- match.arg(grid)
  proposal <- match.arg(proposal)
  resampling <- match.arg(resampling)
  
  # Muutetaan syötedata data.table-muotoona
  y <- as.data.table(y)
  
  # Jos locator_na ja/tai max_plots ei määritetty, käytetään datasta johdettuja vakioarvoja
  if(is.na(locator_n)) {
    locator_n <- length(unique(y$locator_mac))
  }
  if(is.na(max_plots)) {
    max_plots <- t
  }
  if(is.na(t)) {
    t <- length(unique(y$ts))
  }
  
  if(verbose) {
    cat(paste("Running SMC with N: ", N, "; t: ", t, "\nResampling strategy: ", resampling, "; Proposal distribution: ", proposal, "\nPrior grid distribution: ", grid, "; Convergence treshold: ", treshold, "\n", sep=""))
  }
  
  # Luodaan paino-, malli- ja havaintomatriisit
  # arvojen N ja t perusteella
  w <- expand.grid(0:t,1:N)
  # k on aika-askel, i partikkeli
  colnames(w) <- c("k","i")
  x <- w
  
  # Lisätään matrsiisispesifit arvot
  w$w <- rep(NA, nrow(w))
  
  # Alkupainot
  w[w$k==0, "w"] <- rep(log(1/N), N)
  
  # Virhe- ja varianssivektori
  e <- rep(NA, (t-1))
  v <- rep(NA, (t-1))
  
  # Alustetaan algoritmi x_0, tasajakaumasta
  # ei käytetä koko pohjapiirustusta, ainoastaan paikantimien ympäröimää aluetta
  # Deterministinen/satunnainen gridi
  if(grid=="deterministic") {
    lon <- seq(from = min(locators$lon), to = max(locators$lon), length.out = sqrt(N))
    lat <- seq(from = min(locators$lat), to = max(locators$lat), length.out = sqrt(N))
  } else {
    lat <- runif(n = sqrt(N), min = min(locators$lat), max = max(locators$lat))
    lon <- runif(n = sqrt(N), min = min(locators$lon), max = max(locators$lon))
  }
  x_0 <- expand.grid(lat,lon)
  x_0 <- x_0[1:N,]
  colnames(x_0) <- c("lat", "lon")
  # Tallenetaan alkujakauma x-matriisiin
  x[x$k==0, "lat"] <- x_0$lat
  x[x$k==0, "lon"] <- x_0$lon
  
  # Alkuasetelman partikkelikartta
  p_maps <- as.list(rep(NA, (t)))
  if(plotting) {
    locator_coordinates <- locators %>% dplyr::select(lon, lat)
    p_maps[[t]] <- ggplot(x_0, aes(x=lon, y=lat)) + geom_point(shape=15, size=1.9) + annotate("point", x=22.2948889620602, y=60.4481932096263, color="steelblue") + geom_point(data=locator_coordinates, aes(x=lon, y=lat), color="tomato") + ylab("") + xlab(paste("k:", 1, "; N_eff: ", N, sep="")) + theme(axis.title.x = element_text(size = rel(0.7)), axis.text.x=element_blank(), axis.text.y=element_blank())
  }
  
  # Aloitetaan algoritmi
  k <- 1
  for(k_ts in unique(y$ts)) {
    # Tulostetaan aika
    if(verbose) print(paste("Time step ", k, ": ", k_ts, sep=""))
    # Valitaan aika-askeleen havainnot
    y_subset <- y[y$ts==k_ts]
    
    # Iterointi / tärkeytysotanta
    # Lasketaan normalisointivakio, ensin tähän liittyvä vektori
    log_p <- p_yk_xk(y_subset, x[x$k==(k-1),"lat"], x[x$k==(k-1),"lon"], locator_n=locator_n)
    c_k_vector <- w[w$k==(k-1),"w"] + log_p
    # Ja normalisointikerroin vektorin perusteella
    c_k <- logSumExp(c_k_vector)
    # Normalisoidaan partikkelien painot niin että niiden summa on 1
    # (tai tässä tapauksessa, logaritmoituna 0)
    w[w$k==(k-1),"w"] <- c_k_vector-c_k
    
    # Lasketaan effektiivinen otoskoko
    N_eff <- 1/exp(logSumExp(w[w$k==(k-1),"w"]+w[w$k==(k-1),"w"]))
    
    # Luodaan ja tallenetaan partikkelikartat
    if(plotting & (k <= max_plots)) {
      non_zero_weights <- w[w$w!=-Inf & w$k==(k-1),"i"]
      map_df <- cbind(x[x$k==(k-1) & x$i %in% non_zero_weights,3:4],w[w$k==(k-1) & w$i %in% non_zero_weights,"w"])
      colnames(map_df) <- c("lat", "lon", "w")
      # Heitetään pois pienet painot ja luodaan kartta
      p_maps[[k]] <- ggplot(map_df, aes(x=lon, y=lat), alpha=w) + geom_point(size=1.9) + annotate("point", x=22.2948889620602, y=60.4481932096263, color="steelblue", shape=15) + geom_point(data=locator_coordinates, aes(x=lon, y=lat), color="tomato") + ylab("") + xlab(paste("k: ", k, "; N_eff: ", round(N_eff,0), sep="")) + theme(axis.title.x = element_text(size = rel(0.7)), axis.text.x=element_blank(), axis.text.y=element_blank())
    }
    
    # Uudelleenotanta
    if(resampling!="none") {
      if(((resampling=="adaptive") & (N_eff < 2*N/3)) | (resampling=="each-step")) {
        if(verbose) print(paste("Effetive sample size: ", N_eff, ". Resampling.", sep=""))
        resample_index <- sample(x=1:N, size=N, replace=T, prob=exp(w[w$k==(k-1),"w"]))
        x[x$k==(k-1),] <- x[which(x$k==(k-1)),][resample_index,]
        x[x$k==(k-1),"i"] <- 1:N
        w[w$k==(k-1),"w"] <- log(1/N)
      }
    }
    
    # Lasketaan painotettu paikannusvirhe (so. posteriorijakauma ajanhetkenä k)
    # Huom: Tämä on eri asia kuin posteriori 1:k mutta on tarpeeksi hyvä estimaatti.
    # Lasketaan etäisyys jokaisen partikkelin ja fokulaattorin välillä
    e_vector <- c()
    e_vector <- pointDistance(matrix(c(x[x$k==(k-1), "lon"],x[x$k==(k-1), "lat"]), ncol=2),c(22.2948889620602, 60.4481932096263), lonlat = T, allpairs = T)
    # Ja painotetaan partikkelien painoilla
    e[k] <- weighted.mean(e_vector, w[w$k==(k-1),"w"])
    if(verbose) print(paste("Mean error: ", e[k], sep=""))
    if(!is.na(treshold) & k>1) {
      if(abs(e[k]-e[(k-1)])<treshold) {
        if(verbose) print("Algorithm converged.")
        break
      }
    }
    
    # Lasketaan "varianssiestimaatti"
    # Lasketaan partikkelien keskus painotettuna keskiarvona
    c_lon <- weighted.mean(x[x$k==(k-1), "lon"], w[w$k==(k-1),"w"])
    c_lat <- weighted.mean(x[x$k==(k-1), "lat"], w[w$k==(k-1),"w"])
    # Euklidinen normi
    z <- sqrt((x[x$k==(k-1), "lon"]-c_lon)^2+(x[x$k==(k-1), "lat"]-c_lat)^2)
    # Lasketaan "varianssi"
    v[k] <- mean(z)
    
    # Aikapäivitys
    if(k < t) {
      # Prioriotanta (Foculator ei liiku, virhetermi siten nolla)
      x[x$k==k & x$i, "lat"] <- x[x$k==(k-1), "lat"] + rnorm(N,sd=0.000000)
      x[x$k==k & x$i, "lon"] <- x[x$k==(k-1), "lon"] + rnorm(N,sd=0.000000)
      # Päivitetään painot
      if(proposal=="prior") {
        # Jos käytetään prioria, w = w-1 * p(y|x)
        w[w$k==k,"w"] <- w[w$k==(k-1),"w"]  + log_p
      }
      if(proposal=="likelihood") {
        # Jos käytetään uskottavuutta, w = w-1 * p(x|x) = w-1 * I = w-1 (Foculator ei liiku)
        w[w$k==k,"w"] <- w[w$k==(k-1),"w"]
      }
      # Päivitetään aika
      k <- k + 1
    } else {
      break
    }
  }  
  
  return(list(e,x,w,p_maps,v))
}