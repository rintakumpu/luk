# Triangulaatiopaikannus käyttäen ToTal-algoritmia
# Argumentit
# observations: Havainnot, joille paikannus suoritetaan
# max_jitter: Maksimi jitter-arvo, joka sallitaan havainnoille
# min_snr: Minimi signaali-kohinasuhdearvo, joka sallitaan havainnoille
# angle_statistic: Millä metodilla kulmat valitaan 
#   "max_snr": Valitaan kustakin paikantimesta signaali-kohinasuhteeltaan paras kulma
#   "mean": Kunkin paikantimen kulmien aritmeettinen keskiarvo
#   "median": " mediaani
total <- function(observations, max_jitter = NA, min_snr = NA, angle_statistic = "max_snr") {
  
  # Suodatetaan havainnot neliösummilla
  if(!is.na(max_jitter)) {
    observations <- observations %>% filter(ss_jitter < max_jitter)
  }
  if(!is.na(min_snr)) {
    observations <- observations %>% filter(ss_snr > min_snr)
  }
  
  # Katsotaan onko havaintoja tarpeeksi
  if(nrow(observations)<3) { 
    return(list("not enough angles", NA, NA, 0, 0)) } 
  
  # 0. Valmistellaan data
  # Valtaan kulmat keskiarvona tai mediaanina
  if(angle_statistic == "mean") {
    observations_subset <- observations %>% group_by(locator_mac) %>% summarise(lat=lat, lon=lon, azimuth_error=azimuth_error, snr_max=max(ss_snr), mean_azimuth=mean(converted_azimuth), sd=sd(converted_azimuth)) %>% slice(1) %>% ungroup() %>% top_n(3, wt=snr_max)
  }
  if(angle_statistic == "median") {
    observations_subset <- observations %>% group_by(locator_mac) %>% summarise(lat=lat, lon=lon, azimuth_error=azimuth_error, snr_max=max(ss_snr), mean_azimuth=mean(converted_azimuth), sd=sd(converted_azimuth)) %>% slice(1) %>% ungroup() %>% top_n(3, wt=snr_max)
  }
  
  # Tai valitaan jokaisesta paikantimesta SNR-arvoltaan (signaali-kohinasuhteeltaan) paras arvo
  # NOTE: Tämä toimii yleisesti ottaen paremmin kuin keskilukujen käyttö.
  if(angle_statistic == "max_snr") {
    observations_subset <- observations %>% group_by(locator_mac) %>% summarise(max_snr_index = which.max(ss_snr), lat=lat, lon=lon, azimuth_error=azimuth_error, snr_max=max(ss_snr), mean_azimuth=converted_azimuth[max_snr_index], sd=sd(converted_azimuth)) %>% slice(1) %>% ungroup() %>% top_n(3, wt=snr_max)
  }
  
  # Lasketaan keskivirhe (käytetään vain testauksessa, tässä NA)
  mean_azimuth_error <- NA 
  
  # Jos paikantimia <3 palautetaan virhe
  valid_angles <- nrow(observations_subset)
  if(valid_angles<3) { return(list("not enough valid angles", mean_azimuth_error, valid_angles, 0, 0)) }   
  
  mean_snr <- mean(observations_subset$snr_max)
  # Paikantimien koordinaatit
  l1 <- c(observations_subset$lon[1], observations_subset$lat[1])
  l2 <- c(observations_subset$lon[2], observations_subset$lat[2])
  l3 <- c(observations_subset$lon[3], observations_subset$lat[3])
  # Kulmat
  a1 <- as.numeric(observations_subset[1,"mean_azimuth"])
  a2 <- as.numeric(observations_subset[2,"mean_azimuth"])
  a3 <- as.numeric(observations_subset[3,"mean_azimuth"])
  # Jos kaksi kulmaa ovat samat, palautetaan virhe.
  # Tämä olisi mahdollista korjata pudottamalla havainto ja valitsemalla seuraavaksi 
  # paras paikannin. Kulmat menevät samoiksi kuitenkin niin harvoin, että tätä
  # ei toteuteta tässä.
  if((a2==a1)|(a2==a3)|(a1==a3)) {
    return(list("identical angles", mean_azimuth_error, valid_angles, mean_snr, 0)) } 
  
  # Lasketaan vastakkainen kulma (foculator-paikannin), konvertoidaan radiaaneiksi
  a1 <- ((a1+180)%%360)*(pi/180)
  a2 <- ((a2+180)%%360)*(pi/180)
  a3 <- ((a3+180)%%360)*(pi/180)
  
  # Valmistelut tehty,
  # Varsinainen ToTal-algoritmi (http://www.telecom.ulg.ac.be/publi/publications/pierlot/Pierlot2014ANewThree)
  # 1. Lasketaan muunnetut paikanninkoordinaatit
  x_prime_1 <- l1[1]-l2[1]
  y_prime_1 <- l1[2]-l2[2]
  x_prime_3 <- l3[1]-l2[1]
  y_prime_3 <- l3[2]-l2[2]
  
  # 2. Lasketaan kotangentit
  T12 <- cot(a2-a1)
  T23 <- cot(a3-a2)
  T31 <- (1-T12*T23)/(T12+T23)
  
  # 3. Lasketaan muunnetut ympyröiden keskipisteet:
  x_prime_12 <- x_prime_1 + T12*y_prime_1
  y_prime_12 <- y_prime_1 - T12*x_prime_1
  x_prime_23 <- x_prime_3 - T23*y_prime_3
  y_prime_23 <- y_prime_3 + T23*x_prime_3
  x_prime_31 <- (x_prime_3+x_prime_1)+T31*(y_prime_3-y_prime_1)
  y_prime_31 <- (y_prime_3+y_prime_1)-T31*(x_prime_3-x_prime_1)
  
  # 4. Lasketaan k_prime_31
  k_prime_31 <- (x_prime_1*x_prime_3)+(y_prime_1*y_prime_3)+T31*(x_prime_1*y_prime_3-x_prime_3*y_prime_1)
  
  # 5 Lasketaan D (jos D = 0 palautetaan virhe)
  D <- ((x_prime_12-x_prime_23)*(y_prime_23-y_prime_31))-((y_prime_12-y_prime_23)*(x_prime_23-x_prime_31))
  if(D==0){ 
    return(list("D equals zero", mean_azimuth_error, valid_angles, mean_snr, 0)) 
  }
  
  # Palautetaan sijainti + paikannusvirhe
  x_R <- l2[1]+(k_prime_31*(y_prime_12-y_prime_23))/D
  y_R <- l2[2]+(k_prime_31*(x_prime_23-x_prime_12))/D
  return(list(c(x_R,y_R), mean_azimuth_error, valid_angles, mean_snr, D))
}