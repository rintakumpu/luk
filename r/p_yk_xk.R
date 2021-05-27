# Uskottavuusfunktio p(y_k|x_k)
# Argumentit:
# y_input: Havainnot (data.frame), joista uskottavuus lasketaan
# x_input_lat: Partikkelien y-koordinaatit vektorina
# x_input_lon: Partikkelien x-koordinaatit vektorina
# locator_n: Kuinka monta (parasta) paikkaninta käytetään uskottavuuden laskentaan
# log_scale: Palautetaan uskottavuus logaritmisena

p_yk_xk <- function(y_input, x_input_lat, x_input_lon, locator_n=NA, log_scale=T) {
  
  # Partikkelien määrä
  Np <- length(x_input_lon)
  
  # Jos locator_n-muuttujaa ei ole asetettu, käytetään kaikkia
  # datassa olevia paikantimia
  if(is.na(locator_n)) {
    locator_n <- n_distinct(y_input$locator_mac)
  }
  
  # Valitaan vahvimmat havainnot
  y_input_original <- y_input
  y_input <- as_tibble(y_input)
  if(locator_n > n_distinct(y_input$locator_mac)) locator_n <- n_distinct(y_input$locator_mac)
  
  # Käytetään keskiarvoa
  top_locators <- y_input %>% group_by(locator_mac) %>% summarise(snr_max=max(ss_snr), mean_azimuth=median(converted_azimuth), sd=sd(converted_azimuth)) %>% top_n(locator_n, wt=snr_max)
  # Tai vain havaintoja, joilla paras signaali-kohinasuhde
  #top_locators <- y_input %>% group_by(locator_mac) %>% summarise(snr_max_index=which.max(ss_snr), snr_max=max(ss_snr), mean_azimuth=converted_azimuth[snr_max_index], sd=sd(converted_azimuth)) %>% slice(1) %>% ungroup() %>% top_n(locator_n, wt=snr_max)
  
  # Jos varianssia ei voida laskea, asetetaan keskihajonta vakioarvoon (tässä 1)
  top_locators$sd <- replace_na(top_locators$sd, 1)
  
  # Konvertoidaan data.frame tibble-taulukoksi
  y_input <- top_locators %>% inner_join(y_input %>% dplyr::select(lat, lon, locator_mac) %>% filter(locator_mac %in% top_locators$locator_mac) %>% group_by(locator_mac) %>% slice(1), by=c("locator_mac"))
  
  # Luodaan lat-lon-muuttujista yksi merkkijono
  y_input <- y_input %>% unite(latlon_y, c("lat", "lon"), sep=", ", remove = F)
  x_input <- as_tibble(matrix(ncol=2, nrow=Np, dimnames = list(c(), c("lat", "lon"))))
  x_input$lat <- x_input_lat
  x_input$lon <- x_input_lon
  x_input <- x_input %>% unite(latlon_x, c("lat", "lon"), sep=", ", remove = F)
  
  # Luodaan taulukko, jossa on kaikki partikkeli-paikannin-kombinaatiot
  xy_input <- expand_grid(x_input$latlon_x, y_input$latlon_y)
  colnames(xy_input) <- c("latlon_x", "latlon_y")
  xy_input <- xy_input %>% left_join(y_input, by=c("latlon_y")) %>% rename(lat_a = lat, lon_a = lon) %>% separate(latlon_x, c("lat_b", "lon_b"), remove = F, sep = ", ")
  
  # Lasketaan partikkelien ja paikantimien väliset kulmat
  xy_input <- xy_input %>% mutate(angle = ((atan2((as.numeric(lon_b)-as.numeric(lon_a)), (as.numeric(lat_b)-as.numeric(lat_a))) * (180/pi) + 360) %% 360))
  
  # Sekä keskivirhe (laskettu kulma vs. havaintokulma)
  xy_input <- xy_input %>% mutate(diff = pmin(abs(angle-mean_azimuth), 360-abs(angle-mean_azimuth)), se = diff^2)
  
  # Lasketaan uskottavuus
  xy_input <- xy_input %>% mutate(L = exp(-(se/(2*sd^2)))) %>% mutate(l = log(L))
  
  # Uskottavuus logaritmina, joten summataan, jotta saadaan uskottavuuksien tulo
  # havaintojen järjestys katoaa, joten palautetaan se alkuperäisestä taulukosta
  x_order <- x_input[, "latlon_x"]
  l <- (x_order %>% left_join(xy_input %>% group_by(latlon_x) %>% summarise(logSum=sum(l)) %>% ungroup(), by=c("latlon_x")))$logSum
  
  # Jos uskottavuus on 0 (ts. -Inf), tiputetaan yksi paikannin havainnoista
  # ja kutsutaan uskottavuusfunktiota p_yx_xk rekursiivisesti
  if(length(unique(l))==1 & unique(l)[1]==-Inf) {
    ln <- locator_n - 1
    if(ln == 0) {
      # Paikantimet loppuivat, rekursio pohjalla. Palataan.
      return(l)
    } else {
      l <- p_yk_xk(y_input_original, x_input_lat, x_input_lon, locator_n=ln)  
    }
  }
  
  if(!log_scale) {
    L <- exp(l)
    return(L)
  }
  
  return(l)
}