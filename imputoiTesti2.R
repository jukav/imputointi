# R-koodi puuttuvien arvojen imputaation testaamiseen eri matriiseilla


# Eri imputaatiofunktioita

# Funktio matriisin puuttuvien arvojen imputoitiin KNN-menetelmällä
# Etäisyyden mitta on R:n dist()-funktio
# m : Matriisi jota imputoidaan
# k : Lähimpien naapurien määrä
knnImpute <- function(m, k = 5){
  # Palautettava matriisi
  arvot <- m
  # Puuttuvien arvojen indeksit
  puuttuvat <- which(is.na(m), arr.ind = TRUE)
  
  # Matriisi rivien välisistä etäisyyksistä
  etaisyydet.m <- as.matrix(dist(m, method = "euclidian"))
  # Matriisi etäisyyksien järjestyksestä
  etaisyydet <- c()
  
  # Järjestetään etäisyydet
  for(i in 1:nrow(etaisyydet.m)){
    # Poistetaan omaa riviä vastaava indeksi
    etaisyydet <- rbind(etaisyydet, order(etaisyydet.m[i, -i]))
    # Korotetaan oman rivin jälkeisiä indeksejä
    v <- etaisyydet[i,]
    etaisyydet[i, v >= i] <- etaisyydet[i, v >= i] + 1
  }
  
  # Käydään puuttuvat arvot läpi
  for(i in 1:nrow(puuttuvat)){
    rivi = puuttuvat[i,1]
    sarake = puuttuvat[i,2]
    
    # Valitaan lähimmät k riviä
    jarjestys <- etaisyydet[rivi, 1:k]
    
    # Lasketaan keskiarvo ja lisätään matriisiin
    arvot[rivi, sarake] <- mean(m[jarjestys, sarake], na.rm = TRUE)
  }
  
  # Jos jäljellä on NaN-arvoja, korvataan ne nollalla
  arvot[is.na(arvot)] <- 0
  
  # Palautetaan saatu matriisi
  return(arvot)
}

# Funktio matriisin puuttuvien arvojen imputoitiin keskiarvolla
meanImpute <- function(m){
  # Lasketaan keskiarvo
  mean <- mean(as.matrix(m), na.rm = TRUE)
  # Korvataan arvot
  m[is.na(m)] <- mean
  # Jos jäljellä on NaN-arvoja, korvataan ne nollalla
  m[is.na(m)] <- 0
  # Palautetaan matriisi
  return(m)
  
}

# Funktio matriisin puuttuvien arvojen imputoitiin sarakekeskiarvoilla
colMeanImpute <- function(m){
  # Lasketaan sarakekeskiarvot
  means <- colMeans(x = m, na.rm = TRUE);
  # Korvataan arvot
  for(i in 1:ncol(m)){
    m[is.na(m[,i]), i] <- means[i]
  }
  # Jos jäljellä on NaN-arvoja, korvataan ne nollalla
  m[is.na(m)] <- 0
  # Palautetaan matriisi
  return(m)
}



# Funktio imputaation testaamiseen
# Saa syötteenä matriisin
# Valitsee tietyn prosenttimäärän satunnaisia arvoja, jotka poistetaan
# Tallentaa nämä arvot vektoriin indeksiensä kanssa, ja poistaa ne matriisista
# Ajaa sen jälkeen imputaatiofunktion matriisille
# Vertaa imputoituja arvoja alkuperäisiin, ja laskee niistä NMRSE-arvon
# Prosessi toistetaan k kertaa
# Saaduista arvoista lasketaan keskiarvo ja -hajonta
# Tämä toistetaan eri prosenteille
# Saaduista arvoista muodostetaan kuvaaja

testaa <- function(data, toistoja, k, funktio, nimi){
  r <- nrow(data)
  c <- ncol(data)
  
  # Virheet ja niiden keskihajonnat
  virheet <- c()
  hajonnat <- c()
  
  # Valmiiksi puuttuvat arvot
  puuttuvat <- is.na(data)
  
  # Suoritetaan testi eri prosenteilla
  for(i in 1:toistoja){
    # Virheet tietyllä prosentilla
    nrmse.v <- c()
    
    # Toistetaan testiä
    for(j in 1:k){
      # Laske poistettavien arvojen määrä
      n <- i * (1/(toistoja + 1)) * r * c 
      n <- ceiling(n)
      if(n > r*c)
        n <- r*c
      
      # Matriisi satunnaisesti valittuja alkioita
      poista <- vector(length = r*c)
      poista[sample(1:(r*c), n)] <- TRUE
      poista <- matrix(poista, r, c)
      # Hylätään valmiiksi puuttuvat arvot
      poista[puuttuvat] <- FALSE
      
      
      # Poista satunnaiset alkiot
      poistetut <- data[poista]
      data[poista] <- NaN
      
      # Suorita imputaatio
      data <- funktio(data)
      
      # Laske erotukset
      erotus <- poistetut - data[poista]
      # Laske nrmse
      nrmse <- sqrt((sum(erotus^2, na.rm = TRUE)/ length(erotus))) / 
        sqrt(mean(poistetut^2, na.rm = TRUE))
      
      # Lisää virheiden listaan
      nrmse.v <- c(nrmse.v, nrmse)
      
      # Palauta poistetut arvot
      data[poista] <- poistetut
      # Palauta valmiiksi puuttuvat arvot
      data[puuttuvat] <- NaN
    }
    # Lasketaan keskiarvo ja -hajonta tietylle prosentille.
    virheet <- c(virheet, mean(nrmse.v))
    hajonnat <- c(hajonnat, sd(nrmse.v))
  }
  
  # Muunna funktion nimi tekstiksi
  funktio.nimi <- deparse(substitute(funktio))
  
  # Piirrä kuvaaja tuloksista
  x <- 100*(1:toistoja)/(toistoja+1)
  y <- virheet
  
  plot(xlab = paste("Osuudet 1-",t,"% koko matriisin arvoista poistettu ja imputoitu",k,"toistokertaa"), ylab = "NRMSE",
       main = paste("Aineisto:", nimi, "\nFunktio:", funktio.nimi),
       x = x, y = y, xlim = c(0,100), ylim = c(0,1))
  
  # Lisää keskihajonnan viivat
  segments(x, y-hajonnat,x, y+hajonnat)
  epsilon = 0.5
  segments(x-epsilon,y-hajonnat,x+epsilon,y-hajonnat)
  segments(x-epsilon,y+hajonnat,x+epsilon,y+hajonnat)
}

# Testaus

# Luetaan testimatriisit
tst5050 <- readRDS("tst5050.rds")
tst100100 <- readRDS("tst100100.rds")
tst200198 <- readRDS("tst200198.rds")

# Siemennä satunnaislukugeneraattori
set.seed(0)

t <- 20 #poistettujen ja imputoitujen havaintojen prosenttiosuus matriisista
k <- 3 #imputointitoistojen lukumäärä

testaa(tst5050, t, k, meanImpute, "tst5050")
testaa(tst5050, t, k, colMeanImpute, "tst5050")
testaa(tst5050, t, k, knnImpute, "tst5050")

testaa(tst100100, t, k, meanImpute, "tst100100")
testaa(tst100100, t, k, colMeanImpute, "tst100100")
testaa(tst100100, t, k, knnImpute, "tst100100")

testaa(tst200198, t, k ,meanImpute, "tst200198")
testaa(tst200198, t, k, colMeanImpute, "tst200198")
testaa(tst200198, t, k, knnImpute, "tst200198")

