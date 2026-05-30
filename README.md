# VM-gjetteren

En liten Shiny-app som anslår forventet kampresultat i mål mellom to landslag ut fra FIFA-rankingen. Appen bruker FIFAs egen Elo-baserte formel til å regne poengforskjellen om til en forventet vinnersannsynlighet, og oversetter denne til mål ved å modellere hvert lags scoring som en Poisson-fordeling. Resultatet er en full sannsynlighetsfordeling over sluttresultater – med forventede mål, sannsynlighet for seier/uavgjort/tap, fordeling av målforskjellen (Skellam) og mulighet for å simulere enkeltkamper. Laget for nøytral bane (VM 2026), uten hjemmebanefordel. Lagdata leses fra `teams.csv` (kolonnene `Name,Points,Code`).

## Kjøre lokalt

```r
install.packages(c("shiny", "ggplot2", "viridisLite"))
shiny::runApp("app.R")
```
