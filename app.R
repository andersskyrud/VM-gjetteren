# =============================================================================
# Forventet kampresultat fra FIFA-ranking
# -----------------------------------------------------------------------------
# Tar to lags FIFA-ranking og forventet totalmål for en jevn kamp (tau_0), og
# returnerer forventede mål, resultatsannsynligheter, fordelingen for
# målforskjellen (Skellam) samt simulerte enkeltresultater.
#
# Avhengigheter: shiny, ggplot2, viridisLite
#   install.packages(c("shiny", "ggplot2", "viridisLite"))
# Kjør lokalt:  shiny::runApp("app.R")
#
# Lagdata leses fra teams.csv i samme mappe, med kolonnene:
#   Name,Points,Code     (f.eks.  Norge,1620,NO)
# Code er to-bokstavs ISO-landkode og brukes til å lage emoji-flagg.
# =============================================================================

library(shiny)
library(ggplot2)
library(viridisLite)

# ---- Lagdata ----------------------------------------------------------------

# Gjør en to-bokstavs ISO-kode (f.eks. "NO") om til et emoji-flagg.
# Hver bokstav erstattes av tilsvarende "regional indicator"-tegn.
iso_to_flag <- function(code) {
  code <- toupper(trimws(code))
  if (is.na(code) || nchar(code) != 2 || grepl("[^A-Z]", code)) return("")
  pts <- utf8ToInt(code) - utf8ToInt("A") + 0x1F1E6
  intToUtf8(pts)
}

# Les lagdata, eller fall tilbake på et lite innebygd sett om fila mangler.
load_teams <- function(path = "teams.csv") {
  if (file.exists(path)) {
    df <- utils::read.csv(path, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  } else {
    warning("teams.csv ikke funnet - bruker innebygd reserveliste.")
    df <- data.frame(
      Name   = c("Brasil", "Japan", "Norge", "Argentina"),
      Points = c(1761.16, 1660.43, 1550.94, 1874.81),
      Code   = c("BR", "JP", "NO", "AR"),
      stringsAsFactors = FALSE
    )
  }
  df <- df[order(-df$Points), ]
  df$Flag  <- vapply(df$Code, iso_to_flag, character(1))
  # Visningsetikett: flagg + navn + FIFA-poeng. selectizeInput gjengir emoji-flagg.
  df$Label <- trimws(sprintf("%s %s (%d p.)", df$Flag, df$Name, round(df$Points)))
  df
}

teams <- load_teams()

# ---- Modellfunksjoner -------------------------------------------------------

# Resultatmatrise fra to uavhengige Poisson-rater.
# Rad i = mål for lag 1 (0:maxgoals), kolonne j = mål for lag 2.
score_matrix <- function(lambda1, lambda2, maxgoals = 25) {
  outer(dpois(0:maxgoals, lambda1), dpois(0:maxgoals, lambda2))
}

# Forventet poeng for lag 1 (seier = 1, uavgjort = 0,5) gitt resultatmatrisen.
expected_points_1 <- function(M) {
  sum(M[lower.tri(M)]) + 0.5 * sum(diag(M))
}

# FIFAs forventede resultat (vinnersannsynlighet) for lag 1.
# Skala = 600 (FIFA), ingen hjemmebanefordel (nøytral bane).
fifa_we <- function(points1, points2, scale = 600) {
  1 / (10^(-(points1 - points2) / scale) + 1)
}

# Multiplikativ parametrisering
# lambda1 = (tau_match/2)*exp(s), lambda2 = (tau_match/2)*exp(-s),
# der tau_match = tau0 * exp(kappa * |s|) vokser med styrkeforskjellen.
# tau0 er altså totalmål i en JEVN kamp; snittet over alle kamper blir høyere.
solve_supremacy <- function(We, tau0, kappa = 0.7) {
  We <- min(max(We, 1e-6), 1 - 1e-6)
  s_max <- 3.5
  # Største rate ved intervallets endepunkt, for å dimensjonere rutenettet.
  tau_hi   <- tau0 * exp(kappa * s_max)
  lam_hi   <- (tau_hi / 2) * exp(s_max)
  maxgoals <- ceiling(lam_hi + 6 * sqrt(lam_hi) + 5)
  f <- function(s) {
    tau_m <- tau0 * exp(kappa * abs(s))
    l1 <- (tau_m / 2) * exp(s)
    l2 <- (tau_m / 2) * exp(-s)
    expected_points_1(score_matrix(l1, l2, maxgoals)) - We
  }
  uniroot(f, interval = c(-s_max, s_max))$root
}

# Fordeling for målforskjell (lag1 - lag2), kollapset fra resultatmatrisen.
margin_distribution <- function(M) {
  maxg <- nrow(M) - 1
  ds <- -maxg:maxg
  pr <- vapply(ds, function(d) {
    i <- seq_len(nrow(M)); j <- i - d
    keep <- j >= 1 & j <= ncol(M)
    sum(M[cbind(i[keep], j[keep])])
  }, numeric(1))
  data.frame(margin = ds, prob = pr)
}

# Samlet prediksjon.
predict_match <- function(points1, points2, tau0, kappa = 0.7) {
  We    <- fifa_we(points1, points2)
  s     <- solve_supremacy(We, tau0, kappa)
  tau_m <- tau0 * exp(kappa * abs(s))
  l1    <- (tau_m / 2) * exp(s)
  l2    <- (tau_m / 2) * exp(-s)
  lam_hi   <- max(l1, l2)
  maxgoals <- ceiling(lam_hi + 6 * sqrt(lam_hi) + 5)
  M  <- score_matrix(l1, l2, maxgoals)
  list(
    We = We, lambda1 = l1, lambda2 = l2, total = l1 + l2,
    p_win1 = sum(M[lower.tri(M)]),
    p_draw = sum(diag(M)),
    p_win2 = sum(M[upper.tri(M)]),
    matrix = M, margins = margin_distribution(M)
  )
}

# ---- Brukergrensesnitt ------------------------------------------------------

ui <- fluidPage(
  titlePanel("Forventet kampresultat fra FIFA-ranking"),
  sidebarLayout(
    sidebarPanel(
      selectizeInput("lag1", "Lag 1",
                     choices = stats::setNames(teams$Name, teams$Label),
                     selected = "Brasil"),
      selectizeInput("lag2", "Lag 2",
                     choices = stats::setNames(teams$Name, teams$Label),
                     selected = "Japan"),
      tags$hr(),
      sliderInput("tau0", "Forventet antall mål i en jevn kamp (τ₀)",
                  min = 1.5, max = 3.5, value = 2.69, step = 0.01),
      helpText("Standardverdi 2,69 = gjennomsnittlig antall mål per kamp i ",
               "VM 2022 i Qatar."),
      tags$hr(),
      actionButton("simuler", "Simuler kamp", class = "btn-primary"),
      actionButton("simuler1000", "Simuler 1000 kamper"),
      actionButton("nullstill", "Nullstill")
    ),
    mainPanel(
      h4("Simulert enkeltresultat"),
      div(style = "font-size: 1.4em; font-weight: bold; margin-bottom: 6px;",
          textOutput("simresultat")),
      textOutput("simtally"),
      tags$hr(),
      # Grafer til venstre, tabell til høyre.
      fluidRow(
        column(
          width = 8,
          h4("Sannsynlighet for målforskjell"),
          helpText("Søyler = fordeling. Punkter = simulert andel."),
          plotOutput("marginplot", height = "300px"),
          h4("Resultatmatrise (sannsynlighet i %)"),
          plotOutput("matriseplot", height = "340px")
        ),
        column(
          width = 4,
          h4("Forventa resultat"),
          tableOutput("sammendrag")
        )
      ),
      tags$hr(),
      wellPanel(
        h4("Metode – steg for steg"),
        tags$ol(
          tags$li("- Hvert lag har et antall poeng i FIFA-rankingen. Jo større ",
                  "poengforskjell, desto sterkere favoritt."),
          tags$li("- FIFAs egen formel gjør poengforskjellen om til en ",
                  "forventet poengsum for det sterkeste laget: et tall mellom 0 og 1, ",
                  "der man forventer seier ved verdier nære 1 og uavgjort ved 0,5.  ",
                  "Like sterke lag gir 0,5;en stor favoritt nærmer seg 1. (Formelen er ",
                  "Wₑ = 1 / (10^(−Δ/600) + 1), der Δ er poengforskjellen.)"),
          tags$li("- Dette tallet sier hvem som bør vinne, men ikke hvor mange ",
                  "mål det blir. For å komme fra «hvem vinner» til «hvor mange mål» ",
                  "antar vi at mål kommer tilfeldig og følger en ",
                  "en Poisson-fordeling – styrt av ett tall per lag: forventet ",
                  "antall mål (λ)."),
          tags$li("- Vi starter med at begge lag forventes å score like mange mål ",
                  "i en jevn kamp (τ₀ totalt, altså τ₀/2 hver). Deretter ",
                  "vrir vi styrken: favorittens forventning ganges opp og ",
                  "outsiderens ganges ned, helt til modellens forventa poeng ",
                  "stemmer med FIFAs tall fra steg 2. Fordi vi ganger (i stedet for ",
                  "å flytte mål fra det ene laget til det andre), kan det totale ",
                  "målantallet vokse når forskjellen er stor."),
          tags$li("- Nå har hvert lag sitt forventa måltall. Vi regner ut ",
                  "sannsynligheten for hvert mulige sluttresultat (0-0, 1–0, ...), ",
                  "og summerer dem til sannsynlighet for seier, uavgjort og tap, ",
                  "samt fordelingen av målforskjellen.")
        )
      ),
      wellPanel(
        h4("Forutsetninger og antakelser"),
        p("Resultatene hviler på noen bevisste forenklinger som er verdt å ",
          "være klar over:",
          "(1) Mål trekkes for hvert lag fra en Poisson-fordeling – ",
          "antall mål ett lag scorer påvirker ikke det andre",
          "(2) Forventet mål i en jevn kamp (τ₀) settes av brukeren og ",
          "antas likt for alle kamper, uavhengig av hvilke lag som møtes. ",
          "(3) Ingen hjemmebanefordel, siden VM spilles på nøytral bane. ",
          "(4) Målene utledes direkte fra FIFA-poeng. Med andre ord modelleres ikke ",
          "historiske kampdata, som kunne vært et alternativ, og måltalla har derfor ingen ",
          "statistisk usikkerhet – bare selve resultatet varierer fra kamp til ",
          "kamp.",
          "(5) Utregninga antar også at det skåres flere mål i ujevne kamper. Jo større ",
               "rankingforskjell mellom laga - jo flere mål og vice versa. ")
      )
    )
  )
)

# ---- Server -----------------------------------------------------------------

server <- function(input, output, session) {

  # Slå opp poeng og visningsnavn for de to valgte lagene.
  navn1   <- reactive(input$lag1)
  navn2   <- reactive(input$lag2)
  poeng1  <- reactive(teams$Points[match(input$lag1, teams$Name)])
  poeng2  <- reactive(teams$Points[match(input$lag2, teams$Name)])

  pred <- reactive({
    validate(
      need(input$lag1 != input$lag2, "Velg to forskjellige lag.")
    )
    predict_match(poeng1(), poeng2(), input$tau0)
  })

  # Lager for simulerte resultater; nullstilles når inndata endres.
  sims <- reactiveVal(data.frame(g1 = integer(), g2 = integer()))
  observeEvent(list(input$lag1, input$lag2, input$tau0), {
    sims(data.frame(g1 = integer(), g2 = integer()))
  }, ignoreInit = TRUE)

  observeEvent(input$simuler, {
    p <- pred()
    sims(rbind(sims(),
               data.frame(g1 = rpois(1, p$lambda1), g2 = rpois(1, p$lambda2))))
  })
  observeEvent(input$simuler1000, {
    p <- pred()
    sims(rbind(sims(),
               data.frame(g1 = rpois(1000, p$lambda1), g2 = rpois(1000, p$lambda2))))
  })
  observeEvent(input$nullstill, { sims(data.frame(g1 = integer(), g2 = integer())) })

  output$sammendrag <- renderTable({
    p <- pred()
    idx <- which(p$matrix == max(p$matrix), arr.ind = TRUE)[1, ]
    data.frame(
      "Størrelse" = c(
        paste0("Forventa poeng (", navn1(), ")"),
        paste0("Forventa mål (", navn1(), "), λ₁"),
        paste0("Forventa mål (", navn2(), "), λ₂"),
        "Forventa totalmål",
        paste0("P(seier ", navn1(), ")"),
        "P(uavgjort)",
        paste0("P(seier ", navn2(), ")"),
        "Mest sannsynlig resultat"
      ),
      "Verdi" = c(
        sprintf("%.3f", p$We), sprintf("%.2f", p$lambda1), sprintf("%.2f", p$lambda2),
        sprintf("%.2f", p$total),
        sprintf("%.1f %%", 100 * p$p_win1), sprintf("%.1f %%", 100 * p$p_draw),
        sprintf("%.1f %%", 100 * p$p_win2),
        sprintf("%s %d – %d %s", navn1(), idx[1] - 1, idx[2] - 1, navn2())
      ),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  })

  output$simresultat <- renderText({
    s <- sims()
    if (nrow(s) == 0) return("Klikk «Simuler kamp» for å trekke et resultat.")
    last <- s[nrow(s), ]
    sprintf("%s %d – %d %s", navn1(), last$g1, last$g2, navn2())
  })

  output$simtally <- renderText({
    s <- sims(); n <- nrow(s)
    if (n == 0) return("")
    w1 <- sum(s$g1 > s$g2); dr <- sum(s$g1 == s$g2); w2 <- sum(s$g1 < s$g2)
    sprintf("Etter %d simuleringer: %s %.1f %%, uavgjort %.1f %%, %s %.1f %%.",
            n, navn1(), 100 * w1 / n, 100 * dr / n, navn2(), 100 * w2 / n)
  })

  output$marginplot <- renderPlot({
    p <- pred()
    md <- p$margins; md <- md[md$margin >= -5 & md$margin <= 5, ]
    md$utfall <- ifelse(md$margin > 0, navn1(),
                  ifelse(md$margin < 0, navn2(), "Uavgjort"))
    md$utfall <- factor(md$utfall, levels = c(navn2(), "Uavgjort", navn1()))
    vir3 <- viridisLite::viridis(3)   # tre utfall: lag2, uavgjort, lag1
    g <- ggplot(md, aes(factor(margin), prob, fill = utfall)) +
      geom_col() +
      scale_y_continuous(labels = function(x) paste0(round(100 * x), " %")) +
      scale_fill_manual(values = stats::setNames(
        vir3, c(navn2(), "Uavgjort", navn1()))) +
      labs(x = paste0("Målforskjell (", navn1(), " − ", navn2(), ")"),
           y = "Sannsynlighet", fill = "Utfall") +
      theme_minimal(base_size = 13)
    s <- sims()
    if (nrow(s) > 0) {
      emp <- as.data.frame(table(factor(pmax(pmin(s$g1 - s$g2, 5), -5),
                                        levels = -5:5)))
      names(emp) <- c("margin", "n"); emp$prop <- emp$n / sum(emp$n)
      emp <- emp[as.integer(as.character(emp$margin)) %in% -5:5, ]
      g <- g + geom_point(data = emp,
                          aes(factor(margin), prop), inherit.aes = FALSE,
                          colour = "black", fill = "white", shape = 21,
                          size = 2.5, stroke = 1)
    }
    g
  })

  output$matriseplot <- renderPlot({
    p <- pred(); M <- p$matrix; disp <- 0:6
    df <- expand.grid(mal1 = disp, mal2 = disp)
    df$prob <- mapply(function(a, b) M[a + 1, b + 1], df$mal1, df$mal2)
    ggplot(df, aes(factor(mal2), factor(mal1), fill = prob)) +
      geom_tile(colour = "white") +
      geom_text(aes(label = sprintf("%.1f", 100 * prob),
                    colour = prob > max(df$prob) / 2), size = 3, show.legend = FALSE) +
      scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "white")) +
      scale_fill_viridis_c(guide = "none") +
      labs(x = paste0("Mål ", navn2()), y = paste0("Mål ", navn1())) +
      theme_minimal(base_size = 13)
  })
}

shinyApp(ui, server)
