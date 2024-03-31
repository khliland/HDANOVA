#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/


# Modified for use of the Caldana data
# Author: Gooitzen Zwanenburg, April 2018
#
library(shiny)
library(shinydashboard)
library(plyr)
library(stringr)
library(ggplot2)

ui <- dashboardPage(
  dashboardHeader(title = "Data analysis tool",
                  dropdownMenu(type = "notifications",
                               notificationItem(text = "Warning",
                                                icon = icon("exclamation-triangle"),
                                                status = "warning")
                  )
  ),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Uni-variate Analysis", tabName = "univariate", icon = icon("steam-square")),
      menuItem("Model", tabName = "model", icon = icon("image")),
      menuItem("Plots", tabName = "plots", icon = icon("image"))
    )
  ),
  dashboardBody(
    tabItems(

      tabItem(tabName = "univariate",
              fluidRow(
                box(width = 4, status = "primary",
                    selectInput("compound",
                                label = "Compound",
                                choices = list('Alanine' = 1,
                                               'Valine' = 2,
                                               'Leucine' = 3,
                                               'Isoleucine' = 4,
                                               'Proline' = 5,
                                               'Serine' = 6,
                                               'Threonine' = 7,
                                               'beta_alanine' = 8,
                                               'Hydroxyproline' = 9,
                                               'GABA' = 10,
                                               'Aspartate' = 11,
                                               'Asparagine' = 12,
                                               'Methionine' = 13,
                                               'O_acetyl_serine' = 14,
                                               'Glutamate' = 15,
                                               'Phenylalanine' = 16,
                                               'Ornithine' = 17,
                                               'Glutamine' = 18,
                                               'Lysine' = 19,
                                               'Tyrosine' = 20,
                                               'Threonic_acid' = 21,
                                               'Citrulline_Arginine' =22,
                                               'Pyruvic_acid' = 23,
                                               'Citric_acid'  = 24,
                                               'Succinic_acid' = 25,
                                               'Fumaric_acid' = 26,
                                               'Malic_acid' = 27,
                                               'Lactic_acid' = 28,
                                               'Glycolic_acid' = 29,
                                               'Benzoic_acid' = 30,
                                               'Maleic_acid' = 31,
                                               'Nicotinic_acid' = 32,
                                               'Itaconic_acid' = 33,
                                               'Citramalate' = 34,
                                               '4_hydroxy_benzoic_acid' = 35,
                                               'Dehydroascorbic_acid_dimer' = 36,
                                               'Gluconic_acid' = 37,
                                               'Dehydroascorbic_acid' = 38,
                                               'Ascorbic_acid' = 39,
                                               '4_Hydroxycinnamic_acid' = 40,
                                               'Similar_to_Adenine' = 41,
                                               'Shikimate' = 42,
                                               'Erythritol' = 43,
                                               'Arabinose' = 44,
                                               'Arabitol' = 45,
                                               'Fucose' = 46,
                                               'Fructose' = 47,
                                               'Mannitol' = 48,
                                               'Galactose' = 49,
                                               'Glucose' = 50,
                                               'Sucrose' = 51,
                                               'Maltose'= 52,
                                               'Trehalose' = 53,
                                               'Galactinol' = 54,
                                               'myo_inositol' = 55,
                                               'Uracil' = 56,
                                               'Putrescine' = 57,
                                               'Ethanolamine' = 58,
                                               'Glycerol' = 59,
                                               'Indole_3_acetonitrile' = 60,
                                               'Sinapic_acid' = 61,
                                               'Palmitic_acid' = 62,
                                               'Octadecanoic_acid' = 63,
                                               'Docosanoic_acid' = 64,
                                               'Tetracosanoic_acid' = 65,
                                               'Hexacosanoic_acid' =66,
                                               'Octacosanoic_acid' = 67),
                                selected = 1,
                                multiple = FALSE
                    )
                ),
                box(width = 4, status = "primary",
                    radioButtons("scale", "Scaling methods",
                                 c("No scaling" = "noscaling",
                                   "Standardize" = "standardize",
                                   "Square root" = "squareroot",
                                   "Log" = "log")
                    )
                ),
                box(width = 4, status = "primary",
                    checkboxInput("inclinteraction", "Include interaction", FALSE)
                )
              ),
              fluidRow(
                box(width = 12, status = "primary",
                    plotOutput("plot1", height="350px"),
                    plotOutput("plot2", height="350px"),
                    h4("ANOVA table"),
                    verbatimTextOutput("text1")
                )
              )
      ),
      tabItem(tabName = "model",
              fluidRow(
                box(width = 3, status = "primary",
                    selectInput("includefactors",
                                "Include factors",
                                choices = list("1,2",
                                               "1",
                                               "2",
                                               "None"),
                                selected = "None")
                ),
                box(width = 3, status = "primary",
                    selectInput("includeinteraction",
                                "Include interaction",
                                choices = list("None",
                                               "1:2"),
                                selected = "None")
                ),
                box(width = 3, status = "primary",
                    selectInput("combine",
                                "Combine terms",
                                choices = list("None",
                                               "1 + 1:2",
                                               "2 + 1:2"),
                                selected = "None")
                )
              ),
              fluidRow(
                box(width = 9, status = "primary", title = "Model",
                    textOutput("model")
                )
              ),
              fluidRow(
                box(width = 9, status = "primary", title = "Variances",
                    tableOutput("variances")
                )
              )
      ),
      tabItem(tabName = "plots",
              fluidRow(
                uiOutput("toplotUI"),

                box(width = 3, status = "primary",
                    div(style="display: inline-block;vertical-align:top; width: 100px;",
                        selectInput("pc1",
                                    "First PC",
                                    choices = "",
                                    selected = "")
                    ),
                    div(style="display: inline-block;vertical-align:top; width: 20px;",""),
                    div(style="display: inline-block;vertical-align:top; width: 100px;",
                        selectInput("pc2",
                                    "Second PC",
                                    choices = "",
                                    selected = "")
                    )
                )

              ),
              fluidRow(
                box(width = 12, status = "primary",
                    plotOutput("plot")
                )
              )
      )
    )
  )
)

###############################################################################
###############################################################################
#
# Server
#
###############################################################################
###############################################################################

server <- function(input, output, session) {

  toplot                  <- reactiveVal()

  ###############################################################################
  #
  # Functions
  #
  ##############################################################################


  Scaling <- function(data.df) {
    # Scales the data matrix
    # Args:
    #  data.df: (Balanced) dataframe
    # Returns:
    #  data.scaled: scaled dataframe

    # No scaling
    if(input$scale == "noscaling") {
      data.scaled <- data.df
      # Standardize
    } else if (input$scale == "standardize") {
      data.scaled <- as.data.frame(scale(data.df, center = TRUE, scale = TRUE))
      # Square root scaling
    } else if (input$scale == "squareroot") {
      data.scaled <- as.data.frame(sapply(data.df, function(x) sign(x)*sqrt(abs(x))))
      # 10 log scaling
    } else if (input$scale == "log") {
      data.scaled <- as.data.frame(sapply(data.df, function(x) sign(x)*log10(abs(x))))
    }

    return(data.scaled)
  }

  MakeFactorList <- function() {
    # Collects factors to include from input
    # Factors are sorted in ascending order
    # Args:
    #
    # Returns:
    #  factors.to.include: list with factors to include
    # Note: individual factors can be separated by any character except ":"
    #       which is reserved for interactions
    #       a range should be given with a "-": 2-5 is 2,3,4,5

    include.factors <- input$includefactors
    if(include.factors == "1") {
      factors.to.include <- 1
    } else if (include.factors == "2") {
      factors.to.include <- 2
    } else if (include.factors == "1,2") {
      factors.to.include <- c(1,2)
    } else if (include.factors == "None") {
      factors.to.include <- c()
    }

    return(factors.to.include)
  }

  ExtractFactors <- function(st) {
    # Extracts sequences from a character string
    # Args:
    #  st: string with character sequences written as m-n
    # Returns:
    #  factor number or factor sequence

    # Extract number(s)
    tmp.list    <- str_extract_all(st, "[0-9]+")
    from        <- as.numeric(tmp.list[[1]][1])
    to          <- as.numeric(tmp.list[[1]][2])
    if(length(tmp.list[[1]]) == 1) {
      to <- from
    }
    return(as.numeric(seq(from = from, to = to)))
  }

  MakeInteractionList <- function(n.factors) {
    # Collects interactions to include from input
    # Args:
    #  scaled dataframe
    #  number of factors in design matrix
    # Returns:
    #  i.list: list with unique interactions to include
    # Note: interactions are indicated as "n:m"
    if(input$includeinteraction == "1:2") {
      interaction.list <- list(c(1,2))
    } else {
      interaction.list <- list()
    }
    return(interaction.list)
  }

  ExtractInteractions <- function(st, n.factors) {
    # Extracts interacting factors from a character string
    # Factors are sorted in ascending order
    # Args:
    #  st: string with interactions
    # Returns:
    #  interacting.factors: interacting factors

    interacting.factors <- as.numeric(unlist(str_extract_all(st, "[0-9]+")))
    # check for duplicate (eg. 2:2) and out of range factors
    if(any(duplicated(interacting.factors)) | any(interacting.factors > n.factors) ) {
      return(NULL)
    } else {
      return(sort(interacting.factors))
    }
  }

  MakeCombinationList <- function(n.factors) {
    # Get combination of interaction and factors from input field
    # Only two terms may be combined
    # Only a single combination is allowed
    # terms are separated by a "+"
    # Args:
    #   n.factors: number of factors in design matrix
    # Returns:
    #   terms.list: list of factors and interactions to combine

    if(input$combine == "None") {
      terms.list <- list()
    } else if (input$combine == "1 + 1:2") {
      terms.list <- list(1, c(1, 2))
    } else if (input$combine == "2 + 1:2") {
      terms.list <- list(2, c(1, 2))
    }
    return(terms.list)
  }

  MakeCombinationMatrix <- function(design.df, data.df, terms) {
    MakeMatrix <- function(design.df, data.df, term) {

      if(length(term) == 1) {                                          # factor
        matrix.term    <- MakeFactorMatrix(design.df, data.df, as.numeric(term))
      } else if (length(term) > 1) {                                   # interaction
        matrix.term    <- MakeInteractionMatrix(design.df, data.df, as.numeric(term))
      }
      else {
        matrix.term <- NULL
      }
      return(matrix.term)
    }
    matrix <- Reduce("+", lapply(terms,
                                 function(x) MakeMatrix(design.df, data.df, x))
    )
    # For consistency with factor.matrix and interaction.matrix (which are lists)
    # we return a list with a single combination matrix
    combination.matrix <- vector("list", 1)
    combination.matrix[[1]] <- matrix
    return(combination.matrix)
  }

  MakeCombinationText <- function(term, design.df) {
    # Makes text for combinations
    # Args:
    #   term: factor or interaction term
    #   Returns:
    #    text: text for combiation term

    if(length(term) == 1) {   # factor
      text <- colnames(design.df)[term]
    } else {
      text <- colnames(design.df)[term]
    }
    return(text)
  }

  MakeModel <- function(factors.list, interactions.list, combi.list) {

    # check for overlapping interactions
    if(length(combi.list) > 0 & length(interactions.list) > 0) { # combinations and interactions
      # get interaction terms
      combi.interactions           <- Filter(function(x) length(x) > 1, combi.list)
      if(length(combi.interactions) > 0 ) {
        same <- intersect(combi.interactions, interactions.list)
        if(length(same) > 0) {                                   # overlap
          index.same               <- sapply(same, function(y)   # index overlap in interaction.list
            which(sapply(interactions.list,
                         function(x) identical(x, y)))
          )
          for(i in sort(index.same, decreasing = TRUE)) {        # remove overlap from interactions
            interactions.list[[i]] <- NULL
          }
        }
      }
    }
    if(length(combi.list) > 0 & length(factors.list) > 0) { # combinations and factors
      # get factors from combi.list
      combi.factors                <- Filter(function(x) length(x) == 1, combi.list)
      # check for overlapping factors
      if(length(combi.factors) > 0) {
        same                       <- intersect(combi.factors, factors.list)
        if(length(same) > 0) {
          index.same               <- sapply(same, function(x) which(factors.list == x))
          factors.list             <- factors.list[-index.same]
        }
      }
    }
    model                          <- list()
    model$factors                  <- factors.list
    model$interactions             <- interactions.list
    model$combinations             <- combi.list

    return(model)
  }

  MakeModelText <- function(ft, it, ct) {
    ModelText <- "Overall mean"

    if( ft[1] != "" ) {
      ModelText <- paste(ModelText, paste(ft, collapse = " + "), sep = " + ")
    }

    if(length(it) > 0 ) {
      ModelText <- paste(ModelText, paste(it, collapse = " + "), sep = " + ")
    }

    if(ct != "" ) {
      ModelText <- paste(ModelText, paste(ct, collapse = " + "), sep = " + ")
    }

    ModelText <- paste(ModelText, "Residuals", sep = " + ")

    return(ModelText)

  }

  MakeNames <- function(term) {
    # Makes names for factors and interactions
    if(length(term) == 1) {                    # factor
      return(as.character(term))
    } else {                                   # interaction
      return(paste(term, collapse = ":"))
    }
  }

  MakeFactorMatrix <- function(design.df, data.df, factor) {
    # Makes matrix for factor "factor"
    factor.matrix <- MakeAverages(design.df, data.df, factor)[[1]]
    #  overall.means.matrix <- matrix(rep(1, nrow(data.df)))%*%colMeans(data.df)
    #  factor.matrix <- factor.matrix - overall.means.matrix
    return(factor.matrix)
  }


  MakeAverages <- function(design.df, data.df, factors) {
    # Makes column averages over selected rows
    # Args:
    #  design.df: design dataframe
    #  data.df: data dataframe
    #  factors: vector with columns numbers in the design matrix
    # Returns:
    #  cell.averages data frame with cell averages

    GetRows  <- function(design.list, levels) {
      # gets the rows in the design matrix that belong to a level combination
      # Args:
      #  design.list: list with selected columns from the design matrix
      #  levels: level combination of factors
      # Returns:
      #  row.numbers: row numbers for the levels combination

      design.matrix <- sapply(design.list, cbind)
      row.numbers   <- which(sapply(seq(nrow(design.matrix)), function(x)
        identical(design.matrix[x, ], levels))
      )
      return(row.numbers)
    }
    # create empty dataframe
    cell.averages <- as.data.frame(matrix(integer(0),
                                          nrow = dim(data.df)[1],
                                          ncol = dim(data.df)[2])
    )
    design        <- lapply(factors, function(x) as.numeric(design.df[, x]))

    # number of levels for each factor
    num.levels    <- unlist(lapply(design, max))

    # all level combinations
    i.grid        <- as.matrix(expand.grid(lapply(num.levels, seq)))

    # Get matching rows from the design matrix
    matching.rows <- lapply(seq(nrow(i.grid)), function(x) GetRows(design,
                                                                   as.numeric(i.grid[x,])))
    filled.cells  <- matching.rows[sapply(matching.rows, function(x) length(x) > 0)]

    # collect the cell average in selected rows
    for(x in filled.cells) {
      # Note: use rbind to get a matrix out
      cell.averages[x, ] <- as.data.frame(rbind(colMeans(data.df[x, ])))
    }

    # Make vector with levels
    # Note: the length should be the number of rows of desing.

    levels <- rep(0, length(unlist(filled.cells)))
    if(length(levels) != nrow(design.df)) {
      print("Incorrect number of levels")
    }
    for(i in seq(length(filled.cells))) {
      for(j in filled.cells[i]) {
        levels[j] <- i
      }
    }

    return(list(cell.averages, as.factor(levels)))
  }

  MakeInteractionMatrix <- function(design.df, data.df, factors) {
    # Makes interaction matrix from factors
    # Args:
    #  design.df: design matrix
    #  data.df: data matrix
    #  factors: factors that are part of the interaction
    # Returns:
    #  interaction.matrix: matrix of the interaction

    # Make all possible combinations with 'factors':
    # A, B, c, AB, BC, AC, ABC
    combinations       <- lapply(seq(length(factors)), function(x) combn(factors, x))

    # Make matrix with level averages for all combinations
    matrices           <- lapply(combinations, function(x) (lapply(seq(ncol(x)),
                                                                   function(y) MakeAverages(design.df, data.df, x[,y])[[1]]))
    )
    # Make interaction matrix by adding matrices with level averages with the right sign
    # e.g. ABC - AB - AC - BC + A + B + C
    interaction.matrix <- Reduce("+", lapply(seq(from = length(combinations), to = 1, by = -1),
                                             function(x)((-1)^(length(combinations) - x))*Reduce("+",
                                                                                                 matrices[[x]]))
    )
    overall.means.matrix <- matrix(rep(1, nrow(data.df)))%*%colMeans(data.df)

    # Add or subtract overall means from the interaction matrix
    interaction.matrix <- interaction.matrix + ((-1)^(length(factors)))*overall.means.matrix
    return(interaction.matrix)
  }

  MakeResiduals <- function(factors, interactions, combinations, data.df) {
    # Make matrix with residuals
    # Args:
    #  factors: factor matrices
    #  interactions: interaction matrices
    #  combinations: combination matrices
    #  data.df: scaled data matrix
    # Returns:
    #  residuals: matrix with residuals

    overall.means.matrix <- matrix(rep(1, nrow(data.df)))%*%colMeans(data.df)
    residuals   <- data.df
    residuals   <- residuals - overall.means.matrix

    # Factors selected
    if(length(factors) != 0) {
      residuals <- Reduce("-", init = residuals, factors)
    }
    # Interactions selected
    if(length(interactions) != 0 ) {
      residuals <- Reduce("-", init = residuals, interactions)
    }
    # combinations selected
    if(length(combinations) != 0) {
      residuals <- Reduce("-", init = residuals, combinations)
    }
    return(residuals)
  }

  SelectPlotType <- function(plot.type, plot.factor, plot.levels, pc1, pc2,
                             plot.projections = "yes", group.by = 1, x.axis = 0,
                             group.names = "", label.names = "") {
    plot.factor <- as.character(plot.factor)
    if(plot.type == "Scores") {
      p <- PlotScores(pc1, pc2, plot.factor, plot.levels, plot.projections, group.names, label.names)
    } else if (plot.type == "Loadings") {
      p <- PlotLoadings(pc1, pc2, plot.factor)
    } else if (plot.type == "Levels") {
      p <- PlotLevels(pc1, plot.factor, plot.levels, group.by, x.axis, group.names)
    }
    return(p)
  }

  MakePlotParameters <- function(pc1, pc2, pf) {
    # Makes scores, loadings and projections for pc1 and pc2
    # Args:
    #   pc1: first principal component to plot
    #   pc2: second principal component to plot (can be 0)
    #   pf:  factor/interaction/combination to plot
    # Returns: list wwith scores, loadings and projections to plot

    # scores from svd-function are normalized
    u                    <- MakeSession()$svd[[pf]]$u
    d                    <- diag(MakeSession()$svd[[pf]]$d)
    scores.matrix        <- as.data.frame(u %*% d)
    loadings.matrix      <- MakeSession()$svd[[pf]]$v

    if(pc2 == 0) {
      scores               <- as.data.frame(matrix(c(scores.matrix[, pc1],
                                                     rep(0, nrow(scores.matrix))), ncol = 2))
      loadings             <- as.data.frame(matrix(c(loadings.matrix[, pc1],
                                                     rep(0, nrow(loadings.matrix))), ncol = 2))
      projections          <- as.data.frame(as.matrix(MakeSession()$matrices$residuals[[1]])
                                            %*% as.matrix(loadings)) + scores
      percentage.explained <- 100*MakeSession()$svd[[pf]]$d[pc1]^2 /
        sum( (MakeSession()$svd[[pf]]$d)^2 )
    } else {
      scores               <- scores.matrix[, c(pc1, pc2)]
      loadings             <- as.data.frame(MakeSession()$svd[[pf]]$v[, c(pc1, pc2)])
      projections          <- as.data.frame(as.matrix(MakeSession()$matrices$residuals[[1]])
                                            %*% as.matrix(loadings)) + scores
      percentage.explained <- 100*MakeSession()$svd[[pf]]$d[c(pc1, pc2)]^2 /
        sum( (MakeSession()$svd[[pf]]$d)^2 )
    }
    return(list(scores, loadings, projections, percentage.explained))
  }

  PlotScores <- function(pc1, pc2, plot.factor, plot.levels, plot.projections,
                         group.names, label.names) {

    # Makes the score plots
    plot.matrices <-  MakePlotParameters(pc1, pc2, plot.factor)
    scores      <- plot.matrices[[1]]
    projections <- plot.matrices[[3]]
    perc.explnd <- plot.matrices[[4]]

    # Make colors
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    score.colors <- gg_color_hue(length(label.names))
    g <- ggplot(data = scores, aes(x = scores[, 1], y = scores[, 2],
                                   col = plot.levels))
    g <- g + geom_point(size = 6, alpha = 0.7)
    if(plot.projections == "yes") {
      g <- g + geom_point(data = projections, aes(x = projections[, 1],
                                                  y = projections[, 2],
                                                  col = plot.levels))
    }
    if(pc2 == 0) {
      g <- g + labs(title = "Scores",
                    x = paste("PC", as.character(pc1),
                              "(", sprintf("%.2f", perc.explnd[1]),"%", ")",
                              sep = " "),
                    y = ""
      )
    } else {
      g <- g + labs(title = "Scores",
                    x = paste("PC", as.character(pc1),
                              "(", sprintf("%.2f", perc.explnd[1]),"%", ")",
                              sep = " "),
                    y = paste("PC", as.character(pc2),
                              "(", sprintf("%.2f", perc.explnd[2]),"%", ")",
                              sep = " " )
      )
    }
    g <- g + scale_color_discrete(name = group.names, labels = label.names)
    if(length(unique(plot.levels)) > 10)
      g <- g + theme(text = element_text(size = 14),
                     legend.position = "none",
                     legend.text = element_text(size = 14),
                     legend.title = element_text(size = 16),
                     plot.title = element_text(size = 16,
                                               face="bold", hjust = 0.5)
      )
    else {
      g <- g + theme(text = element_text(size = 14),
                     legend.text = element_text(size = 14),
                     legend.title = element_text(size = 16),
                     plot.title = element_text(size = 16,
                                               face="bold", hjust = 0.5)
      )

    }
    return(g)
  }

  PlotLoadings <- function(pc1, pc2, plot.factor) {
    # Makes the loadings plot

    loadings                  <- MakePlotParameters(pc1, pc2, plot.factor)[[2]]
    plot.loadings             <- cbind.data.frame(c(1 : dim(loadings)[1]), loadings)
    if(pc2 == 0) {
      colnames(plot.loadings) <- c("Variable", "L1")
      plot.title              <- paste("Loadings PC ", as.character(pc1))
      p <- ggplot(data = plot.loadings,
                  aes(x = Variable, y = L1)) +
        geom_bar(stat = "identity", fill = "blue") +
        labs(title = plot.title, x = "Variable", y = "Loadings")
    } else {
      colnames(plot.loadings) <- c("Variable", "L1", "L2")
      plot.title              <- paste("Loadings PC ", as.character(pc1), "(blue) and PC ",
                                       as.character(pc2), "(red)")
      p <- ggplot(data = plot.loadings, aes(x = Variable, y = L1)) +
        geom_bar(stat = "identity", fill = "blue") +
        geom_bar(aes(x = Variable, y = L2), stat = "identity",
                 fill = "red", width = 0.5)   +
        labs(title = plot.title, x = "Variable", y = "Loadings")


    }
    return(p)
  }

  PlotLevels <- function(pc, matrix, x.values, group.by, x.axis, group.names) {
    # Plots PC-value against plot.levels and group by group.levels
    scores        <- MakePlotParameters(pc, 0, matrix)[[1]]
    g <- ggplot(data = scores, aes(x = x.values, y = scores[, 1],
                                   col = group.by, group = group.by))
    g <- g + geom_point(size = 4)
    g <- g + geom_line(linewidth = 1)
    g <- g + labs(title = paste("PC ", as.character(pc),
                                " vs levels of factor ",
                                as.character(x.axis)),
                  x = paste("Levels of factor ", as.character(x.axis)),
                  y = paste("PC", as.character(pc), sep = " "))
    g <- g + scale_color_discrete(name = group.names)
    g <- g + theme(text = element_text(size = 14),
                   legend.position = "none",
                   legend.text = element_text(size = 14),
                   legend.title = element_text(size = 16),
                   plot.title = element_text(size = 16,
                                             face="bold", hjust = 0.5)
    )

    return(g)
  }

  SetPCs <- function(n.levels) {
    # Set the values for PC1 an PC2
    # The maximum will be 5, i.e. PC6 and higher cannot be plotted
    if(n.levels == 2) {                           # two levels, single PC
      choices.pc1      <- "1"
      selected.pc1     <- "1"
      choices.pc2      <- "None"
      selected.pc2     <- "None"
    } else if (n.levels > 2) {                    # more than 2 levels
      choices.pc1      <- sapply(seq(n.levels - 1), as.character)
      choices.pc1      <- choices.pc1[as.numeric(choices.pc1) < 6]
      selected.pc1     <- "1"
      choices.pc2      <- sapply(seq(n.levels - 1), as.character)
      choices.pc2      <- choices.pc2[as.numeric(choices.pc2) < 6]
      selected.pc2     <- "2"
    }
    return(list(choices.pc1, selected.pc1, choices.pc2, selected.pc2))
  }

  ########################################################################
  ########################################################################
  #
  #  Reactive expressions
  #
  ########################################################################
  ########################################################################

  ReadInput <- reactive({
    # Reads data- and design file, removes rows with NaN's
    # Returns:
    #  List with datafile, designfile, summary and warnings
    #

    data.df                 <- read.csv("Caldana_data.csv", header=FALSE)
    colnames(data.df)       <- paste("Var", as.character(seq(ncol(data.df))),
                                     sep = "_")
    design.df               <- read.csv("Caldana_F.csv", header = FALSE)
    F                       <- design.df
    colnames(F) <- c("Light", "Time")

    F[F$Light == 1, 1] <- 'Light'
    F[F$Light == 2, 1] <- 'Dark'
    F[F$Light == 3, 1] <- 'Low Light'
    F[F$Light == 4, 1] <- 'High Light'

    F[F$Time == 1, 2] <- '0'
    F[F$Time == 5, 2] <- '40'
    F[F$Time == 2, 2] <- '5'
    F[F$Time == 3, 2] <- '10'
    F[F$Time == 4, 2] <- '20'
    F[F$Time == 6, 2] <- '80'
    F[F$Time == 7, 2] <- '160'

    F$Light <- factor(F$Light, levels = c("Dark", "Low Light", "Light", "High Light"))
    #  F$Light <- factor(F$Light, levels = c("Light", "Dark", "Low Light", "High Light"))
    F$Time  <- factor(F$Time, levels = c("0", "5", "10", "20", "40", "80", "160"))

    colnames(design.df)     <- paste("factor", as.character(seq(ncol(design.df))),
                                     sep = "_")


    number.of.factors       <- dim(design.df)[2]

    design                  <- lapply(seq(ncol(design.df)),
                                      function(x) as.numeric(design.df[, x]))

    # number of levels for each factor
    num.levels              <- unlist(lapply(design, max))
    level.numbers           <- cbind(colnames(design.df), num.levels)
    colnames(level.numbers) <- c("Factor", "Levels")
    return(list(design.df, data.df, level.numbers, F))
  })


  PreTreatData <- reactive({
    # Balances and scales data.
    # Returns:
    #  List with designfile, balanced datafile and scaled datafile

    design.df       <- ReadInput()[[1]]
    data.df         <- ReadInput()[[2]]

    # Scale balanced data
    data.scaled     <- Scaling(data.df)

    return(list(design.df, data.scaled))
  })

  MakeSession <- reactive({

    # Creates session list with model and model matrices
    # Returns:
    #  List: session
    # Get pretreated data
    design.scaled          <- PreTreatData()[[1]]        # design matrix
    data.scaled            <- PreTreatData()[[2]]        # balanced and scaled data matrix

    # Center data
    overall.means.matrix   <- matrix(rep(1, nrow(data.scaled)))%*%colMeans(data.scaled)
    n.factors              <- ncol(design.scaled)        # number of experimental factors
    data.scaled            <- data.scaled - overall.means.matrix
    # Recalculate overal.means.matrix (should be zero now)
    overall.means.matrix <- matrix(rep(1, nrow(data.scaled)))%*%colMeans(data.scaled)

    # Predefine list elements
    factor.matrices                         <- list()
    interaction.matrices                    <- list()
    combination.matrices                    <- list()
    list.svd.factors                        <- list()
    list.svd.interactions                   <- list()
    list.svd.combinations                   <- list()

    # Make session list; contains all matrices
    session                                 <- list()
    session$F                               <- ReadInput()[[4]]
    session$design                          <- design.scaled
    session$scaled                          <- data.scaled
    session$ssq$scaled                      <- sum( (data.scaled)^2 )
    session$ssq$means                       <- sum( (overall.means.matrix)^2 )

    factor.list                             <- MakeFactorList()
    interactions.list                       <- MakeInteractionList(n.factors)
    combinations.list                       <- MakeCombinationList(n.factors)

    model                                   <- MakeModel(factor.list, interactions.list,
                                                         combinations.list)

    # Note: factor.list and interactions.list are no longer valid because MakeModel() may have updated
    # lists. factor.list and interaction.list
    factor.list                             <- model$factors
    interactions.list                       <- model$interactions
    factor.names                            <- lapply(factor.list, function(x) MakeNames(x))
    interaction.names                       <- lapply(interactions.list, function(x) MakeNames(x))
    combination.names                       <- lapply(combinations.list, function(x) MakeNames(x))

    # Factors
    if(length(factor.list) == 0) {                       # no factors
      session$factors                       <- list()    # dummy factor list
      session$matrices$factors              <- list()    # dummy list for matrices for each factor
      session$ssq$factors                   <- list()    # dummy list for ssq of matrix for each factor
      list.svd.factors                      <- list()    # summy list for svd on matrix for each factor
    } else {


      factor.matrices                       <- lapply(factor.list, function(x) MakeFactorMatrix(
        design.scaled, data.scaled, x )
      )
      names(factor.matrices)                <- factor.names
      factor.ssq                            <- lapply(factor.matrices, function(x) sum(x^2))
      list.svd.factors                      <- lapply(factor.matrices, svd)

      session$factors                       <- factor.list           # included factors
      session$matrices$factors              <- vector("list", length(factor.list))
      session$ssq$factors                   <- vector("list", length(factor.list))
      names(session$matrices$factors)       <- factor.names
      names(session$ssq$factors)            <- factor.names
      session$matrices$factors[]            <- factor.matrices
      session$ssq$factors[]                 <- factor.ssq

    }

    # Interactions
    if(length(interactions.list) == 0) {
      session$interactions                  <- list()
      session$interactions.list             <- list()
      session$matrices$interactions         <- list()
      session$ssq$interactions              <- list()
      list.svd.interactions                 <- list()
    } else {
      interaction.matrices                  <- lapply(interactions.list,
                                                      function(x) MakeInteractionMatrix(design.scaled,
                                                                                        data.scaled, x)
      )
      names(interaction.matrices)           <- interaction.names
      interaction.ssq                       <- lapply(interaction.matrices, function(x) sum(x^2))
      list.svd.interactions                 <- lapply(interaction.matrices, svd)

      session$interactions.list             <- interactions.list     # list of the included interactions
      session$interactions                  <- interaction.names     # names of included interactions
      session$matrices$interactions         <- vector("list", length(interactions.list))
      session$ssq$interactions              <- vector("list", length(interactions.list))
      names(session$interactions.list)      <- interaction.names
      names(session$matrices$interactions)  <- interaction.names
      names(session$ssq$interactions)       <- interaction.names
      session$matrices$interactions[]       <- interaction.matrices
      session$ssq$interactions[]            <- interaction.ssq
    }

    # Combinations
    if(length(combinations.list) == 0) {
      session$combination                   <- list()
      session$matrices$combination          <- list()
      session$ssq$combinations              <- list()
      list.svd.combinations                 <- list()
    } else {
      # In this version of the program, this is a single matrix
      combination.matrices                  <- MakeCombinationMatrix(design.scaled,
                                                                     data.scaled, combinations.list)
      combination.ssq                       <- lapply(combination.matrices, function(x) sum(x^2))
      list.svd.combinations                 <- lapply(combination.matrices, svd)
      session$combinations.list             <- combinations.list
      session$combinations                  <- combination.names
      session$matrices$combinations         <- vector("list", length(combination.matrices))
      session$ssq$combinations              <- vector("list", length(combination.matrices))
      if(length(combination.matrices) == 1) {
        names(session$matrices$combinations)  <- "Combination"
        names(session$ssq$combinations)       <- "Combination"
        session$matrices$combinations[]       <- combination.matrices
        session$ssq$combinations[]            <- combination.ssq
      }
    }

    # Residuals
    session$matrices$residuals              <- vector("list", 1)
    session$ssq$residuals                   <- vector("list", 1)
    names(session$matrices$residuals)       <- "Residuals"
    names(session$ssq$residuals)            <- "Residuals"
    session$matrices$residuals[[1]]         <- MakeResiduals(factor.matrices, interaction.matrices,
                                                             combination.matrices, data.scaled)
    session$ssq$residuals[[1]]              <- sum((session$matrices$residual[[1]])^2)
    list.svd.residuals                      <- lapply(session$matrices$residuals, svd)

    # Model components
    if(length(factor.list) >0 ) {
      factor.model.text                     <- paste("factor", as.character(factor.list), sep = "_")
    } else {
      factor.model.text                     <- ""
    }
    interaction.model.text                  <- lapply(interactions.list,
                                                      function(x) paste("interaction",
                                                                        paste(x, collapse = "_"), sep="_")
    )
    if(length(combinations.list) > 0 ) {
      combinations.f                        <- MakeCombinationText(combinations.list[[1]], design.scaled)
      combinations.i                        <- sapply(combinations.list[2],
                                                      function(x) paste("interaction",
                                                                        paste(x, collapse = "_"), sep="_"))
      combination.model.text                <- paste("(", paste(combinations.f, combinations.i, sep = " + "),
                                                     ")")

    } else {
      combination.model.text                <- ""
    }
    model.text                              <- MakeModelText(factor.model.text, interaction.model.text,
                                                             combination.model.text)

    output$model                            <- renderText(model.text)

    # PCA
    list.svd                                <- c(list.svd.factors, list.svd.interactions, list.svd.combinations,
                                                 list.svd.residuals)
    session$svd                             <- vector("list", (length(list.svd.factors) + length(list.svd.interactions) +
                                                                 length(list.svd.combinations) + length(list.svd.residuals))
    )
    names(session$svd)                      <- c(names(session$matrices$factors),
                                                 names(session$matrices$interactions),
                                                 names(session$matrices$combinations), "Residuals")
    session$svd[]                           <- list.svd
    return(session)
  })

  UpdatePlotMenu <- observe({
    # Makes the menu items for the objects (factor, interaction, combination, residual)
    factors                 <- as.character(MakeSession()$factors)
    interactions            <- MakeSession()$interactions
    combinations            <- MakeSession()$combinations



    # Residuals only
    if(length(combinations) == 0 & length(interactions) == 0 & length(factors) == 0) {

      output$toplotUI <-  renderUI({
        tagList(
          box(status = "primary", width = 3,
              radioButtons("toplotres", "Plot",
                           c("Residuals"),
                           selected = "Residuals")),
          box(status = "primary", width = 3,
              radioButtons("plotviewres", "Type",
                           c("Scores", "Loadings"),
                           selected = "Scores")
          ),
          box(status = "primary", width = 3,
              selectInput("reslevel", "Show factor",
                          choices = c("Light", "Time"),
                          selected = "Light")
          )
        )
      })
      updateRadioButtons(session, "toplotres", "Plot",
                         choices = c("Residuals"),
                         selected = "Residuals")
      toplot("Residuals")
      updateRadioButtons(session, "plotviewres", "Type",
                         choices = c("Scores", "Loadings"),
                         selected = input$plotview)
      plot.factor       <- seq(ncol(MakeSession()$design))
    }
    # Factors Only
    else if(length(combinations) == 0 & length(interactions) == 0 & length(factors) != 0) {

      output$toplotUI <-  renderUI({
        tagList(
          box(status = "primary", width = 3,
              radioButtons("toplotfac", "Plot",
                           c("Factors", "Residuals"),
                           selected = "Factors")
          ),
          conditionalPanel(condition = "input.toplotfac == 'Factors'",
                           box(status = "primary", width = 3,
                               selectInput("choose.factor1", "Factor",
                                           choices = factors, selected = factors[1])
                           ),
                           box(status = "primary", width = 3,
                               radioButtons("plotviewfac", "Type",
                                            c("Scores", "Loadings", "Levels"),
                                            selected = "Scores")
                           )
          ),
          conditionalPanel(condition = "input.toplotfac == 'Residuals'",
                           box(status = "primary", width = 3,
                               radioButtons("plotviewres", "Type",
                                            c("Scores", "Loadings"),
                                            selected = "Scores")
                           ),
                           box(status = "primary", width = 3,
                               selectInput("reslevel", "Show factor",
                                           choices = c("Light", "Time"),
                                           selected = "Light")
                           )
          )
        )
      })
      updateRadioButtons(session, "toplotfac", "Plot",
                         choices = c("Factors", "Residuals"),
                         selected = input$toplotfac)
      req(input$toplotfac)
      toplot(input$toplotfac)
      if(input$toplotfac == "Factors") {
        updateSelectInput(session, "choose.factor1", "Factor",
                          choices  = factors,
                          selected = input$choose.factor1)
        req(input$choose.factor1)
        updateRadioButtons(session, "plotviewfac", "Type",
                           choices = c("Scores", "Loadings", "Levels"),
                           selected = input$plotviewfac)
        plot.factor       <- as.numeric(input$choose.factor1)
        updateSelectInput(session, "xvaluesfac", choices = plot.factor,
                          selected = plot.factor[1])
        updateSelectInput(session, "groupfactorfac", choices = "None",
                          selected = "None")
      } else if(input$toplotfac == "Residuals") {
        updateRadioButtons(session, "plotviewres", "Type",
                           choices = c("Scores", "Loadings"),
                           selected = input$plotview)
        plot.factor     <- seq(ncol(MakeSession()$design))
      }
    }

    # Interactions only
    else if(length(combinations) == 0 & length(interactions) != 0 & length(factors) == 0) {
      output$toplotUI <-  renderUI({
        tagList(
          box(status = "primary", width = 3,
              radioButtons("toplotint", "Plot",
                           c("Interaction", "Residuals"),
                           selected = "Interaction")
          ),
          conditionalPanel(condition = "input.toplotint == 'Residuals'",
                           box(status = "primary", width = 3,
                               radioButtons("plotviewres", "Type",
                                            c("Scores", "Loadings"),
                                            selected = "Scores")),
                           box(status = "primary", width = 3,
                               selectInput("reslevel", "Show factor",
                                           choices = c("Light", "Time"),
                                           selected = "Light")
                           )
          ),
          conditionalPanel(condition = "input.toplotint == 'Interaction'",
                           box(status = "primary", width = 3,
                               radioButtons("plotviewint", "Type",
                                            c("Scores", "Loadings", "Levels"),
                                            selected = "Scores")
                           )
          ),
          conditionalPanel(condition = "input.plotviewint == 'Levels'",
                           box(width = 3,
                               selectInput("xvaluesint",
                                           "x-axis: levels of factor",
                                           choices = "",
                                           selected = ""),
                               selectInput("groupfactorint",
                                           "group by",
                                           choices = "",
                                           selected = "")
                           )
          )
        )
      })
      updateRadioButtons(session, "toplotint", "Plot",
                         c("Interaction", "Residuals"),
                         selected = input$toplotint)
      req(input$toplotint)
      toplot(input$toplotint)
      if(input$toplotint == "Interaction") {
        updateRadioButtons(session, "plotviewint", "Type",
                           choices = c("Scores", "Loadings", "Levels"),
                           selected = input$plotviewint)
        plot.factor       <- ExtractInteractions(input$includeinteraction,
                                                 ncol(MakeSession()$design) )
        updateSelectInput(session, "xvaluesint", choices = plot.factor,
                          selected = plot.factor[1])
        updateSelectInput(session, "groupfactorint", choices = plot.factor,
                          selected = plot.factor[2])


      } else if(input$toplotint == "Residuals") {
        updateRadioButtons(session, "plotviewres", "Type",
                           choices = c("Scores", "Loadings"),
                           selected = input$plotviewres)
        plot.factor     <- seq(ncol(MakeSession()$design))
      }
    }

    # Combination only
    else if(length(combinations) != 0 & length(interactions) == 0 & length(factors) == 0) {
      output$toplotUI <-  renderUI({
        tagList(
          box(status = "primary", width = 3,
              radioButtons("toplotcom", "Plot",
                           c("Combinations", "Residuals"),
                           selected = "Combinations")
          ),
          conditionalPanel(condition = "input.toplotcom == 'Combinations'",
                           box(status = "primary", width = 3,
                               radioButtons("plotviewcomb", "Type",
                                            c("Scores", "Loadings", "Levels"),
                                            selected = "Scores")
                           ),
                           conditionalPanel(condition = "input.plotviewcomb == 'Levels'",
                                            box(status = "primary", width = 3,
                                                selectInput("xvaluescomb",
                                                            "x-axis: levels of factor",
                                                            choices = "",
                                                            selected = ""),
                                                selectInput("groupfactorcomb",
                                                            "group by",
                                                            choices = "",
                                                            selected = "")
                                            )
                           )
          ),
          conditionalPanel(condition = "input.toplotcom == 'Residuals'",
                           box(status = "primary", width = 3,
                               radioButtons("plotviewres", "Type",
                                            c("Scores", "Loadings"),
                                            selected = "Scores")),
                           box(status = "primary", width = 3,
                               selectInput("reslevel", "Show factor",
                                           choices = c("Light", "Time"),
                                           selected = "Light")
                           )
          )
        )
      })
      updateRadioButtons(session, "toplotcom", "Plot",
                         c("Combinations", "Residuals"),
                         selected = input$toplotcom)
      req(input$toplotcom)
      toplot(input$toplotcom)
      if(input$toplotcom == "Combinations") {
        updateRadioButtons(session, "plotviewcomb", "Type",
                           choices = c("Scores", "Loadings", "Levels"),
                           selected = input$plotviewcomb)
        plot.factor      <- unique(unlist(MakeCombinationList(ncol(MakeSession()$design))))
        updateSelectInput(session, "xvaluescomb", choices = plot.factor,
                          selected = plot.factor[1])
        updateSelectInput(session, "groupfactorcomb", choices = plot.factor,
                          selected = plot.factor[2])


      } else if(input$toplotcom == "Residuals") {
        updateRadioButtons(session, "plotviewres", "Type",
                           choices = c("Scores", "Loadings"),
                           selected = input$plotviewres)
        plot.factor     <- seq(ncol(MakeSession()$design))
      }
    }

    # Factor plus interaction
    else if(length(combinations) == 0 & length(interactions) != 0 & length(factors) != 0) {
      output$toplotUI <-  renderUI({
        tagList(
          box(status = "primary", width = 3,
              radioButtons("toplotfacint", "Plot",
                           c("Factors", "Interaction", "Residuals"),
                           selected = "Interaction")
          ),
          conditionalPanel(condition = "input.toplotfacint == 'Residuals'",
                           box(status = "primary", width = 3,
                               radioButtons("plotviewres", "Type",
                                            c("Scores", "Loadings"),
                                            selected = "Scores")
                           ),
                           box(status = "primary", width = 3,
                               selectInput("reslevel", "Show factor",
                                           choices = c("Light", "Time"),
                                           selected = "Light")
                           )
          ),
          conditionalPanel(condition = "input.toplotfacint == 'Interaction'",
                           box(status = "primary", width = 3,
                               radioButtons("plotviewint", "Type",
                                            c("Scores", "Loadings", "Levels"),
                                            selected = "Scores")
                           )
          ),
          conditionalPanel(condition = "input.plotviewint == 'Levels'",
                           box(width = 3,
                               selectInput("xvaluesint",
                                           "x-axis: levels of factor",
                                           choices = "",
                                           selected = ""),
                               selectInput("groupfactorint",
                                           "group by",
                                           choices = "",
                                           selected = "")
                           )
          ),
          conditionalPanel(condition = "input.toplotfacint == 'Factors'",
                           box(status = "primary", width = 3,
                               selectInput("choose.factor1", "Factor",
                                           choices = factors, selected = factors[1])
                           ),
                           box(status = "primary", width = 3,
                               radioButtons("plotviewfac", "Type",
                                            c("Scores", "Loadings", "Levels"),
                                            selected = "Scores")
                           )
          )
        )
      })


      updateRadioButtons(session, "toplotfacint", "Plot",
                         choices = c("Factors", "Interaction", "Residuals"),
                         selected = input$toplotfacint)
      req(input$toplotfacint)
      toplot(input$toplotfacint)

      if(input$toplotfacint == "Factors") {                              # Factor
        updateSelectInput(session, "choose.factor1", "Factor",
                          choices  = factors,
                          selected = input$choose.factor1)
        req(input$choose.factor1)
        updateRadioButtons(session, "plotviewfac", "Type",
                           choices = c("Scores", "Loadings", "Levels"),
                           selected = input$plotviewfac)
        plot.factor       <- as.numeric(input$choose.factor1)
        updateSelectInput(session, "xvaluesint", choices = plot.factor,
                          selected = plot.factor[1])
        updateSelectInput(session, "groupfactorint", choices = plot.factor,
                          selected = plot.factor[2])

      } else if(input$toplotfacint == "Interaction") {                   # Interaction
        updateRadioButtons(session, "plotviewint", "Type",
                           choices = c("Scores", "Loadings", "Levels"),
                           selected = input$plotviewint)
        plot.factor       <- ExtractInteractions(input$includeinteraction,
                                                 ncol(MakeSession()$design) )
        updateSelectInput(session, "xvaluesint", choices = plot.factor,
                          selected = plot.factor[1])
        updateSelectInput(session, "groupfactorint", choices = plot.factor,
                          selected = plot.factor[2])


      } else if(input$toplotfacint == "Residuals") {                      # Residual
        updateRadioButtons(session, "plotviewres", "Type",
                           choices = c("Scores", "Loadings"),
                           selected = input$plotviewres)
        plot.factor     <- seq(ncol(MakeSession()$design))

      }
    }
    # Combination plus factor
    else if(length(combinations) != 0 & length(interactions) == 0 & length(factors) != 0) {
      output$toplotUI <-  renderUI({
        tagList(
          box(status = "primary", width = 3,
              radioButtons("toplotfaccom", "Plot",
                           c("Factors", "Combinations", "Residuals"),
                           selected = "Combinations")
          ),
          conditionalPanel(condition = "input.toplotfaccom == 'Residuals'",
                           box(status = "primary", width = 3,
                               radioButtons("plotviewres", "Type",
                                            c("Scores", "Loadings"),
                                            selected = "Scores")
                           ),
                           box(status = "primary", width = 3,
                               selectInput("reslevel", "Show factor",
                                           choices = c("Light", "Time"),
                                           selected = "Light")
                           )
          ),
          conditionalPanel(condition = "input.toplotfaccom == 'Combinations'",
                           box(status = "primary", width = 3,
                               radioButtons("plotviewcomb", "Type",
                                            c("Scores", "Loadings", "Levels"),
                                            selected = "Scores")
                           )
          ),
          conditionalPanel(condition = "input.plotviewcomb == 'Levels'",
                           box(width = 3,
                               selectInput("xvaluescomb",
                                           "x-axis: levels of factor",
                                           choices = "",
                                           selected = ""),
                               selectInput("groupfactorcomb",
                                           "group by",
                                           choices = "",
                                           selected = "")
                           )
          ),
          conditionalPanel(condition = "input.toplotfaccom == 'Factors'",
                           box(status = "primary", width = 3,
                               selectInput("choose.factor1", "Factor",
                                           choices = factors, selected = factors[1])
                           ),
                           box(status = "primary", width = 3,
                               radioButtons("plotviewfac", "Type",
                                            c("Scores", "Loadings", "Levels"),
                                            selected = "Scores")
                           )
          )
        )
      })
      updateRadioButtons(session, "toplotfaccom", "Plot",
                         c("Factors", "Combinations", "Residuals"),
                         selected = input$toplotfaccom)

      req(input$toplotfaccom)
      toplot(input$toplotfaccom)
      if(input$toplotfaccom == "Combinations") {
        updateRadioButtons(session, "plotviewcomb", "Type",
                           choices = c("Scores", "Loadings", "Levels"),
                           selected = input$plotviewcomb)

        plot.factor      <- unique(unlist(MakeCombinationList(ncol(MakeSession()$design))))

        updateSelectInput(session, "xvaluescomb", choices = plot.factor,
                          selected = plot.factor[1])
        updateSelectInput(session, "groupfactorcomb", choices = plot.factor,
                          selected = plot.factor[2])
      } else if (input$toplotfaccom == "Factors")  {
        updateSelectInput(session, "choose.factor1", "Factor",
                          choices  = factors,
                          selected = factors[1])
        req(input$choose.factor1)
        updateRadioButtons(session, "plotviewfac", "Type",
                           choices = c("Scores", "Loadings", "Levels"),
                           selected = input$plotviewfac)
        plot.factor       <- as.numeric(input$choose.factor1)
      } else if(input$toplotfaccom == "Residuals") {
        updateRadioButtons(session, "plotviewres", "Type",
                           choices = c("Scores", "Loadings"),
                           selected = input$plotviewres)
        plot.factor     <- seq(ncol(MakeSession()$design))
      }
    }

    plot.levels        <- MakeAverages(MakeSession()$design,
                                       MakeSession()$scaled, plot.factor)[[2]]
    n.levels           <- min(max(as.numeric(plot.levels)),
                              ncol(MakeSession()$scaled),
                              nrow(MakeSession()$scaled))
    pcs                <- SetPCs(n.levels)
    updateSelectInput(session, "pc1", choices = pcs[[1]],
                      selected = pcs[[2]])
    updateSelectInput(session, "pc2", choices = pcs[[3]],
                      selected = pcs[[4]])
  })

  DoPlot <- reactive({

    plotFig <- toplot()
    req(plotFig)
    if(plotFig == "Factors") {
      req(input$choose.factor1, input$plotviewfac, input$pc1, input$pc2)
      type.to.plot   <- input$plotviewfac
      factor.to.plot <- as.numeric(input$choose.factor1)
      plot.levels <- MakeSession()$F[, factor.to.plot]
      #     plot.levels        <- MakeAverages(MakeSession()$design,
      #                                         MakeSession()$scaled, factor.to.plot)[[2]]
      pc1            <- as.numeric(input$pc1)
      pc2            <- input$pc2
      if(pc2 == "None") {
        pc2 <- 0
      } else {
        pc2 <- as.numeric(input$pc2)
      }
      group.names    <- colnames(MakeSession()$F)[factor.to.plot]
      label.names    <- levels(MakeSession()$F[, factor.to.plot])

      if(type.to.plot == "Levels") {
        index       <- as.numeric(gsub("factor_", "", as.character(factor.to.plot)))
        group.index <- as.numeric(gsub("factor_", "", as.character(factor.to.plot)))
        plot.levels <- MakeSession()$F[, index]
        group.by    <- as.factor(rep(1, nrow(MakeSession()$F)))
        p <- SelectPlotType(type.to.plot, input$choose.factor1, plot.levels, pc1, pc2,
                            group.by = group.by, x.axis = colnames(MakeSession()$F)[index],
                            group.names = colnames(MakeSession()$F)[group.index])
      } else {

        p <- SelectPlotType(type.to.plot, factor.to.plot, plot.levels, pc1, pc2,
                            group.names = group.names, label.names = label.names)
      }

    } else if (plotFig == "Interaction") {
      req(input$includeinteraction == "1:2", input$plotviewint, input$pc1, input$pc2)
      type.to.plot       <- input$plotviewint
      factors            <- ExtractInteractions(input$includeinteraction,
                                                ncol(MakeSession()$design) )
      interaction.name   <- MakeNames(factors)
      levels             <- MakeAverages(MakeSession()$design,
                                         MakeSession()$scaled, factors)
      plot.levels        <- levels[[2]]
      pc1            <- as.numeric(input$pc1)
      pc2            <- input$pc2
      if(pc2 == "None") {
        pc2 <- 0
      } else {
        pc2 <- as.numeric(input$pc2)
      }
      if(type.to.plot == "Levels") {
        req(input$xvaluesint, input$groupfactorint)
        index       <- as.numeric(gsub("factor_", "", input$xvaluesint))
        group.index <- as.numeric(gsub("factor_", "", input$groupfactorint))
        plot.levels <- MakeSession()$F[, index]
        group.by    <- MakeSession()$F[, group.index]
        p <- SelectPlotType(type.to.plot, interaction.name, plot.levels, pc1, pc2,
                            group.by = group.by, x.axis = colnames(MakeSession()$F)[index],
                            group.names = colnames(MakeSession()$F)[group.index])
      } else {
        print(interaction.name)
        p <- SelectPlotType(type.to.plot, interaction.name, plot.levels, pc1, pc2)
      }
      return(p)

    } else if (plotFig == "Combinations") {
      req(input$combine, input$plotviewcomb, input$pc1, input$pc2)
      type.to.plot   <- input$plotviewcomb
      factors        <- unique(unlist(MakeCombinationList(ncol(MakeSession()$design))))
      levels         <- MakeAverages(MakeSession()$design,
                                     MakeSession()$scaled, factors)
      plot.levels    <- levels[[2]]
      pc1            <- as.numeric(input$pc1)
      pc2            <- input$pc2
      if(pc2 == "None") {
        pc2 <- 0
      } else {
        pc2 <- as.numeric(input$pc2)
      }
      if(type.to.plot == "Levels") {
        req(input$xvaluescomb, input$groupfactorcomb)
        index       <- as.numeric(gsub("factor_", "", input$xvaluescomb))
        plot.levels <- MakeSession()$F[, index]
        group.index <- as.numeric(gsub("factor_", "", input$groupfactorcomb))
        group.by    <- MakeSession()$F[, group.index]
        p <- SelectPlotType(type.to.plot, "Combination", plot.levels, pc1, pc2,
                            group.by = group.by,
                            x.axis = colnames(MakeSession()$F)[index],
                            group.names = colnames(MakeSession()$F)[group.index])
      } else {
        p <- SelectPlotType(type.to.plot, "Combination", plot.levels, pc1, pc2)
      }
      return(p)

    } else if (plotFig == "Residuals") {

      req(input$plotviewres, input$pc1, input$pc2)
      type.to.plot   <- input$plotviewres
      plot.levels    <- as.factor(rep(1, nrow(MakeSession()$design)))
      if(input$reslevel == "Light") {
        factors <- 1
      } else if (input$reslevel == "Time") {
        factors <- 2
      }

      n.levels       <- min(ncol(MakeSession()$scaled), nrow(MakeSession()$scaled))
      levels         <- MakeAverages(MakeSession()$design,
                                     MakeSession()$scaled, factors)
      plot.levels    <- levels[[2]]
      group.names    <- colnames(MakeSession()$F)[factors]
      label.names    <- levels(MakeSession()$F[, factors])
      pcs            <- SetPCs(n.levels)
      pc1            <- as.numeric(input$pc1)
      pc2            <- input$pc2
      if(pc2 == "None") {
        pc2 <- 0
      } else {
        pc2 <- as.numeric(input$pc2)
      }
      p <- SelectPlotType(type.to.plot, "Residuals", plot.levels, pc1, pc2, plot.projections = "no",
                          group.names = group.names, label.names = label.names)
      return(p)

    }
    return(p)
  })

  SsqTable <- reactive ({
    if(length(MakeSession()$factors) > 0) {
      factors           <- unlist(MakeSession()$factors)
      ssq.factors       <- unlist(MakeSession()$ssq$factors)
    } else {
      factors           <- NULL
      ssq.factors       <- NULL
    }
    if(length(MakeSession()$interactions) > 0) {
      interactions      <- unlist(MakeSession()$interactions)
      ssq.interactions  <- unlist(MakeSession()$ssq$interactions)
    } else {
      interactions      <- NULL
      ssq.interactions  <- NULL
    }
    if(length(MakeSession()$combinations) > 0) {
      combinations      <- "Combination"
      ssq.combinations  <- unlist(MakeSession()$ssq$combinations)
    } else {
      combinations      <- NULL
      ssq.combinations  <- NULL
    }
    ssq.names           <- c("Overall Means", "Data", factors, interactions, combinations, "Residuals")
    ssq.values          <- c(MakeSession()$ssq$means, MakeSession()$ssq$scaled, ssq.factors,
                             ssq.interactions, ssq.combinations,
                             MakeSession()$ssq$residuals[[1]])
    ssq.percentages     <- (ssq.values / MakeSession()$ssq$scaled)*100
    ssq.table           <- cbind.data.frame(ssq.names, ssq.values, ssq.percentages)
    colnames(ssq.table) <- c("Source", "Sum of Squares", "Percentage of variation")
    return(ssq.table)

  })

  #######################################Caldana data############################

  doUnivariate <- reactive({

    # Modified specifically for Caldana data

    just_data <- PreTreatData()[[2]]
    colnames(just_data) <- c('Alanine',
                             'Valine',
                             'Leucine',
                             'Isoleucine',
                             'Proline',
                             'Serine',
                             'Threonine',
                             'beta_alanine',
                             'Hydroxyproline',
                             'GABA',
                             'Aspartate',
                             'Asparagine',
                             'Methionine',
                             'O_acetyl_serine',
                             'Glutamate',
                             'Phenylalanine',
                             'Ornithine',
                             'Glutamine',
                             'Lysine',
                             'Tyrosine',
                             'Threonic_acid',
                             'Citrulline_Arginine',
                             'Pyruvic_acid',
                             'Citric_acid',
                             'Succinic_acid',
                             'Fumaric_acid',
                             'Malic_acid',
                             'Lactic_acid',
                             'Glycolic_acid',
                             'Benzoic_acid',
                             'Maleic_acid',
                             'Nicotinic_acid',
                             'Itaconic_acid',
                             'Citramalate',
                             '4_hydroxy_benzoic_acid',
                             'Dehydroascorbic_acid_dimer',
                             'Gluconic_acid',
                             'Dehydroascorbic_acid',
                             'Ascorbic_acid',
                             '4_Hydroxycinnamic_acid',
                             'Similar_to_Adenine',
                             'Shikimate',
                             'Erythritol',
                             'Arabinose',
                             'Arabitol',
                             'Fucose',
                             'Fructose',
                             'Mannitol',
                             'Galactose',
                             'Glucose',
                             'Sucrose',
                             'Maltose',
                             'Trehalose',
                             'Galactinol',
                             'myo_inositol',
                             'Uracil',
                             'Putrescine',
                             'Ethanolamine',
                             'Glycerol',
                             'Indole_3_acetonitrile',
                             'Sinapic_acid',
                             'Palmitic_acid',
                             'Octadecanoic_acid',
                             'Docosanoic_acid',
                             'Tetracosanoic_acid',
                             'Hexacosanoic_acid',
                             'Octacosanoic_acid')
    F <- PreTreatData()[[1]]
    colnames(F) <- c("Light", "Time")

    F[F$Light == 1, 1] <- 'Light'
    F[F$Light == 2, 1] <- 'Dark'
    F[F$Light == 3, 1] <- 'Low Light'
    F[F$Light == 4, 1] <- 'High Light'

    F[F$Time == 1, 2] <- '0'
    F[F$Time == 5, 2] <- '40'
    F[F$Time == 2, 2] <- '5'
    F[F$Time == 3, 2] <- '10'
    F[F$Time == 4, 2] <- '20'
    F[F$Time == 6, 2] <- '80'
    F[F$Time == 7, 2] <- '160'

    caldana_data <- cbind(just_data, F)

    # We need to make the 'Light' column a factor (is now numeric)
    caldana_data$Light  <- as.factor(caldana_data$Light)

    # We need to make the 'Time' column a factor (is now numeric)
    caldana_data$Time  <- as.factor(caldana_data$Time)

    ######################  Anova  ##############################


    compound_nr     <- as.numeric(input$compound)
    metabolite      <- names(just_data[compound_nr])
    metabolite_data <- caldana_data[, metabolite]

    # Make anova data frame (columns 68:69 are "Light" and "Time")
    anova_data <- cbind( 'metabolite' = metabolite_data, caldana_data[, 68:69])

    # Set Time and Light levels in the right order.
    anova_data$Time <- factor(anova_data$Time, levels = c("0", "5", "10", "20", "40", "80", "160"))
    anova_data$Light <- factor(anova_data$Light, levels = c("Dark", "Low Light", "Light", "High Light"))

    # Do anova with the  aov function
    if (input$inclinteraction) {
      aov_met <- aov(metabolite ~ Time*Light, data = anova_data)
    }
    else {
      aov_met <- aov(metabolite ~ Time + Light, data = anova_data)
    }


    # Calculate averages per Time and per Light level
    metabolite_avg <- ddply(anova_data, .(Time, Light), summarize, value = mean(metabolite))

    # Set order of light levels (again)
    metabolite_avg$Light <- factor(metabolite_avg$Light, levels = c("Dark", "Low Light", "Light", "High Light"))

    # Use ggplot to make a plot of the average data.
    p <- ggplot(metabolite_avg, aes(x = Time, y = value, col = Light)) +
      geom_point(size = 3.5) +
      geom_line(data = metabolite_avg, aes(x = Time, y = value, group = Light, col = Light))  +
      geom_point(data = anova_data, aes(x = Time, y = anova_data[, 1] ), size = 1.6, alpha = 0.8) +
      labs(title = paste("Measured", metabolite),  x = 'Time', y = 'Concentration', color = "Light levels") +
      theme(text = element_text(size = 14),
            legend.text = element_text(size = 14),
            legend.title = element_text(size = 16),
            plot.title = element_text(size = 16,
                                      face="bold", hjust = 0.5)
      )


    # ### Estimated data

    fitted_data = cbind(aov_met$fitted.values, caldana_data[,68:69])
    colnames(fitted_data) = c("value", "Light", "Time")

    # Set Time and Light levels in the right order.
    fitted_data$Time <- factor(anova_data$Time, levels = c("0", "5", "10", "20", "40", "80", "160"))
    fitted_data$Light <- factor(anova_data$Light, levels = c("Dark", "Low Light", "Light", "High Light"))

    # Use ggplot to make a plot anova estimates.
    q <- ggplot(data = fitted_data, aes(x = Time, y = value, col = Light)) +
      geom_point(size = 3) +
      geom_line( aes(x = Time, y = value, group = Light, col = Light))  +
      geom_point(data = metabolite_avg, aes(x = Time, y = value, group = Light, col = Light),
                 size=3 , shape=24) +
      labs(title = paste("Estimated", metabolite),  x = 'Time', y = 'Concentration', color= "Light levels") +
      theme(text = element_text(size = 14),
            legend.text = element_text(size = 14),
            legend.title = element_text(size = 16),
            plot.title = element_text(size = 16,
                                      face="bold", hjust = 0.5)
      )

    aov_table         <- summary(aov_met)

    return(list(p, q, aov_table))
  })



  output$plot1         <- renderPlot(doUnivariate()[[1]])

  output$plot2         <- renderPlot(doUnivariate()[[2]])
  output$text1         <- renderPrint(doUnivariate()[[3]])

  output$plot          <- renderPlot(print(DoPlot()))
  output$variances     <- renderTable(SsqTable())

  #  output$upc            <- renderTable(UpdatePcs())


}




shinyApp(ui, server)
