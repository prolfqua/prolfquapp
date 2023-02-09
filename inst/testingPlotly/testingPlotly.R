
library(plotly)
library(crosstalk)
REPORTDATA <- readRDS("/REPORTDATA.rds")
datax <- REPORTDATA$cse$get_contrasts()


xx <- datax |>
  group_by(Bait) |> nest()


makeshared <- function(x){
  SharedData$new(x, ~Prey, group = "Choose protein")
}
xx <- mutate(xx, shared = purrr::map( data,  makeshared  ))


res <- list()

xd <- purrr::map( xx$shared, function(x){
  p <- plot_ly(x, x = ~log2_EFCs, y = ~ I(-log10(BFDR)), type = "scatter", mode = "markers" , color = I("black"), text = ~ Prey, showlegend = FALSE) |>
    highlight("plotly_click", off = "plotly_doubleclick" )
  return(p)
}
)

DT::datatable(shared_prot)



xd |>  subplot(nrows = 1, shareX = TRUE, shareY = TRUE)
