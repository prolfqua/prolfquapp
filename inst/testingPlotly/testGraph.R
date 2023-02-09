library(plotly)
library(igraph)
library(dplyr)

set.seed(42)
A <- data.frame(Category=rep("A",9))
B <- data.frame(Authors=c("John Snow, 2016 (0, 0)",
                          "Daenerys Targaryen, 2016 (0, 0)",
                          "Arya Stark, 2016 (0, 0)",
                          "Cersei Lannister, 2016 (0, 0)","Tyrion Lannister, 2016 (1, 1)","Brienne of Tarth, 2016 (0, 0)","Theon Greyjoy, 2016 (1, 0)","Khal Drogo, 2015 (16, 0)","Bran Stark, 2015 (3, 3)"))

edgelist <- bind_cols(A,B)
C <- data.frame(Category=rep("B",5))
edgelistB <- bind_cols(C, B[sample(1:nrow(B),5),, drop=F])

edgelist <- rbind(edgelist,edgelistB)

graph <- graph_from_data_frame(edgelist)
plot(graph)
L <- layout_nicely(graph)
vs <- V(graph)
es <- as.data.frame(get.edgelist(graph))
Ne <- length(es[1]$V1)
Xn <- L[,1]
Yn <- L[,2]

network <- plot_ly(type = "scatter",
                   x = Xn,
                   y = Yn,
                   mode = "markers+text",
                   text = names(vs),
)

edge_shapes <- list()
for(i in 1:Ne) {
  v0 <- es[i,]$V1
  v1 <- es[i,]$V2

  edge_shape = list(
    type = "line",
    line = list(color = "#030303", width = 0.3),
    x0 = Xn[match(v0,names(vs))],
    y0 = Yn[match(v0,names(vs))],
    x1 = Xn[match(v1,names(vs))],
    y1 = Yn[match(v1,names(vs))],
    opacity = 1
  )

  edge_shapes[[i]] <- edge_shape}

network <- layout(
  network,
  shapes = edge_shapes,
  xaxis = list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
  yaxis = list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
)

network
