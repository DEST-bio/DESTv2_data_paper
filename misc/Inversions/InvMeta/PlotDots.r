legend.col <- function(col, lev) {
    opar <- par
    n <- length(col)
    bx <- par("usr")
    box.cx <- c(
        bx[2] + (bx[2] - bx[1]) / 1000,
        bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50
    )
    box.cy <- c(bx[3], bx[3])
    box.sy <- (bx[4] - bx[3]) / n
    xx <- rep(box.cx, each = 2)
    par(xpd = TRUE)
    for (i in 1:n) {
        yy <- c(
            box.cy[1] + (box.sy * (i - 1)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i - 1))
        )
        polygon(xx, yy, col = col[i], border = col[i])
    }
    par(new = TRUE)
    plot(0, 0,
        type = "n",
        ylim = c(min(lev), max(lev)),
        yaxt = "n", ylab = "",
        xaxt = "n", xlab = "",
        frame.plot = FALSE
    )
    map.axes(side = 4, las = 2, tick = FALSE, line = .25)
    par <- opar
}

library(rworldmap)
library(kriging)
library(tidyverse)
library(maps)
library(viridis)
dataset <- read.table("/media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/InvMeta/FullInvDestv2.txt", header = T, na.string = "NA")
# dataset=read.table("/media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/InvMeta/EurInv.txt",header=T,na.string="NA")

dataset <- na.omit(dataset)
attach(dataset)

data_long <- gather(dataset, Inversion, Frequency, In.2L.t:In.3R.Payne, factor_key = TRUE)

data_long$DEST[data_long$DEST == 21] <- "previous data"
data_long$DEST[data_long$DEST == 22] <- "DEST dataset"

panel.grid.major = element_blank(), panel.grid.minor = element_blank()

data_long$DataSource <- as.factor(data_long$DEST)
world_map <- map_data("world")

world_plot <- ggplot() +
    geom_polygon(
        data = world_map, aes(x = long, y = lat, group = group),
        fill = "lightgray", color = "black"
    ) +
    geom_point(data = data_long, aes(x = Longitude, y = Latitude, color = Frequency, pch = DataSource), size = 2) +
    scale_color_gradientn(
        name = "Inversion Frequency",
        colours = c("white", "yellow", "red", "purple", "blue")
    ) + # You can change the color palette as needed
    facet_wrap(. ~ Inversion, nrow = 2, ncol = 2) +
    xlab("Longitude") +
    ylab("Latitude") +
    theme_bw() +
    theme(panel.background = element_rect(fill = "lightblue"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())
world_plot

ggsave("/media/inter/mkapun/projects/DESTv2_data_paper/misc/Inversions/InvMeta/PlotDots.pdf",
    world_plot,
    width=12,
    height=6)