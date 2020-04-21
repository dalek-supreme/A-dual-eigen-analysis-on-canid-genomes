# see if the 'main positions' intersect
u.cutoff <- 0.8 -> v.cutoff
num_layers <- 10
main.pos.x <- list()
main.pos.y <- list()
for (layer in seq(num_layers)) {
    u <- filter.loadings(snp.svd$u[,layer],method='2',norm.cutoff=u.cutoff)
    v <- filter.loadings(snp.svd$v[,layer],method='2',norm.cutoff=v.cutoff)
    main.pos.x[[layer]] <- nonzero.pos(u)
    main.pos.y[[layer]] <- nonzero.pos(v)
}

# main.points <- data.frame(
#     x <- numeric(0),
#     y <- numeric(0),
#     layer <- numeric(0)
# )

# for (layer in seq(num_layers)){
#     for (i in main.pos.x[[layer]]){
#         for (j in main.pos.y[[layer]]){
#             main.points <- rbind(main.points,c(i,j,layer))
#         }
#     }
# }
# colnames(main.points)<-c('row','col','layer')

main.points.x <- data.frame(
    x <- numeric(0),
    layer <- numeric(0)
)

for (layer in seq(num_layers)){
    for (i in main.pos.x[[layer]]){
        main.points.x <- rbind(main.points.x,c(i,layer))
    }
}
colnames(main.points.x)<-c('row','layer')

library(ggplot2)
main.plot.x <- ggplot(main.points.x, aes(x = row, y = layer, colour = as.factor(layer))) + geom_point(alpha = .5)
pdf('main-row-plot.pdf')
main.plot.x
dev.off()

main.points.y <- data.frame(
    y <- numeric(0),
    layer <- numeric(0)
)

for (layer in 2:num_layers){
    for (i in main.pos.y[[layer]]){
        main.points.y <- rbind(main.points.y,c(i,layer))
    }
    print(layer)
}
colnames(main.points.y)<-c('col','layer')

library(ggplot2)
main.plot.y <- ggplot(main.points.y, aes(x = col, y = layer, colour = as.factor(layer))) + geom_point(alpha = .5)
pdf('main-col-plot.pdf')
main.plot.y
dev.off()
