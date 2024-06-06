plot_bar2 <- function (physeq, x = "Sample", y = "Abundance", fill = NULL, 
                       title = NULL, facet_grid = NULL) 
{
  mdf = psmelt(physeq)
  #mdf$Sample = factor(mdf$Sample, levels = c('10-20','30-40','50-100'))
  mdf$Sample = factor(mdf$Sample, levels = c('10','20','30','40','50','60','70','80','90','100'))
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack")
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}