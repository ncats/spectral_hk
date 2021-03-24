library(ggplot2)
library(dplyr)

X <- read.csv('spectral_h1_chembl_pIC50_2.csv')

groups = list(
  c('target_class_l4', 'TK protein kinase group',
    'TK protein kinase', 'spectral_h1_tk_pIC50.png'),
  c('target_class_l3', 'Protein Kinase', 'Kinase',
    'spectral_h1_kinase_pIC50.png'),
  c('target_class_l2', 'Family A G protein-coupled receptor',
    'GCPR', 'spectral_h1_gcpr_pIC50.png'),
  c('target_class_l2', 'Voltage-gated ion channel',
    'Ion channel', 'spectral_h1_ic_pIC50.png') 
)

for (g in groups) {
    df <-X[X[g[1]]==g[2],]
    df$act <- as.numeric(df$avg_act)
    df$sd <- as.numeric(df$std_act)

    hash <- new.env(hash = TRUE)
    labels <- c()
    for (k in df$spectral_h1) {
       x <- 1
       if (exists(k, env = hash)) {
          x <- hash[[k]] + 1
       }
       hash[[k]] <- x
       labels <- c(labels, paste(k, x, sep='_'))
    }
    df$labels <- labels
    
    # only keep unique hash key + target pairs
    #df <- df %>%
    #   group_by(spectral_h1) %>%
    #   filter(n() == 1)
    
    # order by target 
    #df$spectral_h1 <- factor(df$spectral_h1,
    #                    levels=df$spectral_h1[order(df$target)])
    df$labels <- factor(df$labels, levels=df$labels[order(df$target)])
    
    # show the first 50
    data <- df[1:50,]
    ggplot(data) +
      geom_bar(aes(x=labels,y=act), stat='identity', alpha=0.3) +
      geom_text(data=data, aes(x=labels, y=0.1, label=target),
                           hjust=0, size=2, color='blue') +
      geom_errorbar(aes(x=labels,ymin=act-sd,ymax=act+sd)) +
      coord_flip() +
      ggtitle(paste(g[3], '(ChEMBL 28)', sep=' ')) +
      xlab("Spectral key (h1)") +
      ylab("Activity (pIC50)") +
      theme(axis.text.y=element_text(family="mono")) + 
      ggsave(g[4])
}
