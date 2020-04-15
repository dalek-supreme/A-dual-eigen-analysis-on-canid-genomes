load('snp.svd.R')
pdf('SVplot(from2).pdf')
plot(sort(snp.svd$d[2:length(snp.svd$d)]))
dev.off()

load('snp.rsvd.R')
pdf('RSVplot(from2).pdf')
plot(sort(snp.rsvd$d[2:length(snp.rsvd$d)]))
dev.off()
