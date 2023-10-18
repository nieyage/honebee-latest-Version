RNAfold -j6 --noconv --auto-id NRT_OR_TES_upstream50bp.fa > NRT_OR_TES_upstream50bp.out
RNAfold -j6 --noconv --auto-id RT_OR_TES_upstream50bp.fa > RT_OR_TES_upstream50bp.out


sort -t ">" -k 2,2 -o sorted_sequences.fasta NRT_OR_TES_upstream50bp.fa

grep ">" RT_OR_TES_upstream50bp.fa>RT_ID
grep -E -o "\(.+?\)" RT_OR_TES_upstream50bp.out >RTresult
paste RT_ID RTresult |sort -k1,1 


R

data<- read.csv("/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/Fig7/dotplot_feature_TES_type_hairpin2.csv")
# The number of hairpin in OR TES region 
my_comparisons <- list( c("NRT", "RT") )
library(ggpubr)
pdf("./00_Figure/Fig7/Fig7M_hairpin_number_NRT_vs_RT.pdf",width=4,height=4)
ggboxplot(data, x="type", y="Hairpin_number", 
  bxp.errorbar = T,width=0.6, notch = F)+
stat_compare_means(comparisons = my_comparisons)+theme(legend.position="none")+ylab("# of hairpin")
dev.off()

RNAfold --noconv -p < LOC100578045.fa > LOC100578045.res

/md01/nieyg/software/ViennaRNA-2.6.4/src/Utils/relplot.pl LOC100578045_ss.ps LOC100578045_dp.ps > LOC100578045_rss.ps
magick convert -density 300 LOC100578045_rss.ps LOC100578045_rss.pdf # 生成二级结构pdf图






