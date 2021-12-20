generatePDFreport = function() 
  {
  pdf(file=pdf.report.file, width = 6.5, height = 6.5)

  #cover
  grid.newpage()
  cover <- textGrob("LRGASP Challenge 1\ncomparison report",
                    gp=gpar(fontface="italic", fontsize=40, col="orangered"))
  grid.draw(cover)

  
  # Plot Table 1 
  #grid.arrange(t1.1)
  print(pt1.1)
  # Plot Table 1.2 (UJC collapsed)
  #grid.arrange(t1.2)
  print(pt1.2)
  print("Plotting gene characterization...") 
  # Gene Characterization
  s <- textGrob("Gene Characterization", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
  grid.arrange(s)
  print(p1.1)
  print(p1.2)
  print(p2)
  
  print("Plotting distances to TSS/TTS..")
  # Distances to TSS and TTS of FSM and ISM
  s <- textGrob("Distance to annotated TSS and TTS", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
  grid.arrange(s)
  print(p4.1)
  print(p4.2)
  print(p5.1)
  print(p5.2)
  print(p6.1)
  print(p6.2)
  print(p8)
  
  print("Plotting Jaccard...")
  # Jaccard Index
  s <- textGrob("Presence/absence analysis \nof UJC across pipelines", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
  grid.arrange(s)
  corrplot(jac_index_mat, order="hclust", is.corr = F, method = "circle", tl.col = "black",
           title="Pair-wise Jaccard index between submissions", mar = c(2,1.5,1.5,1.5),
           col = COL2("RdBu", 100 ))
  print(p12)
  for (i in 1:length(p12.1)) {
    print(p12.1[[i]])
  }
  print(p12.2)
  
  
  # Intersections
  if (opt$upset){
    print("Plotting Upset plots...")
    s <- textGrob("Intersection regarding to UJC", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
    grid.arrange(s)
    print(p10)
    print(p10.1)
    print(p11)
   # for (p in 1:length(names(p12))){
   #   print(p12[[p]])
   # }
  }
  
  # Strandard Deviation of TSS and TTS
  print("Plotting SD...")
  s <- textGrob("Standard Deviation of TSS and TTS", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
  grid.arrange(s)
  print(p13)
  print(p14)

  print("Plotting metrics...")
  s <- textGrob("LRGASP metrics Comparison", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
  grid.arrange(s)
  # SIRVs
  #grid.arrange(t3)
  for (i in 1:length(pt3)) {
    print(pt3[[i]])
  }
  # FSM
  #grid.arrange(t4.1)
  for (i in 1:length(pt4.1)) {
    print(pt4.1[[i]])
  }
  #grid.arrange(t4.2)
  for (i in 2:(length(pt4.2)-1)) {
    print(pt4.2[[i]])
  }
  # ISM
  #grid.arrange(t5.1)
  for (i in 1:length(pt5.1)) {
    print(pt5.1[[i]])
  }
  #grid.arrange(t5.2)
  for (i in 2:(length(pt5.2)-1)) {
    print(pt5.2[[i]])
  }
  # NIC
  #grid.arrange(t6.1)
  for (i in 1:length(pt6.1)) {
    print(pt6.1[[i]])
  }
  #grid.arrange(t6.2)
  for (i in 2:length(pt6.2)) {
    print(pt6.2[[i]])
  }
  # NNC
  #grid.arrange(t7.1)
  for (i in 1:length(pt7.1)) {
    print(pt7.1[[i]])
  }
  #grid.arrange(t7.2)
  for (i in 2:length(pt7.2)) {
    print(pt7.2[[i]])
  }
  
  dev.off()
  
  print("SQANTI3 report successfully generated!")
}
