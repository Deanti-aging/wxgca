WXGCAcode = function()
{#约束传播只对于以约束作为先验的情况有用，对于标签作为先验没用因为直接排列就行
  #此外，我这种情况是不允许用标签调参的，要调参只能用约束符合数来调参——聚类结果排列成必连勿连，然后一个矩阵上必连对应1，勿连对应0，我的已有约束也对应一个这样的矩阵，然后统计两个矩阵元素相同数量
  #一般安装过一次后这部分要注释掉，否则每运行都要install，其实只要libaray
  #约束数也不是越多越好，不仅慢很多，对SSMDK反而会变差
  # b=a[a$dataset_name=="FLAG"& a$clustering_method=="SSMDK"&a$Link_num==10 ,4]
  # install.packages("dplyr") #交并补运算
  # install.packages("ggplot2")
  # install.packages("tidytext")
  # install.packages("pacman") 
  # install.packages("ggthemes")
  # install.packages("readr")
  # install.packages("showtext")
  # install.packages("patchwork")
  # install.packages("cluster") #计算轮廓系数的函数silhouette()，返回聚类的平均轮廓宽度
  # install.packages("kmed") #fastkmed
  # install.packages("conclust") #ccls lcvqe ckmeans
  # install.packages("clv") 
  # install.packages("SSLR") #wine(dataset)
  # install.packages("dbscan") #dbscan hdbscan
  # install.packages("plot3D") #scatter3D函数 画3D图 
  #install.packages("flexclust") #Rand Index
  # install.packages("export")#保存图片为svg并且字体为中文GB1 ，ggsave和svg无法处理中文字体为白体
  # library(export)
  # install.packages("wesanderson") 
  # library(wesanderson)#wes_palette函数 画图
  library(flexclust)
  library(plot3D)
  library(dplyr)
  library(cluster)
  library(ggplot2)
  library(tidytext)
  library(patchwork)
  library(kmed)
  library(conclust)
  library(clv)
  library(SSLR)
  library(pacman)
  library(ggthemes)
  library(readr)
  library(showtext)
  library(dbscan)
  datalist = data_initialize()
  tround = 5
  loop = 1
  dataset_name = Link_num = mustLink_num = 
    cantLink_num = clustering_method = 
    RI = ARI = SIL = sd = 
    minptsd = minptsh = eps = c()
  
  
  for(i in c(1:length(datalist))){#1:length(datalist)
    dataset = datalist[[i]]
    ncluster =length(unique(dataset$datalabel))
    #普通的半监督聚类算法和聚类算法只能跑纯数值的，所以分两个data
    if(length(dataset$idcat)==0 && length(dataset$idbin)==0){
      data1 = data2 = data.frame(dataset$data) 
    }else{
      data1 = data.frame(dataset$data[, dataset$idnum]) #因为TAE只有一个num属性，结果提取出来就是向量，下面ncol(data1)直接报错
      data2 = data.frame(dataset$data)
    }
    #data1需要我给他们0-1标准化，而data2在混合变量函数里会标准化,注意标准化是一列列地来，如果属性量纲一致就不标准化而损失信息
    #好像不用归一化，这是测试算法为目的，不是建模为目的
    # for (j in c(1:ncol(data1))) {
    #   max_c = max(data1[,j])
    #   min_c = min(data1[,j])
    #   for (p in c(1:nrow(data1))) {
    #     data1[p,j] = (data1[p,j]-min_c)/(max_c-min_c)
    #   }
    # }
    distdata1 = as.matrix(dist(data1))
    
    datalabel = dataset$datalabel
    
    # if(dataset$dataname == "NOISE"){
    #   ncluster = length(unique(dataset$datalabel))-1 #因为噪声簇实际只有两簇
    # }else{
    #   ncluster = length(unique(dataset$datalabel))#实际应用中要分的簇数可能是未知的，这里我懒得加这个调参选项,直接设置为实际簇数
    # }
    
    
    
    idnum = dataset$idnum
    idbin = dataset$idbin
    idcat = dataset$idcat
    #用data2和gower distance生成的distdata2来找密度最高点产生高质量约束，因为他们算法无所谓，但是对我影像很大
    # if(length(dataset$idcat)==0 && length(dataset$idbin)==0){
    #   distdata2 = distdata1
    # }else{
    #   distdata2 = distmix(data2, method = "gower", idnum = idnum, idbin = idbin, idcat = idcat)
    # }
    # high_density_points = sort(rowSums(distdata2),decreasing = T)[c(1:(0.8*length(datalabel)))]#一半高质量簇
    # high_density_points = as.numeric(names(high_density_points))
    consnum_list = round(seq(from = 10, to =1000, length.out = 5)) 
    comb = t(combn(c(1:length(datalabel)),2))#high_density_points 
    
    
    for (j in c(1:length(consnum_list))){#
      constrainedLink_num = consnum_list[j]
      
      result = list(lcvqe_ri=c(), lcvqe_ari=c(), lcvqe_sil=c(),
                    ckmeans_ri=c(), ckmeans_ari=c(), ckmeans_sil=c(),
                    ccls_ri=c(), ccls_ari=c(), ccls_sil=c(),
                    # fastkmed_ri=c(), fastkmed_ari=c(), fastkmed_sil=c(),
                    # fastkmedo_ri=c(), fastkmedo_ari=c(), fastkmedo_sil=c(),
                    # dbscan_ri=c(), dbscan_ari=c(), dbscan_sil=c(),
                    # hdbscan_ri=c(), hdbscan_ari=c(), hdbscan_sil=c(),
                    ssmdk_ri=c(), ssmdk_ari=c(),ssmdk_sil=c(),
                    ssmdka_ri=c(), ssmdka_ari=c(),ssmdka_sil=c(),
                    ssmdd_ri=c(), ssmdd_ari=c(),ssmdd_sil=c(),
                    ssmdh_ri=c(), ssmdh_ari=c(),ssmdh_sil=c()
      )
      
      for (p in c(1:tround)) {#每个数据集-每个约束数-每个算法跑10波不同的伪随机抽样的约束对算均值
        #生成must-link和cannot-link并且保证必连和勿连约束至少都存在一个
        mustLink = NULL
        cantLink = NULL
       
        
        while(is.null(mustLink) || is.null(cantLink)){
          set.seed(loop) 
          loop = loop + 1 #为了每波约束对不同
          labelpair_idx = sample(c(1:nrow(comb)), size = constrainedLink_num, replace = FALSE)#
          allink = comb[labelpair_idx,]
          for(k in c(1:nrow(allink))){
            if(datalabel[allink[k,1]]==datalabel[allink[k,2]]){mustLink = rbind(mustLink, allink[k,])}
            if(datalabel[allink[k,1]]!=datalabel[allink[k,2]]){cantLink = rbind(cantLink, allink[k,])}
          }
        }
        #LCVQE
        if (length(idnum)==1){
          result$lcvqe_ri = 0
        }else{
          lcvqe_result = lcvqe(data1, k=ncluster, mustLink, cantLink)#迭代数默认10就10吧
          result$lcvqe_ri = c(result$lcvqe_ri, indicator(datalabel, lcvqe_result, distdata1)$ri)#计算指标给我返回来4给数？
          # result$lcvqe_ari = c(result$lcvqe_ari, indicator(datalabel, lcvqe_result, as.matrix(dist(data1)))$ari)
          # result$lcvqe_sil = c(result$lcvqe_sil, indicator(datalabel, lcvqe_result, as.matrix(dist(data1)))$sil)
        }
        cat("Processed ",i,consnum_list[j],p," LCVQE", " .\n")
        # COP-kmeans 通过我复制其源码和用其例子但是k=1测试发现，返回0时表示出现约束违反
        ckmeans_result = ckmeans(data1, k=ncluster, mustLink, cantLink)
        if (ckmeans_result[1] == 0) {
          result$ckmeans_ri = c(result$ckmeans_ri, 0)
          # result$ckmeans_ari = c(result$ckmeans_ari, -1)
          # result$ckmeans_sil = c(result$ckmeans_sil, -1)
        } else {
          result$ckmeans_ri = c(result$ckmeans_ri, indicator(datalabel, ckmeans_result, distdata1)$ri)
          # result$ckmeans_ari = c(result$ckmeans_ari, indicator(datalabel, ckmeans_result, as.matrix(dist(data1)))$ari)
          # result$ckmeans_sil = c(result$ckmeans_sil, indicator(datalabel, ckmeans_result, as.matrix(dist(data1)))$sil)
        }
        # tryCatch({
        #   ckmeans_result = ckmeans(data1, k=ncluster, mustLink, cantLink)#默认100
        # }, error = function(e){
        #   result$ckmeans_ri = c(result$ckmeans_ri, -1)
        #   result$ckmeans_ari = c(result$ckmeans_ari, -1)
        #   result$ckmeans_sil = c(result$ckmeans_sil, -1)
        # }, finally = {
        #   result$ckmeans_ri = c(result$ckmeans_ri, indicator(datalabel, ckmeans_result, as.matrix(dist(data1)))$ri)
        #   result$ckmeans_ari = c(result$ckmeans_ari, indicator(datalabel, ckmeans_result, as.matrix(dist(data1)))$ari)
        #   result$ckmeans_sil = c(result$ckmeans_sil, indicator(datalabel, ckmeans_result, as.matrix(dist(data1)))$sil)
        # })
        cat("Processed ",i,consnum_list[j],p," COP-kmeans", " .\n")
        
        #CCLS 如果数据集进行过标准化，这个算法就几乎不会随着约束增加而改进了
        ccls_result = ccls(data1, k=ncluster, mustLink, cantLink)#默认1
        result$ccls_ri = c(result$ccls_ri, indicator(datalabel, ccls_result, distdata1)$ri)
        # result$ccls_ari = c(result$ccls_ari, indicator(datalabel, ccls_result, as.matrix(dist(data1)))$ari)
        # result$ccls_sil = c(result$ccls_sil, indicator(datalabel, ccls_result, as.matrix(dist(data1)))$sil)
        cat("Processed ",i,consnum_list[j],p," CCLS", " .\n")
        #WXCGA 
        wxgca_result = wxgca(data2,datalabel,ncluster,idnum,idbin,idcat,mustLink,cantLink)
        result$ssmdk_ri = c(result$ssmdk_ri, wxgca_result$ssmdk)
        # result$ssmdk_ari = c(result$ssmdk_ari, indicator(datalabel, wxgca_result$ssmdk, wxgca_result$distdata)$ari)
        # result$ssmdk_sil = c(result$ssmdk_sil, indicator(datalabel, wxgca_result$ssmdk, wxgca_result$distdata)$sil)
        # result$ssmdka_ri = c(result$ssmdka_ri, wxgca_result$ssmdka)
        # result$ssmdko_ari = c(result$ssmdko_ari, indicator(datalabel, wxgca_result$ssmdk, wxgca_result$distdata)$ari)
        # result$ssmdko_sil = c(result$ssmdko_sil, indicator(datalabel, wxgca_result$ssmdk, wxgca_result$distdata)$sil)
        result$ssmdd_ri = c(result$ssmdd_ri, wxgca_result$ssmdd)
        # result$ssmdd_ari = c(result$ssmdd_ari, indicator(datalabel, wxgca_result$ssmdd, wxgca_result$distdata)$ari)
        # result$ssmdd_sil = c(result$ssmdd_sil, indicator(datalabel, wxgca_result$ssmdd, wxgca_result$distdata)$sil)
        result$ssmdh_ri = c(result$ssmdh_ri, wxgca_result$ssmdh)
        # result$ssmdh_ari = c(result$ssmdh_ari, indicator(datalabel, wxgca_result$ssmdh, wxgca_result$distdata)$ari)
        # result$ssmdh_sil = c(result$ssmdh_sil, indicator(datalabel, wxgca_result$ssmdh, wxgca_result$distdata)$sil)
        
        # minptsd = c(minptsd,wxgca_result$minptsd)
        # minptsh = c(minptsh,wxgca_result$minptsh)
        # eps = c(eps,wxgca_result$eps)
        cat("Processed ", i,consnum_list[j],p," wxgca"," .\n")
      }
      
      dataset_name = c(dataset_name, names(datalist[i]))
      Link_num = c(Link_num, constrainedLink_num)
      # mustLink_num = c(mustLink_num, nrow(mustLink))
      # cantLink_num = c(cantLink_num, nrow(cantLink))
      clustering_method = c(clustering_method, "LCVQE")
      RI = c(RI, mean(result$lcvqe_ri) )
      # ARI = c(ARI, mean(result$lcvqe_ari) )
      # SIL = c(SIL, mean(result$lcvqe_sil) )
      sd = c(sd, sd(result$lcvqe_ri))
      
      dataset_name = c(dataset_name, names(datalist[i]))
      Link_num = c(Link_num, constrainedLink_num)
      # mustLink_num = c(mustLink_num, nrow(mustLink))
      # cantLink_num = c(cantLink_num, nrow(cantLink))
      clustering_method = c(clustering_method, "COP-kmeans")
      RI = c(RI, mean(result$ckmeans_ri) )
      # ARI = c(ARI, mean(result$ckmeans_ari) )
      # SIL = c(SIL, mean(result$ckmeans_sil) )
      sd = c(sd, sd(result$ckmeans_ri))
      
      dataset_name = c(dataset_name, names(datalist[i]))
      Link_num = c(Link_num, constrainedLink_num)
      # mustLink_num = c(mustLink_num, nrow(mustLink))
      # cantLink_num = c(cantLink_num, nrow(cantLink))
      clustering_method = c(clustering_method, "CCLS")
      RI = c(RI, mean(result$ccls_ri) )
      # ARI = c(ARI, mean(result$ccls_ari) )
      # SIL = c(SIL, mean(result$ccls_sil) )
      sd = c(sd, sd(result$ccls_ri))
      
      dataset_name = c(dataset_name, names(datalist[i]))
      Link_num = c(Link_num, constrainedLink_num)
      # mustLink_num = c(mustLink_num, nrow(mustLink))
      # cantLink_num = c(cantLink_num, nrow(cantLink))
      clustering_method = c(clustering_method, "SSMDK")
      RI = c(RI, mean(result$ssmdk_ri) )
      # ARI = c(ARI, mean(result$ssmdk_ari) )
      # SIL = c(SIL, mean(result$ssmdk_sil) )
      sd = c(sd, sd(result$ssmdk_ri))
      
      # dataset_name = c(dataset_name, names(datalist[i]))
      # Link_num = c(Link_num, constrainedLink_num)
      # # mustLink_num = c(mustLink_num, nrow(mustLink))
      # # cantLink_num = c(cantLink_num, nrow(cantLink))
      # clustering_method = c(clustering_method, "SSMDKA")
      # RI = c(RI, mean(result$ssmdka_ri) )
      # # ARI = c(ARI, mean(result$ssmdko_ari) )
      # # SIL = c(SIL, mean(result$ssmdko_sil) )
      # sd = c(sd, sd(result$ssmdka_ri))
      
      dataset_name = c(dataset_name, names(datalist[i]))
      Link_num = c(Link_num, constrainedLink_num)
      # mustLink_num = c(mustLink_num, nrow(mustLink))
      # cantLink_num = c(cantLink_num, nrow(cantLink))
      clustering_method = c(clustering_method, "SSMDD")
      RI = c(RI, mean(result$ssmdd_ri) )
      # ARI = c(ARI, mean(result$ssmdd_ari) )
      # SIL = c(SIL, mean(result$ssmdd_sil) )
      sd = c(sd, sd(result$ssmdd_ri))
      
      dataset_name = c(dataset_name, names(datalist[i]))
      Link_num = c(Link_num, constrainedLink_num)
      # mustLink_num = c(mustLink_num, nrow(mustLink))
      # cantLink_num = c(cantLink_num, nrow(cantLink))
      clustering_method = c(clustering_method, "SSMDH")
      RI = c(RI, mean(result$ssmdh_ri) )
      # ARI = c(ARI, mean(result$ssmdh_ari) )
      # SIL = c(SIL, mean(result$ssmdh_sil) )
      sd = c(sd, sd(result$ssmdh_ri))
      
      #在跳出不同约束但同约束量和同一数据集的这里加上3给聚类算法的实验
      #注意临时变量用p而不是i和j，否则后面用i取数据集名称就出错了
      # densitymaxlist = rowSums(distdata1)
      # densitypoint = c()
      # for (p in c(1:ncluster)) {
      #   densitypoint = c(densitypoint, which(densitymaxlist == max(densitymaxlist))[1])
      #   densitymaxlist[which(densitymaxlist == max(densitymaxlist))[1]] = min(densitymaxlist)
      # }
      n = nrow(distdata1)
      # k_list = seq(from = 2, to = n/2)
      # max_result = -2 #ARI和SIL范围-1到1，RI是0-1，
      # for (p in k_list){
      #   fastkmed_ri = tryCatch(
      #     {
      #       fastkmed_ri = indicator(datalabel, fastkmed(distdata1, p)$cluster, distdata1)$ri
      #     },
      #     error = function(e){0}
      #   )
      #   if(fastkmed_ri > max_result){
      #     max_result = fastkmed_ri
      #   }
      # }
      # fastkmed_ri = max_result
      
      fastkmed_ri = indicator(datalabel, fastkmed(distdata1, ncluster)$cluster, distdata1)$ri
      
      # fastkmeda_ri = indicator(datalabel, fastkmed(distdata1, ncluster, init = densitypoint )$cluster, distdata1)$ri #默认10迭代
      #我说怎么原算法表现这么好，这里用标签调参了
      
      
      minpts_list = seq(from = 2, to = n/ncluster)
      max_result = -2 #ARI和SIL范围-1到1，RI是0-1，
      for (p in minpts_list){
        hdbscan_ri = indicator(datalabel, hdbscan(distdata1, p)$cluster, distdata1)$ri
        if(hdbscan_ri > max_result){
          max_result = hdbscan_ri
          # minpts_h = p
        }
      }
      hdbscan_ri = max_result
      # hdbscan_ri = indicator(datalabel, hdbscan(distdata1, 3)$cluster, distdata1)$ri
      
      max_value = max(distdata1)
      distdata0 = distdata1
      distdata0[which(distdata0==0)] = max_value
      min_value = min(distdata0)
      minpts_list = seq(from = 2, to = n/ncluster)
      eps_list = seq(from = min_value, to = max_value,  length.out = 100)#
      max_result = -2
      for (p in eps_list){#不用调参mintps后整个算法都快了有10倍
        dbscan_ri = indicator(datalabel, dbscan(distdata1, p, minPts = 5, weights = NULL, borderPoints = TRUE)$cluster, distdata1)$ri
        if(dbscan_ri > max_result){
          max_result = dbscan_ri
        }
      }
      dbscan_ri = max_result
      
      # minpts_list = seq(from = 2, to = n/2)
      # max_result = -2
      # for (p in minpts_list){#不用调参mintps后整个算法都快了有10倍
      #   k_dist = c()
      #   for (k in c(1:n)) {
      #     k_dist[k] = sum(sort(distdata[k,])[1:p]) #每个点k近邻之和并从排序
      #   }
      #   k_dist = sort(k_dist)
      #   k_dist_thr = max(k_dist[2:n]-k_dist[1:(n-1)])
      #   dbscan_ri = indicator(datalabel, dbscan(distdata, k_dist_thr, p)$cluster, distdata)$ri
      #   if(dbscan_ri > max_result){
      #     max_result = dbscan_ri
      #   }
      # }
      # dbscan_ri = max_result
      
      # minpts = 2*(length(idbin)+length(idcat)+length(idnum))-1
      # k_dist = c()
      # for (p in c(1:n)) {
      #   sum_n = distdata1[p,]
      #   sum_n = sort(sum_n)
      #   sum_n = sum_n[1:minpts]
      #   k_dist[p] = sum(sum_n) #每个点k近邻之和并从排序
      # }
      # k_dist = sort(k_dist)
      # k_dist_thr = max(k_dist[2:n]-k_dist[1:(n-1)])
      # dbscan_ri = indicator(datalabel, dbscan(distdata1, k_dist_thr, minPts = minpts)$cluster, distdata1)$ri

      
      dataset_name = c(dataset_name, rep(names(datalist[i]),3) )
      Link_num = c(Link_num, rep(constrainedLink_num,3) )
      clustering_method = c(clustering_method,"FASTKMED","HDBSCAN","DBSCAN")
      RI = c(RI,fastkmed_ri,hdbscan_ri,dbscan_ri)
      sd = c(sd,0,0,0)#一个单独值算不出sd，会返回一个NA，导致作图自动删除相关行
    }
    #picture
    out = data.frame(
      dataset_name,
      clustering_method,
      Link_num,
      # mustLink_num,
      # cantLink_num,
      RI,
      # ARI,
      # SIL,
      sd
      # minptsh,
      # minptsd,
      # eps
    )
    #which里的逻辑号用|和&表示或和且
    out_now = out[which(out$dataset_name==names(datalist[i]) & out$clustering_method!="FASTKMED" & out$clustering_method!="HDBSCAN" & out$clustering_method!="DBSCAN"),]#因为out累计了所有数据
    #没有sd显示的话，全实线无点的折线图就行，不过要自定义配色并且把坐标轴题目等正过来，
    print(ggplot(out_now, aes(x=Link_num, y=RI,linetype = clustering_method,  color= clustering_method, group = clustering_method))#, shape = clustering_method
          + geom_line(linewidth = 1)
          # + geom_point(size = 3)
          + labs(x="Constraints Number",#约束数
                 y = "RI",
                 title = names(datalist[i]))
          + theme_bw() 
          + theme(legend.position = "bottom",
                  legend.direction = "horizontal",
                  legend.title = element_blank(),
                  panel.grid = element_blank(),
                  plot.title = element_text(hjust = 0.5))
          + theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
          + scale_x_continuous(breaks = consnum_list)
          + scale_color_manual(values = c("#16048A","#6201A9","#9E189F","#CC4A76","#EA7854","#FDB330","#f8e620")) #自定义颜色
          + ylim(0,1))
    scriptpath = rstudioapi::getSourceEditorContext()$path
    scriptpathn = nchar(scriptpath)
    suppath = substr(scriptpath, 1, scriptpathn-12)
    datasetpath = paste(suppath,"/dataset/",sep="")
    ggsave(paste(datasetpath,names(datalist[i]),"_results.pdf", sep=""),device = cairo_pdf,width =5, height =3.75)
    # ggsave(paste(datasetpath,names(datalist[i]),"_results硕士毕业.pdf", sep=""),device = cairo_pdf,width =5, height =3.75,family="GB1")

    
    if (TRUE) {
      out_now = out[which(out$dataset_name==names(datalist[i]) & (out$clustering_method=="FASTKMED" | out$clustering_method=="SSMDK" )),]
      #没有sd显示的话，全实线无点的折线图就行，不过要自定义配色并且把坐标轴题目等正过来，
      print(ggplot(out_now, aes(x=Link_num, y=RI,linetype = clustering_method,  color= clustering_method, group = clustering_method))#, shape = clustering_method
            + geom_line(linewidth = 1)
            # + geom_point(size = 3)
            + labs(x="Constraints Number",#Constraints Number
                   y = "RI",
                   title = names(datalist[i]))
            + theme_bw() 
            + theme(legend.position = "bottom",
                    legend.direction = "horizontal",
                    legend.title = element_blank(),
                    panel.grid = element_blank(),
                    plot.title = element_text(hjust = 0.5))
            + theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
            + scale_x_continuous(breaks = consnum_list)
            + scale_color_manual(values = c("#16048A","#6201A9","#9E189F","#CC4A76","#EA7854","#FDB330","#f8e620")) #自定义颜色
            + ylim(0,1))
      scriptpath = rstudioapi::getSourceEditorContext()$path
      scriptpathn = nchar(scriptpath)
      suppath = substr(scriptpath, 1, scriptpathn-12)
      datasetpath = paste(suppath,"/dataset/",sep="")
      ggsave(paste(datasetpath,names(datalist[i]),"_FASTKMED_results.pdf", sep=""),device = cairo_pdf,width =5, height =3.75)
      
      
      out_now = out[which(out$dataset_name==names(datalist[i]) & (out$clustering_method=="DBSCAN" | out$clustering_method=="SSMDD" )),]
      #没有sd显示的话，全实线无点的折线图就行，不过要自定义配色并且把坐标轴题目等正过来，
      print(ggplot(out_now, aes(x=Link_num, y=RI,linetype = clustering_method,  color= clustering_method, group = clustering_method))#, shape = clustering_method
            + geom_line(linewidth = 1)
            # + geom_point(size = 3)
            + labs(x="Constraints Number",#Constraints Number
                   y = "RI",
                   title = names(datalist[i]))
            + theme_bw() 
            + theme(legend.position = "bottom",
                    legend.direction = "horizontal",
                    legend.title = element_blank(),
                    panel.grid = element_blank(),
                    plot.title = element_text(hjust = 0.5))
            + theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
            + scale_x_continuous(breaks = consnum_list)
            + scale_color_manual(values = c("#16048A","#6201A9","#9E189F","#CC4A76","#EA7854","#FDB330","#f8e620"))#自定义颜色
            + ylim(0,1))
      scriptpath = rstudioapi::getSourceEditorContext()$path
      scriptpathn = nchar(scriptpath)
      suppath = substr(scriptpath, 1, scriptpathn-12)
      datasetpath = paste(suppath,"/dataset/",sep="")
      ggsave(paste(datasetpath,names(datalist[i]),"_DBSCAN_results.pdf", sep=""),device = cairo_pdf,width =5, height =3.75)
      
      out_now = out[which(out$dataset_name==names(datalist[i]) & (out$clustering_method=="HDBSCAN" | out$clustering_method=="SSMDH" )),]
      #没有sd显示的话，全实线无点的折线图就行，不过要自定义配色并且把坐标轴题目等正过来，
      print(ggplot(out_now, aes(x=Link_num, y=RI,linetype = clustering_method,  color= clustering_method, group = clustering_method))#, shape = clustering_method
            + geom_line(linewidth = 1)
            # + geom_point(size = 3)
            + labs(x="Constraints Number",#Constraints Number
                   y = "RI",
                   title = names(datalist[i]))
            + theme_bw() 
            + theme(legend.position = "bottom",
                    legend.direction = "horizontal",
                    legend.title = element_blank(),
                    panel.grid = element_blank(),
                    plot.title = element_text(hjust = 0.5))
            + theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
            + scale_x_continuous(breaks = consnum_list)
            + scale_color_manual(values = c("#16048A","#6201A9","#9E189F","#CC4A76","#EA7854","#FDB330","#f8e620")) #自定义颜色
            + ylim(0,1))
      scriptpath = rstudioapi::getSourceEditorContext()$path
      scriptpathn = nchar(scriptpath)
      suppath = substr(scriptpath, 1, scriptpathn-12)
      datasetpath = paste(suppath,"/dataset/",sep="")
      ggsave(paste(datasetpath,names(datalist[i]),"_HDBSCAN_results.pdf", sep=""),device = cairo_pdf,width =5, height =3.75)
      
    }
  }
  
  # save(real_out, file=paste(suppath,names(datalist[i])))
  save(out, file=paste(datasetpath, "wxgca_results.RData", sep=""))
  #对于rdata文件，读取时是get(load(xxx))
  return(out)
}
wxgca = function(data, datalabel, ncluster,idnum=NULL,idbin=NULL,idcat=NULL,mustLink=NULL, cantLink=NULL){
  #获得ncluster个密度最小的点
  if(is.null(idbin) && is.null(idcat)){
    distdata = as.matrix(dist(data, method = "euclidean"))
  } else {
    distdata = distmix(data, method = "gower", idnum = idnum, idbin = idbin, idcat = idcat)#还好，里面对数值变量自带归一化
  }
  
  original_distdata = distdata
  # densitymaxlist = rowSums(distdata)
  # densitypoint = c()
  # for (i in c(1:ncluster)) {
  #   densitypoint = c(densitypoint, which(densitymaxlist == max(densitymaxlist))[1])
  #   densitymaxlist[which(densitymaxlist == max(densitymaxlist))[1]] = min(densitymaxlist)
  # }
  
  
  #pro part
  # if(FALSE){
  #   n = nrow(distdata)
  #   m = nrow(mustLink)
  #   c = nrow(cantLink)
  #   promatrix = matrix(0,n,n)
  #   for (i in c(1:m)) {promatrix[mustLink[i,1], mustLink[i,2]] = promatrix[mustLink[i,2], mustLink[i,1]] = 1}
  #   for (i in c(1:c)) {promatrix[cantLink[i,1], cantLink[i,2]] = promatrix[cantLink[i,2], cantLink[i,1]] = 2}
  #   mustpoint = unique(c(mustLink))
  #   mustLink = matrix(,0,2)
  #   cantLink = matrix(,0,2)
  #   while (length(mustpoint) > 0) {
  #     pro = used = subi = c(mustpoint[1]) 
  #     mustpoint = mustpoint[-1]
  #     cant = c()
  #     while (length(pro) > 0) {
  #       for (i in c(1:n)) {
  #         if (promatrix[pro[1],i]==1 && !i%in%used) {
  #           pro = c(pro,i); used = c(used,i); subi = c(subi,i); mustpoint=mustpoint[-which(mustpoint==i)]
  #         }else if (promatrix[pro[1],i]==2) {
  #           cant = c(cant,i)
  #         }
  #       }
  #       pro = pro[-1]
  #       if (length(pro)==0 && length(subi)>0){
  #         tmustLink = t(combn(subi, 2))
  #         mustLink = rbind(mustLink, tmustLink)
  #       }
  #       if (length(pro)==0 && length(subi)>0 && length(cant)>0){
  #         tcantLink = matrix(unlist(expand.grid(subi,cant)), ncol=2)
  #         cantLink = rbind(cantLink, tcantLink)
  #       }
  #     }
  #   }
  # }
  
  #就是我那个示意图作为测试集,总的过程类似一个必连约束的多图的深度遍历，先遍历完一个子图，然后子图的结点两两间都是必连，将新生成的必连放入必连集中
  #找子图结点相关的勿连，将对应的勿连点和该子图的所有结点生成新的勿连对加入勿连集合
  #检查是否遍历完所有必连点了，没有的话开始遍历下一个子图
  # mustLink=matrix(c(1,2,5,4,2,3,6,5),nrow=4) #测试用 6g
  # cantLink=matrix(c(4,2,3,3,7,7),nrow=3) #9g
  mustLink_vector = c(mustLink)#必连向量，代替必连集
  mustLink_vector_length = length(mustLink_vector)
  cantLink_vector = c(cantLink)
  cantLink_vector_length = length(cantLink_vector)
  unvisited_queue = unique(mustLink_vector)
  visiting_queue = c()
  o_cantLink= cantLink #我说怎么约束效果不佳，合着我的cl只用了新生成的，旧的全仍了
  mustLink= matrix(NA,0,2)#不能和原总配对集做并集，必须更新
  cantLink= matrix(NA,0,2)
  #主要数据结构有：未访问队列，待访问队列，当前访问点，当前子图，当前子图的必连集
  while (length(unvisited_queue) > 0) {#因为多图
    current_subgraph_points = c() #新建/清空子图存储器
    visiting_queue = c(visiting_queue,unvisited_queue[1]) #将未访问点的第一个插入待访问队列
    while (length(visiting_queue) > 0) {#当前子图的深度优先遍历，将当前子图遍历后i.e.待访问队列清空才停止，注意不能用!is.null(visiting_queue)，因为删除完是empty，不算null
      visiting_point = visiting_queue[1] #将待访问队列的第一个元素作为当前访问点
      unvisited_queue = setdiff(unvisited_queue,visiting_queue) #将待访问队列从未访问队列中删除
      current_subgraph_points = c(current_subgraph_points,visiting_queue[1]) #将待访问队列的第一个元素放入当前子图
      visiting_queue = visiting_queue[-1]#将待访问队列的第一个元素删除
     
     
      visiting_point_position = which(mustLink_vector == visiting_point) #找到当前访问点在必连向量中的全部位置
      # updata_points_position = visiting_point_position
      if (length(visiting_point_position)>0) {
        corresponding_points = c()#当前访问点的对应点集
        setdiff_mustlink_vector = visiting_point_position
        for (i in c(1:length(visiting_point_position))) {#遍历当前访问点的所有必连约束
          point_position = visiting_point_position[i]
          if (point_position <= mustLink_vector_length/2) {#用必连向量从位置找到对应点
            corresponding_point_postion = point_position+mustLink_vector_length/2
            corresponding_point = mustLink_vector[corresponding_point_postion]
            corresponding_points = c(corresponding_points,corresponding_point)
            setdiff_mustlink_vector = c(setdiff_mustlink_vector,corresponding_point_postion)
            # updata_points_position = c(updata_points_position, visiting_point_position[i]+mustLink_vector_length/2)
          } else if (point_position > mustLink_vector_length/2) {
            corresponding_point_postion = point_position-mustLink_vector_length/2
            corresponding_point = mustLink_vector[corresponding_point_postion]
            corresponding_points = c(corresponding_points,corresponding_point)
            setdiff_mustlink_vector = c(setdiff_mustlink_vector,corresponding_point_postion)
            # updata_points_position = c(updata_points_position, visiting_point_position[i]-mustLink_vector_length/2)
          }
          # current_subgraph_points = c(current_subgraph_points,corresponding_point)#将对应点们加入当前子图-不需要，遍历中将当前访问点放入当前子图就能完成这一步
          # visiting_queue = unique(visiting_queue)
          # unvisited_queue = setdiff(unvisited_queue,corresponding_point)#将对应点们从未访问队列中删除-不需要，同上在遍历中会从未访问队列中删除
        }
        # corresponding_points = intersect(corresponding_points, unvisited_queue)#因为必连向量不更新，所以要防止遍历过的点作为当前访问点的必连约束的对应点又插入待访问序列，但是这个没用，因为未访问点只包含当前访问点集合，在当前访问序列里的没有排除。比如第一波1对应点345那么未=2345在序列=345子图1，第二波，在点=3然后3对应45导致在序列=34545，这不仅导致多了很多次循环找点还会一个个多放入子图
        # visiting_queue = c(visiting_queue,setdiff(corresponding_points,visiting_queue) )#将对应点们插入待访问队列
        
        visiting_queue = c(visiting_queue,corresponding_points)
        mustLink_vector = mustLink_vector[-setdiff_mustlink_vector]
        mustLink_vector_length = length(mustLink_vector)
      }
      
      
    }
    #当前子图结点通过两两配对构成完全图，从完全图中提取所有边作为当前子图必连集，加入总必连集
    current_subgraph_mustlink_set = t(combn(current_subgraph_points, 2))
    mustLink = rbind(mustLink, current_subgraph_mustlink_set)
    #经验：有进有出或者进出不是一个个地固定地的遍历用while和一个queue，一个个地固定地只出不进的遍历用for从1:n就行
    co_points = intersect(cantLink_vector, current_subgraph_points)#通过交集运算找到勿连向量与当前子图结点的共同点集（如果有）
    corresponding_cantlink_points = c() #新建/清空对应勿连点集
    if (length(co_points)>0) {#判断是否有对应勿连点子集，两个空集交集为NULL，如果没有交集为numeric(0)，但都长度=0
      for (i in c(1:length(co_points))) {#遍历共同点集用勿连向量找到每个共同点的对应勿连点（，因此可能多个对应勿连点，）
        visiting_point = co_points[i]
        visiting_point_position = which(cantLink_vector == visiting_point) #因为是基于勿连的共同点，至少存在一个对应勿连点
        for (j in c(1:length(visiting_point_position))) {#每个共同点可能有多个勿连约束，因此可能有多个对应勿连点
          if (visiting_point_position[j] <= cantLink_vector_length/2) {
            corresponding_point = cantLink_vector[ visiting_point_position[j]+cantLink_vector_length/2 ]
          } else if (visiting_point_position[j] > cantLink_vector_length/2) {
            corresponding_point = cantLink_vector[ visiting_point_position[j]-cantLink_vector_length/2 ]
          }
          corresponding_cantlink_points = union(corresponding_cantlink_points, corresponding_point)#共同点的对应勿连点还可能有重复的，所以要以并集方式将对应勿连点放入对应勿连点集，好在不是必连那里因为必连-必连约束要深度优先遍历，因此不用搞个标记队列
        }
      }
      #对应勿连点集中每个对应勿连点都和当前子图结点构成勿连约束，并将之加入总勿连集
      for (i in c(1:length(corresponding_cantlink_points))) {
        for (j in c(1:length(current_subgraph_points))) {
          cantLink = rbind(cantLink, matrix(c(corresponding_cantlink_points[i],current_subgraph_points[j]),nrow=1))
        }
      }
    }
  }
  cantLink = rbind(o_cantLink, cantLink)
  #注意排除重复的
  for (i in c(1:nrow(cantLink))) {
    if (cantLink[i,1] > cantLink[i,2]) {
      tem = cantLink[i,1]
      cantLink[i,1] = cantLink[i,2] 
      cantLink[i,2] = tem
    }
  }
  tem_cantlink = data.frame(cantLink)
  tem_cantlink = unique(tem_cantlink)
  cantLink = as.matrix(tem_cantlink, ncol=2)
  #我首先发现约束后距离矩阵值过大，然后在实际数据中发现mustlinkcantlink有大量重复，看看这个重复出现的位置是抽样还是传播这里
  #检查的代码用mustLink[which(mustLink==100,arr.ind = T)[,1],1]和mustLink[which(mustLink==100,arr.ind = T)[,1],2]，如果一一对应没有重复就没问题
  
  #测试
  # distdata = as.matrix(dist(matrix(c(0, 1, 1, 0, 0, 0, 1, 1), nrow = 4)))
  # mustLink = matrix(c(1, 2), nrow = 1)
  # cantLink = matrix(c(1, 4), nrow = 1)
  #constrained part
  m = nrow(mustLink)
  c = nrow(cantLink)
  n = nrow(distdata) 
  distdata_old = distdata_new = distdata
  # maxdist = distdata_new[which(distdata==0 ,arr.ind = TRUE)] = max(distdata_new)
  # mindist = min(distdata_new)
  if(m>0){
    for(i in c(1:m)){
      distdata[mustLink[i,2], mustLink[i,1]] = distdata[mustLink[i,1], mustLink[i,2]] = distdata[mustLink[i,1], mustLink[i,2]]/2#mindist/n
    }
  }
  
  if(c>0){
    for(i in c(1:c)){
      distdata[cantLink[i,2], cantLink[i,1]] = distdata[cantLink[i,1], cantLink[i,2]] = distdata[cantLink[i,1], cantLink[i,2]]*2#maxdist*n
    }
  }
  
  #clustering 调参这里只能以一个评价指标为准，所以结果呈现一般只出一个评价指标
  # ssmdk = indicator(datalabel, fastkmed(distdata, ncluster)$cluster, distdata)$ri
  # ssmdka = indicator(datalabel, fastkmed(distdata, ncluster, init = densitypoint )$cluster, distdata)$ri #默认10迭代
  
  # k_list = seq(from = 2, to = n/2)
  # max_result = -2 #ARI和SIL范围-1到1，RI是0-1，
  # for (i in k_list){
  #   ssmdk = tryCatch(
  #     {
  #       ssmdk = indicator(datalabel, fastkmed(distdata, i)$cluster, distdata)$ri
  #     },
  #     error = function(e){0}
  #   )
  #   if(ssmdk > max_result){
  #     max_result = ssmdk
  #   }
  # }
  # ssmdk = max_result
  
  ssmdk = tryCatch(
    {
      ssmdk = indicator(datalabel, fastkmed(distdata, ncluster)$cluster, distdata)$ri
    },
  error = function(e){0}
  )
  
  # ssmdka = tryCatch(
  #   {
  #     ssmdka = indicator(datalabel, fastkmed(distdata, ncluster, init = densitypoint)$cluster, distdata)$ri
  #   }, 
  #   error = function(e){0}
  # )
  
  #如果算法输出是乱序输出的，那就无法调参
  minpts_list = seq(from = 2, to = n/ncluster)
  max_result = -2 #ARI和SIL范围-1到1，RI是0-1，
  for (i in minpts_list){
    #过度拉伸的距离反而导致算法报错和结果出错，甚至在ssmdd中出现评估全是噪声点，所以拉大肯定要控制，拉近其实也要控制
    ssmdh = hdbscan(distdata, i)$cluster
    ssmdh = indicator(datalabel, ssmdh, distdata)$ri
    if(ssmdh > max_result){
      max_result = ssmdh
      minpts_h = i
    }
  }
  ssmdh = max_result

  # minpts_list = seq(from = 2, to = n/2)
  # max_result = -2 #ARI和SIL范围-1到1，RI是0-1，
  # for (i in minpts_list){
  #   #过度拉伸的距离反而导致算法报错和结果出错，甚至在ssmdd中出现评估全是噪声点，所以拉大肯定要控制，拉近其实也要控制
  #   ssmdh = hdbscan(distdata, i)$cluster
  #   
  #   if(ssmdh > max_result){
  #     max_result = ssmdh
  #     minpts_h = i
  #   }
  # }
  # ssmdh = indicator(datalabel, ssmdh, distdata)$ri
  
  # ssmdh = hdbscan(distdata, 3)$cluster#这是不调参的
  # ssmdh = indicator(datalabel, ssmdh, distdata)$ri
  
  # ssmdh = hdbscan(distdata, minpts_h)$cluster
  # ssmdh = indicator(datalabel, ssmdh, distdata)$ri
  
  #增加minpts的调参对结果几乎没改变，反而增加了好几倍的整体运行时间，minpts应该包括了核心点本身
  #其实我的半监督的可以用已知约束对算法进行调参
  
  
  
  # for (i in minpts_list){
  #   for (j in eps_list){
  #     ssmdd = indicator(datalabel, dbscan(distdata, j, i, weights = NULL, borderPoints = TRUE)$cluster, distdata)$ri
  #     if(ssmdd > max_result){
  #       max_result = ssmdd
  #       minpts_d = i
  #       eps = j
  #     }
  #   }
  # }
  
  
  max_value = max(distdata)
  distdata0 = distdata
  distdata0[which(distdata0==0)] = max_value
  min_value = min(distdata0)
  minpts_list = seq(from = 2, to = n/ncluster)
  eps_list = seq(from = min_value, to = max_value, length.out = 100)#
  max_result = -2
  for (j in eps_list){#不用调参mintps后整个算法都快了有10倍
    ssmdd = indicator(datalabel, dbscan(distdata, j, minPts = 5)$cluster, distdata)$ri
    if(ssmdd > max_result){
      max_result = ssmdd
    }
  }
  ssmdd = max_result
  
  # #找k-dist value的threshold
  # minpts_list = seq(from = 2, to = n/2)
  # max_result = -2
  # for (j in minpts_list){#不用调参mintps后整个算法都快了有10倍
  #   k_dist = c()
  #   for (i in c(1:n)) {
  #     k_dist[i] = sum(sort(distdata[i,])[1:j]) #每个点k近邻之和并从排序
  #   }
  #   k_dist = sort(k_dist)
  #   k_dist_thr = max(k_dist[2:n]-k_dist[1:(n-1)])
  #   ssmdd = indicator(datalabel, dbscan(distdata, k_dist_thr, j)$cluster, distdata)$ri
  #   if(ssmdd > max_result){
  #     max_result = ssmdd
  #   }
  # }
  # ssmdd = max_result
  
  # ssmdd = indicator(datalabel, dbscan(distdata, k_dist_thr, minPts = minpts$cluster, distdata)$ri
  
  # ssmdd = indicator(datalabel, dbscan(distdata, eps, minPts = 5, weights = NULL, borderPoints = TRUE)$cluster, distdata)$ri
  
  
  wxgca_result = list(
    ssmdk = ssmdk,
    # ssmdka = ssmdka,
    ssmdd = ssmdd,
    ssmdh = ssmdh
    # minpts_h = minpts_h,
    # minpts_d = minpts_d,
    # eps = eps
  )
  
  return(wxgca_result)
}

indicator = function(datalabel, resultlabel, datadist){
  #DBSCAN和Hdbscan输出可能乱序，sil测得不对，此外，k变体不用调参直接设置几个簇就行，因此也不涉及用对的标签调参，我需要看看DBSCAN之类的聚类算法怎么设计实验的
  tab = table(datalabel, resultlabel)
  ri = randIndex(tab, correct=TRUE, original=TRUE)[2]
  ari = randIndex(tab, correct=TRUE, original=TRUE)[1]
  #注意要算sil得区分是单纯数据data1产生的distdata1还是混合变量data2产生的distdata2
  sil = mean(silhouette(resultlabel, datadist))
  res = list(ri=ri,ari=ari,sil=sil)
  return(res)
}

data_initialize = function(){
  #注意标签值不能为负或者0因为比如dbscan里0表示噪声
  scriptpath = rstudioapi::getSourceEditorContext()$path #"D:/科研/混合密度聚类算法小论文/文档/MDPI模板/WXGCA.R"
  scriptpathn = nchar(scriptpath)
  suppath = substr(scriptpath, 1, scriptpathn-12) #"D:/科研/混合密度聚类算法小论文/文档/MDPI模板"
  datasetpath = paste(suppath,"/dataset/",sep="") #"D:/科研/混合密度聚类算法小论文/文档/MDPI模板/dataset/"
  datasetmaking()
  load(file=paste(datasetpath, "overlapping_clusters.RData", sep=""))
  load(file=paste(datasetpath, "noise20_100_clusters.RData", sep=""))
  load(file=paste(datasetpath, "halfring_clusters.RData", sep=""))
  load(file=paste(datasetpath, "ring_clusters.RData", sep=""))
  
  data = overlapping_clusters
  data[which(data=="cluster 1", arr.ind = TRUE)] = 1
  data[which(data=="cluster 2", arr.ind = TRUE)] = 2
  data[,4] = as.numeric(data[,4])
  datalabel = data[,4]
  data = data[,1:3]
  idnum = c(1,2)
  idbin = c(3)
  idcat = NULL
  OVERLAPPING = list(
    data = data,
    datalabel = datalabel,
    idnum = idnum,
    idbin = idbin,
    idcat = idcat,
    dataname = c("OVERLAPPING")
  )
  
  data = noise20_100_clusters
  data[which(data=="cluster", arr.ind = TRUE)] = 1
  # data[which(data=="cluster 2", arr.ind = TRUE)] = 2
  data[which(data=="noise", arr.ind = TRUE)] = 2
  data[,3] = as.numeric(data[,3])
  datalabel = data[,3]
  data = data[,1:2]
  idnum = c(1,2)
  idbin = NULL
  idcat = NULL
  NOISE = list(
    data = data,
    datalabel = datalabel,
    idnum = idnum,
    idbin = idbin,
    idcat = idcat,
    dataname = c("NOISE")
  )
  
  data = ring_clusters
  data[which(data=="cluster 1", arr.ind = TRUE)] = 1
  data[which(data=="cluster 2", arr.ind = TRUE)] = 2
  data[,3] = as.numeric(data[,3])
  datalabel = data[,3]
  data = data[,1:2]
  idnum = c(1,2)
  idbin = NULL
  idcat = NULL
  RING = list(
    data = data,
    datalabel = datalabel,
    idnum = idnum,
    idbin = idbin,
    idcat = idcat,
    dataname = c("RING")
  )
  
  data = halfring_clusters
  data[which(data=="cluster 1", arr.ind = TRUE)] = 1
  data[which(data=="cluster 2", arr.ind = TRUE)] = 2
  data[,3] = as.numeric(data[,3])
  datalabel = data[,3]
  data = data[,1:2]
  idnum = c(1,2)
  idbin = NULL
  idcat = NULL
  HALFRING = list(
    data = data,
    datalabel = datalabel,
    idnum = idnum,
    idbin = idbin,
    idcat = idcat,
    dataname = c("HALFRING")
  )
  
  data = datasets::iris
  data$Species = as.numeric(data$Species)
  datalabel = data$Species
  data = data[,1:4]
  idnum = c(1:4)
  idbin = NULL
  idcat = NULL
  IRIS = list(
    data = data,
    datalabel = datalabel,
    idnum = idnum,
    idbin = idbin,
    idcat = idcat,
    dataname = c("IRIS")
  )
  #http://archive.ics.uci.edu/dataset/53/iris
  
  data = wine
  datalabel = as.numeric(data$Wine) 
  data = data[,1:13]
  idnum = c(1:13)
  idbin = NULL
  idcat = NULL
  WINE = list(
    data = data,
    datalabel = datalabel,
    idnum = idnum,
    idbin = idbin,
    idcat = idcat,
    dataname = c("WINE")
  )
  #http://archive.ics.uci.edu/dataset/109/wine
  
  fertility_Diagnosis = read.csv(paste(datasetpath, "fertility_Diagnosis.txt", sep=""), header=F)
  data = fertility_Diagnosis
  data[which(data=="N", arr.ind = TRUE)] = 1
  data[which(data=="O", arr.ind = TRUE)] = 2
  datalabel = as.numeric(data[,10])
  data = data[,1:9]
  idnum = c(2,9)
  idbin = c(3,4,5)
  idcat = c(1,6,7,8)
  FERTILITY = list(
    data = data,
    datalabel = datalabel,
    idnum = idnum,
    idbin = idbin,
    idcat = idcat,
    dataname = c("FERTILITY")
  )
  #http://archive.ics.uci.edu/dataset/244/fertility
  
  zoo = read.csv(paste(datasetpath, "zoo.data", sep=""), header=F)
  data = zoo[,2:18]
  datalabel = as.numeric(data[,17])
  data = data[,1:16]
  idnum = NULL
  idbin = c(1:16)
  idcat = NULL
  ZOO = list(
    data = data,
    datalabel = datalabel,
    idnum = idnum,
    idbin = idbin,
    idcat = idcat,
    dataname = c("ZOO")
  )
  #http://archive.ics.uci.edu/dataset/111/zoo
  
  tae = read.csv(paste(datasetpath, "tae.data", sep=""), header=F)
  data = tae
  datalabel = as.numeric(data[,6])
  data = data[,1:5]
  idnum = c(5)
  idbin = c(1,4)
  idcat = c(2,3)
  TAE = list(
    data = data,
    datalabel = datalabel,
    idnum = idnum,
    idbin = idbin,
    idcat = idcat,
    dataname = c("TAE")
  )
  #http://archive.ics.uci.edu/dataset/100/teaching+assistant+evaluation
  
  flag = read.csv(paste(datasetpath, "flag.data", sep=""), header=F)
  data = flag[,c(2:ncol(flag))]
  type18 = unique(flag[,18])
  type29 = unique(flag[,29])
  type30 = unique(flag[,30])
  #可以发现颜色通用，所以用一个循环就行
  for(i in c(1:length(type18))){
    data[which(data==type18[i], arr.ind = TRUE)] = i
  }
  # for(i in c(1:length(type29))){
  #   data[which(data==type29[i], arr.ind = TRUE)] = i
  # }
  # for(i in c(1:length(type30))){
  #   data[which(data==type30[i], arr.ind = TRUE)] = i
  # }
  datalabel = data[,6] + 1
  for(i in c(1:29)){data[,i] = as.numeric(data[,i])}
  idnum = c(3,4,7,8,9,18,19,20,21,22)
  idbin = c(10:16,23:27)
  idcat = c(1,2,5,17,28,29)#6是宗教当作预测标签了
  FLAG = list(
    data = data,
    datalabel = datalabel,
    idnum = idnum,
    idbin = idbin,
    idcat = idcat,
    dataname = c("FLAG")
  )
  
  datalist = list(
    RING=RING,
    HALFRING=HALFRING,
    
    OVERLAPPING=OVERLAPPING,
    NOISE=NOISE,
    #在纯数值数据集比如iris和wine上表现就不如那些半监督聚类算法，在混合变量数据集上表现又比他们好了
    WINE=WINE,
    IRIS=IRIS,
    
    # ZOO=ZOO, 因为只有二值变量，没混合公式算法没结果的，注意跳过
    #我的算法表现不稳定是一个大问题
    # FLAG=FLAG,
    
    FERTILITY=FERTILITY,
    TAE=TAE
  )
  return(datalist)
}

datasetmaking = function() {
  #实际上，那种边界分明的数据集半监督反而没效果
  scriptpath = rstudioapi::getSourceEditorContext()$path
  scriptpathn = nchar(scriptpath)
  suppath = substr(scriptpath, 1, scriptpathn-12)#注意这个脚本名字的字符数改变，这里也要变，否则会导致程序报错 
  datasetpath = paste(suppath,"/dataset/",sep="")
  #rnorm(20, mean = 465, sd = 10)是正态分布，这种样本集对于均质密度聚类DBSCAN是灾难对HBDSCAN也不好，而且点很大且固定比例，点数量多了根本判断不出是否均密度
  #要用runif函数生成在某个区间内的均匀随机数，比如原本是20均值10sd，那么就该为10-30均匀，而且均匀范围越小点数越多越好避免密度不均
  #runif也没有均质问题,加个round也无大勇，对于簇间重合问题，一个簇加个0.5就行以及扩大范围,用sample问题是不均质
  #改成均值+边界模糊确实原算法差而改进算法好
  # Generate data set with normal distribution and overlap
  if (TRUE) {
    set.seed(1)
    # x1 = x2 = y1 = y2 = c()
    # for (i in c(1:10)) {
    #   x1 = c(x1,rep(i,10))
    #   y1 = c(y1,c(1:10))
    #   x2 = c(x2,rep(i+2.5,10))
    #   y2 = c(y2,c(3.5:12.5))
    # }

    # x1 = rnorm(50, mean = 60, sd = 1)
    # y1 = rnorm(50, mean = 60, sd = 1)
    x1 = round(runif(50, 200, 700))
    y1 = round(runif(50, 200, 700))
    # # x = sample(c(200:700),100)
    # # y = sample(c(200:700),100)
    # # x1 = sample(x,50)
    # # y1 = sample(y,50)
    # # x2 = intersect(x,x1)
    # # y2 = intersect(y,y1)
    
    # x2 = rnorm(50, mean = 60, sd = 1)
    # y2 = rnorm(50, mean = 60, sd = 1)
    x2 = round(runif(50, 200, 700))+0.5
    y2 = round(runif(50, 200, 700))+0.5
    z1 = rep(100,length(x1))
    z2 = rep(1,length(x2))
    cluster1 = data.frame(x = x1, y = y1, z = z1, label = "cluster 1")
    cluster2 = data.frame(x = x2, y = y2, z = z2, label = "cluster 2")
    overlapping_clusters <- rbind(cluster1, cluster2)
    #search overlap points
    equalidx = c()
    distdata = as.matrix(dist(overlapping_clusters[,c(1:2)]))
    equalist = which(distdata==0, arr.ind = TRUE)
    for(i in c(1:nrow(equalist))){
      if(equalist[i,1] != equalist[i,2]){
        equalidx = c(equalidx, equalist[i,1], equalist[i,2])
      }
    }
    equalidx
    save(overlapping_clusters, file=paste(datasetpath, "overlapping_clusters.RData", sep=""))
    
    #3D plot
    svg(paste(datasetpath, "overlapping_clusters_3D.svg", sep=""))
    # pdf(paste(datasetpath, "overlapping_clusters_3D.pdf", sep=""))
    # svg(paste(paste("D:/research/master/graduation oral examination/", "overlapping_clusters_3D.svg", sep=""), sep=""))
    pmar <- par(mar = c(1, 1, 1, 1))
    colors = c(rep("#16048A", length(x1)), rep("#9E189F", length(x2)))
    with(overlapping_clusters, scatter3D(x = x, y = y, z = z,
                                         pch = 21, cex = 1.5,col="black",bg=colors,
                                         xaxt = F,
                                         xlab = NA, ylab = NA,
                                         zlab = NA,#"X" "Y" "Z"
                                         ticktype = "simple",bty = "f",box = TRUE,
                                         #panel.first = panelfirst,
                                         theta = 140, phi = 20, d=3,
                                         colkey = FALSE)#list(length = 0.5, width = 0.5, cex.clab = 0.75))
    )
    # legend("bottom",title =  "",legend=c("cluster 1", "cluster 2"),pch=21,
    #        cex=1,y.intersp=1.5,pt.bg = colors,bg="white",bty="n")
    # pdf.options(reset = TRUE)
    dev.off()#把这个删掉就不显示图片了
   
    
    #2D plot
    # print(ggplot(data = overlapping_clusters, mapping = aes(x = x, y = y, shape = label, color = label)) + geom_point(size = 2) + labs(title = "OVERLAPPING", x ="X", y = "Y") + theme_bw() + theme(legend.position = "bottom",axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 45),axis.text.y = element_text(vjust = 0.5, hjust = 0.5, angle = 45), panel.grid=element_blank(), plot.title = element_text(hjust = 0.5), legend.direction = "horizontal",legend.title = element_blank()))
    print(ggplot(data = overlapping_clusters, mapping = aes(x = x, y = y,fill = as.factor(label))) 
          + geom_point(colour = "white",size = 3,shape = 21) 
          + labs(title = "OVERLAPPING", x ="X", y = "Y") #
          + theme_bw() 
          + theme(legend.position = "bottom",axis.text.x = element_text(vjust = 0.5, hjust = 0.5),axis.text.y = element_text(vjust = 0.5, hjust = 0.5))
          + theme(panel.grid=element_blank(), plot.title = element_text(hjust = 0.5), legend.direction = "horizontal",legend.title = element_blank())
          + scale_fill_manual(values = c("#16048A","#9E189F"))
          # + coord_fixed(ratio = 0.6)
          + theme(axis.title.y = element_text(angle = 0, vjust = 0.5)))
    ggsave(paste(datasetpath, "overlapping_clusters.pdf", sep=""),width =5, height =3.75)
    # ggsave(file = paste("D:/research/master/graduation oral examination/", "overlapping_clusters.svg", sep=""),font="Times New Roman")
    # svglite(filename  = paste("D:/research/master/graduation oral examination/", "overlapping_clusters.svg", sep=""),system_fonts ="SimSun")
    # graph2svg(file = paste("D:/research/master/graduation oral examination/", "overlapping_clusters.svg", sep=""),font="Times New Roman")
  }
  
  # Generate data set with normal distribution and noise 噪声点要多边界要不分明-我算法才凸显作用
  if (TRUE) {
    set.seed(1)
    # generateRingShapePointsbyRnom = function(xmean,ymean,num,r,sd,label){
    #   cluster = vector()
    #   for (i in 1:num) {
    #     angle = i * 18
    #     rest_angle = angle%%180
    #     if ( (rest_angle>=0 && rest_angle<90)) {
    #       r = r + 10
    #     } else {
    #       r = r - 5
    #     }
    #     # x = rnorm(1, mean = r * cos(angle) + xmean, sd = sd)
    #     # y = rnorm(1, mean = r * sin(angle) + ymean, sd = sd)
    #     x = runif(1, r * cos(angle) + xmean - sd, r * cos(angle) + xmean + sd)
    #     y = runif(1, r * sin(angle) + ymean - sd, r * sin(angle) + ymean + sd)
    #     cluster = rbind(cluster, data.frame(x,y,label))
    #   }
    #   return(cluster)
    # }
    # noise <- generateRingShapePointsbyRnom(500,500,20,50,1,'noise')
    # # x1 <- rnorm(20, mean = 465, sd = 10)
    # # y1 <- rnorm(20, mean = 510, sd = 10)
    # # x2 <- rnorm(20, mean = 515, sd = 10)
    # # y2 <- rnorm(20, mean = 520, sd = 10)
    # x1 = runif(40, 440, 490)
    # y1 = runif(40, 480, 530)
    # x2 = runif(40, 500, 550)
    # y2 = runif(40, 480, 530)
    
    # cluster1 <- data.frame(x = x1, y = y1, label = 'cluster 1')
    # cluster2 <- data.frame(x = x2, y = y2, label = 'cluster 2')
    # noise20_100_clusters <- rbind(cluster1, cluster2, noise)
    
    generateRingShapePointsbyRnom = function(xmean,ymean,num,r,sd,label){
      cluster = vector()
      for (i in 1:num) {
        angle = i * 18
            rest_angle = angle%%180
            if ( (rest_angle>=0 && rest_angle<90)) {
              r = r + 10
            } else {
              r = r - 5
            }
        x = round(runif(1, r * cos(angle) + xmean - sd, r * cos(angle) + xmean + sd))
        y = round(runif(1, r * sin(angle) + ymean - sd, r * sin(angle) + ymean + sd))
        cluster = rbind(cluster, data.frame(x,y,label))
      }
      return(cluster)
    }
    noise <- generateRingShapePointsbyRnom(500,500,20,100,20,'noise')
    x1 = round(runif(50, 460, 550))
    y1 = round(runif(50, 460, 550))
    cluster <- data.frame(x = x1, y = y1, label = 'cluster')
    noise20_100_clusters <- rbind(cluster, noise)
    
    equalidx = c()
    distdata = as.matrix(dist(noise20_100_clusters[,c(1:2)]))
    equalist = which(distdata==0, arr.ind = TRUE)
    for(i in c(1:nrow(equalist))){
      if(equalist[i,1] != equalist[i,2]){
        equalidx = c(equalidx, equalist[i,1], equalist[i,2])
      }
    }
    equalidx
    save(noise20_100_clusters, file=paste(datasetpath, "noise20_100_clusters.RData", sep=""))
    
    print(ggplot(data = noise20_100_clusters, mapping = aes(x = x, y = y,fill = as.factor(label))) 
          + geom_point(colour = "white",size = 3,shape = 21) 
          + labs(title = "NOISE", x ="X", y = "Y") 
          + theme_bw() 
          + theme(legend.position = "bottom",axis.text.x = element_text(vjust = 0.5, hjust = 0.5),axis.text.y = element_text(vjust = 0.5, hjust = 0.5))
          + theme(panel.grid=element_blank(), plot.title = element_text(hjust = 0.5), legend.direction = "horizontal",legend.title = element_blank())
          + scale_fill_manual(values = c("#16048A","#9E189F","#FDB330"))
          + theme(axis.title.y = element_text(angle = 0, vjust = 0.5)))
    ggsave(paste(datasetpath, "noise20_100_clusters.pdf", sep=""),device = cairo_pdf,width =5, height =3.75)
  }
  
  # Generate Ring data set
  if (TRUE) {
    set.seed(1)
    # generatePointsByRnom <- function(xmean, ymean, sd, num, label) {
    #   x <- rnorm(num, mean = xmean, sd = sd)
    #   y <- rnorm(num, mean = ymean, sd = sd)
    #   data.frame(x, y, label)
    # }
    # generateRingShapePointsbyRnom <- function(r, class) {
    #   cluster = vector()
    #   for (i in 1:50) {
    #     angle = i * 18
    #     x = r * cos(angle)
    #     y = r * sin(angle)
    #     cluster <- rbind(cluster, generatePointsByRnom(x+50, y+50, 1, 1, class))
    #   }
    #   cluster
    # }
    # cluster1 <- generateRingShapePointsbyRnom(10, 'cluster 1')
    # cluster2 <- generateRingShapePointsbyRnom(20, 'cluster 2')
    # ring_clusters <- rbind(cluster1, cluster2)
    
    generateRingShapePointsbyRnom = function(xmean,ymean,num,r,sd,label){
      cluster = vector()
      for (i in 1:num) {
        angle = i * 7
        # x = rnorm(1, mean = r * cos(angle) + xmean, sd = sd)
        # y = rnorm(1, mean = r * sin(angle) + ymean, sd = sd)
        x = runif(1, r * cos(angle) + xmean - sd, r * cos(angle) + xmean + sd)
        y = runif(1, r * sin(angle) + ymean - sd, r * sin(angle) + ymean + sd)
        cluster = rbind(cluster, data.frame(x,y,label))
      }
      return(cluster) 
    }
    cluster1 <- generateRingShapePointsbyRnom(50,50,100,5,1,'cluster 1')
    cluster2 <- generateRingShapePointsbyRnom(50,50,100,15,1,'cluster 2')
    ring_clusters <- rbind(cluster1, cluster2)
    # cluster1$x = cluster1$x
    # cluster1$y = cluster1$y
    # cluster2$x = cluster2$x
    # cluster2$y = cluster2$y
    equalidx = c()
    distdata = as.matrix(dist(noise20_100_clusters[,c(1:2)]))
    equalist = which(distdata==0, arr.ind = TRUE)
    for(i in c(1:nrow(equalist))){
      if(equalist[i,1] != equalist[i,2]){
        equalidx = c(equalidx, equalist[i,1], equalist[i,2])
      }
    }
    equalidx
    save(ring_clusters, file=paste(datasetpath, "ring_clusters.RData", sep=""))
    
    print(ggplot(data = ring_clusters, mapping = aes(x = x, y = y,fill = as.factor(label))) 
          + geom_point(colour = "white",size = 3,shape = 21) 
          + labs(title = "RING", x ="X", y = "Y") 
          + theme_bw() 
          + theme(legend.position = "bottom",axis.text.x = element_text(vjust = 0.5, hjust = 0.5),axis.text.y = element_text(vjust = 0.5, hjust = 0.5))
          + theme(panel.grid=element_blank(), plot.title = element_text(hjust = 0.5), legend.direction = "horizontal",legend.title = element_blank())
          + scale_fill_manual(values = c("#16048A","#9E189F"))
          + theme(axis.title.y = element_text(angle = 0, vjust = 0.5)))
    ggsave(paste(datasetpath, "ring_clusters.pdf", sep=""),device = cairo_pdf,width =5, height =3.75)
  }
  
  # Generate Half Ring data set
  # 对HDBSCAN，如果两个簇密度一样反而不好区别，还必须是密度不一样才行
  if (T) {
    set.seed(1)
    generatePointsByUniform <- function(x, y, num, label) {
      x <- round(runif(num,x-15,x+15))
      y <- round(runif(num, y-15,y+15))
      data.frame(x, y, label)
    }
    generateRingShapePointsbyUniform <- function(r,xoffset,yoffset,divangle, class) {
      cluster = vector()
      for (i in 1:60) {
        angle = i /divangle
        x = round(r * cos(angle)+xoffset)
        y = round(r * sin(angle)+yoffset)
        cluster <- rbind(cluster, generatePointsByUniform(x+100,y+100,1,class))
      }
      cluster
    }
    cluster1 <- generateRingShapePointsbyUniform(100,0,0,20,'cluster 1')
    cluster2 <- generateRingShapePointsbyUniform(100,100,50,-20,'cluster 2')
    halfring_clusters <- rbind(cluster1, cluster2)
    equalidx = c()
    distdata = as.matrix(dist(halfring_clusters[,c(1:2)]))
    equalist = which(distdata==0, arr.ind = TRUE)
    for(i in c(1:nrow(equalist))){
      if(equalist[i,1] != equalist[i,2]){
        equalidx = c(equalidx, equalist[i,1], equalist[i,2])
      }
    }
    equalidx
    save(halfring_clusters, file=paste(datasetpath, "halfring_clusters.RData", sep=""))
    
    print(ggplot(data = halfring_clusters, mapping = aes(x = x, y = y,fill = as.factor(label))) 
          + geom_point(colour = "white",size = 3,shape = 21) 
          + labs(title = "HALFRING", x ="X", y = "Y") 
          + theme_bw() 
          + theme(legend.position = "bottom",axis.text.x = element_text(vjust = 0.5, hjust = 0.5),axis.text.y = element_text(vjust = 0.5, hjust = 0.5))
          + theme(panel.grid=element_blank(), plot.title = element_text(hjust = 0.5), legend.direction = "horizontal",legend.title = element_blank())
          + scale_fill_manual(values = c("#16048A","#9E189F"))
          + theme(axis.title.y = element_text(angle = 0, vjust = 0.5)))
    ggsave(paste(datasetpath, "halfring_clusters.pdf", sep=""),device = cairo_pdf,width =5, height =3.75)
  }
  
}