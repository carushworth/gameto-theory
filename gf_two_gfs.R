library(tidyverse, lib.loc="/home/crusher/R/x86_64-pc-linux-gnu-library/3.6")

# Meiosis

# Meiosis for the gf complex
meiosisGFcomplex <- function(this.dip, tmp.haplo.names, r12,r23, this.order = "AMF") {
  orig.names <- names(this.dip)
  names(this.dip) <- str_remove(string = names(this.dip), "x|y")
  this.order <- strsplit(this.order,"")[[1]]
  reorder <- order(this.order)[c(1,3,2)]
  # oder in dataframe order: "A.mat" "M.mat" "F.mat" "A.pat" "M.pat" "F.pat" "pop" ...
  # physical order AMF
  # genet
  meiosis.haps <-
    zapsmall(.5 * (1 - r12)* (1 - r23) *(as.numeric(tmp.haplo.names == paste(this.dip[paste(this.order,c("mat","mat","mat"),sep=".")][reorder], collapse = ""))) + # no_rec
               .5 * (1 - r12)* (1 - r23) *(as.numeric(tmp.haplo.names == paste(this.dip[paste(this.order,c("pat","pat","pat"),sep=".")][reorder], collapse = ""))) + # no_rec
               #
               .5 * (r12)    * (1 - r23) *(as.numeric(tmp.haplo.names == paste(this.dip[paste(this.order,c("mat","pat","pat"),sep=".")][reorder], collapse = ""))) + # 12 rec
               .5 * (r12)    * (1 - r23) *(as.numeric(tmp.haplo.names == paste(this.dip[paste(this.order,c("pat","mat","mat"),sep=".")][reorder], collapse = ""))) + # 12 rec
               #
               .5 * (1 - r12)* (    r23) *(as.numeric(tmp.haplo.names == paste(this.dip[paste(this.order,c("mat","mat","pat"),sep=".")][reorder],collapse = ""))) + # 23_rec
               .5 * (1 - r12)* (    r23) *(as.numeric(tmp.haplo.names == paste(this.dip[paste(this.order,c("pat","pat","mat"),sep=".")][reorder], collapse = ""))) + # 23_rec
               #
               .5 * (    r12)* (    r23) *(as.numeric(tmp.haplo.names == paste(this.dip[paste(this.order,c("mat","pat","mat"),sep=".")][reorder],collapse = ""))) + # both_rec
               .5 * (    r12)* (    r23) *(as.numeric(tmp.haplo.names == paste(this.dip[paste(this.order,c("pat","mat","pat"),sep=".")][reorder],collapse = ""))))   # both_rec
  names(this.dip) <- orig.names
  return(meiosis.haps)  # helps with numeric precision
}

# meiosis for unlinked loci
unlinkedMeiosis <- function(tmp.diplo, tmp.haplo.names, haplo.vals){
  tmp.geno <-  abs(1-(tmp.diplo[grep("mat",names(tmp.diplo))] + tmp.diplo[grep("pat",names(tmp.diplo))] )/2)
  apply(abs(tmp.geno - t(haplo.vals)) ,2,prod)
}



# whom mates with whome
matingCombos <- function(my.haps, my.dips, discrim, meiotic.prod){
  how.random.x  <- 1 - ((rowSums(my.dips[,c("Fx.mat","Fx.pat")]) != 0) * discrim[1]) 
  how.random.y  <- 1 - ((rowSums(my.dips[,c("Fy.mat","Fy.pat")]) != 0) * discrim[2]) 
  # first take F alleles for mat and pat, if there's a F present in either the sum will be >0
  # for those that have a sum >0, multiply by preference intensity ("discrim") and subtract that amt from 1
  # so this says if you have a 1 at either F locus you'll be discriminated against
  mating.freqs <- matrix(1,ncol = length(my.haps), nrow = length(how.random.x))
  # these are the frequencies of matings if you're a certain mom/dad combo
  # but right now it's just a matrix of all 1s
  colnames(mating.freqs) <- names(my.haps) # dad names as column names
  rownames(mating.freqs) <- rownames(meiotic.prod) # mom names as row names, taken from meiotic.prod
  mating.freqs           <- t(t(mating.freqs) * my.haps)
  # transposing 2x makes values in columns the rows, wait this just gives us back 8 cols and 64 rows
  mating.freqs[,grep(pattern = ".0.", substring(names(my.haps),first = 1, last=3))]  <- mating.freqs[,grep(pattern = ".0.", substring(names(my.haps),first = 1, last=3))]  * how.random.x
  mating.freqs[,grep(pattern = ".0.", substring(names(my.haps),first = 4, last=6))]  <- mating.freqs[,grep(pattern = ".0.", substring(names(my.haps),first = 4, last=6))]  * how.random.y
  # these are grabbing all the little m's cause we're AMF and multiplying by how permissive each mom is/how random the mating is
  denom <- rowSums(mating.freqs)
  mating.freqs / ifelse( denom  == 0,1,denom)# standardizing so that discriminating moms aren't producing fewer kids
}


diploidSel <- function(tmp.freq, my.dips, s, focal.pop){
  s_linked <- s_unlinked <- s[paste("s",focal.pop,sep="")] # Currently selection on linked and unlinked loci is eaual... this can be changed
  #n.maladapt_unlinked <- rowSums(my.dips[,grep("U",colnames(my.dips))] != focal.pop)
  n.maladapt_linked <- rowSums(my.dips[,grep("A",colnames(my.dips))] != focal.pop)
  w          <- zapsmall( (1-s_linked)^n.maladapt_linked )
  wbar       <- sum((zapsmall(tmp.freq) *  zapsmall(w))) # mean fitness = the sum of t
  wrel       <- zapsmall(w / wbar)  #relative fitness
  newfreq    <- zapsmall(wrel  * tmp.freq) 
  return(newfreq/ (sum(newfreq)) ) # again for precision
}

migrateMateReproduceSelect <- function(diplos, focal.pop, prop.replaced, meiotic.prod, discrim, s, delta_hap_components, hap.ids){
  dhap <- rep(NA, 4*length(hap.ids))
  home.haps  <- zapsmall(colSums(meiotic.prod * diplos$freqs[diplos$pop == focal.pop] ) * (1-prop.replaced))
  away.haps  <- zapsmall(colSums(meiotic.prod * diplos$freqs[diplos$pop != focal.pop] ) * (prop.replaced)) # prop.replaced is migration rate
  my.haps    <- home.haps + away.haps# pollen pool
  #if(delta_hap_components){
  #  init.hap <- colSums(meiotic.prod * diplos$freqs[diplos$pop == focal.pop])
    # change in hap freq by migrations
    # change in pollen       no change in ovules   average
  #  dp.migration.pollen <- my.haps - init.hap
  #}
  #recover()
  my.dips    <- diplos[diplos$pop == focal.pop,]                                                                     # females
  pat.haps   <- matingCombos(my.haps = my.haps , my.dips =  my.dips, discrim =discrim, meiotic.prod = meiotic.prod)  # mating rules
  mat.haps   <- meiotic.prod * my.dips$freqs / sum( my.dips$freqs)                                                                         # meiosis in feamles
  new.genos  <- t(pat.haps) %*% (mat.haps) # rows are paternal, colums are maternal
  p.migrant  <- away.haps  /  my.haps # proportion of each pollen haplotype that is a migrant.. this is useful for downstream accounting.
  pop.reinforce <- 1 - sum(rowSums(new.genos) * p.migrant ,na.rm=TRUE) /prop.replaced   # quantifying reinforcement as 1 - prob mating with migrant
  # first thing in pop0 and pop1
  before.sel <- c(t(new.genos)) # flip rows and columns and then flatten to make  new diploid freqsQ
  after.sel  <- diploidSel(tmp.freq = before.sel, my.dips = my.dips, s = s, focal.pop = focal.pop)
  print(c(focal.pop, sum(pat.haps)/4096,sum(mat.haps),sum(new.genos) , sum(before.sel), sum(after.sel)))
  return(c(pop.reinforce, after.sel))
}


run2GFsim <-function(n.gen = Inf, r12 = 1e-4, r23 =  1e-4,
                    init.freqs = c(fAx_0 = 0, fMx_0 = 0, fFx_0 = 0,  
                                   fAy_0 = 0, fMy_0 = 1, fFy_0 = .01, 
                                   fAx_1 = 1, fMx_1 = 1, fFx_1 = .01,
                                   fAy_1 = 1, fMy_1 = 0, fFy_1 = 0),
                    discrim = c(1,1), s = NULL, tol = 1e-5, max.gen = 2e6, min.gen = 1e3,
                    prop.replaced0 = .1,
                    prop.replaced1 = .1,
                    this.order = "AMF",
                    n.unlinked =0,
                    get.blank = FALSE,
                    delta_hap_components = FALSE,
                    return.reinforce.only = TRUE,
                    s0 = NULL,
                    s1 = NULL
){
  print(init.freqs)
  if(is.null(s) & is.null(s0) & is.null(s1) ){
    warning("no s specified")
    return(":(")
  }
  if(is.null(s) & !(is.null(s0) & is.null(s1))){
    s <- c(s0 = s0, s1 = s1)
  }
  if(is.na(tol)){  max.gen  = n.gen; tol <- Inf}
  if(length(s)==2){my.s = -99}
  if(length(s)==1){
    my.s <- s
    s    <- c(s[1],s[1])
    names(s) <- c("s0","s1")
  }
  # SETUP (this makes all the different diploid genos)
  print(sprintf("s = %s, m = %s, r12 = %s, r23 = %s, init_freq = %s, s0 = %s, s1 = %s", my.s, prop.replaced0, r12, r23, init.freqs["fF_1"],s["s0"], s["s1"]))
  diplos <- expand.grid(data.frame(rbind(numeric(length = 2*(6+n.unlinked)+1),1)))
  tmp.names <- c("Ax","Mx","Fx","Ay","My","Fy",paste("U", 0:n.unlinked,sep ="")[-1])
  colnames(diplos) <- c(paste(rep(tmp.names,times = 2), rep(c("mat","pat"), each = length(tmp.names)), sep="."),"pop")
  print(init.freqs)
  diplos$freqs <-  apply(diplos, 1, function(X){
    this.pop        <- X[["pop"]]
    freqs           <- grep( this.pop, names(init.freqs))
    tmp.big         <- c(init.freqs[grep( this.pop, names(init.freqs))],rep(this.pop,n.unlinked))
    print(length(tmp.big))
    print(length(c(names(tmp.big)[1:6],paste(paste(rep("U",n.unlinked), (0:n.unlinked) , sep=""),this.pop,sep="_"  )[-1])))
    print(init.freqs)
    print(init.freqs[grep( this.pop, names(init.freqs))])
    print(rep(this.pop,n.unlinked))
    print(grep(this.pop, names(init.freqs)))
    print(this.pop)
    names(tmp.big)  <- c(names(tmp.big)[1:6],paste(paste(rep("U",n.unlinked), (0:n.unlinked) , sep=""),this.pop,sep="_"  )[-1]) # assumes we begin as locally adapted at all loci
    f.big           <- rep(tmp.big, 2)
    prod(as.numeric((1 - X[-length(X)]) == 0) * f.big  + as.numeric((X[-length(X)]) == 0) * (1-f.big))
  })
  haplos <- expand.grid(data.frame(rbind(numeric(length = (6+n.unlinked)+1),1))) # makes all the different haplotypes
  names(haplos) <- c("Ax","Mx","Fx","Ay","My","Fy",paste("U", 0:n.unlinked,sep ="")[-1],"pop")
  haplo.names  <-  apply(unique(haplos[,-ncol(haplos)]), 1, paste, collapse = "")
  #
  ### MEIOSIS
  # First worry about our linked complex ones
  gf.complex <- unique(diplos[,grep("Ax|Mx|Fx", colnames(diplos))])
  tmp.gf     <- apply(diplos[diplos$pop == 0,grep("Ax|Mx|Fx", colnames(diplos))],1, paste,collapse="")
  three.loc.meiosis <- t(apply(gf.complex,1, meiosisGFcomplex, tmp.haplo.names = unique(str_sub(haplo.names,1,3)),
                               r12 = r12, r23 = r23,  this.order = this.order) )
  three.loc.meiosis  <- (three.loc.meiosis / rowSums(three.loc.meiosis))  # fix for numerical precision
  rownames(three.loc.meiosis) <- apply(gf.complex,1, paste,collapse="")
  colnames(three.loc.meiosis) <- unique(str_sub(haplo.names,1,3))
#  three.loc.meiosis        <- three.loc.meiosis[ tmp.gf,str_sub(haplo.names,1,3)]
#  colnames(three.loc.meiosis) <- haplo.names
   #
  two.x.three.loc.meiosis <- do.call(rbind,lapply(1:nrow(three.loc.meiosis), function(ROW){
    #print(ROW)
    xxx <- do.call(cbind,lapply(1:ncol(three.loc.meiosis), function(COL){
      #print(c(ROW,COL))
      tmp.meios <- three.loc.meiosis * three.loc.meiosis[ROW,COL]
      colnames( tmp.meios ) <- paste(colnames(tmp.meios),
                                   colnames(tmp.meios)[COL],
                                   sep = "")
      tmp.meios
     }))
    rownames(xxx) <- paste(str_sub(rownames(three.loc.meiosis),1,3),   # so that we list them as mat / pat.. rather than  matx patx mat paty
                           str_sub(rownames(three.loc.meiosis)[ROW],1,3),
                           str_sub(rownames(three.loc.meiosis),4,6),
                           str_sub(rownames(three.loc.meiosis)[ROW],4,6),
                           sep="")
    return(xxx)
  }))
  rm(three.loc.meiosis)


  if(n.unlinked == 0){
    meiotic.prod <- two.x.three.loc.meiosis
    rm(two.x.three.loc.meiosis)
  }
  # Now meioisis at unlinked loci
  if(n.unlinked > 0){
    unlinked.complex <- unique(diplos[,grep("U", colnames(diplos))])
    tmp.complex      <- apply(diplos[diplos$pop == 0,grep("U", colnames(diplos))],1, paste,collapse="")
    unlinked.haps    <- unique(str_sub(haplo.names,start = 4))
    unlinked.meiosis <- t(apply(X = unlinked.complex,
                                MARGIN = 1,
                                FUN = unlinkedMeiosis,
                                tmp.haplo.names = unlinked.haps,
                                haplo.vals = apply(do.call(rbind,strsplit(unlinked.haps,"")), 2,as.numeric)
    ))
    rownames(unlinked.meiosis) <- apply(unlinked.complex,1, paste,collapse="")
    colnames(unlinked.meiosis) <- unlinked.haps
    unlinked.meiosis           <- unlinked.meiosis[tmp.complex ,str_sub(haplo.names,4)]
    meiotic.prod               <- two.x.three.loc.meiosis * unlinked.meiosis # jam them together
    rownames( meiotic.prod)    <- apply(diplos[diplos$pop == 0,1:(ncol(diplos)-2)], 1, paste, collapse = "")
    rm(unlinked.meiosis, two.x.three.loc.meiosis,tmp.complex,unlinked.haps)
  }
  print(head(diplos))
  aaaaa <- diplos %>% dplyr::filter(pop == 0) 
  zzz   <- apply(aaaaa[,1:12], 1, paste, collapse ="") %>% as.character()
  meiotic.prod <- meiotic.prod[zzz,]
  ### assigning names to everything
  names.haplo.time       <- c(apply(unique(diplos[,grep("U|q|pat",names(diplos),invert = T)]), 1, paste, collapse = ""))
  haplo.time             <- matrix(ncol = 3 + length(names.haplo.time), nrow = max.gen)
  colnames(haplo.time)   <- c("gen","reinf_0","reinf_1", names.haplo.time)
  #meanUs <- geno.time
  #colnames(meanUs)[c(2:3)] <- c("U_0","U_1")
  ##
  hap.ids <- unique(str_sub(haplo.names,1,3))
  dhap_components <- matrix(ncol = 1 + 2 * length(hap.ids) * 4, nrow = max.gen)
  colnames(dhap_components) <- c("gen",paste(rep(c("p0","p1"), each=32),rep(paste(hap.ids,rep(c("migrationPollen", "matingPollen", "selPat", "selMat"),each = 8),sep="_"),2),sep="_"))
  # this loop can give a blank output
  if(get.blank){
    geno.time <- data.frame(geno.time)
    meanUs    <- data.frame(meanUs)
    geno.time[,"gen"] <- 1: max.gen
    meanUs[,"gen"] <- 1: max.gen
    return(list(geno.time = geno.time, meanUs = meanUs))
  }
  time.to.stop = FALSE
  g <- 0
  while(! time.to.stop){
    g <- g + 1
    pop0 <- migrateMateReproduceSelect(diplos = diplos,
                                       focal.pop = 0,
                                       prop.replaced = prop.replaced0,
                                       meiotic.prod = meiotic.prod ,
                                       discrim = discrim ,
                                       s = s,
                                       delta_hap_components = FALSE,
                                       hap.ids = hap.ids)
    pop1 <- migrateMateReproduceSelect(diplos = diplos,
                                       focal.pop = 1,
                                       prop.replaced = prop.replaced1,
                                       meiotic.prod = meiotic.prod ,
                                       discrim = discrim,
                                       s = s,
                                       delta_hap_components = FALSE,
                                       hap.ids = hap.ids)
    
#    print(c(g, "pop0=",sum(pop0[-1]),  "pop1=" ,sum(pop1[-1])))
#    dhap_components[g,]<- c(g,pop0[(length(pop0)-31):length(pop0)],pop1[(length(pop1)-31):length(pop1)])
    #pop0 <- pop0[-((length(pop0)-31):length(pop0))]
    #pop1 <- pop1[-((length(pop1)-31):length(pop1))]
    #if(g%%100 == 0){print(sprintf("generation %s of %s",g, n.gen))}
    new.freqs     <- c(pop0[-1], pop1[-1])
    diff.freqs    <- sum(abs(new.freqs  - diplos$freqs ))
    diplos$freqs  <- new.freqs   # put together pop0 and pop1
    #
    mat_haps      <- tibble(diplos) %>% 
      group_by(Ax.mat, Mx.mat, Fx.mat , Ay.mat, My.mat, Fy.mat, pop) %>% 
      summarise(freqs  = sum(freqs),.groups = "drop") %>%
      mutate(hap = paste(Ax.mat , Mx.mat, Fx.mat , Ay.mat ,My.mat, Fy.mat,pop, sep=""),
             parent = "mat") %>%
      dplyr::select(hap,parent, freqs)
    #
    pat_haps     <-tibble(diplos) %>% 
      group_by(Ax.pat, Mx.pat, Fx.pat , Ay.pat, My.pat, Fy.pat, pop) %>% 
      summarise(freqs  = sum(freqs),.groups = "drop") %>%
      mutate(hap = paste(Ax.pat , Mx.pat, Fx.pat , Ay.pat ,My.pat, Fy.pat,pop, sep=""),
             parent = "pat")%>%
      dplyr::select(hap,parent, freqs)
    #
    hap_freqs <- bind_rows(mat_haps, pat_haps) %>%#mean haps across parent of origin
      group_by(hap) %>%
      summarise(freq = mean(freqs),.groups = "drop") %>%
      mutate(hap = factor(hap, levels = colnames(haplo.time)[-c(1:3)])) %>% # matching order
      arrange() %>%
      pull(freq)
    haplo.time[g,]<-c(g,pop0[1],pop1[1],hap_freqs)%>% as.numeric() 
  #  recover()

    #if(n.unlinked == 0){ geno.time[g,]  <- as.numeric(c(g, pop0[1], pop1[1],  diplos$freqs))  }
    # for each gen of geno.time we have gen, reinforcement in each pop, and the diplo frequencies
    #if(n.unlinked > 0){
    #  tmp.diplos <- diplos %>%
    #    group_by(A.mat, M.mat, F.mat, A.pat ,M.pat, F.pat, pop) %>%
    #    summarise(freqs = sum(freqs)) %>%
    #    ungroup() %>%
    #    mutate(geno = paste(A.mat, M.mat, F.mat, A.pat, M.pat, F.pat,   pop, sep = ""))
      # FIX THIS
    #  geno.time[g,]  <- c(g, pop0[1], pop1[1],
    #                      unlist(tmp.diplos %>% select(freqs, geno ) %>%
    #                               spread(key = geno, value = freqs))[names.geno.time]   )
    #  tmp.diplos2 <- diplos %>%
    #    mutate(mean_U = diplos %>% select(starts_with("U")) %>% rowMeans() )      %>%
    #    group_by(A.mat, M.mat, F.mat, A.pat ,M.pat, F.pat, pop) %>%
    #    summarise(U = sum(freqs * mean_U), freqs = sum(freqs))       %>%
    #    ungroup()

    #  meanUs[g,] <- c(g,
    #                  tmp.diplos2 %>% group_by(pop) %>% summarise(x = sum(U)) %>%pull(),
    #                  unlist(tmp.diplos2 %>%
    #                           mutate(geno = paste(A.mat, M.mat, F.mat, A.pat, M.pat, F.pat,   pop, sep = ""),
    #                                  meanGenoU    = U / freqs) %>%
    #                           select(meanGenoU, geno ) %>%
    #                           spread(key = geno, value = meanGenoU))[names.geno.time])
    #}
    time.to.stop <- ((g+1) > min.gen) & (diff.freqs < tol | (g+1) > n.gen | (g+1) > max.gen)
    #print(paste("gen = ",g, "   & diff.freqs = ",diff.freqs))
    print(g)
  }
  haplo.time <- data.frame(haplo.time[1:g,])  %>%
    mutate(teo_reinforce = (reinf_1 - reinf_1[1])  / (1- reinf_1[1])  )
  #meanUs    <- meanUs[1:g,]
#  dhap_components <-  dhap_components[1:g,]
  max_reinforce_teo <- max(haplo.time[,"teo_reinforce"],na.rm=TRUE)
  #max_freq_M0   <- data.frame(geno.time) %>% select(matches("X.1")) %>% select(ends_with("0")) %>% rowSums() %>%max()
  #max_freq_F1   <- data.frame(geno.time) %>% select(matches("X..1")) %>% select(ends_with("1")) %>% rowSums() %>%max()
  #final_adapt_diff_unlinked <- diff(rev(meanUs[g,c("U_0","U_1")])) # how much adaptive divergence do we have at the end (we start with complete)
  #start_reinforce <- min(which(geno.time[,"teo_reinforce"]>0))
  #end_reinforce <- max(which(geno.time[,"teo_reinforce"]>0))
  #names(final_adapt_diff_unlinked) <- NULL
  if(return.reinforce.only)  {
    to.keep  <- which(haplo.time$gen %% floor(nrow(haplo.time) / 4000) ==1 |  c(1, rowSums(abs(haplo.time[-1,-c(1:3)] - haplo.time[-nrow(haplo.time),-c(1:3)]  ))) > 5e-4)
    geno.time <- slice(haplo.time, to.keep)  %>% dplyr::select(gen, reinf_1)
    return(list(geno.time =  geno.time,
                params = c(this.order = this.order, n.gen = n.gen, r12 =r12, r23 =r23, s=my.s, n.unlinked = n.unlinked, prop.replaced0 = prop.replaced0 , prop.replaced1 = prop.replaced1, discrim = discrim, s0 = s["s0"],s1 = s["s1"])
    ))
  }
  return(list(haplo.time =  haplo.time,
              #meanUs = data.frame(meanUs),
#              dhaps  = data.frame(dhap_components),
              #summary.stats = c(max_reinforce_teo = max_reinforce_teo, min_diff_A = min_diff_A, max_freq_M0 = max_freq_M0, max_freq_F1 = max_freq_F1, final_adapt_diff_unlinked = final_adapt_diff_unlinked, start_reinforce,end_reinforce,tot.gen = g),
              params = c(this.order = this.order, n.gen = n.gen, r12 =r12, r23 =r23, s=my.s, n.unlinked = n.unlinked, prop.replaced0 = prop.replaced0 , prop.replaced1 = prop.replaced1, discrim = discrim, s0 = s["s0"],s1 = s["s1"])
  ))
}
