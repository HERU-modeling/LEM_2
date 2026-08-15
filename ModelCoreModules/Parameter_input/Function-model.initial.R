## function used to generate the initial population for all groups at different states 
## output: a 42 * 19 matrix, initial population for 42 groups at 19 different states

model.initial = function(par, diag18, gp.gn){
  with(as.list(c(par, diag18 = diag18, gp.gn)), {
    # diag18: # of hiv diagnosed in 18 groups, 
  	# created by "targetNYC.R", from NYC surveillance reports

    # population in 18 groups, stratified by gender, ethnicity, and risk ("pop")
    pop.total  = pop.male + pop.female
    pop.msm    = pop.male  *prop.msm.among.male
    pop.pwid.f = pop.total *prop.pwid.among.total *(1-prop.male.among.pwid)
    pop.pwid.m = pop.total *prop.pwid.among.total *prop.male.among.pwid
  
    # # midu -> mwid
    # pop.mwidP  = pop.pwid.m * prop.mpwid.among.pwid # MSM/PWID from M/PWID
    # pop.mwidM  = pop.msm    * prop.mpwid.among.msm # MSM/PWID from MSM
    # pop.mwid   = (pop.mwidP + pop.mwidM) / 2  # MSM/PWID average of 2 above
    # midu -> mwid, Azzam removed the averaging method, only using mpwid from male pwid
    pop.mwidP = pop.pwid.m * prop.mpwid.among.pwid
    pop.mwid = pop.mwidP
    
    
    pop.pwid.m = pop.pwid.m - pop.mwid #PWID
    pop.msm    = pop.msm    - pop.mwid
    
    pop.het.m  = pop.male   - pop.msm - pop.pwid.m - pop.mwid
    pop.het.f  = pop.female - pop.pwid.f
    pop        = round(c(pop.msm, pop.mwid, pop.pwid.m, pop.pwid.f, pop.het.m, pop.het.f))
#browser()  
#    if (any(prop.curr.art <= viral_sup)) {warning(paste0("prop.curr.art <= viral_sup detected; values adjusted in ", CITY))}
#    if (any(prop.ever.art <= prop.curr.art)) {warning(paste0("prop.ever.art <= prop.curr.art detected", CITY))}
    
    # Ensure prop.curr.art>viral_sup
    for (i in seq_len(length(prop.curr.art))) {
    		if (prop.curr.art[i] - viral_sup[i] <= 0) {
#    			warning(paste0("prop.curr.art = ", prop.curr.art[i],", viral_sup = ", viral_sup[i]," in Juris ", CITY,", state ", names18[i]))
#    			browser()
    			prop.curr.art[i] = viral_sup[i] + 0.01
    	}
    }  
    # Ensure prop.ever.art>prop.curr.art
    for (i in seq_len(length(prop.ever.art))) {
    	if (prop.ever.art[i] - prop.curr.art[i] <= 0) {
#    		warning(paste0("prop.curr.art = ", prop.curr.art[i],", viral_sup = ", viral_sup[i]," in Juris ", CITY,", state ", names18[i]))
#    		browser()
    		prop.ever.art[i] = prop.curr.art[i] + 0.01
    	}
    } 
    
    # infected but not diag 
    # D/(I+D)=prop.HIV.aware => I=(D-prop*D)/prop
    I18 = (diag18 - prop.HIV.aware *diag18) /prop.HIV.aware # increase prop.aware will decrease no. of I18
#browser()    
############
#Re-adjust I18 based on ndiag target, Should not impact Prevalence since we don't touch diag 18     
    #Level 1: msm, het,pwid prop among ndiag
    ndiag.msm.prop <- sum(calib.target$ndiag.obs[1,c("Young","Old")])/sum(calib.target$ndiag.obs[1,c("W","B","H")])
    ndiag.het.prop <- sum(calib.target$ndiag.obs[1,c("M.HET","F.HET")])/sum(calib.target$ndiag.obs[1,c("W","B","H")])
    ndiag.pwid.prop <- sum(calib.target$ndiag.obs[1,c("PWID")])/sum(calib.target$ndiag.obs[1,c("W","B","H")])
    
#Level 2: Y, O among msm; H,L among HEt
    ndiag.Y.msm.prop <- calib.target$ndiag.obs[1,c("Young")] /sum(calib.target$ndiag.obs[1,c("Young","Old")]) #Not use for I18
    ndiag.O.msm.prop <- calib.target$ndiag.obs[1,c("Old")]   /sum(calib.target$ndiag.obs[1,c("Young","Old")]) #Not use for I18
    ndiag.M.het.prop <- calib.target$ndiag.obs[1,c("M.HET")] /sum(calib.target$ndiag.obs[1,c("M.HET","F.HET")])
    ndiag.F.het.prop <- calib.target$ndiag.obs[1,c("F.HET")] /sum(calib.target$ndiag.obs[1,c("M.HET","F.HET")])
    
#Level 3: W,B,H among all
    ndiag.W.prop <- sum(calib.target$ndiag.obs[1,c("W")])/sum(calib.target$ndiag.obs[1,c("W","B","H")])
    ndiag.B.prop <- sum(calib.target$ndiag.obs[1,c("B")])/sum(calib.target$ndiag.obs[1,c("W","B","H")])
    ndiag.H.prop <- sum(calib.target$ndiag.obs[1,c("H")])/sum(calib.target$ndiag.obs[1,c("W","B","H")])
    
#COmbine, prop for I18
    ndiag.W.Y.msm.prop <- ndiag.W.prop * ndiag.Y.msm.prop * ndiag.msm.prop
    ndiag.B.Y.msm.prop <- ndiag.B.prop * ndiag.Y.msm.prop * ndiag.msm.prop
    ndiag.H.Y.msm.prop <- ndiag.H.prop * ndiag.Y.msm.prop * ndiag.msm.prop #Sum these 3 = ndiag.Y.msm.prop * ndiag.msm.prop
    ndiag.W.O.msm.prop <- ndiag.W.prop * ndiag.O.msm.prop * ndiag.msm.prop
    ndiag.B.O.msm.prop <- ndiag.B.prop * ndiag.O.msm.prop * ndiag.msm.prop
    ndiag.H.O.msm.prop <- ndiag.H.prop * ndiag.O.msm.prop * ndiag.msm.prop #Sum these 3 = ndiag.O.msm.prop * ndiag.msm.prop
    
    ndiag.W.M.het.prop <- ndiag.W.prop * ndiag.M.het.prop * ndiag.het.prop
    ndiag.B.M.het.prop <- ndiag.B.prop * ndiag.M.het.prop * ndiag.het.prop
    ndiag.H.M.het.prop <- ndiag.H.prop * ndiag.M.het.prop * ndiag.het.prop #Sum these 3 = ndiag.M.het.prop * ndiag.het.prop
    ndiag.W.F.het.prop <- ndiag.W.prop * ndiag.F.het.prop * ndiag.het.prop
    ndiag.B.F.het.prop <- ndiag.B.prop * ndiag.F.het.prop * ndiag.het.prop
    ndiag.H.F.het.prop <- ndiag.H.prop * ndiag.F.het.prop * ndiag.het.prop #Sum these 3 = ndiag.M.het.prop * ndiag.het.prop
    
    ndiag.W.pwid.prop  <- ndiag.W.prop * ndiag.pwid.prop
    ndiag.B.pwid.prop  <- ndiag.B.prop * ndiag.pwid.prop
    ndiag.H.pwid.prop  <- ndiag.H.prop * ndiag.pwid.prop
# 1= ndiag.W.Y.msm.prop+ndiag.B.Y.msm.prop+ndiag.H.Y.msm.prop+ndiag.W.O.msm.prop+ndiag.B.O.msm.prop+ndiag.H.O.msm.prop+ndiag.W.M.het.prop+ndiag.B.M.het.prop+ndiag.H.M.het.prop+ndiag.W.F.het.prop+ndiag.B.F.het.prop+ndiag.H.F.het.prop+ndiag.W.pwid.prop+ndiag.B.pwid.prop+ndiag.H.pwid.prop    
    sum.I18 <- sum(I18)
    I18[names18 == "MSM/white"]      <- sum.I18 * ndiag.W.prop*ndiag.msm.prop #W.msm: ndiag.W.Y.msm.prop+ndiag.W.O.msm.prop
    I18[names18 == "MSM/black"]      <- sum.I18 * ndiag.B.prop*ndiag.msm.prop #B.msm: ndiag.B.Y.msm.prop+ndiag.B.O.msm.prop
    I18[names18 == "MSM/hispanic"]   <- sum.I18 * ndiag.H.prop*ndiag.msm.prop #H.msm: ndiag.H.Y.msm.prop+ndiag.H.O.msm.prop
# ndiag.W.prop*ndiag.msm.prop+ndiag.B.prop*ndiag.msm.prop+ndiag.H.prop*ndiag.msm.prop =  ndiag.W.Y.msm.prop+ndiag.B.Y.msm.prop+ndiag.H.Y.msm.prop+ndiag.W.O.msm.prop+ndiag.B.O.msm.prop+ndiag.H.O.msm.prop
    pwid.sum.W <- sum(I18[names18 %in% c("MSM/PWID/white", "PWID/male/white", "PWID/female/white")])
    I18[names18 == "MSM/PWID/white"] <- sum.I18 * ndiag.W.pwid.prop          *(I18[names18 == "MSM/PWID/white"]/pwid.sum.W)
    I18[names18 =="PWID/male/white"] <- sum.I18 * ndiag.W.pwid.prop          *(I18[names18 =="PWID/male/white"]/pwid.sum.W)
    I18[names18 =="PWID/female/white"]<- sum.I18 * ndiag.W.pwid.prop         *(I18[names18 =="PWID/female/white"]/pwid.sum.W)
    
    pwid.sum.B <- sum(I18[names18 %in% c("MSM/PWID/black", "PWID/male/black", "PWID/female/black")])
    I18[names18 == "MSM/PWID/black"] <- sum.I18 * ndiag.B.pwid.prop          *I18[names18 == "MSM/PWID/black"]/pwid.sum.B
    I18[names18 =="PWID/male/black"] <- sum.I18 * ndiag.B.pwid.prop          *I18[names18 =="PWID/male/black"]/pwid.sum.B
    I18[names18 =="PWID/female/black"]<- sum.I18 * ndiag.B.pwid.prop         *I18[names18 =="PWID/female/black"]/pwid.sum.B
    pwid.sum.H <- sum(I18[names18 %in% c("MSM/PWID/hispanic", "PWID/male/hispanic", "PWID/female/hispanic")])
    I18[names18 == "MSM/PWID/hispanic"] <- sum.I18 * ndiag.H.pwid.prop          *I18[names18 == "MSM/PWID/hispanic"]/pwid.sum.H
    I18[names18 =="PWID/male/hispanic"] <- sum.I18 * ndiag.H.pwid.prop          *I18[names18 =="PWID/male/hispanic"]/pwid.sum.H
    I18[names18 =="PWID/female/hispanic"]<- sum.I18 * ndiag.H.pwid.prop         *I18[names18 =="PWID/female/hispanic"]/pwid.sum.H
    I18[names18 == "HET/male/white"]     <- sum.I18 * ndiag.W.M.het.prop
    I18[names18 == "HET/male/black"]     <- sum.I18 * ndiag.B.M.het.prop
    I18[names18 == "HET/male/hispanic"]  <- sum.I18 * ndiag.H.M.het.prop
    I18[names18 == "HET/female/white"]   <- sum.I18 * ndiag.W.F.het.prop
    I18[names18 == "HET/female/black"]   <- sum.I18 * ndiag.B.F.het.prop
    I18[names18 == "HET/female/hispanic"]<- sum.I18 * ndiag.H.F.het.prop
    
    # cause I18 do not include young,old -dist them later; prop.young.I ES alos use ndiag from atlasplus
    prop.young.I[msm.gn$Y.w]   <- ndiag.Y.msm.prop
    prop.young.I[msm.gn$Y.b]   <- ndiag.Y.msm.prop
    prop.young.I[msm.gn$Y.h]   <- ndiag.Y.msm.prop
#browser()    
 ############    
    
    # susceptible
    S18 = pop - diag18 - I18                                # diag18 = Da+D1+D2+D3+T1+T2+T3+O1+O2+O3
 
    # ever on ART
    T.ever18 = diag18 *prop.ever.art # prop.ever.art = (T1+T2+T3+O1+O2+O3)/diag18
    # on ART among ever on ART
    T18 = T.ever18 * prop.curr.art/prop.ever.art  # (T1+T2+T3) /(T1+T2+T3+O1+O2+O3) = prop.curr.art/prop.ever.art & T.ever18 = (T1+T2+T3+O1+O2+O3)
    T18_3  = (prop.curr.art - viral_sup) * diag18 # prop.curr.art = (T1+T2+T3) /diag18 & viral_sup = (T1+T2)/diag18 
    T18_12 = T18 - T18_3
    # drop off ART
    O18 = T.ever18 - T18
    # diagnosed but not on ART
    D18 = diag18 - T18 - O18
 # browser()   
 #  init18            = cbind(S18, I18, D18, T18, O18)
    init18            = cbind(S18, I18, D18, T18_12,T18_3, O18)
    row.names(init18) = names18
   
# Check diag.allmsm: sum(init18[gp18.gn$all.msm ,c("D18","T18_12","T18_3","O18")]) 
# Check diag.het.m:  sum(init18[gp18.gn$het.m ,c("D18","T18_12","T18_3","O18")])
# Check diag.het.f:  sum(init18[gp18.gn$het.f ,c("D18","T18_12","T18_3","O18")])
# Check diag.pwid:   sum(init18[gp18.gn$pwid  ,c("D18","T18_12","T18_3","O18")]) 
# Check diag.W: 		 sum(init18[gp18.gn$white ,c("D18","T18_12","T18_3","O18")])
# Check diag.B: 		 sum(init18[gp18.gn$black ,c("D18","T18_12","T18_3","O18")])
# Check diag.H: 		 sum(init18[gp18.gn$hisp  ,c("D18","T18_12","T18_3","O18")]) 

    # Temporary Fix - Check for negative value, turn it positive
    for (i in seq_len(nrow(init18))) {
    	for (j in seq_len(ncol(init18))) {
    		if (init18[i, j] < 0) {
 #   			warning(paste0("Model.initial: value in Juris ", CITY," strata ", rownames(init18)[i], ", state ", colnames(init18)[j], " was negative: pop =  ", pop[i], ", diag18= ", diag18[i],  ", prop.HIV.aware= ", prop.HIV.aware[i]))
 #   			browser()
    			init18[i, j] <- -init18[i, j]
    		}
    	}
    }
#browser()    
#    if (any(init18< 0)) {warning(paste0("Negative in init18 ", CITY))}   

    n.gp  = length(names.gp)
    # change n.gp group names without OAT, low, high
    rname = gsub(paste(c("/low", "/high", "/old", "/young"), collapse="|"), "", names.gp)

    #init has n.gp in row, populated with same n if same rname
#   init  = matrix(0, n.gp, 5) 
    init  = matrix(0, n.gp, 6) 
    row.names(init) = names.gp
#   colnames(init)  = c("S", "I", "D", "T", "O")
    colnames(init)  = c("S", "I", "D", "T12", "T3", "O")
    prop.high.sus.all =numeric(n.gp)
    prop.high.inf.all =numeric(n.gp)
    prop.young.sus.all =numeric(n.gp)
    prop.young.inf.all =numeric(n.gp)
    for (i in 1:18){
      ind = which(rname %in% names18[i])
      prop.high.sus.all[ind] = prop.high.sus[i]
      prop.high.inf.all[ind] = prop.high.inf[i]
      prop.young.sus.all[ind] = prop.young.sus[i]
      prop.young.inf.all[ind] = prop.young.inf[i]
      for (j in 1:length(ind)){
        init[ind[j], ] = init18[i, ]
      }
    }
  


    # low/high for MSM, MSM/PWID and HET
     init[high,"S"]                              = init[high, "S"]  * as.numeric(prop.high.sus.all[high])
     init[young,"S"]                             = init[young, "S"] * as.numeric(prop.young.sus.all[young])

     init[low,"S"]                    = init[low, "S"]  *(1 - as.numeric(prop.high.sus.all[high]))
     init[old,"S"]                    = init[old, "S"]  *(1 - as.numeric(prop.young.sus.all[young]))

#     init[high , c("I","D","T","O")]  = init[high , c("I","D","T","O")] *as.numeric(prop.high.inf.all[high])
#     init[young, c("I","D","T","O")]  = init[young, c("I","D","T","O")] *as.numeric(prop.young.inf.all[young])
     init[high , c("I", "D", "T12", "T3", "O")]  = init[high , c("I", "D", "T12", "T3", "O")] *as.numeric(prop.high.inf.all[high])
     init[young, c("I", "D", "T12", "T3", "O")]  = init[young, c("I", "D", "T12", "T3", "O")] *as.numeric(prop.young.inf.all[young])
     
#     init[low , c("I","D","T","O")]   = init[low,  c("I","D","T","O")] *(1-as.numeric(prop.high.inf.all[high]))
#     init[old , c("I","D","T","O")]   = init[old,  c("I","D","T","O")] *(1-as.numeric(prop.young.inf.all[young]))
     init[low , c("I", "D", "T12", "T3", "O")]   = init[low,  c("I", "D", "T12", "T3", "O")] *(1-as.numeric(prop.high.inf.all[high]))
     init[old , c("I", "D", "T12", "T3", "O")]   = init[old,  c("I", "D", "T12", "T3", "O")] *(1-as.numeric(prop.young.inf.all[young]))
#browser() #check init
    # Re-distribute number of diagnosed PLHIV young/MSM to Infected young/MSM - Accounted for Lower awareness of HIV among young/MSM

      sum.all.msm.I             = (init[young, c("I")] + init[old, c("I")])
      init[young, c("I")]       = sum.all.msm.I  * as.numeric(prop.young.I[msm.gn$young])
      init[old  , c("I")]       = sum.all.msm.I  * (1-as.numeric(prop.young.I[msm.gn$young]))

    S2 = prop.S2 *init[ ,"S"]
    S1 = init[ ,"S"] - S2
#   I  = init[ ,"I"] *cbind(prop.Ia, prop.I1, prop.I2, prop.I3)    
# browser()
    # Re-do I dist CD4: account for I3.prop = p.ldiag by races, rest keep there prop between them (1-prop.I3)
    race_mult <- numeric(nrow(init))
    race_mult[gp.gn$white] <- calib.target$late_diag.obs[1,"W"]*I3_multi[gp18.gn$white][1]
    race_mult[gp.gn$black] <- calib.target$late_diag.obs[1,"B"]*I3_multi[gp18.gn$black][1]
    race_mult[gp.gn$hisp ] <- calib.target$late_diag.obs[1,"H"]*I3_multi[gp18.gn$hisp ][1]
    I3    <- init[ , "I"] * race_mult
    Ia12  <- (init[ , "I"]-I3) *cbind(prop.Ia/(1-prop.I3), prop.I1/(1-prop.I3), prop.I2/(1-prop.I3))
    I  = cbind(Ia12,I3) #Check: rowSums(I) = init[,"I"], and I3/rowSums(I) = calib.target$late_diag.obs[1,]
    
    D  = init[ ,"D"] *cbind(prop.Da, prop.D1, prop.D2, prop.D3)
#   T  = init[ ,"T"] *cbind(prop.T1, prop.T2, prop.T3)
    prop.T <- prop.T1/(prop.T1+prop.T2) #New par: Proportion of Prop_T1 among T1+T2, can replace current prop.T, rep(0.6,30)
    T12  = init[ ,"T12"] *cbind(prop.T, (1- prop.T))  
    T3 = init[ ,"T3"]
    O  = init[ ,"O"] *cbind(prop.O1, prop.O2, prop.O3)

#    init2 = cbind(S1, S2, I, D, T, O) 
    init2 = cbind(S1, S2, I, D, T12,T3, O) 
    
    colnames(init2) = vparameters$state.name[c(1:16)]
#browser() 
# Check viral_sup: rowSums(init2[,c("T1","T2")])/rowSums(init2[,c(7:16)]) = calib.target$viral_sup.obs[1,]
# Check diag.allmsm: sum(init2[gp.gn$all.msm ,c(7:16)]) 
# Check diag.het.m:  sum(init2[gp.gn$het.m,c(7:16)])
# Check diag.het.f:  sum(init2[gp.gn$het.f,c(7:16)])
# Check diag.pwid:   sum(init2[gp.gn$pwid ,c(7:16)]) 
# Check diag.W: 		 sum(init2[gp.gn$white,c(7:16)])
# Check diag.B: 		 sum(init2[gp.gn$black,c(7:16)])
# Check diag.H: 		 sum(init2[gp.gn$hisp ,c(7:16)]) 
    return(init2) #Convert init 30x5 to 30x16
  })
}

  