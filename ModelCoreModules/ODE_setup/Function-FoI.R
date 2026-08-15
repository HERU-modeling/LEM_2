################################################################################
## Force of infection function
## Transmission matrix lambda
## Details see supplement of P12
## update: Oct 14, 2022
################################################################################

## FoI parameters ##
# y: number of individuals in each compartment (n.gp*19 matrix)
# no: number of opposite sexual partners for each compartment (n.gp*19 matrix)
# uoC: condom use adjustment term for opposite sex (n.gp*19 matrix)
# ns: number of same sexual partners for each compartment (18*19 matrix)
# usC: condom use adjustment term for same sex (18*19 matrix)
# eO: assortative coefficient for heterosexual mixing between ethnic groups (vector of 3 for white, black, hispanic) different by risk level
# eS: assortative coefficient for homosexual mixing between ethnic groups   (vector of 3 for white, black, hispanic) different by risk level
# sigmaFM: probability of transmission by F->M for 16 HIV+ states
# sigmaMF: probability of transmission by M->F for 16 HIV+ states
# simgaM: probability of transmission by same sex for 16 HIV+ states
# tau: probability of transmission by injection for 16 HIV+ states
# eff.prep: % reduction in risk of infection for PrEP
# d: number of injection
# s: proportion of injections that are shared (9 groups)
# eff.oat: % reduction in # of shared injections reduced due to OAT
# cov.ssp: coverage of syringe services program
# eff.ssp: Percentage reduction for shared injections due to SSP
# s.multi: Reduction rate of sharing injection probability after diagnosis

# outputs: matrix of c(lambda, gamma, beta_o,  beta_s),
# each of them with dim (n.gp,2) with row of each group, second column for PrEP
# lambda- over all transition, gamma - needle sharing, beta_o - heterosexual, beta_s - same sex

# sufficient contact rate

#' Generate transmission matrix, lambda, gamma, beta_o, beta_s 
#'
#' @param y number of individuals in each compartment 
#' @param no number of opposite sexual partners for each compartment
#' @param uoC condom use adjustment term for opposite sex
#' @param ns number of same sexual partners for each compartment
#' @param usC condom use adjustment term for same sex
#' @param eO assortative coefficient for heterosexual mixing between ethnic groups (vector of 3 for white, black, hispanic) different by risk level
#' @param eS assortative coefficient for homosexual mixing between ethnic groups (vector of 3 for white, black, hispanic) different by risk level
#' @param sigmaFM Probability of transmission by F->M for 16 HIV+ states
#' @param sigmaMF Probability of transmission by M->F for 16 HIV+ states
#' @param sigmaM Probability of transmission by same sex for 16 HIV+ states
#' @param tau Probability of transmission by injection for 16 HIV+ states
#' @param eff.prep Percentage reduction in risk of infection for PrEP
#' @param d Number of injections
#' @param s Proportion of the injections are shared (9 dim)
#' @param eff.oat Percent reduction in number of shared injections reduced to OAT
#' @param cov.ssp Coverage of syringe services programs
#' @param eff.ssp Percentage reduction for shared injections due to SSP
#' @param s.multi Reduction rate of sharing injection probability after diagnosis
#' @param bal Balance the number of sexual partnership between males & females
#'
#' @return a list(lambda = lambda, gamma = gamma, beta_o = beta_o, beta_s = beta_s)
#' @export
#'
#' @examples
FoI <- function(y, no, uoC, ns, usC, eO, eS, sigmaFM, sigmaMF, sigmaM, tau, eff.prep, d, s,
                eff.oat, prop.ssp.t, eff.ssp = eff.ssp, s.multi, bal = FALSE, prop.oat.t, adherence.prep, mix.prob.HL, mix.prob.YO) {
  with(as.list(gp.gn), {
    #### indicator for gender, ethnicity, risk groups ####
    names.gp <- row.names(y)
    n.gp <- length(names.gp)

    # set the max e, need to revisit
    eO[eO > 1] <- 1 # 3 low w b h and 3 high w b h
    eS[eS > 1] <- 1

    # total number of sexual partners, 1 number
    #' Calculate total number of sexual partners
    #'
    #' \code{xnSum} xx
    #'
    #' Description
    #' @param e assortativity
    #' @param ind ind is a vector indicating the locations the interesting group among 42 groups (should be fixed via indexing)
    #' @param n number of sexual partners for specific group
    #' @param uC condom adjustment
    #' @return
    #' A single number for all sexual partnerships (?)
    #' @examples
    #' g
    #' @export
    xnSum <- function(e, ind, n, uC) {
#      sum((1 - e) * y[ind, ] * n[ind, ] * uC[ind, ])
    	sum((1 - e) * y[ind, ] * n[ind, ])
    }

    ########### mixing matrix mij ##########
    #' Calculate mixing matrix
    #'
    #' \code{mixMatrix} xx
    #'
    #' Description
    #' @param ind ind is a vector indicating the locations the interesting group among 42 groups (should be fixed via indexing); for male partner it has a length of 30 (42 excludes 12 female related strata)
    #' @param ind.low
    #' @param ass.e assortativity coefficient
    #' @param n number of sexual partners for specific group
    #' @param uC condom adjustment
    #' @param indM logical (TRUE if male)
    #' @param indF logical (TRUE if female)
    #' @return
    #' A matrix of dimension x
    #' @examples
    #' generate
    #' @export
#    mixMatrix <- function(ind, ind.low, ass.e, n, uC, indM = FALSE, indF = FALSE) {
    mixMatrix <- function(ind, ind.low, ass.e, n, indM = FALSE, indF = FALSE) {
      if (indM == TRUE) {
#        ass.e_pt <- c(rep(ass.e[4:6], 8), ass.e)     # assortativeness for   male partner: (Originally 30: 6 MSM (high/low), 12 MSM-PWID (high/low/OAT), 6 PWID (high/low/OAT), 3 HET/low, 3 HET/high
      	ass.e_pt <- c(rep(ass.e[4:6], 5), ass.e)     # assortativeness for   male partner: 21 include 6 MSM (young/old), 6 MSM-PWID (young/old), 3 PWID, 3 HET/low, 3 HET/high
      } else if (indF == TRUE) {
#        ass.e_pt <- c(ass.e[4:6], ass.e[4:6], ass.e) # assortativeness for female partner: (Originally 12: 6 PWID, 3 HET/low, 3 HET/high)
      	ass.e_pt <- c(ass.e[4:6], ass.e) # assortativeness for female partner: 9 include 3 PWID, 3 HET/low, 3 HET/high #Combine w/ 21 from M to 30 (gp.gn)
      } else {
        ass.e_pt <- ass.e
      }
#browser() # For manually check xnSum: sum((1 - ass.e_pt) * y[ind, ] * n[ind, ] * uC[ind, ])
      # total number of sexual partnerships for a given group, before distributing them across mixing matrix entries.
#      xn.total <- xnSum(e = ass.e_pt, ind = ind, n = n, uC = uC)
      xn.total <- xnSum(e = ass.e_pt, ind = ind, n = n)

      # n: number of sexual partnership for specific group,
      # ind is a vector indicating the locations the interesting group among 42 groups

      # sum(sus & undiagosed states)* condom adjusted sexual partnership, length(ind) * 1
      # 6 states: S1,S2,Ia,I1,I2,I3
#      xn.undiag <- rowSums(y[ind, 1:6]) * n[ind, 1] * uC[ind, 1]
      xn.undiag <- rowSums(y[ind, 1:6]) * n[ind, 1]

      # sum(diagnosed states) * condom adjusted sexual partnership, length(ind) * 1
      # 10 states: Da,D1,D2,D3,T1,T2,T3,O1,O2,O3
#      xn.diag <- rowSums(y[ind, 7:16]) * n[ind, 7] * uC[ind, 7]
      xn.diag <- rowSums(y[ind, 7:16]) * n[ind, 7]
      xn.j <- matrix((1 - ass.e_pt) * c(xn.undiag, xn.diag), 6, length(ind) * 2, byrow = TRUE) # 1-ass.e_pt random partnership, 6 by 2*length(ind)
#browser()
      l.loc <- match(ind.low, ind) # index location of het.low group in ind
      h.loc <- setdiff(seq_along(ind), l.loc) # index location of not het.low group among ind, all be high

      # not include the cols in l.loc, high risk, 3 by (length(h.loc)*2)
      # row 3 low risk by race, high-low risk indivdual partnership
      high.prop <- xn.j[1:3, c(-l.loc, -(length(ind) + l.loc))] / rowSums(xn.j[1:3, c(-l.loc, -(length(ind) + l.loc))])
      # not include the cols in h.loc, low risk, 3 by (length(l.loc)*2)
      # row 3 high risk by race
      low.prop <- xn.j[4:6, c(-h.loc, -(length(ind) + h.loc))] / rowSums(xn.j[4:6, c(-h.loc, -(length(ind) + h.loc))])
#browser()
      # To avoid overmix, standardized mixing prob within same group to 99% & between dif groups to 1%- This step is to clear out the corresponded cells in r.mix & to be populated later
      xn.j[1:3, c(-l.loc, -(length(ind) + l.loc))] <- 0 
      xn.j[4:6, c(-h.loc, -(length(ind) + h.loc))] <- 0

      # random mixing, 6 by 2*length(ind)
      r.mix <- xn.j / rowSums(xn.j)

      if (xn.total == 0) r.mix <- 0

#### Above are the random mixing (r.mix: left part of equa 1), below are Kronecker delta (k.delta.prop) - Apply for assortative mixing ( mixing within the same group, within which would be fully random)
      # xn.j without ass.e, similar to xn.j with ass.e=0 or mixing at random
      xn.j.wo.e <- matrix(c(xn.undiag, xn.diag), 6, length(ind) * 2, byrow = TRUE)

      # several diagonal matrices, horizontally aligned
      # low risk, 3 by (length(ind)*2) with l.loc cols and 2 *l.loc cols to be diag(1, 3, 3)
      k.delta.l <- matrix(0, 3, length(ind))

      k.delta.l[, l.loc] <- diag(1, 3, 3)
      k.delta.l <- cbind(k.delta.l, k.delta.l)

      # high risk, 3 by (length(ind)*2) with h.loc cols and 2 *h.loc cols to be diag(1, 3, 3)
      k.delta.h <- matrix(0, 3, length(ind))

      k.delta.h[, h.loc] <- diag(1, 3, 3)
      k.delta.h <- cbind(k.delta.h, k.delta.h)

      # combine low/high
      k.delta <- rbind(k.delta.l, k.delta.h) # 6 x length(ind)*2

      # among the same ethnicity, mixing proportion depends on the number of partners
      k.delta.prop <- (xn.j.wo.e * k.delta) / rowSums(xn.j.wo.e * k.delta)
####
      # compute the mixing matrix
      mmr <- ass.e * k.delta.prop + (1 - ass.e) * r.mix # 6 x (length(ind)*2), ass.e length 6

     if (indM == TRUE || indF == TRUE) {
      mmr[1:3, c(-l.loc, -(length(ind) + l.loc))] <- mix.prob.HL * high.prop # HET/high-low risk partnership for mf & mO
      mmr[4:6, c(-h.loc, -(length(ind) + h.loc))] <- mix.prob.HL * low.prop # HET/high-low risk partnership
#browser()      
     } else {
     	mmr[1:3, c(-l.loc, -(length(ind) + l.loc))] <- mix.prob.YO * high.prop # MSM/young-old partnership for mS
     	mmr[4:6, c(-h.loc, -(length(ind) + h.loc))] <- mix.prob.YO * low.prop # MSM/young-old partnership
#browser()      	
     }
      mmr <- mmr / rowSums(mmr) # 6 x (length(ind)*2), normalize to be probability
    
      return(mmr)
    }

    # mij for female (index) to male (partner):  proportion of total heterosexual partnerships that females in each stratified group (by race/ethnicity and risk level) form with each male group.   (6 rows; length(m)*2 columns)
#    mf <- mixMatrix(male, L.het.m, eO, no, uoC, indM = TRUE)
    mf <- mixMatrix(male, L.het.m, eO, no, indM = TRUE)
    # ind = m; ind.low = L.het.m; ass.e = eO; n = no; uC = uoC; indM = T
    # only "HET/male/low/white"    "HET/male/low/black"    "HET/male/low/hispanic", low

    # write.table(mf, file = file.path("Data_files", "Intermediate", paste0("Matrixmf-", CITY, ".csv")), append = TRUE, row.names = FALSE, sep = ",")
#browser()    

    # mij for female to male  opposite sex (6 rows; length(f)*2 columns)
#    mO <- mixMatrix(female, L.het.f, eO, no, uoC, indF = TRUE)
    mO <- mixMatrix(female, L.het.f, eO, no, indF = TRUE)

    # write.table(mO, file = file.path("Data_files", "Intermediate", paste0("MatrixmO-", CITY, ".csv")), append = TRUE,  row.names = F, sep = ",")


    # mij for male same sex (3 rows; length(all.msm)*2 columns)
#    mS <- mixMatrix(all.msm, old, eS, ns, usC)
    mS <- mixMatrix(all.msm, old, eS, ns)
    # ind = all.msm; ind.low = L.allmsm; ass.e = eS; n = ns; uC = usC; indM = F; indF = F


    ########### balance the sexual partnership ##########
    ## need to revisit later to see whether there is a better way to do it

    if (bal == TRUE) {
      # number of total sexual partners for female (X*n*mij),
      # stratified by ethnicity and risk behaviour (HET low vs. others)
      xnmf <- matrix(0, 6, 6)
      xnmf.tot <- numeric(6)
#browser()     
      # total partners for female het low
#      xnmf.tot[1] <- xnSum(0, intersect(white, L.het.f), no, uoC) # partners for white females
#      xnmf.tot[2] <- xnSum(0, intersect(black, L.het.f), no, uoC) # partners for black females
#      xnmf.tot[3] <- xnSum(0, intersect(hisp, L.het.f), no, uoC) # partners for hisp females
      xnmf.tot[1] <- xnSum(0, intersect(white, L.het.f), no) # partners for white females
      xnmf.tot[2] <- xnSum(0, intersect(black, L.het.f), no) # partners for black females
      xnmf.tot[3] <- xnSum(0, intersect(hisp, L.het.f), no) # partners for hisp females
      # total partners for female others
#      xnmf.tot[4] <- xnSum(0, intersect(white, high.f), no, uoC) # partners for white females
#      xnmf.tot[5] <- xnSum(0, intersect(black, high.f), no, uoC) # partners for black females
#      xnmf.tot[6] <- xnSum(0, intersect(hisp, high.f), no, uoC) # partners for hisp females
      xnmf.tot[4] <- xnSum(0, intersect(white, high.f), no) # partners for white females
      xnmf.tot[5] <- xnSum(0, intersect(black, high.f), no) # partners for black females
      xnmf.tot[6] <- xnSum(0, intersect(hisp, high.f), no) # partners for hisp females

      # number of partners by ethnicity and risk, mf from male to female
      l.loc.m <- (rep(male, 2) %in% L.het.m) # true if het male low
      np.l <- sum(l.loc.m) / 3

      other.loc.m <- (rep(male, 2) %in% high.m) # true if male others
      np.h <- sum(other.loc.m) / 3
      for (i in 1:6) {
        xnmf[i, 1] <- sum(mf[i, l.loc.m][(1:np.l) * 3 - 2] * xnmf.tot[i]) # white het male low partners
        xnmf[i, 2] <- sum(mf[i, l.loc.m][(1:np.l) * 3 - 1] * xnmf.tot[i]) # black het male low partners
        xnmf[i, 3] <- sum(mf[i, l.loc.m][(1:np.l) * 3] * xnmf.tot[i]) # hisp het male low partners
        xnmf[i, 4] <- sum(mf[i, other.loc.m][(1:np.h) * 3 - 2] * xnmf.tot[i]) # white male other partners
        xnmf[i, 5] <- sum(mf[i, other.loc.m][(1:np.h) * 3 - 1] * xnmf.tot[i]) # black male other partners
        xnmf[i, 6] <- sum(mf[i, other.loc.m][(1:np.h) * 3] * xnmf.tot[i]) # hisp male other partners
      }
      
#browser()   # Check unadjusted xnmf, sum(xnmf) (due to unadjsuted mf)

      # number of total sexual partners for male,
      # stratified by ethnic groups (X*n*mij)
      xnmm <- matrix(0, 6, 6)
      xnmm.tot <- numeric(6)
      # total partners for male het low
#      xnmm.tot[1] <- xnSum(0, intersect(white, L.het.m), no, uoC) # partners for white males
#      xnmm.tot[2] <- xnSum(0, intersect(black, L.het.m), no, uoC) # partners for black males
#      xnmm.tot[3] <- xnSum(0, intersect(hisp, L.het.m), no, uoC) # partners for hisp males
      xnmm.tot[1] <- xnSum(0, intersect(white, L.het.m), no) # partners for white males
      xnmm.tot[2] <- xnSum(0, intersect(black, L.het.m), no) # partners for black males
      xnmm.tot[3] <- xnSum(0, intersect(hisp, L.het.m), no) # partners for hisp males
      # total partners for male others
#      xnmm.tot[4] <- xnSum(0, intersect(white, high.m), no, uoC) # partners for white males
#      xnmm.tot[5] <- xnSum(0, intersect(black, high.m), no, uoC) # partners for black males
#      xnmm.tot[6] <- xnSum(0, intersect(hisp, high.m), no, uoC) # partners for hisp males
      xnmm.tot[4] <- xnSum(0, intersect(white, high.m), no) # partners for white males
      xnmm.tot[5] <- xnSum(0, intersect(black, high.m), no) # partners for black males
      xnmm.tot[6] <- xnSum(0, intersect(hisp, high.m), no) # partners for hisp males

      # number of partners by ethnicity and risk, mO from female to male
      l.loc.f <- (rep(female, 2) %in% L.het.f) # true if het female low
      np.l <- sum(l.loc.f) / 3
      other.loc.f <- (rep(female, 2) %in% high.f) # true if female others
      np.h <- sum(other.loc.f) / 3

      for (i in 1:6) {
        xnmm[i, 1] <- sum(mO[i, l.loc.f][(1:np.l) * 3 - 2] * xnmm.tot[i]) # white female partners,low
        xnmm[i, 2] <- sum(mO[i, l.loc.f][(1:np.l) * 3 - 1] * xnmm.tot[i]) # black female partners,low
        xnmm[i, 3] <- sum(mO[i, l.loc.f][(1:np.l) * 3] * xnmm.tot[i]) # hisp female partners, low
        xnmm[i, 4] <- sum(mO[i, other.loc.f][(1:np.h) * 3 - 2] * xnmm.tot[i]) # white female partners, high
        xnmm[i, 5] <- sum(mO[i, other.loc.f][(1:np.h) * 3 - 1] * xnmm.tot[i]) # black female partners , high
        xnmm[i, 6] <- sum(mO[i, other.loc.f][(1:np.h) * 3] * xnmm.tot[i]) # hisp female partners , high
      }
#browser()  # Check unadjusted xnmm, sum(xnmm) (due to unadjusted mO)
      
      # need to revisit the following two lines
      # imbalance. xnmm f to m, xnmf m to f.
      v <- xnmm / t(xnmf)
      v[is.na(v)] <- 0 # for fully assortative mixing

      # save the imbalance matrix V
      # write.table(v, file = file.path("Data_files", "Intermediate", paste0("MatrixV-", CITY, ".csv")), append = TRUE, row.names = F, sep = ",")

      # degree of compromise between the two sexes.
      theta <- 0.5 # set to 0.5 to be equal between the two groups

      # adjust the imbalance, male to female
      imb.adj.m <- matrix(0, 6, length(male) * 2)
      imb.adj.m[, l.loc.m] <- v[, 1:3]^theta
      imb.adj.m[, other.loc.m] <- v[, 4:6]^theta
#browser() # check unadjsuted mf, mO, v, t(xnmf)      - Everything above is use to calculate V then, V^theta. mf does not change before this point 
      mf <- mf * imb.adj.m
#browser() # Check new mf
      # write.table(mf, file = file.path("Data_files", "Intermediate", paste0("Matrixbalmf-", CITY, ".csv")), append = TRUE, row.names = F,sep = ",")

      # female to male
      imb.adj.f <- matrix(0, 6, length(female) * 2)
      imb.adj.f[, l.loc.f] <- v[, 1:3]^(theta - 1)
      imb.adj.f[, other.loc.f] <- v[, 4:6]^(theta - 1)
      imb.adj.f[mapply(is.infinite, imb.adj.f)] <- 0 # could it be infinite?
      
      mO <- mO * imb.adj.f
#browser() # Check new mO
      # write.table(mO, file = file.path("Data_files", "Intermediate", paste0("MatrixbalmO-", CITY, ".csv")), append = TRUE, row.names = F, sep = ",")
    }


### Tri 7-21- Added code to check to see if balance adjustment result in v=1

#browser()
# number of total sexual partners for female (X*n*mij),
# stratified by ethnicity and risk behaviour (HET low vs. others)
xnmf <- matrix(0, 6, 6)
xnmf.tot <- numeric(6)

# total partners for female het low
xnmf.tot[1] <- xnSum(0, intersect(white, L.het.f), no) # partners for white females
xnmf.tot[2] <- xnSum(0, intersect(black, L.het.f), no) # partners for black females
xnmf.tot[3] <- xnSum(0, intersect(hisp, L.het.f), no) # partners for hisp females
# total partners for female others
xnmf.tot[4] <- xnSum(0, intersect(white, high.f), no) # partners for white females
xnmf.tot[5] <- xnSum(0, intersect(black, high.f), no) # partners for black females
xnmf.tot[6] <- xnSum(0, intersect(hisp, high.f), no) # partners for hisp females

# number of partners by ethnicity and risk, mf from male to female
l.loc.m <- (rep(male, 2) %in% L.het.m) # true if het male low
np.l <- sum(l.loc.m) / 3

other.loc.m <- (rep(male, 2) %in% high.m) # true if male others
np.h <- sum(other.loc.m) / 3
for (i in 1:6) {
	xnmf[i, 1] <- sum(mf[i, l.loc.m][(1:np.l) * 3 - 2] * xnmf.tot[i]) # white het male low partners
	xnmf[i, 2] <- sum(mf[i, l.loc.m][(1:np.l) * 3 - 1] * xnmf.tot[i]) # black het male low partners
	xnmf[i, 3] <- sum(mf[i, l.loc.m][(1:np.l) * 3] * xnmf.tot[i]) # hisp het male low partners
	xnmf[i, 4] <- sum(mf[i, other.loc.m][(1:np.h) * 3 - 2] * xnmf.tot[i]) # white male other partners
	xnmf[i, 5] <- sum(mf[i, other.loc.m][(1:np.h) * 3 - 1] * xnmf.tot[i]) # black male other partners
	xnmf[i, 6] <- sum(mf[i, other.loc.m][(1:np.h) * 3] * xnmf.tot[i]) # hisp male other partners
}

#browser() #New adjusted mf, sum(xnmf) should still be the same, but individual cell values should have changed

# number of total sexual partners for male,
# stratified by ethnic groups (X*n*mij)
xnmm <- matrix(0, 6, 6)
xnmm.tot <- numeric(6)
# total partners for male het low
xnmm.tot[1] <- xnSum(0, intersect(white, L.het.m), no) # partners for white males
xnmm.tot[2] <- xnSum(0, intersect(black, L.het.m), no) # partners for black males
xnmm.tot[3] <- xnSum(0, intersect(hisp, L.het.m), no) # partners for hisp males
# total partners for male others
xnmm.tot[4] <- xnSum(0, intersect(white, high.m), no) # partners for white males
xnmm.tot[5] <- xnSum(0, intersect(black, high.m), no) # partners for black males
xnmm.tot[6] <- xnSum(0, intersect(hisp, high.m), no) # partners for hisp males

# number of partners by ethnicity and risk, mO from female to male
l.loc.f <- (rep(female, 2) %in% L.het.f) # true if het female low
np.l <- sum(l.loc.f) / 3
other.loc.f <- (rep(female, 2) %in% high.f) # true if female others
np.h <- sum(other.loc.f) / 3

for (i in 1:6) {
	xnmm[i, 1] <- sum(mO[i, l.loc.f][(1:np.l) * 3 - 2] * xnmm.tot[i]) # white female partners,low
	xnmm[i, 2] <- sum(mO[i, l.loc.f][(1:np.l) * 3 - 1] * xnmm.tot[i]) # black female partners,low
	xnmm[i, 3] <- sum(mO[i, l.loc.f][(1:np.l) * 3] * xnmm.tot[i]) # hisp female partners, low
	xnmm[i, 4] <- sum(mO[i, other.loc.f][(1:np.h) * 3 - 2] * xnmm.tot[i]) # white female partners, high
	xnmm[i, 5] <- sum(mO[i, other.loc.f][(1:np.h) * 3 - 1] * xnmm.tot[i]) # black female partners , high
	xnmm[i, 6] <- sum(mO[i, other.loc.f][(1:np.h) * 3] * xnmm.tot[i]) # hisp female partners , high
}

#browser() #New adjusted mO, sum(xnmm) should still be the same, but individual cell values should have changed

# need to revisit the following two lines
# imbalance. xnmm f to m, xnmf m to f.
v <- xnmm / t(xnmf)

#browser() # Check v (after adjustment) - Should be 1
###

    #' Sum of sufficient contact rates among susceptible i with HIV-infected j
    #'
    #' \code{betaI} xx
    #'
    #' Description
    #' @param i susceptible index (male - susceptible (S states) HET/M, female - sus HET/F, all.msm - sus MSM)
    #' @param ind ind is a vector indicating the locations the interesting group among 42 groups (should be fixed via indexing)
    #' @param mij passed in row by row of mf,mO, mS
    #' @param sigma Probability of transmission by F->M or M->F or M->M for 16 HIV+ states
    #' @return
    #' A matrix of dimension x
    #' @examples
    #' generate
    #' @export
    betaI <- function(i, ind, mij, sigma, homo = FALSE) {
      # den_beta is the number of potential partners among 6 undiagnosed (S1,S2,Ia,I1,I2,I3) AND diag(D,T,O states) 
    	# => (length(ind)*2) since length(ind) for either undiag (1:6) OR diag 
#browser(expr = {i >= 16})    	#Note the i only match up to 15, after that 16-21 corresponds 19-24 male, 
      den_beta <- c(rowSums(y[ind, 1:6]), rowSums(y[ind, 7:16]))
      # c.prob: Monthly PROB of a single person from any strata select a partnership 
      # Probability of mixing (mij) divided by number of potential partners (den_beta)
      c.prob <- mij / den_beta # 6 by (length(ind)*2)
      # c.prob2: make the c.prob into matrix of dim=c(length(ind),14) where 14 is the HIV+ states (I,D,T,O states)
      # 4 infected but undiagnosed (Ia,I1,I2,I3)
      # 10 diagnosed (Da,D1,D2,D3,T1,T2,T3,O1,O2,O3)
      
      c.prob2 <- matrix(
        c(
          rep(c.prob[seq_along(ind)], 4),
          rep(c.prob[(length(ind) + 1):(length(ind) * 2)], 10)
        ),
        length(ind), 14
      )

      # Monthly probability selecting j (infected) compartment, length(ind) * 14, 14 hiv+ states
      #prob <- y[ind, 4:19] * c.prob2
       prob <- y[ind, 3:16] * c.prob2
       
# browser()
if (any(is.nan(-log(1 - prob * rep(sigma, each = length(ind)))))) {
stop("NaNs in log argument")
}

      # Rate of monthly HIV infection - not on PreP, length(ind) * 14 infected states
      beta <- -log(1 - prob * rep(sigma, each = length(ind)))

      # Rate of monthly HIV infection - if on PreP, the transmission prob. will be reduced by eff.prep
      betaP <- -log(1 - prob * rep(sigma * (1 - eff.prep * adherence.prep) , each = length(ind)))

      if (homo == TRUE) {
        # susceptible states number of partners
        nuC <- ns[i, 1] * usC[i, 1]
      } else {
        nuC <- no[i, 1] * uoC[i, 1]
      }
 #browser(expr = {i >= 16})     
      return(c(
        nuC * sum(beta), # no PrEP - 1st Col of beta_o or beta_s
        nuC * sum(betaP) # PrEP    - 2nd Col of beta_o or beta_s
      )) 
    }

    ####### sufficient contact rate for heterosex #######
    beta_o <- matrix(0, n.gp, 2) # on PrEP and not on PrEP

    for (i in male) { # heterosexual for male susceptible
      if (i %in% L.het) {
        if (i %in% white) beta_o[i, ] <- betaI(i = i, ind = female, mij = mO[1, ], sigma = sigmaFM)
        if (i %in% black) beta_o[i, ] <- betaI(i = i, ind = female, mij = mO[2, ], sigma = sigmaFM)
        if (i %in% hisp) beta_o[i, ] <- betaI(i = i, ind = female, mij = mO[3, ], sigma = sigmaFM)
      } else {
        if (i %in% white) beta_o[i, ] <- betaI(i = i, ind = female, mij = mO[4, ], sigma = sigmaFM)
        if (i %in% black) beta_o[i, ] <- betaI(i = i, ind = female, mij = mO[5, ], sigma = sigmaFM)
        if (i %in% hisp) beta_o[i, ] <- betaI(i = i, ind = female, mij = mO[6, ], sigma = sigmaFM)
      }
    }
    for (i in female) { # heterosexual for female susceptible
      if (i %in% L.het) {
        if (i %in% white) beta_o[i, ] <- betaI(i = i, ind = male, mij = mf[1, ], sigma = sigmaMF)
        if (i %in% black) beta_o[i, ] <- betaI(i = i, ind = male, mij = mf[2, ], sigma = sigmaMF)
        if (i %in% hisp) beta_o[i, ] <- betaI(i = i, ind = male, mij = mf[3, ], sigma = sigmaMF)
      } else {
        if (i %in% white) beta_o[i, ] <- betaI(i = i, ind = male, mij = mf[4, ], sigma = sigmaMF)
        if (i %in% black) beta_o[i, ] <- betaI(i = i, ind = male, mij = mf[5, ], sigma = sigmaMF)
        if (i %in% hisp) beta_o[i, ] <- betaI(i = i, ind = male, mij = mf[6, ], sigma = sigmaMF)
      }
    }
#browser()
    ###### sufficient contact rate for homosex #######
    beta_s <- matrix(0, n.gp, 2)

    for (i in all.msm) {
      if (i %in% old) {
        if (i %in% white) beta_s[i, ] <- betaI(i = i, ind = all.msm, mij = mS[1, ], sigma = sigmaM, homo = TRUE)
        if (i %in% black) beta_s[i, ] <- betaI(i = i, ind = all.msm, mij = mS[2, ], sigma = sigmaM, homo = TRUE)
        if (i %in% hisp) beta_s[i, ] <- betaI(i = i, ind = all.msm, mij = mS[3, ], sigma = sigmaM, homo = TRUE)
      } else {
        if (i %in% white) beta_s[i, ] <- betaI(i = i, ind = all.msm, mij = mS[4, ], sigma = sigmaM, homo = TRUE)
        if (i %in% black) beta_s[i, ] <- betaI(i = i, ind = all.msm, mij = mS[5, ], sigma = sigmaM, homo = TRUE)
        if (i %in% hisp) beta_s[i, ] <- betaI(i = i, ind = all.msm, mij = mS[6, ], sigma = sigmaM, homo = TRUE)
      }
    }
#browser()
    ###### sufficient contact rate for shared needle injection #######
    gamma <- matrix(0, n.gp, 2)
    ds.all <- numeric(n.gp) # 0 if not PWID

    ds <- d * s #T:d from Master, s from master but also convert in ode_model
    # Weighted average impact of OAT: Reduce needle shared
#     ds.all[all.pwid] <- ds[all.pwid] * (1 - cov.ssp[all.pwid] * (1 - eff.ssp)) * prop.oat[all.pwid] * (1-eff.oat) + 
#    	                  ds[all.pwid] * (1 - cov.ssp[all.pwid] * (1 - eff.ssp)) * (1-prop.oat[all.pwid]) 
     ds.all[all.pwid] <- ds[all.pwid] * (1 - prop.ssp.t[all.pwid] * (1 - eff.ssp)) * prop.oat.t[all.pwid] * (1-eff.oat) + 
        	               ds[all.pwid] * (1 - prop.ssp.t[all.pwid] * (1 - eff.ssp)) * (1-prop.oat.t[all.pwid]) 
     
    # if diagnosed, sharing prob gets decreased by s.multi
    ds.all2 <- matrix(ds.all, n.gp, 16) # T:convert ds.all into a mtrix 30 groups by 16 states
    ds.all2[, 7:16] <- ds.all * (1 - s.multi) # T:Updates the col (states) 7 (Da)-> 16 (O3) with value ds.all multiply by (1-s.multi) instead of ds.all
    den_gamma <- sum(y * ds.all2) # denominator of selecting a partner

    # not on PreP
    gamma.sum <- sum(-log(1 - (y[, -(1:2)] * ds.all2[, -(1:2)] / den_gamma) * rep(tau, each = n.gp)))
    # on PreP
    gamma.sum.p <- sum(-log(1 - (y[, -(1:2)] * ds.all2[, -(1:2)] / den_gamma) * rep(tau, each = n.gp) * (1 - eff.prep * adherence.prep) ))

    if (den_gamma == 0) {
      gamma.sum <- 0
      gamma.sum.p <- 0
    }
    gamma <- cbind(
      ds.all * gamma.sum, # no PrEP
      ds.all * gamma.sum.p #   PrEP
    ) 
#browser()
    ##### overall sufficient contact rate #####
    lambda <- gamma + beta_o + beta_s # each of them n.gp by 2 - for on/off prep

    return(list(lambda = lambda, gamma = gamma, beta_o = beta_o, beta_s = beta_s))
  })
}
