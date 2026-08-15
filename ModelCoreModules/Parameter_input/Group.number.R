################################################################################
#### Functions to derive group number within different sub domains          ####
################################################################################

# subgroups under pwid
pwid.gn.fun = function (names.pwid){
  female = grep("female",  names.pwid)
  male   = grep("female",  names.pwid, invert =T)
  white  = grep("white",   names.pwid)
  black  = grep("black",   names.pwid)
  hisp   = grep("hispanic",names.pwid)
  F.w    = intersect(female, white)
  F.b    = intersect(female, black)
  F.h    = intersect(female, hisp)
  M.w    = intersect(male, white)
  M.b    = intersect(male, black)
  M.h    = intersect(male, hisp)
  return(list(female=female, male=male, white=white, black=black, hisp=hisp, F.w=F.w, F.b=F.b, F.h=F.h, M.w=M.w, M.b=M.b, M.h=M.h))
}

# subgroups of msm
msm.gn.fun = function (names.msm){
  white  = grep("white",   names.msm)
  black  = grep("black",   names.msm)
  hisp   = grep("hispanic",names.msm)
  young  = grep("young",   names.msm)
  old    = grep("old"  ,   names.msm)
  mwid   = grep("MSM/PWID",names.msm)
  all.msm= grep("MSM",     names.msm)
  msm    = setdiff(all.msm, mwid)
  Y.w    = intersect(young, white)
  Y.b    = intersect(young, black)
  Y.h    = intersect(young, hisp)
  O.w    = intersect(old  , white)
  O.b    = intersect(old  , black)
  O.h    = intersect(old  , hisp)
  return(list(white=white, black=black, hisp=hisp, young=young, old=old, mwid=mwid, msm=msm, Y.w=Y.w, Y.b=Y.b, Y.h=Y.h, O.w=O.w, O.b=O.b, O.h=O.h))
}

# no risk level specified, subgroups under 18 groups(gender*risk group*ethnicity)
gp18.gn.fun = function (names.18){
  white    = grep("white",   names.18)
  black    = grep("black",   names.18)
  hisp     = grep("hispanic",names.18)
  female   = grep("female",  names.18)
  male     = grep("female",  names.18, invert =T)
  all.msm  = grep("MSM",     names.18)
  all.pwid = grep("PWID",    names.18)
  mwid     = grep("MSM/PWID",names.18)
  het      = grep("HET",     names.18)
  het.m    = intersect(het, male)
  het.f    = intersect(het, female)
  msm      = setdiff(all.msm, all.pwid)
  pwid     = setdiff(all.pwid, all.msm)
  return(list(white=white, black=black, hisp=hisp, female=female, male=male, all.msm=all.msm, all.pwid=all.pwid, msm=msm, pwid=pwid, mwid=mwid, het=het, het.m=het.m, het.f=het.f))
}

# subgroups among all 30 population groups
gp.gn.fun = function (names.gp){
  white    = grep("white",   names.gp)  # white: 14 groups
  black    = grep("black",   names.gp)  # black: 14 groups
  hisp     = grep("hispanic",names.gp)  # hispanic: 14 groups
  female   = grep("female",  names.gp)  # female: 12 groups
  male     = grep("female",  names.gp, invert =T) # males: 30 groups
  old      = grep("old",     names.gp) # NEW
  young    = grep("young",   names.gp) # NEW
  low      = grep("low",     names.gp) # Low-risk: 15 groups, no PWID
  high     = grep("high",    names.gp) # High-risk: 15 groups, no PWID
  all.msm  = grep("MSM",     names.gp) # All MSM, including MWID: 18 groups
  all.pwid = grep("PWID",    names.gp) # All PWID, including MWID: 24 groups
  msm      = setdiff(all.msm, all.pwid) # MSM only, excluding PWID: 6 groups
  pwid     = setdiff(all.pwid, all.msm)  # PWID only, excluding MSM: 12 groups
  mwid     = grep("MSM/PWID",names.gp)  # MWID: 12 groups
  het      = grep("HET",     names.gp) # HET: 12 groups
# oat      = grep("OAT",     names.gp) # OAT: 12 groups
# off.oat  = setdiff(all.pwid, oat)  # All PWID (including MWID) not on OAT: 12 groups
  Y.msm    = intersect(young, msm)  # Originally H.allmsm
  O.msm    = intersect(old ,  msm)  # Orignally L.allmsm
  Y.mwid   = intersect(young, mwid)  # Originally H.allmsm
  O.mwid   = intersect(old ,  mwid)  # Orignally L.allmsm
  het.m    = intersect(het, male)  # HET male: 6 groups
  het.f    = intersect(het, female) # HET female: 6 groups
  pwid.m   = intersect(pwid, male) # male PWID, excluding MWID: 6 groups
  pwid.f   = intersect(pwid, female)  # female PWID: 6 groups
  H.het.m  = intersect(high, het.m)
  H.het.f  = intersect(high, het.f)
  L.het    = intersect(low, het) # low-risk HET: 6 groups
  L.het.m  = intersect(low, het.m)  # low-risk male HET: 3 groups
  L.het.f  = intersect(low, het.f) # low-risk female HET: 3 groups
  H.m      = union(intersect(young, male),intersect(high, male)) # Originally, intersect(high, male)
  H.f      = intersect(high, female)
  L.m      = union(intersect(old, male),intersect(low, male))  # Originally, intersect(low, male)
  L.f      = intersect(low, female)

  ## used in mixing matrix and balance, definition of low/high group (diff from sexual intensity level)
  high.m   = setdiff(male,  L.het.m)   #high-risk male: Including PWID, high-risk MSM, high-risk MWID, high-risk HET
  high.f   = setdiff(female, L.het.f)  #high-risk female: Including PWID, high-risk HET
  
  
  return(list(white=white, black=black, hisp=hisp, female=female, male=male, 
              low=low, high=high, old=old, young=young, all.msm=all.msm, all.pwid=all.pwid, 
              msm=msm, pwid=pwid, mwid=mwid, het=het, het.m=het.m, het.f=het.f,
              pwid.m=pwid.m, pwid.f=pwid.f, Y.msm=Y.msm, O.msm=O.msm, Y.mwid=Y.mwid, O.mwid=O.mwid,
              H.het.m=H.het.m, H.het.f=H.het.f, L.het = L.het, L.het.m=L.het.m, L.het.f=L.het.f,
              H.m=H.m, H.f=H.f, L.m=L.m,  L.f=L.f, high.m = high.m, high.f = high.f))
}