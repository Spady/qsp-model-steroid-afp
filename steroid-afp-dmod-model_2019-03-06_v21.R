#dMod-based Model of Steroid-Binding Fusion Protein
#Model as submitted with thesis 2019-03-06
#Minor edits for clarity made 2024-01-03

library(dMod)
library(ggplot2)
library(dplyr)
setwd("C:/Users/espad/Documents/Pharmacokinetics Model/R Model Work/2024-01-03_2019modeltest")

#eqnlist$volumes is a structure that needs to have each species name and a matching value
#Reference volume suffixes must be specific to chemical species in setup because of grepl

#Volume and Surface area constants
OfftoOnRatio = 2     #Ratio of off-target cells to on-target cells
VolB = 18          #B = blood
VolC1 = 5.9e-02     #C = cytoplasm
VolE = 3.3e-03    #E = endosome
VolC2 = VolC1*OfftoOnRatio     #Off-target cytoplasm
SAendo = 3.9e+03     #Endosome membrane surface area
SAplasma1 = 1.3e+04    #Plasma membrane surface area
SAplasma2 = SAplasma1*OfftoOnRatio     #Off-target plasma membrane surface area

refvolumes <- numeric(4)
names(refvolumes) <- c("_B", "_C1", "_E", "_C2")
refvolumes["_B"] = VolB    #B = blood
refvolumes["_C1"] = VolC1     #C = cytoplasm
refvolumes["_E"] = VolE    #E = endosome
refvolumes["_C2"] = VolC2    #Off-target cytoplasm+endosome

#Parameter value assignment - starting conc for species, values for rate constants

initamts = numeric()
initamts["abgbxster_B"] = 1e-9
initamts["ster_C1"] = 0
initamts["ster_B"] = 0
initamts["rcep_B"] = 3.99362e-9
initamts["rcep_E"] = 6.53502e-6
initamts["glcr_C1"] = 8.99459E-08
initamts["glcr_C2"] = 8.99459E-08


#Names for these variables must match exactly with names of elements in the ODEsystem setup
#Solver may get grumpy when initamts are 0 and something processes away from that. (I think)
#But it might have just been a function of typos in the addReaction setup. 
#See 2019-02-12 Pharmacokinetics in R Package Notes for how I finally handled the receptor 
#internalization and externalization rates correctly.
#Also it seems that for compartment transfer you need to multiply rate by vol of destination
#compartment. This is an error in how the setup is handled, if I am correct. 
#This destination volume corr ends up included in the receptor transfer already

rateconst <- numeric()
rateconst["k01"] = 1e+06     #antibody-receptor association rate
rateconst["k21"] = 1e+06     #antibody-receptor association rate in endosome
rateconst["k02"] = 3.86E-04     #antibody-receptor dissociation rate
rateconst["k22"] = 3.86E-04     #antibody-receptor dissociation rate in endosome
rateconst["k03"] = 3.98E+04     #glucocorticoid binder-steroid association rate
rateconst["k23"] = 3.98E+04     #glucocorticoid binder-steroid association rate in endosome
rateconst["k04"] = 2.20E-05     #glucocorticoid binder-steroid dissociation rate
rateconst["k24"] = 2.20E-05     #glucocorticoid binder-steroid dissociation rate in endosome
rateconst["k05"] = 5.5e-3 #5.00E-04 is slow;  5.5E-3 is fast   #receptor internalization rate with or without ab (1st order)
rateconst["k25"] = 5.5e-3*(SAplasma1/SAendo)#5.5e-3*(SAplasma1/SAendo) #CHANGE WITH k05!     #Receptor externalization rate (1st order)
rateconst["k06"] = 2.01E-06     #antibody degradation in blood
rateconst["k07"] = 3.08E-05     #steroid degradation in blood
rateconst["k08"] = 2.10E-9*SAplasma1*VolC1    #steroid diffusion from blood to cytoplasm
rateconst["k28"] = 2.10E-9*SAplasma1*VolB     #steroid diffusion from cytoplasm to blood
rateconst["k38"] = 2.10E-9*SAplasma2*VolC2    #steroid diffusion from blood to off-target cell cytoplasm
rateconst["k48"] = 2.10E-9*SAplasma2*VolB     #steroid diffusion from off-target cell cytoplasm to blood
rateconst["k09"] = 2.10E-9*SAendo*VolE      #steroid diffusion from cytoplasm to endosome
rateconst["k29"] = 2.10E-9*SAendo*VolC1     #steroid diffusion from endosome to cytoplasm
rateconst["k10"] = 3.70E-04     #free antibody degradation in endosome
rateconst["k11"] = 1.65e-4 #2.26E-05 is slow; 1.65e-4 is fast     #receptor degradation in endosome; synthesis of new receptor is included
rateconst["k12"] = 3.98E+04     #glucocorticoid receptor on rate
rateconst["k13"] = 2.20E-05     #glucocorticoid receptor off rate


#abgb = antibody-glucocorticoid binding fusion protein
#ster = steroid
#rcep = on-target receptor
#glcr = glucocorticoid receptor
#complex denoted with x
#abgbxsterxrcep is the order

ODEsys <- eqnlist()

#Receptor and steroid binding in blood
ODEsys <- addReaction(ODEsys, "abgbxster_B + rcep_B", "abgbxsterxrcep_B", "k01*abgbxster_B*rcep_B", "RecepAssnB1")
ODEsys <- addReaction(ODEsys, "abgbxsterxrcep_B", "abgbxster_B + rcep_B", "k02*abgbxsterxrcep_B", "RecepDssnB1")
ODEsys <- addReaction(ODEsys, "abgb_B + rcep_B", "abgbxrcep_B", "k01*abgb_B*rcep_B", "RecepAssnB2")
ODEsys <- addReaction(ODEsys, "abgbxrcep_B", "abgb_B + rcep_B", "k02*abgbxrcep_B", "RecepDssnB2")

ODEsys <- addReaction(ODEsys, "abgb_B + ster_B", "abgbxster_B", "k03*abgb_B*ster_B", "SteroidAssnB1")
ODEsys <- addReaction(ODEsys, "abgbxster_B", "abgb_B + ster_B", "k04*abgbxster_B", "SteroidDssnB1")
ODEsys <- addReaction(ODEsys, "abgbxrcep_B + ster_B", "abgbxsterxrcep_B", "k03*abgbxrcep_B*ster_B", "SteroidAssnB2")
ODEsys <- addReaction(ODEsys, "abgbxsterxrcep_B", "abgbxrcep_B + ster_B", "k04*abgbxsterxrcep_B", "SteroidDssnB2")

#Receptor interalization
ODEsys <- addReaction(ODEsys, "rcep_B", "rcep_E", "k05*rcep_B", "RecepIntern1")
ODEsys <- addReaction(ODEsys, "rcep_E", "rcep_B", "k25*rcep_E", "RecepExtern1")
ODEsys <- addReaction(ODEsys, "abgbxrcep_B", "abgbxrcep_E", "k05*abgbxrcep_B", "RecepIntern2")
ODEsys <- addReaction(ODEsys, "abgbxrcep_E", "abgbxrcep_B", "k25*abgbxrcep_E", "RecepExtern2")
ODEsys <- addReaction(ODEsys, "abgbxsterxrcep_B", "abgbxsterxrcep_E", "k05*abgbxsterxrcep_B", "RecepIntern3")
ODEsys <- addReaction(ODEsys, "abgbxsterxrcep_E", "abgbxsterxrcep_B", "k25*abgbxsterxrcep_E", "RecepExtern3")

#Receptor and steroid binding in endosome
ODEsys <- addReaction(ODEsys, "abgbxster_E + rcep_E", "abgbxsterxrcep_E", "k21*abgbxster_E*rcep_E", "RecepAssnE1")
ODEsys <- addReaction(ODEsys, "abgbxsterxrcep_E", "abgbxster_E + rcep_E", "k22*abgbxsterxrcep_E", "RecepDssnE1")
ODEsys <- addReaction(ODEsys, "abgb_E + rcep_E", "abgbxrcep_E", "k21*abgb_E*rcep_E", "RecepAssnE2")
ODEsys <- addReaction(ODEsys, "abgbxrcep_E", "abgb_E + rcep_E", "k22*abgbxrcep_E", "RecepDssnE2")

ODEsys <- addReaction(ODEsys, "abgb_E + ster_E", "abgbxster_E", "k23*abgb_E*ster_E", "SteroidAssnE1")
ODEsys <- addReaction(ODEsys, "abgbxster_E", "abgb_E + ster_E", "k24*abgbxster_E", "SteroidDssnE1")
ODEsys <- addReaction(ODEsys, "abgbxrcep_E + ster_E", "abgbxsterxrcep_E", "k23*abgbxrcep_E*ster_E", "SteroidAssnE2")
ODEsys <- addReaction(ODEsys, "abgbxsterxrcep_E", "abgbxrcep_E + ster_E", "k24*abgbxsterxrcep_E", "SteroidDssnE2")

#Antibody and steroid degradation in blood
ODEsys <- addReaction(ODEsys, "abgbxster_B", "ster_B", "k06*abgbxster_B", "AntibodyDegB1")
ODEsys <- addReaction(ODEsys, "abgb_B", "", "k06*abgb_B", "AntibodyDegB2")
ODEsys <- addReaction(ODEsys, "ster_B", "", "k07*ster_B", "SterDegB1")

#Steroid diffusion between compartments
ODEsys <- addReaction(ODEsys, "ster_B", "ster_C1", "k08*ster_B", "SterDiffBloodtoCyto")
ODEsys <- addReaction(ODEsys, "ster_C1", "ster_B", "k28*ster_C1", "SterDiffCytotoBlood")
ODEsys <- addReaction(ODEsys, "ster_B", "ster_C2", "k38*ster_B", "SterDiffBloodtoOffCyto2")
ODEsys <- addReaction(ODEsys, "ster_C2", "ster_B", "k48*ster_C2", "SterDiffOffCyto2toBlood")
ODEsys <- addReaction(ODEsys, "ster_C1", "ster_E", "k09*ster_C1", "SterDiffCytotoEndo")
ODEsys <- addReaction(ODEsys, "ster_E", "ster_C1", "k29*ster_E", "SterDiffEndotoCyto")

#Receptor and antibody degradation in endosome including receptor synthesis
ODEsys <- addReaction(ODEsys, "abgb_E", "", "k10*abgb_E", "FreeAntibodyDegE1")
ODEsys <- addReaction(ODEsys, "abgbxster_E", "ster_C1", "k10*abgbxster_E", "FreeAntibodyDegE2")
ODEsys <- addReaction(ODEsys, "abgbxrcep_E", "rcep_E", "k11*abgbxrcep_E", "RecepDegE1")
ODEsys <- addReaction(ODEsys, "abgbxsterxrcep_E", "ster_C1+rcep_E", "k11*abgbxsterxrcep_E", "RecepDegE2")

#Steroid binding to glucocorticoid receptor
ODEsys <- addReaction(ODEsys, "glcr_C1+ster_C1", "glcrxster_C1", "k12*glcr_C1*ster_C1", "GlucrecepSteroidBind")
ODEsys <- addReaction(ODEsys, "glcrxster_C1", "glcr_C1+ster_C1", "k13*glcrxster_C1", "GlucrecepSteroidDiss")
ODEsys <- addReaction(ODEsys, "glcr_C2+ster_C2", "glcrxster_C2", "k12*glcr_C2*ster_C2", "GlucrecepSteroidBindOffTarg")
ODEsys <- addReaction(ODEsys, "glcrxster_C2", "glcr_C2+ster_C2", "k13*glcrxster_C2", "GlucrecepSteroidDissOffTarg")

#Matches volumes to each chemical species by compartment suffix
fullvolumes <- numeric( length(colnames(ODEsys$smatrix)) )
names(fullvolumes) = colnames(ODEsys$smatrix)


for (i in c(1:length(fullvolumes)) )  {
  for (j in c(1:length(refvolumes)) )  {
    if ( grepl( names(refvolumes)[j], names(fullvolumes)[i] ) == TRUE )  
    { refvolumes[j] -> fullvolumes[i] }
  }
}
ODEsys$volumes <- fullvolumes

#Solving the system of ODEs
ODEnext <- odemodel(ODEsys, modelname = "ODEsys1", compile =  TRUE)  
ODEpredfunc <- Xs(ODEnext)

#Parameter assignment - run if you change parameters up top
paramnames <- getParameters(ODEpredfunc)
paramfull <- numeric(length(paramnames))
names(paramfull) <- paramnames

for (i in c(1:length(paramfull)) )  {
  for (j in c(1:length(initamts)) )  {
    if (names(initamts)[j] == names(paramfull)[i] )  
    { initamts[j] -> paramfull[i] }
  }
}

for (i in c(1:length(paramfull)) )  {
  for (j in c(1:length(rateconst)) )  {
    if (names(rateconst)[j] == names(paramfull)[i] )  
    { rateconst[j] -> paramfull[i] }
  }
}

#Now choose times and run prediction
times <- seq(0, 5e5, len = 1000)
prediction <- ODEpredfunc(times, paramfull)
plot(prediction)

#Generating plots and predictions with altered parameters
#Active receptor endocytosis is required
times <- seq(0, 5e5, len = 1000)

paramfull_pas <- paramfull
paramfull_pas["k05"] = 5.00E-04 #receptor internalization rate with or without ab (1st order)
paramfull_pas["k25"] = 5.00e-4*(SAplasma1/SAendo) #Receptor externalization rate (1st order)
paramfull_pas["k11"] = 2.26E-05     #receptor degradation in endosome; synthesis of new receptor is included
PassivePred <- ODEpredfunc(times,paramfull_pas)

paramfull_act <- paramfull
paramfull_act["k05"] = 5.5e-3 #receptor internalization rate with or without ab (1st order)
paramfull_act["k25"] = 5.5e-3*(SAplasma1/SAendo)  #Receptor externalization rate (1st order)
paramfull_act["k11"] = 1.65e-4  #receptor degradation in endosome; synthesis of new receptor is included
ActivePred <- ODEpredfunc(times,paramfull_act)

GraphSetA <- ActivePred[[1]]
GraphSetB <- PassivePred[[1]]
plot(GraphSetA[,1], GraphSetA[,"glcrxster_C1"]/paramfull["glcr_C1"],
     main="",
     ylab="Fraction GR Occupied",
     xlab="Time (s)",
     type="l",
     col="#006699" ,
     lty = 1)
lines(GraphSetA[,1], GraphSetA[,"glcrxster_C2"]/paramfull["glcr_C1"], 
      col="#006699" , 
      lty = 2)
lines(GraphSetB[,1], GraphSetB[,"glcrxster_C1"]/paramfull["glcr_C1"], 
      col="#CC0000", 
      lty = 1)
lines(GraphSetB[,1], GraphSetB[,"glcrxster_C2"]/paramfull["glcr_C1"], 
      col="#CC0000", 
      lty = 2)
par(xpd=TRUE)
legend(2e5,1,
       c("On-Target Cells, Fast Int.",
         "Off-Target Cells, Fast Int.", 
         "On-Target Cells, Slow Int.",
         "Off-Target Cells, Slow Int."),
       lty=c(1,2,1,2), 
       col=c("#006699","#006699","#CC0000","#CC0000")
)


#Endosomal release of steroid from fusion protein is not required

releaseparam = paramfull
releaseparam["k03"] = 3.98E+04     #glucocorticoid binder-steroid association rate
releaseparam["k23"] = releaseparam["k03"]    #glucocorticoid binder-steroid association rate in endosome
releaseparam["k04"] = 2.20E-05     #glucocorticoid binder-steroid dissociation rate

relfactors = c(1,5e3,5e4)
times <- seq(0, 5e5, len = 1000)
PredOut = list()

for (i in c(1:length(relfactors))) {
  releaseparam["k24"] = releaseparam["k04"]*relfactors[i]     #glucocorticoid binder-steroid dissociation rate in endosome
  PredOut[[i]] <- ODEpredfunc(times,releaseparam)[[1]]
}

GraphSetA <- PredOut[[1]]
GraphSetB <- PredOut[[2]]
GraphSetC <- PredOut[[3]]
plot(GraphSetA[,1], GraphSetA[,"glcrxster_C1"]/paramfull["glcr_C1"],
     main="",
     ylab="Fraction GR Occupied",
     xlab="Time (s)",
     type="l",
     col="#006699" ,
     lty = 1)
lines(GraphSetA[,1], GraphSetA[,"glcrxster_C2"]/paramfull["glcr_C1"], 
      col="#006699" , 
      lty = 2)
lines(GraphSetB[,1], GraphSetB[,"glcrxster_C1"]/paramfull["glcr_C1"], 
      col="#CC0000", 
      lty = 1)
lines(GraphSetB[,1], GraphSetB[,"glcrxster_C2"]/paramfull["glcr_C1"], 
      col="#CC0000", 
      lty = 2)
lines(GraphSetC[,1], GraphSetC[,"glcrxster_C1"]/paramfull["glcr_C1"], 
      col="#CC6600", 
      lty = 1)
lines(GraphSetC[,1], GraphSetC[,"glcrxster_C2"]/paramfull["glcr_C1"], 
      col="#CC6600", 
      lty = 2)
par(xpd=TRUE)
legend(2e5,1.1,
       c("On-Target, 1x Dissociation",
         "Off-Target, 1x Dissociation", 
         "On-Target, 5,000x Dissociation",
         "Off-Target, 5,000x Dissociation",
         "On-Target, 50,000x Dissociation",
         "Off-Target, 50,000x Dissociation"),
       lty=c(1,2,1,2,1,2), 
       col=c("#006699","#006699","#CC0000","#CC0000", "#CC6600", "#CC6600")
)


#Endogenous cortisol added to model overwhelms receptors (for supplement)

cortparam <- paramfull
cortparam["abgbxster_B"] = 0
cortparam["ster_B"] = 4.138e-8
cortparam["k12"] = 3.76e4 #3.98E+04     #glucocorticoid receptor on rate
cortparam["k13"] = 7.71e-4 #2.20E-05*10     #glucocorticoid receptor off rate
cortparam["k07"] = 0#3.08E-05     #steroid degradation in blood


cortisols = c(4.138e-8,6.483e-8,1.793e-8)
times <- seq(0, 4e4, len = 1000)
PredOut2 = list()

for (i in c(1:length(cortisols))) {
  cortparam["ster_B"] = cortisols[i]     #cortisol concentration
  PredOut2[[i]] <- ODEpredfunc(times,cortparam)[[1]]
}

GraphSetA <- PredOut2[[1]]
GraphSetB <- PredOut2[[2]]
GraphSetC <- PredOut2[[3]]
plot(GraphSetA[,1], GraphSetB[,"glcrxster_C1"]/paramfull["glcr_C1"],
     main="",
     ylab="Fraction GR Occupied",
     xlab="Time (s)",
     type="l",
     col="#CC0000" ,
     lty = 1)
lines(GraphSetB[,1], GraphSetA[,"glcrxster_C1"]/paramfull["glcr_C1"], 
      col="#006699", 
      lty = 1)
lines(GraphSetC[,1], GraphSetC[,"glcrxster_C1"]/paramfull["glcr_C1"], 
      col="#CC6600", 
      lty = 1)
par(xpd=TRUE)
legend(2e4,0.3,
       c("On-Target, Mean Cortisol",
         "On-Target, +1 SD Cortisol",
         "On-Target, -1 SD Cortisol"),
       lty=c(1,1,1), 
       col=c("#006699","#CC0000", "#CC6600")
)


#Antibody drug could target endogenous cortisol

targparam <- paramfull
targparam["abgbxster_B"] = 0
targparam["ster_B"] = 4e-9
targparam["abgb_B"] = 4e-9
targparam["k12"] = 3.76e4 #3.98E+04     #glucocorticoid receptor on rate
targparam["k13"] = 7.71e-4 #2.20E-05*10     #glucocorticoid receptor off rate
targparam["k07"] = 0#3.08E-05     #steroid degradation in blood

abgbamts = c(8e-9,4e-9,2e-9)
times <- seq(0, 2.5e5, len = 1000)
PredOut3 = list()

for (i in c(1:length(abgbamts))) {
  targparam["abgb_B"] = abgbamts[i]
  PredOut3[[i]] <- ODEpredfunc(times,targparam)[[1]]
}

GraphSetA <- PredOut3[[1]]
GraphSetB <- PredOut3[[2]]
GraphSetC <- PredOut3[[3]]
plot(GraphSetA[,1], GraphSetA[,"glcrxster_C1"]/paramfull["glcr_C1"],
     main="",
     ylab="Fraction GR Occupied",
     xlab="Time (s)",
     type="l",
     col="#CC0000" ,
     lty = 1)
lines(GraphSetA[,1], GraphSetA[,"glcrxster_C2"]/paramfull["glcr_C1"], 
      col="#CC0000" , 
      lty = 2)
lines(GraphSetB[,1], GraphSetB[,"glcrxster_C1"]/paramfull["glcr_C1"], 
      col="#006699", 
      lty = 1)
lines(GraphSetB[,1], GraphSetB[,"glcrxster_C2"]/paramfull["glcr_C1"], 
      col="#006699", 
      lty = 2)
lines(GraphSetC[,1], GraphSetC[,"glcrxster_C1"]/paramfull["glcr_C1"], 
      col="#CC6600", 
      lty = 1)
lines(GraphSetC[,1], GraphSetC[,"glcrxster_C2"]/paramfull["glcr_C1"], 
      col="#CC6600", 
      lty = 2)
par(xpd=TRUE)
legend(1.15e5,0.8,
       c("On-Target, 2 ab : 1 cortisol",
         "Off-Target, 2 ab : 1 cortisol", 
         "On-Target, 1 ab : 1 cortisol",
         "Off-Target, 1 ab : 1 cortisol",
         "On-Target, 1 ab : 2 cortisol",
         "Off-Target, 1 ab : 2 cortisol"),
       lty=c(1,2,1,2,1,2), 
       col=c("#CC0000","#CC0000","#006699","#006699", "#CC6600", "#CC6600")
)


#Steroid-antibody and steroid-receptor affinity can decrease and still allow targeting
affinityparam <- paramfull
affinityparam["k04"] = 2.20E-05     #glucocorticoid binder-steroid dissociation rate
affinityparam["k24"] = 2.20E-05     #glucocorticoid binder-steroid dissociation rate in endosome
affinityparam["k13"] = 2.20E-05     #glucocorticoid receptor off rate

dissfactor = c(1,2,4)
times <- seq(0, 5e5, len = 1000)
PredOut4 = list()

for (i in c(1:length(dissfactor))) {
  affinityparam["k04"] = paramfull["k04"]*dissfactor[i]
  affinityparam["k24"] = paramfull["k24"]*dissfactor[i]
  affinityparam["k13"] = paramfull["k13"]*dissfactor[i]
  PredOut4[[i]] <- ODEpredfunc(times,affinityparam)[[1]]
}

GraphSetA <- PredOut4[[1]]
GraphSetB <- PredOut4[[2]]
GraphSetC <- PredOut4[[3]]
plot(GraphSetA[,1], GraphSetA[,"glcrxster_C1"]/paramfull["glcr_C1"],
     main="",
     ylab="Fraction GR Occupied",
     xlab="Time (s)",
     type="l",
     col="#006699" ,
     lty = 1)
lines(GraphSetA[,1], GraphSetA[,"glcrxster_C2"]/paramfull["glcr_C1"], 
      col="#006699" , 
      lty = 2)
lines(GraphSetB[,1], GraphSetB[,"glcrxster_C1"]/paramfull["glcr_C1"], 
      col="#CC0000", 
      lty = 1)
lines(GraphSetB[,1], GraphSetB[,"glcrxster_C2"]/paramfull["glcr_C1"], 
      col="#CC0000", 
      lty = 2)
lines(GraphSetC[,1], GraphSetC[,"glcrxster_C1"]/paramfull["glcr_C1"], 
      col="#CC6600", 
      lty = 1)
lines(GraphSetC[,1], GraphSetC[,"glcrxster_C2"]/paramfull["glcr_C1"], 
      col="#CC6600", 
      lty = 2)
par(xpd=TRUE)
legend(2.5e5,1.1,
       c("On-Target, 1x Affinity",
         "Off-Target, 1x Affinity", 
         "On-Target, 0.5x Affinity",
         "Off-Target, 0.5x Affinity",
         "On-Target, 0.25x Affinity",
         "Off-Target, 0.25x Affinity"),
       lty=c(1,2,1,2,1,2), 
       col=c("#006699","#006699","#CC0000","#CC0000", "#CC6600", "#CC6600")
)


#Permeability coefficient effects are very strong

permeabparam <- paramfull
permeabparam["k08"] = 2.10E-9*SAplasma1*VolC1    #steroid diffusion from blood to cytoplasm
permeabparam["k28"] = 2.10E-9*SAplasma1*VolB     #steroid diffusion from cytoplasm to blood
permeabparam["k38"] = 2.10E-9*SAplasma2*VolC2    #steroid diffusion from blood to off-target cell cytoplasm
permeabparam["k48"] = 2.10E-9*SAplasma2*VolB     #steroid diffusion from off-target cell cytoplasm to blood
permeabparam["k09"] = 2.10E-9*SAendo*VolE      #steroid diffusion from cytoplasm to endosome
permeabparam["k29"] = 2.10E-9*SAendo*VolC1     #steroid diffusion from endosome to cytoplasm

permfactor = c(0.1,0.5,1)
times <- seq(0, 5e5, len = 1000)
PredOut5 = list()

for (i in c(1:length(permfactor))) {
  permeabparam["k08"] = paramfull["k08"]*permfactor[i]
  permeabparam["k28"] = paramfull["k28"]*permfactor[i]
  permeabparam["k38"] = paramfull["k38"]*permfactor[i]
  permeabparam["k48"] = paramfull["k48"]*permfactor[i]
  permeabparam["k09"] = paramfull["k09"]*permfactor[i]
  permeabparam["k29"] = paramfull["k29"]*permfactor[i]
  PredOut5[[i]] <- ODEpredfunc(times,permeabparam)[[1]]
}

GraphSetA <- PredOut5[[1]]
GraphSetB <- PredOut5[[2]]
GraphSetC <- PredOut5[[3]]
plot(GraphSetA[,1], GraphSetA[,"glcrxster_C1"]/paramfull["glcr_C1"],
     main="",
     ylab="Fraction GR Occupied",
     xlab="Time (s)",
     type="l",
     col="#CC6600",
     lty = 1)
lines(GraphSetA[,1], GraphSetA[,"glcrxster_C2"]/paramfull["glcr_C1"], 
      col="#CC6600", 
      lty = 2)
lines(GraphSetB[,1], GraphSetB[,"glcrxster_C1"]/paramfull["glcr_C1"], 
      col="#CC0000", 
      lty = 1)
lines(GraphSetB[,1], GraphSetB[,"glcrxster_C2"]/paramfull["glcr_C1"], 
      col="#CC0000", 
      lty = 2)
lines(GraphSetC[,1], GraphSetC[,"glcrxster_C1"]/paramfull["glcr_C1"], 
      col="#006699", 
      lty = 1)
lines(GraphSetC[,1], GraphSetC[,"glcrxster_C2"]/paramfull["glcr_C1"], 
      col="#006699", 
      lty = 2)
par(xpd=TRUE)
legend(2.5e5,1.35,
       c("On-Target, 1x Permeab.",
         "Off-Target, 1x Permeab.", 
         "On-Target, 0.5x Permeab.",
         "Off-Target, 0.5x Permeab.",
         "On-Target, 0.1x Permeab.",
         "Off-Target, 0.1x Permeab."),
       lty=c(1,2,1,2,1,2), 
       col=c("#006699","#006699","#CC0000","#CC0000", "#CC6600", "#CC6600")
)


