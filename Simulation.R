source('group_pmom_new.R')
source('SSS-group new.R')
source('functions.R')

############################################################
################ Design 1 (Baseline Design) ################ 
############################################################

########## g = 50, t = 3 ##########

#######~~~~~~~~~~~~~~~~~~~~~~~~~~#######
##~~~~~~~~~~~~~# case 1 #~~~~~~~~~~~~~##
#######~~~~~~~~~~~~~~~~~~~~~~~~~~#######
library(MASS)
#***************************** setting 1 #***************************** 
set.seed(20211212)
data.1.1.1 <- gen_data(n = 100, g = 50, g.size = 4, t = 3, setting = 1, case = 1
                       #, seeds = 20211212
                       )

########### SSS-GROUP ############
SSS.dt.1.1.1 = data.1.1.1
g = 50
g.size = 4
colnames(SSS.dt.1.1.1$X.sd) = rep(1:g, each=g.size)

#set.seed(484)
fit.sss.1.1.1=SSS(X=SSS.dt.1.1.1$X.sd, y=SSS.dt.1.1.1$Y.sd
                  , ind_fun=group_pmom_log_post, delta=0.01
                  ,r=2, alpha_1=0.01, alpha_2=0.01
                  ,N=3,C0=1,verbose=TRUE)

sss.rt.1.1.1 = result.evl.SSS(fit.sss.1.1.1, SSS.dt.1.1.1, g=g, t=3)
sss.rt.1.1.1

########## group lasso ###########
print("group lasso")
library(grpreg)
gplasso.rt.1.1.1 = result.evl.gplasso(data=data.1.1.1, g, g.size, t=3
                                      , penalty="grLasso")
gplasso.rt.1.1.1

########## group SCAD ###########
print("group SCAD")
gpSCAD.rt.1.1.1 = result.evl.gplasso(data=data.1.1.1, g, g.size, t=3
                                     , penalty="grSCAD")
gpSCAD.rt.1.1.1

########## group MCP ###########
print("group MCP")
gpMCP.rt.1.1.1 = result.evl.gplasso(data=data.1.1.1, g, g.size, t=3
                                     , penalty="grMCP")
gpMCP.rt.1.1.1

########## group gel ###########
print("group gel")
gpgel.rt.1.1.1 = result.evl.gplasso(data=data.1.1.1, g, g.size, t=3
                                    , penalty="gel")
gpgel.rt.1.1.1

### spike and slab group lasso ###
library(devtools)
#install_github(repo = "jantonelli111/SSGL")
library(SSGL)
print("spike and slab group lasso")
SSGL.rt.1.1.1 = result.evl.SSGL(data=data.1.1.1, g=g, g.size, t=3
                                , lambda0seq=seq(21, 100, by=1))
SSGL.rt.1.1.1




#***************************** setting 2 #***************************** 
set.seed(151501212)
data.1.1.2 <- gen_data(n = 100, g = 50, g.size = 4, t = 3, setting = 2, case = 1
                       #, seeds = 15151212
                       )

########### SSS-GROUP ############
SSS.dt.1.1.2 = data.1.1.2
g = 50
g.size = 4
colnames(SSS.dt.1.1.2$X.sd) = rep(1:g, each=g.size)

#set.seed(2390)
fit.sss.1.1.2=SSS(X=SSS.dt.1.1.2$X.sd, y=SSS.dt.1.1.2$Y.sd
                  , ind_fun=group_pmom_log_post, delta=0.01
                  ,r=2, alpha_1=0.01, alpha_2=0.01
                  ,N=3,C0=1,verbose=TRUE)

sss.rt.1.1.2 = result.evl.SSS(fit.sss.1.1.2, SSS.dt.1.1.2, g=g, t=3)
sss.rt.1.1.2

########## group lasso ###########
print("group lasso")
gplasso.rt.1.1.2 = result.evl.gplasso(data=data.1.1.2, g, g.size, t=3
                                      , penalty="grLasso")
gplasso.rt.1.1.2

########## group SCAD ###########
print("group SCAD")
gpSCAD.rt.1.1.2 = result.evl.gplasso(data=data.1.1.2, g, g.size, t=3
                                     , penalty="grSCAD")
gpSCAD.rt.1.1.2

########## group MCP ###########
print("group MCP")
gpMCP.rt.1.1.2 = result.evl.gplasso(data=data.1.1.2, g, g.size, t=3
                                    , penalty="grMCP")
gpMCP.rt.1.1.2

########## group gel ###########
print("group gel")
gpgel.rt.1.1.2 = result.evl.gplasso(data=data.1.1.2, g, g.size, t=3
                                    , penalty="gel")
gpgel.rt.1.1.2

### spike and slab group lasso ###
print("spike and slab group lasso")
SSGL.rt.1.1.2 = result.evl.SSGL(data=data.1.1.2, g=g, g.size, t=3
                                , lambda0seq=seq(1, 100, by=2))
SSGL.rt.1.1.2



#***************************** setting 3 #***************************** 

set.seed(20241212)
data.1.1.3 <- gen_data(n = 100, g = 50, g.size = 4, t = 3, setting = 3, case = 1
                       #, seeds = 20221212
                       )

########### SSS-GROUP ############
SSS.dt.1.1.3 = data.1.1.3
g = 50
g.size = 4
colnames(SSS.dt.1.1.3$X.sd) = rep(1:g, each=g.size)

#set.seed(394)
fit.sss.1.1.3=SSS(X=SSS.dt.1.1.3$X.sd, y=SSS.dt.1.1.3$Y.sd
                  , ind_fun=group_pmom_log_post, delta=0.01
                  ,r=2, alpha_1=0.01, alpha_2=0.01
                  ,N=3,C0=1,verbose=TRUE)

sss.rt.1.1.3 = result.evl.SSS(fit.sss.1.1.3, SSS.dt.1.1.3, g=g, t=3)
sss.rt.1.1.3

########## group lasso ###########
print("group lasso")
gplasso.rt.1.1.3 = result.evl.gplasso(data=data.1.1.3, g=g, g.size, t=3
                                      , penalty="grLasso")
gplasso.rt.1.1.3

########## group SCAD ###########
print("group SCAD")
gpSCAD.rt.1.1.3 = result.evl.gplasso(data=data.1.1.3, g, g.size, t=3
                                     , penalty="grSCAD")
gpSCAD.rt.1.1.3

########## group MCP ###########
print("group MCP")
gpMCP.rt.1.1.3 = result.evl.gplasso(data=data.1.1.3, g, g.size, t=3
                                    , penalty="grMCP")
gpMCP.rt.1.1.3

########## group gel ###########
print("group gel")
gpgel.rt.1.1.3 = result.evl.gplasso(data=data.1.1.3, g, g.size, t=3
                                    , penalty="gel")
gpgel.rt.1.1.3

### spike and slab group lasso ###
print("spike and slab group lasso")
SSGL.rt.1.1.3 = result.evl.SSGL(data=data.1.1.3, g=g, g.size, t=3
                                , lambda0seq=seq(1, 100, by=2))
SSGL.rt.1.1.3



#***************************** setting 4 #***************************** 
set.seed(191212)
data.1.1.4 <- gen_data(n = 100, g = 50, g.size = 4, t = 3, setting = 4, case = 1
                       #, seeds = 19991212
                       )

########### SSS-GROUP ############
SSS.dt.1.1.4 = data.1.1.4
g = 50
g.size = 4
colnames(SSS.dt.1.1.4$X.sd) = rep(1:g, each=g.size)

#set.seed(8659)
fit.sss.1.1.4=SSS(X=SSS.dt.1.1.4$X.sd, y=SSS.dt.1.1.4$Y.sd
                  , ind_fun=group_pmom_log_post, delta=0.02
                  ,r=2, alpha_1=0.01, alpha_2=0.01
                  ,N=3,C0=1,verbose=TRUE)

sss.rt.1.1.4 = result.evl.SSS(fit.sss.1.1.4, SSS.dt.1.1.4, g=g, t=3)
sss.rt.1.1.4

########## group lasso ###########
print("group lasso")
gplasso.rt.1.1.4 = result.evl.gplasso(data=data.1.1.4, g=g, g.size, t=3
                                      , penalty="grLasso")
gplasso.rt.1.1.4

########## group SCAD ###########
print("group SCAD")
gpSCAD.rt.1.1.4 = result.evl.gplasso(data=data.1.1.4, g, g.size, t=3
                                     , penalty="grSCAD")
gpSCAD.rt.1.1.4

########## group MCP ###########
print("group MCP")
gpMCP.rt.1.1.4 = result.evl.gplasso(data=data.1.1.4, g, g.size, t=3
                                    , penalty="grMCP")
gpMCP.rt.1.1.4

########## group gel ###########
print("group gel")
gpgel.rt.1.1.4 = result.evl.gplasso(data=data.1.1.4, g, g.size, t=3
                                    , penalty="gel")
gpgel.rt.1.1.4

### spike and slab group lasso ###
print("spike and slab group lasso")
SSGL.rt.1.1.4 = result.evl.SSGL(data=data.1.1.4, g=g, g.size, t=3
                                , lambda0seq=seq(1, 100, by=2))

SSGL.rt.1.1.4