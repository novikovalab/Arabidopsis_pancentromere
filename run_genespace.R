# Based on Nikita Tikhomirov's script
wd='Pref' 
genomeIDs=c('NT1_v2', 'Pref')
#ncores
n=60


library(GENESPACE)
path2mcscanx <- '/opt/share/software/scs/appStore/stretchApps/synteny/MCScanX/va8443a9/bin'

# If rerunning from scratch, remove orthofinder results
# unlink(paste0(wd, 'orthofinder'), recursive=T)

# Check input data
gpar <- init_genespace(wd=wd, 
                       genomeIDs=genomeIDs,
                       ploidy=rep(1, 2),
                       path2mcscanx=path2mcscanx,
                       nCores=n
)

# Run synteny analysis
out <- run_genespace(gpar, overwrite = T)

# Plot
ripd <- plot_riparian(gsParam=out, useRegions=F, syntenyWeight = 1)
