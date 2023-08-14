
# Load the analysis script that has lots of functions in it
source("analyzeNNK.R");

# Setup for the specific library requires sequence of library positions ONLY and position IDs in order. 
setupNNKAnalysis(theSeq="VSTQLLIRSEQILNNAKIIIVTVKSIRIGPGQAFYYTFAQSSGGDLEITTHSIINMWQRAGQAMYTRDGGKDNNVNETFRPGGSDMRDNWRS",thePositions=c(219:224,236:250,267:282,323:337,385:396,417:443)); # define seq and pos CA

# Unsorted Cells
cat("example_unsort.data\n"); # print CA
A1 = loadFile("example_unsort.data",posBase="",numnt_min=-10,numnt_max=10); # load csv into dat subset CA

# Sorted Cells
cat("exmaple_sorted.data\n"); # print CA
A2 = loadFile("example_sorted.data",posBase="",numnt_min=-10,numnt_max=10); # load csv into dat subset CA

# Compute frequencies and propensities
Abs1sort.unsort.prop  = heatmapPlotAA.merge.propensity2(A2,A1,noplot=TRUE,posBase=""); # filters data and funds propensities CAA

# Using Abs1sort.unsort.prop you can do lots of things, including plotting data (more on that in a second)

#  function heatmapPlotAA.merge.propensity2 should be renamed at some point, but it filters data (zero counts, low counts, etc) and collects into 4 dataframes.
#  Abs1sort.unsort.prop is a list of 4 dataframes:  list(prop,probHH,probHH2,countsHH,countsHH2)

#  These datasets have values for every amino acid at every position

#  prop      = propensity enrichment     (probHH/probHH2)
#  probHH    = prob in sorted            (freq: counts of given AA,position/total at that position)
#  probHH2   = prob in unsorted(ref)     (freq: counts of given AA,position/total at that position)
#  countsHH  = counts in sorted          (counts: each AA at each position)
#  countsHH2 = counts in unsorted(ref)   (counts: each AA at each position)


# Plot AAs as follows:
#plotPieDistributionAA(c("A1", "A2"),pos=223,landscape=TRUE);

# Get most enriched mutations:
#getMostEnriched(Abs1sort.unsort.prop[[1]],cutoff=0.3);

# Create all tables and figures
#createTablesAndFigures(c("Abs1sort.unsort.prop"),"Abs1sort_analysis"); # this did not work for me CA
#createTablesAndFigures(Abs1sort.unsort.prop,"Abs1sort_analysis"); # this did work for me CA

# Condense mutations across sorts if you have them
#b = condenseTables(tables=c("V033I1.unsort.prop","V033I1r1.unsort.prop","V033I1r2.unsort.prop"),minNumTables=3,tables.index=1,lower_cutoff=0);


cat("Done Load\n");