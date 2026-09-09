# Freeze software/session information --------------------------------------
dir.create("results", showWarnings = FALSE)
sink(file.path("results", "sessionInfo.txt"))
print(sessionInfo())
sink()
message("Wrote results/sessionInfo.txt")
