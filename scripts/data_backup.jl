# use Robocopy to back up all data to the shock drive
backuploc = raw"S:\users\shocktube\Alex\data\vortex_ring\thesis"

run(`robocopy $(datadir()) $(backuploc) /MIR /E /R:3 /W:3 /NP /NDL /NJH /NJS /TEE /LOG:$(datadir("backup.log"))`)
