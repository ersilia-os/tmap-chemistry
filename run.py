from src.tmap import Tmap

infile = "data.csv"
outfolder = "output"

tm = Tmap(infile, outfolder)
tm.run()
