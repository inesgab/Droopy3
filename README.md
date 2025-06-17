# Experimental Data analysis

## Dropsignal
1. From millidrop computer, export “EXPXXXXXXX” data folder to monod computer
2. access monod from local computer, add `template.csv` and `protocol.json` to new EXPXXXXXXX folders (`$ cp .../template.csv ...`)
3. `$ ma` and `$ cd 1machine` (or 2machine)
4. `$ drps` (dropsignal) and open dropsignal on navigator
5. **Edit protocol** to check template and drop coalescence
6. **Run analysis**
7. **Open results** and eventually manage coalescence errors by checking single run scatter and editing protocol
8. Once the data is ok, import on local computer (`$ scp -r …`) and analyse with droopy

## Droopy
0. Launch app.py
1. **Load parameters**: Enter the path of the folder containing the EXPXXXXX folders
2. **Load parameters**: Parameter fit in (a, b, electronic precision)
3. **Check the data**: Check the data by label button
4. **Check the data**: : If good, start fitting, and wait until treatment is over (takes up to an hour approx.)
5. **Check outliers**: *Start check outliers* : Manually check the generated curves. If one seems off, sort it out of the batch. Continue until all are done
6. **reFit non outliers**: start fitting the purified data: also a bit long, and recheck the fits. Manually add the droplet number of an eventual wrong fit in `EXPXXXX/analysis/outliersByEye.csv`


## Plotting: 

If antibio: 
`plotLagTimeAntibio.py`, after renaming the main folder (it needs to contain “antib”)
Ex: GFP_**antib**_nal/EXP20250611_1103/…
