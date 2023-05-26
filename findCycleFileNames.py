# Goal: Reading the CSV File, Find ecDNA in the CSV File => list of cycle files name

# Imports
import pandas as pd

# Variables
filePath = "~/Downloads/results/"
fileName = "aggregated_results.csv"
suffixOfCycleFileName = "_annotated_cycles.txt"
cycleFileNames = []

# Filters rows that have the classification 'ecDNA'
dataTable = pd.read_csv(filePath + fileName)
dataTable = dataTable[dataTable['Classification'].notna()]
dataTable = dataTable[dataTable['Classification'].str.contains('ecDNA')]

dataTable = dataTable['Feature ID']

# Makes a list of cycle file names
for id in dataTable:
    index = id.index('amplicon')
    cycleFileNames.append(id[:index+9] + suffixOfCycleFileName)