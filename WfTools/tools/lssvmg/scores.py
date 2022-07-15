import numpy as np
import csv

targets = []
with open("data/targets.csv") as fh:
    reader = csv.reader(fh)
    for row in reader:
        targets.append(row)

targets = np.array(targets, np.float)
train_targets = targets[:644]
test_targets = targets[644:]

results = []
with open("result.csv") as fh:
    reader = csv.reader(fh)
    for row in reader:
        results.append(row)
results = np.array(results, np.float)
train_results = results[:644]
test_results = results[644:]

train_score = np.sqrt(np.mean((train_results-train_targets)**2))
test_score  = np.sqrt(np.mean((test_results-test_targets)**2))

print ("train rmse: %.5f, test rmse: %.5f" %(train_score, test_score))
