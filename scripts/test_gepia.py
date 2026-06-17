#!/usr/bin/env python3

import gepia

c = gepia.correlation()

c.setParams({
  "dataset": ["BRCA_Tumor", "BLCA_Tumor"],
    "methodoption": "pearson",
    "signature1": ["HMBS"],
    "signature1_norm": "",
    "signature2": ["UROS"],
    "signature2_norm": ""
})

print("Parameters:")
c.showParams()

print("\nRunning query...")
result = c.query()

print("\nResult:")
print(result)
