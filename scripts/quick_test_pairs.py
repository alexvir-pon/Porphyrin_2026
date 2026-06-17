import gepia

pairs = [
    ("H2AFZ", "HMBS"),
    ("ATAD2", "HMBS"),
    ("MKI67", "HMBS"),
    ("PCNA", "HMBS"),
]

dataset = ["GBM_Tumor", "BLCA_Tumor", "ESCA_Tumor", "HNSC_Tumor", "OV_Tumor", "STAD_Tumor"]

for g1, g2 in pairs:
    print(f"Running {g1} vs {g2}")
    c = gepia.correlation()
    c.setParams({
        "dataset": dataset,
        "methodoption": "spearman",
        "signature1": [g1],
        "signature2": [g2],
        "signature1_norm": "",
        "signature2_norm": "",
    })
    c.query()
