import requests
API_KEY = '8612d3021dfae8ed745525754ada890b'

# test something random
params = {"name": ["HIP 12345", "HIP 56789", "HIP 113044"]}
entry = requests.get("https://hypatiacatalog.com/hypatia/api/v2/star", auth=(API_KEY, "api_token"), params=params)
print('retrieved from catalogue:', entry.json())

# retrieve star composition
star = '2MASS 00380401+3442416'
els = ['mg', 'si', 'fe', 'ca', 'al']
params = {"name": [star] * len(els), "element": els, "solarnorm": ["lod09"] * len(els)}
entry = requests.get("https://hypatiacatalog.com/hypatia/api/v2/composition", auth=(API_KEY, "api_token"),
                     params=params)
print('retrieved from catalogue:', entry.json())


# params = {"hip": "98355", "element": "ca", "solarnorm": "asp05"}
# entry = requests.get("https://hypatiacatalog.com/hypatia/api/v2/composition",auth=(API_KEY,"api_token"), params=params)

# params = {"solarnorm": "lod09", "filter1_1": "mg", }
# entry = requests.get("https://hypatiacatalog.com/hypatia/api/v2/data",auth=(API_KEY,"api_token"), params=params)
#
# print(entry)
# print(entry.json())
