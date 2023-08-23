import requests

API_KEY = '8612d3021dfae8ed745525754ada890b'

# test something random from https://www.hypatiacatalog.com/api
params = {"name": ["HIP 12345", "HIP 56789", "HIP 113044"]}
entry = requests.get("https://hypatiacatalog.com/hypatia/api/v2/star", auth=(API_KEY, "api_token"), params=params)
print('retrieved from catalogue:\n', entry.json())
# works fine

# retrieve star composition
star = '2MASS 00380401+3442416'
els = ['mg', 'si', 'fe', 'ca', 'al']
params = {"name": [star] * len(els), "element": els, "solarnorm": ["lod09"] * len(els)}
entry = requests.get("https://hypatiacatalog.com/hypatia/api/v2/composition", auth=(API_KEY, "api_token"),
                     params=params)
try:
    print('retrieved from catalogue:\n', entry.json())
except Exception as e:
    print(e)
    # simplejson.errors.JSONDecodeError: Expecting value: line 1 column 1 (char 0)

try:
    print(entry)
except Exception as e:
    print(e)
    # <Response [500]>


# # try from official instructions
# print('trying official demo')
# params = {"name": "HIP 98355", "element": "ca", "solarnorm": "asp05"}
# entry = requests.get("https://hypatiacatalog.com/hypatia/api/v2/composition", auth=(API_KEY, "api_token"), params=params)
# print(entry.json())

