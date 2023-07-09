import pickle

merges, passives = {}, {}
with open("Possible Combos.txt") as combos:
    key, values = '', []
    for i in combos.readlines():
        if i.startswith("="):
            key = i.split('=')[-1].strip()
            print(key)
        elif len(i.strip()) and not i.startswith("("):
            values = list(map(lambda x : x.strip().upper(), i.split(',')))
        else:
            merges[key] = values
            key, values = '', []

with open("passive edges.txt") as pas:
    key, values = '', []
    for i in pas.readlines():
        if i.startswith("="):
            key = i.split('=')[-1].strip()
            print(key)
        elif len(i.strip()) and not i.startswith("("):
            values = list(map(lambda x : x.strip().upper(), i.split(',')))
        else:
            passives[key] = values
            key, values = '', []

print(len(merges))
with open("merged_edges", 'wb') as pickle_edge:
    pickle.dump(merges, pickle_edge)
with open("passive_edges", 'wb') as passive_edge:
    pickle.dump(passives, passive_edge)