import xml.etree.ElementTree as ET
import json
tree = ET.parse('unimod.xml')
root = tree.getroot()
print(len(root[1]))
tag = '{http://www.unimod.org/xmlns/schema/unimod_2}'
unimod_db = dict()
for mod in root[1]:
    name = mod.attrib['full_name']
    mass = float(mod.findall('{http://www.unimod.org/xmlns/schema/unimod_2}delta')[0].attrib['mono_mass'])
    for spec in mod.findall('{http://www.unimod.org/xmlns/schema/unimod_2}specificity'):
        site, cla = spec.attrib['site'], spec.attrib['classification']
        unimod_db.setdefault(site, [])
        unimod_db[site].append([mass, name, cla])
for k, v in unimod_db.items():
    v.sort()
json_object = json.dumps(unimod_db, indent=4)
with open('db_unimod', 'w') as f:
    f.write(json_object)