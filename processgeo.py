import pandas as pd
import GEOparse
import os
import json

studies = [
#"GSE131320"
#"GSE134068"
#"GSE120299"
#"GSE101207"
]




def make_meta_file(studies: list, organism):
    for study in studies:     
        gse = GEOparse.get_GEO(geo = study, destdir = f'./static/data/{organism}/{study}', silent=True)

        if organism == 'human':
            geo_organism = 'Homo sapiens'
        
        samples = gse.metadata['sample_id']
        meta_dict = {}
        meta_write_path = f'./static/data/{organism}/{study}/{study}_Metadata.json'
        for samp in samples:
            source = gse.gsms[samp].metadata['source_name_ch1'][0]
            orga = gse.gsms[samp].metadata['organism_ch1'][0]
            conditions =  gse.gsms[samp].metadata['characteristics_ch1']
            print(samp, source, ":",conditions)
            condition = ""
            for cond in conditions:
                key, value = cond.split(":")[0].strip(), cond.split(":")[1].strip()
                if 'condition' in key or 'state' in key or 'treatment' in key :
                    condition += value
                elif 'time' in key:
                    condition =  condition + " " + value
            if source not in meta_dict and orga == geo_organism:
                meta_dict[source] = {condition: [samp]}
            elif orga == geo_organism:
                if condition in meta_dict[source]:
                    meta_dict[source][condition].append(samp)
                else:
                    meta_dict[source][condition] = [samp]
        with open(meta_write_path, 'w') as f:
            json.dump(meta_dict, f)


#make_meta_file(studies, 'mouse')






def make_meta_file2(studies: list, organism):
    for study in studies:     
        gse = GEOparse.get_GEO(geo = study, destdir = f'./static/data/{organism}/{study}', silent=True)

        if organism == 'human':
            geo_organism = 'Homo sapiens'
        
        samples = gse.metadata['sample_id']
        meta_dict = {}
        meta_write_path = f'./static/data/{organism}/{study}/{study}_Metadata.json'
        for samp in samples:
            source = gse.gsms[samp].metadata['source_name_ch1'][0]
            orga = gse.gsms[samp].metadata['organism_ch1'][0]
            conditions =  gse.gsms[samp].metadata['characteristics_ch1']
            print(samp, source, ":",conditions)
            condition = ""
            for cond in conditions:
                key, value = cond.split(":")[0].strip(), cond.split(":")[1].strip()
                if 'quality' in key or 'state' in key:
                    condition += value
                elif 'time' in key:
                    condition =  condition + " " + value
            if source not in meta_dict and orga == geo_organism:
                meta_dict[source] = {condition: [samp]}
            elif orga == geo_organism:
                if condition in meta_dict[source]:
                    meta_dict[source][condition].append(samp)
                else:
                    meta_dict[source][condition] = [samp]
        with open(meta_write_path, 'w') as f:
            json.dump(meta_dict, f)


#make_meta_file2(studies, 'human')




def get_genes(organism):
    
    for study in os.listdir('static/data/' + organism):
        if study[0:3] == 'GSE':
            path = 'static/data/%s/%s/%s_Expression.txt' %  (organism, study, study)
            table = pd.read_csv(path, delimiter='\t', index_col=0)
            with open('static/data/%s/%s/genes.json' %  (organism, study), 'w') as f:
                json.dump(list(table.index.values), f, ensure_ascii=False, indent=4)

get_genes('mouse')
get_genes('human')








