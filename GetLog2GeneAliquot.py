import requests
import pandas as pd
import json
import os 
from DatabaseGenericInteraction import * 
import math 


# Funzione per scaricare gli Ensembl Id di tutti i singoli geni, questo è necesario poichè la query utilizzata sotto ritorna solo i nomi dei geni
# Funzione necessaria perchè non è certo che un determinato gene indicato sotto sia già stato salvato all'interno del database e quindi non si può fare una traduzione attraverso il database direttamente 
def getEnsemblId(gene_name):
    # Spostamento nella cartella di Load e apertura del file usato come cache 
    os.chdir("PATH\\load")
    with open("Ensembl_Gene_Translation.txt",'a',encoding='utf-8') as f:
        # Costruisco la request
        server = "https://rest.ensembl.org"
        ext = "/xrefs/symbol/homo_sapiens/"+ gene_name +"?"
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
        if not r.ok:
                r.raise_for_status()
        decoded = r.json()
        if len(decoded) == 0:
            return None, None
        f.write(gene_name + ":" + decoded[0]['id']+'\n')
    return gene_name, decoded[0]['id']


def query_pdc(pdc_study_id, data_type):

    # Graphql query passata per l'ottenimento delle informazioni
    query = '''
{ 
    quantDataMatrix(
    pdc_study_id: "''' + pdc_study_id +'''" data_type: "''' + data_type + '''" acceptDUA: true
    )
}'''
    # PDC API url
    url = 'https://pdc.cancer.gov/graphql'

    # Invia una richiesta POST attraverso una graphql query
    pdc_response = requests.post(url, json={'query': query})

    if pdc_response.ok:
        return pdc_response.json()
    else:
        return pdc_response.raise_for_status()

# Per ogni studio possiamo ottenere i log2_ratio e gli unshared_log2_ratio
def getLog2RatioInfo(program_pdc, project_id, cursor, conn):
    data_type = [ 'log2_ratio', 'unshared_log2_ratio' ]
    geni_trans = dict()
    for type in data_type:
        decoded = query_pdc(program_pdc, type)

        if 'errors' in decoded:
            continue
        matrix = decoded['data']['quantDataMatrix']
        # pd.DataFrame crea una struttura di 2 dimensioni, avendo così una sorta di tabella da poter navigare nell'estrazione dei dati 
        ga = pd.DataFrame(matrix[1:], columns=matrix[0])
        # lista di liste di gene_name
        gene_name = ga.iloc[:,:1].values

        # Utilizziamo una sorta di cache per la traduzione gene name : Ensembl gene id perchè non tutti i geni potrebbero essere stati campionati
        try:
            f = ( open("PATH\\load\\Ensembl_Gene_Translation.txt", "r") ).read()
            rows = f.split('\n')
            for row in rows:
                r = row.split(':')
                if len(r)>1:
                    geni_trans[r[0]] = r[1]
        except:
            pass
    
        # Verifico che esiste il file di cache con le traduzioni dei gene_name
        if len(geni_trans) == 0:
            os.chdir("PATH\load")
            with open("Ensembl_Gene_Translation.txt",'w',encoding='utf-8') as f:
                # Effettuo la conversione dei gene name in Ensembl Id che è lo standard usato in tutto il progetto 
                server = "https://rest.ensembl.org"
                for nome in gene_name:
                    ext = "/xrefs/symbol/homo_sapiens/"+nome[0]+"?"
                    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
                    if not r.ok:
                        r.raise_for_status()
                    decoded = r.json()
                    geni_trans[nome[0]] = decoded[0]['id']
                    f.write(nome[0] + ":" + decoded[0]['id']+'\n')

        # Controllo che ho già i corrispondenti Ensembl Id dei geni ottenuti
        # Se manca qualche gene allora lo scarico, salvo nel file usato come cache e lo aggiungo al dizionario usato durante la corrente esecuzione 
        gen_trans_key = list(geni_trans.keys())
        for gen in gene_name:
            if ':' in gen[0]:
                geni_trans[gen[0]] = None
            elif  gen[0] not in gen_trans_key:
                print("Nuovo Gene Trovato " + gen[0])
                gen_name , gene_id = getEnsemblId(gen[0])
                if gen_name != None:
                    geni_trans[gen_name] = gene_id
                
        if not 'Gene/Aliquot' in ga:
            continue

        ga = ga.set_index('Gene/Aliquot')

        # ottengo le colonne 
        oldnames = list(ga.columns)
        newnames = [ x.split(':')[1] for x in oldnames ]

        ga.rename(columns=dict(zip(oldnames, newnames)), inplace=True)
        ga = ga.sort_index(axis=1)

        aliquot = list(ga)
        # numero colonne -> len(aliquot)
        # numero righe -> len(ga)
        gen_key = list(geni_trans.keys())
        # Bisogna recuperaer i gene_id e gene_name 
        for x in range(len(ga)):
            for aliq in aliquot:

                # Bisogna controllare che geni_trans[gene_name[x][0]] esista perchè qualche volta viene restituito il valore None dalla traduzione
                if gene_name[x][0] not in gen_key or geni_trans[gene_name[x][0]] == None:
                    continue

                # Effettuiamo un controllo se tutti i valori che dovranno essere usati come primary key sono presenti nel database, in modo da evitare errori nell'esecuzione delle insert o update 
                if checkExistGene(geni_trans[gene_name[x][0]], cursor) and checkExistProject(project_id, cursor)  and checkExistAliquote(aliq, cursor):

                    # Effettuiamo il check se i valori sono "NaN"
                    try:
                        # Se è NaN lo sostituisco con 0
                        if math.isnan( float( (ga.iloc[x])[aliq] ) ) :
                                log_val = 0
                        else:
                                # Se NON è NaN lo converto in stringa e lo utilizzo
                                log_val = str((ga.iloc[x])[aliq])
                    except:
                        log_val = 0

                    # controllo quale colonna è da aggiornare nella tabella in base al valore di type 
                    if type == 'log2_ratio':                        
                        if checkExistProtein_PDC(geni_trans[gene_name[x][0]], project_id, aliq,cursor)  :
                            if getLog2Ratio(geni_trans[gene_name[x][0]], project_id, aliq,cursor) == 0 or  getLog2Ratio(geni_trans[gene_name[x][0]], project_id, aliq,cursor) == None:
                                query = "UPDATE  public.protein_PDC SET log2_ratio = {} WHERE gene_id = '{}' and aliquot = '{}' and project_id = '{}' ;".format(log_val, geni_trans[gene_name[x][0]], aliq, project_id )
                                cursor.execute(query)
                                conn.commit()
                            else:
                                query = None
                        else:
                            query = "INSERT INTO public.protein_PDC (gene_id, log2_ratio, aliquot, project_id) VALUES ('{}', {}, '{}', '{}') ON CONFLICT (gene_id,project_id,aliquot) DO NOTHING;".format(geni_trans[gene_name[x][0]], log_val, aliq, project_id )
                            cursor.execute(query)
                            conn.commit()
                    else:
                        if checkExistProtein_PDC(geni_trans[gene_name[x][0]], project_id, aliq,cursor):
                            if getUnsharedLog2Ratio(geni_trans[gene_name[x][0]], project_id, aliq,cursor) == 0 or getUnsharedLog2Ratio(geni_trans[gene_name[x][0]], project_id, aliq,cursor) == None:
                                query = "UPDATE  public.protein_PDC SET unshared_log2_ratio = {} WHERE gene_id = '{}' and aliquot = '{}' and project_id = '{}' ;".format(log_val, geni_trans[gene_name[x][0]], aliq, project_id )
                                cursor.execute(query)
                                conn.commit()
                        else:
                            query = "INSERT INTO public.protein_PDC (gene_id, unshared_log2_ratio, aliquot, project_id) VALUES ('{}', '{}', '{}', '{}') ON CONFLICT (gene_id,project_id,aliquot) DO NOTHING;".format(geni_trans[gene_name[x][0]], log_val, aliq, project_id )
                            cursor.execute(query)
                            conn.commit()
                else:
                    pass

                        

