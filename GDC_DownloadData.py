from DatabaseManager import *
from DatabaseGenericInteraction import *
import params 
import sys
import requests
import psycopg2
import json
import csv
import pandas as pd
from io import StringIO

# Funzione per scaricare e processare i dati di espressione da GDC
'''
Gli passo params in modo da potervi accedere 
'''
def download_and_process_expression_data(param):

    try:

        # Creo connessione con il database 
        cursor, conn = databaseConnection()

        if cursor == None or conn == None:
            print("Impossibile connettersi al database")
            return 
        else: 
            print("Connessione riuscita")

        # URL dell'API GDC per scaricare i file corrispondenti alle analisi
        gdc_api_url = "https://api.gdc.cancer.gov/files"

        '''
 # Filtro riguardante il sito primario della malattia che vogliamo analizzare
                {
                    "op": "=",
                    "content": {
                        "field": "cases.primary_site",
                        "value": "testis"
                    }
                },
'''

        # Download dei dati
        filters = {
            "op": "and",
            "content": [                

                # Filtro riguardante il tipo di dati che vogliamo analizzare
                {
                    "op": "in",
                    "content": {
                        "field": "data_type",
                        #"value": ["Gene Expression Quantification"]   # clsProtein Expression Quantification
                        "value": ["Protein Expression Quantification"]   # clsProtein Expression Quantification
                    }
                },

                # Filtro riguardante l'accesso dei dati
                {
                    "op": "=",
                    "content": {
                        "field": "access",
                        "value": "open"
                    }
                },

                # Filtro riguardante il formato del file su cui è riportata l'analisi
                {
                    "op": "=",
                    "content": {
                        "field": "data_format",
                        "value": "TSV"
                    }
                }
            ]
        }

        # Parametri usati per efettuare la get 
        params = {
            "filters": json.dumps(filters),
            # Altri campi per ottenere maggiori informazioni sui file scaricati 
            "fields": "file_name,file_size,created_datetime,updated_datetime,data_type,experimental_strategy,data_category,cases.project.project_id,cases.case_id,cases.submitter_id,associated_entities.entity_submitter_id",
            "format": "JSON",
            "size": "12000000",  # Numero massimo di file da scaricare per richiesta
            "pretty": "true"  # pretty indica che la response viene formattata con spazi aggiuntivi per migliorare la leggibilità
        }
        
        try:
            response = requests.get(gdc_api_url, params=params)
        except:
            print("Errore connessione server")
            return 

        # Elaborazione dei dati e inserimento nel database
        for file_info in json.loads(response.content.decode("utf-8"))["data"]["hits"]:
            if "id" in file_info:
                file_id = file_info["id"]
            else:
                break
            print(file_id)

            # Check vari per evitare al massimo errori durante l'esecuzione del codice 
            if "cases" not in file_info:
                break 
            elif "submitter_id" not in file_info["cases"][0]:
                break
            elif "case_id" not in file_info["cases"][0]:
                break

            if "project" not in file_info["cases"][0]:
                break 
            elif "project_id" not in file_info["cases"][0]["project"]:
                break
            

            # Gestione dei progetti 
            project_id = file_info["cases"][0]["project"]["project_id"]
            # Controllo se il progetto esiste già nel database 
            if not checkExistProject(project_id, cursor): 
                # searchProject restituisce None in caso di errore
                project_name = searchProject(project_id)

                # Se ho recuperato informazioni sul progetto allora faccio l'insert nel database 
                if project_name != None:
                    insertNewProject(project_id, project_name, cursor, conn)
                    conn.commit()
                else:
                    break

            # Gestione dei case
            if not checkExistCase(file_info["cases"][0]["submitter_id"], cursor):
                # searchCase restituisce None in caso di errore 
                case_data = searchCase(file_info["cases"][0]["case_id"])
                if case_data != None:
                    insertNewCase(case_data, project_id, cursor, conn)
                    conn.commit()
                else:
                    break


            type_id = None

            # Gestione dei file 
            if not checkExistFile(file_id, cursor):
                type_id = insertNewAnalysis(file_id, file_info, project_id, cursor, conn)
                conn.commit()
                expression_data = downloadFile(file_id, type_id)
                if expression_data == None:
                    return 

                if "associated_entities"  in file_info:
                     
                    for entity in file_info["associated_entities"]: 
                            if "entity_submitter_id" in entity:
                                if checkExistBiospecimen(entity["entity_submitter_id"], cursor):
                                    cursor.execute("INSERT INTO analysis_entity VALUES (%s, %s)", (file_id, entity["entity_submitter_id"]))
                    conn.commit()
                    if type_id == 1:
                            for data_row in expression_data:
                                # Inserimento dei dati di espressioni GENICHE nel database
                                if "gene_type"  in data_row:
                                    gene_type_id = getGeneType(data_row["gene_type"], cursor, conn)
                                else:
                                    continue

                                if "gene_id" not in data_row or "gene_name" not in data_row:
                                    continue

                                gene_id = data_row["gene_id"]
                                gene_name = data_row["gene_name"]
                                # Inserisci il gene nel database
                                cursor.execute("INSERT INTO gene VALUES (%s, %s, %s) ON CONFLICT (gene_id) DO NOTHING;", (gene_id, gene_name , gene_type_id))

                                if "stranded_first" not in data_row or "stranded_second" not in data_row or "tpm_unstranded" not in data_row or "fpkm_unstranded" not in data_row or "fpkm_uq_unstranded" not in data_row or "unstranded" not in data_row:
                                    continue

                                stranded_first = data_row["stranded_first"]
                                stranded_second = data_row["stranded_second"]
                                if stranded_first != 0 and stranded_second != 0: 
                                    cursor.execute("INSERT INTO gene_expression_file VALUES (%s, %s, %s, %s, %s, %s, %s, %s)", (file_id, gene_id, data_row["tpm_unstranded"], data_row["fpkm_unstranded"], data_row["fpkm_uq_unstranded"], data_row["unstranded"], stranded_first, stranded_second))
                            conn.commit()
                    elif type_id == 2:
                            # Inserimento dei dati di espressioni PROTEICHE nel database
                            for data_row in expression_data:
                                #print(data_row)
                                #print("\n")
                                if "AGID" not in data_row:
                                    continue

                                agid = data_row["AGID"]


                                if "lab_id" not in data_row or "catalog_number" not in data_row or "set_id" not in data_row or "peptide_target" not in data_row:
                                    continue
                                # Inserimento della proteina nel database
                                cursor.execute("INSERT INTO protein_GDC VALUES (%s, %s, %s, %s, %s) ON CONFLICT (agid) DO NOTHING;", (agid, data_row["lab_id"], data_row["catalog_number"], data_row["set_id"], data_row["peptide_target"]))

                                if  "protein_expression" in data_row:
                                    if data_row["protein_expression"] != "NaN":
                                        cursor.execute("INSERT INTO protein_expression_file VALUES (%s, %s, %s)", (file_id, agid, data_row["protein_expression"]))

                            conn.commit()
                    print("File inserito nel database")
            else: 
                #print("Il file è gia presente nel database")
                pass
            conn.commit()

        print(f"Download, elaborazione e inserimento dei dati completati.")
        

    except Exception as error:
        # Gestione generica degli errori
        conn.rollback()
        print(f"Errore sconosciuto: {error}")

    finally:
        # Ripristina l'autocommit
        conn.autocommit = True
        # Chiudi la connessione
        cursor.close()
        conn.close()


'''
Funzione che ricerca i dati di un singolo progetto
'''
def searchProject(id):

    try:
        project_url = "https://api.gdc.cancer.gov/projects/" + id
        params = {
            "fields": "name",
            "format": "JSON",
            "pretty": "true"
        }

        response = requests.get(project_url, params=params)

    except Exception as error:
        # Gestione generica degli errori
        print(f"Errore connessione server: {error}")
        return None 


    if response.status_code == 200:

        js = json.loads(response.content.decode("utf-8"))

        if "data" not in js:
            print("Errore nella response per la funzione searchProject")
            return None
        elif "name" not in js["data"]:
            print("Errore nella response per la funzione searchProject")
            return None

        data = js["data"]

        return data["name"]
    else:
        print(f"Errore durante il download del progetto: {response.status_code}")
        return None

'''
Funzione che ricerca i dati di un singolo case
'''
def searchCase(id):
    
    try:
        cases_url = "https://api.gdc.cancer.gov/cases/" + id

        params = {
            "fields": "submitter_id,demographic.ethnicity,demographic.gender,demographic.race,demographic.vital_status,primary_site,disease_type,samples.submitter_id,samples.sample_type,samples.sample_type_id,samples.tumor_code,samples.tumor_code_id,samples.tumor_descriptor,samples.portions.submitter_id,samples.portions.analytes.submitter_id,samples.portions.analytes.concentration,samples.portions.analytes.aliquots.submitter_id,samples.portions.analytes.aliquots.concentration", #samples.portions.slides.submitter_id
            "format": "JSON",
            "pretty": "true"
        }
        response = requests.get(cases_url, params=params)
    except Exception as error:
        # Gestione generica degli errori
        print(f"Errore connessione server: {error}")
        return None 

    if response.status_code == 200:

        js = json.loads(response.content.decode("utf-8"))
        if "data" not in js:
            print("Errore nella response per la funzione searchCase")
            return None
        return js["data"]
    
    print(f"Errore durante il download del case: {response.status_code}")
    return None 

'''
Funzione che inserisce i dati di un singolo nuovo case all'interno del database
'''
def insertNewCase(data, project_id, cursor, conn):
    #search_site_malattia = "SELECT site_id, disease_id FROM primary_site, disease WHERE LOWER(site) = LOWER('{}') AND LOWER(type) = LOWER('{}')"
    insert_case = "INSERT INTO public.case VALUES ('{}', '{}', '{}', '{}', '{}', '{}', '{}', '{}');"

    if "primary_site" not in data or "disease_type" not in data or "submitter_id" not in data or "demographic" not in data or "samples" not in data:
        return 
    elif "ethnicity" not in data["demographic"] or "gender" not in data["demographic"] or "race" not in data["demographic"] or "vital_status" not in data["demographic"]:
        return 
    
    primary_site = getPrimarySite(data["primary_site"], cursor, conn)
    disease = getDisease(data["disease_type"], cursor,conn)

    if "disease_type" in data:        
        cursor.execute(insert_case.format(data["submitter_id"], data["demographic"]["ethnicity"], data["demographic"]["gender"], data["demographic"]["race"], data["demographic"]["vital_status"], project_id, primary_site , disease))
        insertNewSamples(data["samples"], data["submitter_id"], cursor)

'''
Funzione che inserisce i dati dei nuovi campioni all'interno del database
'''
def insertNewSamples(samples, case_id, cursor):
    inserisci_biospecie = "INSERT INTO biospecimen VALUES (%s, %s, %s)"
    inserisci_tumore = "INSERT INTO tumor VALUES (%s, %s, %s) ON CONFLICT (tumor_code_id) DO NOTHING;"
    inserisci_tipo_campione = "INSERT INTO sample_type VALUES (%s, %s) ON CONFLICT (type_id) DO NOTHING;"
    inserisci_campione = "INSERT INTO sample VALUES (%s, %s, %s)"
    inserisci_porzione = "INSERT INTO portion VALUES (%s, %s)"
    inserisci_analita = "INSERT INTO analyte VALUES (%s, %s, %s)"
    inserisci_aliquota = "INSERT INTO aliquote VALUES (%s, %s, %s)"
    #print(samples)
    for sample in samples:

        if "submitter_id" not in sample:
            continue
        sample_id = sample["submitter_id"]
        if "tumor_code_id" in sample and sample["tumor_code_id"] != None and "tumor_code" in sample and "tumor_descriptor" in sample: 
            tumor_code = sample["tumor_code_id"]
            cursor.execute(inserisci_tumore, (tumor_code, sample["tumor_code"], sample["tumor_descriptor"]))
        else: tumor_code = None
        if "sample_type_id" in sample and sample["sample_type_id"] != None and "sample_type" in sample: 
            type_id = sample["sample_type_id"]
            cursor.execute(inserisci_tipo_campione, (type_id, sample["sample_type"]))
        else: type_id = None   


        cursor.execute(inserisci_biospecie, (sample_id, case_id, 1))
        cursor.execute(inserisci_campione, (sample_id, type_id, tumor_code))
        if "portions" in sample:
            for portion in sample["portions"]:
                try:

                    if "submitter_id" not in portion:
                        continue
                    portion_id = portion["submitter_id"]

                    cursor.execute(inserisci_biospecie, (portion_id, case_id, 2))
                    cursor.execute(inserisci_porzione, (portion_id, sample_id))
                    if "analytes" in portion:
                        for analyte in portion["analytes"]:

                            if "submitter_id" not in analyte:
                                continue
                            analyte_id = analyte["submitter_id"]
                            cursor.execute(inserisci_biospecie, (analyte_id, case_id, 3))

                            if "concentration" not in analyte:
                                continue
                            cursor.execute(inserisci_analita, (analyte_id, portion_id, analyte["concentration"]))

                            if "aliquots" in analyte:
                                for aliquote in analyte["aliquots"]:
                                    if "submitter_id" not in aliquote:
                                        continue
                                    aliquote_id = aliquote["submitter_id"]
                                    cursor.execute(inserisci_biospecie, (aliquote_id, case_id, 4))

                                    if "concentration" not in aliquote:
                                        continue
                                    cursor.execute(inserisci_aliquota, (aliquote_id, analyte_id, aliquote["concentration"]))
                except:
                    pass

'''
Funzione che inserisce i dati di un singolo nuovo file all'interno del database
'''
def insertNewAnalysis(file_id, file_info, project_id, cursor, conn):

    if "data_type" not in file_info or "data_category" not in file_info or "experimental_strategy" not in file_info or "file_name" not in file_info or "file_size" not in file_info or "created_datetime" not in file_info or "updated_datetime" not in file_info:
        return None
    
    dataType = getDataType(file_info["data_type"], cursor, conn)
    dataCategory = getDataCategory(file_info["data_category"], cursor, conn)
    expStrategy = getExperimentalStrategy(file_info["experimental_strategy"], cursor, conn)

    # Inserisci i dettagli del file nel database
    cursor.execute("INSERT INTO analysis VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)", (file_id, file_info["file_name"], file_info["file_size"], file_info["created_datetime"], file_info["updated_datetime"], project_id, dataType, dataCategory, expStrategy))
    return dataType

'''
Funzione per il download di un singolo file 
'''
def downloadFile(file_id, type_id):
    try:
        file_url = "https://api.gdc.cancer.gov/data/" + file_id
        response = requests.get(file_url)
    except Exception as error:
        # Gestione generica degli errori
        print(f"Errore connessione server: {error}")
        return None

    if response.status_code == 200:
        # Elabora i dati dal file scaricato
        if type_id == 1: 
            data = pd.read_csv(StringIO(response.text), sep="\t", comment="#", skiprows=[2,3,4,5])
        elif type_id == 2: 
            data = pd.read_csv(StringIO(response.text), sep="\t", comment="#")
    
        # Trasforma i dati in una lista di dizionari
        expression_data = data.to_dict(orient="records")

        return expression_data
    else:
        print(f"Errore durante il download del file: {response.status_code}")
        return None

'''
cursor,conn = databaseConnection()
if cursor == None:
    print("Errore nello stabilire la connessione col server")
    sys.exit()

print("CONNESSIONE AL DATABASE STABILITA")
'''


'''
if __name__ == '__main__':
    download_and_process_expression_data(params)
'''