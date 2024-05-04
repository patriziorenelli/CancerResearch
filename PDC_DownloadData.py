from DatabaseManager import *
import requests
import json
from DatabaseGenericInteraction import * 
from GetGeneProtInformation import *
from GetLog2GeneAliquot import *


#  I dati su PDC girano intorno ai singoli studi, organizzando i dati in studi basati su su caratteristiche come proteoma o PTM e tipi sperimentali come senza etichetta o TMT 
# In genere gli studi sono specifici per la malattia e per primary site 
class Program:
    def __init__(self, program_id, program_submitter_id, program_name, projects):
        self.program_id = program_id
        self.program_submitter_id = program_submitter_id
        self.program_name = program_name
        self.projects = projects    # è un array

class Project:
   def __init__(self,  project_id, project_submitter_id, project_name, studies):
    self. project_id =  project_id
    self.project_submitter_id = project_submitter_id
    self.project_name = project_name
    self.studies = studies  # è un array
    

class Study:
    def __init__(self, study_id, study_submitter_id, study_name, disease_types, primary_sites, pdc_study):
        self.study_id = study_id
        self.study_submitter_id = study_submitter_id
        self.study_name = study_name
        self.disease_types = disease_types  # è un array
        self.primary_sites = primary_sites  # è un array
        self.pdc_study = pdc_study

'''
Funzione che ottiene i dati di tutti gli studi, programmi e progetti disponibili, necessari per ottenere tutti gli altri dati disponibili -> Superflui vista la struttura del nostro db
'''
def getPrograms2():
    try:
        gdc_api_url = "https://pdc.cancer.gov/graphql?query={allPrograms (acceptDUA: true)  {program_id  program_submitter_id  name projects  {project_id  project_submitter_id  name  studies  {pdc_study_id study_id study_submitter_id submitter_id_name analytical_fraction study_name disease_types primary_sites embargo_date experiment_type acquisition_type} }}}"
        response = requests.get(gdc_api_url)        
        programs = []
        prog = json.loads(response.content.decode("utf-8"))

        if "data" not in prog and "allPrograms" not in prog["data"]:
            return 

        # Da qui possiamo prendere i programmi == project e tutti i studi associati ad ogni programma 
        for program in prog["data"]["allPrograms"]:
            studies = []
            projects = []
            if "projects" not in program:
                continue
            for pr in program["projects"]:
                if "studies" not in pr:
                    continue 
                for st in pr["studies"]:
                    # Qui si puo' applicare il filtro sul primary_site che è relativo ai study 
                    if "study_id" not in st or "study_submitter_id" not in st or "study_name" not in st or "disease_types" not in st or "primary_sites" not in st:
                        continue 
                    studies.append(Study(st["study_id"], st["study_submitter_id"], st["study_name"], st["disease_types"],st["primary_sites"],  st["pdc_study_id"]))
                if "project_id" not in pr or "project_submitter_id" not in pr or "name" not in pr:                
                    continue
                projects.append(Project(pr["project_id"], pr["project_submitter_id"], pr["name"], studies))
                if "program_id" not in program or "program_submitter_id" not in program or "name" not in program:
                    continue
            programs.append( Program(program["program_id"], program["program_submitter_id"],program["name"], projects ) )

        return programs
    except:
        print("ERRORE FATALE!")
        return None

'''
La funzione che usiamno realmente per ottenere i progetti 
'''
def getProgram():

        gdc_api_url = "https://pdc.cancer.gov/graphql?query={ studyCatalog (acceptDUA: true) { pdc_study_id versions { study_id study_submitter_id submitter_id_name study_shortname study_version is_latest_version } } }"
        response = requests.get(gdc_api_url)
        studies = []
        st = json.loads(response.content.decode("utf-8"))
        
        for study in st["data"]["studyCatalog"]:
            studies.append(Study(study["versions"][0]["study_id"], study["versions"][0]["study_submitter_id"], study["versions"][0]["submitter_id_name"], None, None, study["pdc_study_id"] ))
        return studies

'''
Funzione che ottiene i case relativi ad un singolo programma
'''
def getCases(study, study_id, cursor, conn):
        # Otteniamo il numero di cases relativi allo studio
        gdc_api_url = 'https://pdc.cancer.gov/graphql?query={paginatedCaseDemographicsPerStudy (study_id: "' + study +  '" offset: 0 limit: 1 acceptDUA: true) { total caseDemographicsPerStudy { case_id case_submitter_id disease_type primary_site demographics { demographic_id ethnicity gender demographic_submitter_id race cause_of_death days_to_birth days_to_death vital_status year_of_birth year_of_death age_at_index premature_at_birth weeks_gestation_at_birth age_is_obfuscated cause_of_death_source occupation_duration_years country_of_residence_at_enrollment} } pagination { count sort from page total pages size } }}'
        response = requests.get(gdc_api_url)
        # Controlliamo se l'esito della query all'API è valido
        if response.status_code == 204 or response == None or response.content == None:
            return

        cases_number = json.loads(response.content.decode("utf-8"))
        if "data" not in cases_number or "paginatedCaseDemographicsPerStudy" not in cases_number["data"] or "total" not in cases_number["data"]["paginatedCaseDemographicsPerStudy"]:
            return 
        cases_number = str( json.loads(response.content.decode("utf-8"))["data"]["paginatedCaseDemographicsPerStudy"]["total"] )

        # Otteniamo i cases relativi allo studio utilizzando il giusto valore per il limit 
        gdc_api_url = 'https://pdc.cancer.gov/graphql?query={paginatedCaseDemographicsPerStudy (study_id: "' + study +  '" offset: 0 limit: ' + cases_number + ' acceptDUA: true) { total caseDemographicsPerStudy { case_id case_submitter_id disease_type primary_site demographics { demographic_id ethnicity gender demographic_submitter_id race cause_of_death days_to_birth days_to_death vital_status year_of_birth year_of_death age_at_index premature_at_birth weeks_gestation_at_birth age_is_obfuscated cause_of_death_source occupation_duration_years country_of_residence_at_enrollment} } pagination { count sort from page total pages size } }}'
        response = requests.get(gdc_api_url)
        re = json.loads(response.content.decode("utf-8"))
        if "data" not in re or "paginatedCaseDemographicsPerStudy" not in re["data"]:
            return

        # Otteniamo i valori necessari per inizializzare un nuovo case (paziente), prendendo i valori sulla sua "demografia"
        for case in json.loads(response.content.decode("utf-8"))["data"]["paginatedCaseDemographicsPerStudy"]["caseDemographicsPerStudy"]:
                try:
                    # Controlliamo che un case non sia già presente all'interno del database e otteniamo i valori necessari già registrati sul database oppure li aggiungiamo (come da definizione delle funzioni getDisease, getPrimarySite)
                    if not checkExistCase(case['case_submitter_id'],cursor):
                        if 'demographics' in case:
                            disease = getDisease(case['disease_type'], cursor, conn)
                            primary_site = getPrimarySite( case['primary_site'], cursor, conn)
                            sub_id = case['case_submitter_id']
                            case = case['demographics'][0]
                            insertNewCase(sub_id, case['ethnicity'], case['gender'], case['race'], case['vital_status'], study_id, primary_site, disease, cursor, conn)
                except:
                    print("AIUTO CHECK ")
                    conn.rollback()
                    continue

'''
Funzione che ottiene tutti i campioni biologici di un singolo studio
'''
def getSample(study, cursor, conn):

    # Otteniamo il numero dei sample associati ad uno studio
    query = 'https://pdc.cancer.gov/graphql?query={ paginatedCasesSamplesAliquots(pdc_study_id:"' + study +  '"  offset:0 limit: 1 ) { total casesSamplesAliquots { case_id case_submitter_id days_to_lost_to_followup disease_type index_date lost_to_followup primary_site samples { sample_id sample_submitter_id sample_type sample_type_id gdc_sample_id gdc_project_id biospecimen_anatomic_site composition current_weight days_to_collection days_to_sample_procurement diagnosis_pathologically_confirmed freezing_method initial_weight intermediate_dimension longest_dimension method_of_sample_procurement pathology_report_uuid preservation_method sample_type_id shortest_dimension time_between_clamping_and_freezing time_between_excision_and_freezing tissue_type tumor_code tumor_code_id tumor_descriptor diagnoses{ diagnosis_id diagnosis_submitter_id annotation} aliquots { aliquot_id aliquot_submitter_id analyte_type aliquot_run_metadata { aliquot_run_metadata_id label experiment_number fraction replicate_number date alias analyte} } } } pagination { count sort from page total pages size } } }'
    gdc_api_url = query
    
    response = requests.get(gdc_api_url)
    # Controlliamo che il response sia valido 
    if json.loads(response.content.decode("utf-8"))["data"]["paginatedCasesSamplesAliquots"] == None:
         return
    samples_number = str( json.loads(response.content.decode("utf-8"))["data"]["paginatedCasesSamplesAliquots"]["total"] )

    # Ottengo i valori sui semple e gli altri campioni biologici associati
    gdc_api_url = 'https://pdc.cancer.gov/graphql?query={ paginatedCasesSamplesAliquots(pdc_study_id:"' + study +  '"  offset:0 limit: ' + samples_number + ' ) { total casesSamplesAliquots { case_id case_submitter_id days_to_lost_to_followup disease_type index_date lost_to_followup primary_site samples { sample_id sample_submitter_id sample_type sample_type_id gdc_sample_id gdc_project_id biospecimen_anatomic_site composition current_weight days_to_collection days_to_sample_procurement diagnosis_pathologically_confirmed freezing_method initial_weight intermediate_dimension longest_dimension method_of_sample_procurement pathology_report_uuid preservation_method sample_type_id shortest_dimension time_between_clamping_and_freezing time_between_excision_and_freezing tissue_type tumor_code tumor_code_id tumor_descriptor diagnoses{ diagnosis_id diagnosis_submitter_id annotation} aliquots { aliquot_id aliquot_submitter_id analyte_type concentration aliquot_run_metadata { aliquot_run_metadata_id label experiment_number fraction replicate_number date alias analyte} } } } pagination { count sort from page total pages size } } }'
    response = requests.get(gdc_api_url)
    #print("sample number: " + str( samples_number)  )
    
    samples =  json.loads(response.content.decode("utf-8"))["data"]["paginatedCasesSamplesAliquots"]["casesSamplesAliquots"]

    for s in samples: 
        case_submitter_id = s['case_submitter_id']
        for sample in s["samples"]:
                    # Ottengo i valori sul tumore analizzato
                    tumor_code = sample['tumor_code']
                    tumor_code_id =  sample['tumor_code_id']
                    tumor_description = sample['tumor_descriptor']
                    # Se i valori dei tumori sono validi e il tumore non è già presente nel database lo aggiungo
                    if  len(str(tumor_code_id)) != 0 and len(str(tumor_code)) != 0 and len(str(tumor_description)) != 0 and  str(tumor_code_id) != 'None' and str(tumor_code_id) != 'null' and not checkExistTumor(tumor_code_id, cursor):
                            insertNewTumor(tumor_code_id, tumor_code, tumor_description, cursor, conn)
                    else:
                        tumor_code_id = None
                    # Ottengo i valori sul sample 
                    sample_type_id = sample['sample_type_id']
                    sample_type = sample['sample_type']
                    
                    if sample_type_id == None:
                        sample_type_id = searchSampleTypeId(sample_type, cursor)
                        if sample_type_id == None:
                                sample_type_id = None
                        else:
                            sample_type_id = sample_type_id[0] 


                    # Se un sample type non è presente lo salvo all'interno del database
                    elif not searchSampleTypeId(sample_type, cursor):
                        if len(sample_type_id) > 0:
                            sample_type_id = int(sample_type_id.replace('"',''))
                            insertNewSampleType(sample_type_id, sample_type, cursor, conn)
                            sample_type_id =  searchSampleTypeId(sample_type, cursor)
                        else:
                            # imposto un sample_type_id = None perchè impostare un qualunque altro valore potrebbe causare poi conflitti con altri sample_type con type_id correttamente impostato
                            sample_type_id =  None

                    else:
                        # Effettuo comunque la conversione del sample_type_id nel sample_type_id ottenuto dal database ottenuti attraverso il sample_type
                        # Questa azione è necessaria poichè potrebbero esserci sample_type_id che sono collegati ad un sample_type che è salvato nel database con un altro sample_type_id causando errori nel momento dell'insermento nel database 
                        sample_type_id =  searchSampleTypeId(sample_type, cursor)


                    # sample['sample_submitter_id'] può essere un'array quindi bisogna fare uno split su , e fare l'insert per i singoli samp in sample['sample_submitter_id]
                    for samp in sample['sample_submitter_id'].split(','):
                        sample_id =  samp

                        # 2 IF DA RIATTIVARE -> SI POSSONO LASCIARE ANCHE COSI' TANTO LE QUERY HANNO LE CLAUSOLO ON CONFLICT 
                        # TANTO COSI' EVITEREMMO TENTATIVI DI INSERT, MA FAREMMO COMUNQUE 2 SELECT SEMPRE + EVENTUALMENTE 5 INSERT IN CASO DI ESITO POSITIVO DEI CHECK 
                        #if not checkExistBiospecimen(sample_id, cursor):
                        #    if checkExistCase(case_submitter_id, cursor):
                                
                        insertNewBiospecimen(sample_id, case_submitter_id, 1, cursor, conn)
                        insertNewSample(sample_id, sample_type_id, tumor_code_id, cursor, conn)
                    for aliquote in sample['aliquots']:
                                    concentration = aliquote['concentration']
                                    aliquote_id = aliquote['aliquot_submitter_id']
                                    # In PDC non abbiamo questi dati separati quindi l'unico modo per rispettare i vincoli imposti dal database per la parte di GDC è quella di creare portion e analyte fittizi se no ci basterebbe aliquote
                                    insertNewPortion(sample_id, sample_id, cursor, conn)
                                    insertNewAnalyte(sample_id, sample_id, concentration, cursor,conn)
                                    # SI E' TOLTA LA FK TRA ALIQUOTE E ANALYTE & ALIQUOTE E BIOSPECIMEN ->  PER EFFETTUARE QUESTA MODIFICA, PRIMA SI USAVA SAMPLE_ID ANCHE IN ALIQUOTE AL POSTO DI ALIQUOTE_ID -> bisogna vedere se la funzione di scraping web funziona così
                                    insertNewAliquote(aliquote_id, sample_id, concentration, cursor, conn)

'''
Funzione che ottiene i dati dei geni disponibili 
'''
def getGenes(cursor,conn):
    # Otteniamo tutti i geni con le proteine e informazioni associate 
    gdc_api_url = 'https://pdc.cancer.gov/graphql?query={ getPaginatedGenes(offset: 0 limit: 20000 acceptDUA: true) { total genesProper { gene_id gene_name NCBI_gene_id authority description organism chromosome locus proteins assays } pagination { count sort from page total pages size } } }'
    response = requests.get(gdc_api_url)
    re = json.loads(response.content.decode("utf-8"))
    # Controllo la response
    if "data" not in re or "getPaginatedGenes" not in re["data"] or "genesProper" not in re["data"]["getPaginatedGenes"]:
        return 
    
    geni_found = dict()
    all_program = dict()

    
    for gene in  json.loads(response.content.decode("utf-8"))["data"]["getPaginatedGenes"]["genesProper"]:
        gene_name = gene["gene_name"]
        # Effetto la request per ottenere i dati del singolo gene 
        gdc_api_url= 'https://pdc.cancer.gov/graphql?query={geneSpectralCount (gene_name: "' + gene_name+ '" acceptDUA: true){gene_id gene_name NCBI_gene_id authority description organism chromosome locus proteins assays spectral_counts { project_submitter_id plex spectral_count distinct_peptide unshared_peptide study_submitter_id pdc_study_id} }}'
        response = requests.get(gdc_api_url)

        if response == None or response.content == None:
            continue
        try:
            re = json.loads(response.content.decode("utf-8"))
        except:
             continue
        
        if "data" not in re or "geneSpectralCount" not in re["data"]:
            break
        
        geneSpectralCount = re["data"]["geneSpectralCount"]


        for genSpec in geneSpectralCount:
            if "spectral_counts" not in genSpec or "proteins" not in genSpec or "gene_name" not in genSpec and "project_submitter_id" in genSpec["spectral_counts"] :
                break
            # PROTEINE ASSOCIATE AI GENE 
            proteins = (genSpec["proteins"]).split(";")

            ncbi_gene_id = str( genSpec["NCBI_gene_id"])
            gene_name = genSpec["gene_name"] 
            study = genSpec["spectral_counts"]

            studies = []

            # Ottengo i submitter_id che identificano i vari studi all'interno del database
            # Necessario poichè la query per i geni ottiene soltanto i pdc_study_id
            gdc_api_url = "https://pdc.cancer.gov/graphql?query={allPrograms (acceptDUA: true)  {program_id  program_submitter_id  name projects  {project_id  project_submitter_id  name  studies  {pdc_study_id study_id study_submitter_id submitter_id_name analytical_fraction study_name disease_types primary_sites embargo_date experiment_type acquisition_type} }}}"
            response = requests.get(gdc_api_url)
            if response == None or response.content == None:
                #print(response, response.content)
                continue
            try:
                prog = json.loads(response.content.decode("utf-8"))
            except:
                 continue

            if "data" not in prog and "allPrograms" not in prog["data"]:
                continue 

            # Da qui possiamo prendere i programmi == project e tutti i studi associati ad ogni programma 
            for program in prog["data"]["allPrograms"]:
                studies = []
                if "projects" not in program:
                    break
                for pr in program["projects"]:
                    if "studies" not in pr:
                        break 
                    for st in pr["studies"]:
                        # Qui si puo' applicare il filtro sul primary_site che è relativo ai study 
                        if "study_submitter_id" not in st or "pdc_study_id" not in st :
                            break 
                        all_program[st["pdc_study_id"]] = st["study_submitter_id"]


            # Array degli studi ottenuti dalla query dei geni da convertire                
            for st in study:
                studies.append(st['pdc_study_id'])

            headers = {
                'Accept': 'application/json',
                'Content-Type': 'application/json',
            }
            json_data = {
                'gene_ids': [
                ncbi_gene_id
                ]

            }

            # Effettuo la traduzione da gene_id NCBI in gede_id STANDARD Ensembl ed ottengo il gene_type non presente in PDC, per farlo usiamo l'API offerta da NCBI
            gene_tr = json.loads( requests.post('https://api.ncbi.nlm.nih.gov/datasets/v2alpha/gene', headers=headers, json=json_data).content.decode("utf-8") )

            if "reports" in gene_tr and len(gene_tr["reports"]) >=1 and "gene" in gene_tr["reports"][0] and "ensembl_gene_ids" in gene_tr["reports"][0]["gene"] and len(gene_tr["reports"][0]["gene"]["ensembl_gene_ids"]) >= 1 and ( (gene_tr["reports"][0]["gene"]["ensembl_gene_ids"])[0]).startswith('ENSG') and "type" in  gene_tr["reports"][0]["gene"]:
                gene_id = (gene_tr["reports"][0]["gene"]["ensembl_gene_ids"])[0]
                #print(gene_id, gene_name, gene_tr["reports"][0]["gene"]["type"], getGeneType(gene_tr["reports"][0]["gene"]["type"], cursor, conn))
                insertNewGene(gene_id, gene_name, getGeneType(gene_tr["reports"][0]["gene"]["type"], cursor, conn), cursor, conn)
                geni_found['gene_id'] = gene_name
                #Nella tabella gene potremo avere comunque situazioni del genere:
                #ENSG00000109576.14   AADAT 1  -> Causato dalle transcription che sfrutta GDC 
                #ENSG00000109576    AADAT 1

                # Si effettua l'inserimento delle proteine associate a ciascun gene 
                for prot in proteins:
                    for st in studies:
                        if st != None:
                            # Vado ad associare ad ogni gene le proteine individuate e in quale studio il gene e le proteine sono state individuate 
                            if checkExistProject(all_program[st],cursor):
                                insertNewGeneProteinStudy(gene_id, prot, all_program[st], cursor, conn)
                 
    '''
    #   Con la funzione che sfrutta selenium                  
    for gen in list(geni_found.keys()):
        GetGeneProInformation(geni_found[gen], gen, cursor, conn)

    # Con la funzione per la raccolta dei log2
    for st in range(0, len(all_program) ):
            key = list(all_program.keys())[st]
            getLog2RatioInfo(key, all_program[key], cursor, conn) 
    '''
        
    pass

'''
Funzione utilizata per ottenere maggiori infomazioni sui geni e sui loro valori
'''
def getProteinInfo(cursor, conn):
    geni_found = dict()
    # Seleziono tutti i geni dal database 
    cursor.execute("SELECT * FROM gene ")
    result = cursor.fetchall()
    for x in result:
         geni_found[x[0]] = x[1]

    for gen in list(geni_found.keys()):
        getGeneProInformation(geni_found[gen], gen, cursor, conn)
        pass
    all_program = dict()
    # Ottengo i valori di riferimento di tutti i study disponibili in PDC
    gdc_api_url = "https://pdc.cancer.gov/graphql?query={allPrograms (acceptDUA: true)  {program_id  program_submitter_id  name projects  {project_id  project_submitter_id  name  studies  {pdc_study_id study_id study_submitter_id submitter_id_name analytical_fraction study_name disease_types primary_sites embargo_date experiment_type acquisition_type} }}}"
    response = requests.get(gdc_api_url)
    prog = json.loads(response.content.decode("utf-8"))


    # Da qui possiamo prendere i programmi == project e tutti i studi associati ad ogni programma 
    for program in prog["data"]["allPrograms"]:
                if "projects" not in program:
                    break
                for pr in program["projects"]:
                    if "studies" not in pr:
                        break 
                    for st in pr["studies"]:
                        # Qui si puo' applicare il filtro sul primary_site che è relativo ai study 
                        if "study_submitter_id" not in st or "pdc_study_id" not in st :
                            break 
                        all_program[st["pdc_study_id"]] = st["study_submitter_id"]

    for st in range(0, len(all_program) ):
            key = list(all_program.keys())[st]
            getLog2RatioInfo(key, all_program[key], cursor, conn) 
         

def sviluppo(cursor, conn):           
    programmi = []
    # getProgram è la variante sopra indicata che consente di raccogliere i dai anche su progetti e programmi, se si utilizza questa funzione sarà necessario eseguire il blocco di codice commentato sotto
    #programmi = getPrograms2()
    programmi = getProgram()

    '''
    Creo un array studies in modo da evitare l'iterazione dei 3 cicli annidati per le funzioni di ottenimento case e sample
    studies = []
    for programma in programmi:
        for progetto in programma.projects:
            for studio in progetto.studies: 
                studies.append(studio)
                # Effettuiamo l'insert degli study che sono l'equivalente dei project all'interno di PDC 
                # RIATTIVARE QUESTO IF
                #if not checkExistProject(studio.study_submitter_id, cursor):
                insertNewProject(studio.study_submitter_id, studio.study_name, cursor, conn)
                print(progetto.project_name)
            #print("finita fase di inserimento progetti")
    '''
    # Insertiamo i nuovi progeti all'interno del database
    for studio in programmi: 
        insertNewProject(studio.study_submitter_id, studio.study_name, cursor, conn)
        

    x = 0
    # La transformazione in set e in list è fondamentale quanto si utilizza la funzione getPrograms2() poichè elimino i study duplicati 
    # portandoli così da 942 a 144 come indicato sul sito PDC
    studies_set = set(programmi)
    studies = list(studies_set)

    x = 0
    for st in range(len(studies), len(studies)):
                x+=1
                getCases(studies[st].study_id, studies[st].study_submitter_id, cursor, conn)
                getSample(studies[st].pdc_study, cursor, conn)
    getProteinInfo(cursor, conn)

# Creo connessione con il database 
cursor, conn = databaseConnection()
sviluppo(cursor, conn)
