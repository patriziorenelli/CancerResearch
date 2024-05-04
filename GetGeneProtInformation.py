from selenium import webdriver
from time import sleep
from bs4 import BeautifulSoup
from selenium.webdriver.common.by import By
import random
from DatabaseGenericInteraction import *

# Con questa funzione otteniamo i valori dei peptide per ogni aliquot associata al gene, fa utilizzo di Selenium per effettuare il web scraping 
def getGeneProInformation(geneName, gene_id, cursor, conn):
      # Il chromedriver deve essere di versione compatibili con la versione di chrome (applicazione) installata
      p = r"PATH\chromedriver.exe"
      driver = webdriver.Chrome(p)
      # Costruzione link per ottenere informazioni su un gene 
      link = "https://pdc.cancer.gov/pdc/gene/" + geneName
      # Il Chromedriver aprirà il link richiesto
      driver.get(link)
      #Attesa di 5 secondi per consentire il completamento del caricamento della pagina
      sleep(5) 
      # Individua il bottone contenente il testo 'Continue' nel popup di disclaimer presente nella pagina e simula il click di un utente
      driver.find_element_by_xpath ("//button[contains( text( ), 'Continue')]").click()
      # Attesa di 3 secondi 
      sleep(5)
      cambio = True
      pagina = 0
      while cambio:
            pagina+=1
            # Ottengo il codice html della pagina web 
            html_source = driver.page_source
            #driver.quit()
            # Trovo la posizione della stringa indicata
            pos = html_source.index('Biospecimens/Samples in Which the Gene Product Was Detected')
            # Escludo una parte del codice html trovato
            valid = html_source[pos:]
            
            # BeautifulSoup ci consente di effettuare il parsing del codice html, crea un albero di analisi per il codice della pagina web scaricata
            tr = BeautifulSoup(valid)
            # La funzione get_text() ci consente di recuperare il contenuto dei tag html, utilizzando il carattere ' ' = spazio come divisore degli elementi
            str_value = tr.get_text(' ')

            valori = (str_value.split('Aliquot')[1]).split("Plex (Dataset alias) Label Study Experiment Type Spectral Counts Distinct Peptides Unshared Peptides Precursor Area Unshared Area Log2 Ratio Unshared Log2 Ratio")
            # Array contenente le Aliquot
            aliquot = [elem for elem in valori[0].split(' ') if elem.strip()]
            # Otteniamo i valori da salvare per ogni aliquote 
            str_value = tr.get_text('\n')
            valori = [element.replace('\n', '*') for element in str_value.split('\n')][25:]
            # Array che contiene le etichette dei reagenti 
            lab_c = ["label_free", 
            "itraq_113",
            "itraq_114",
            "itraq_115",
            "itraq_116",
            "itraq_117",
            "itraq_118",
            "itraq_119",
            "itraq_121",
            "tmt_126",
            "tmt_127n",
            "tmt_127c",
            "tmt_128n",
            "tmt_128c",
            "tmt_129n",
            "tmt_129c",
            "tmt_130n",
            "tmt_130c",
            "tmt_131",
            "tmt_131c",
            "tmt_132n",
            "tmt_132c",
            "tmt_133n",
            "tmt_133c",
            "tmt_134n",
            "tmt_134c",
            "tmt_135n",]

            # La tabella in ogni pagina contiene al massimo 10 righe 
            for y in range(0, 10):
                  # Ci consente di individuare la posizione corretta degli elementi all'interno della tabella fare un for che trova la posizione corretta da usare se y*1 y*2 y*3 y*4 y*5 y*6 y*7 y*8 y*9 y*10 y*11
                  for x in range(0, 12):
                        # Controlliamo il valore di label per capire se mancano dei dati nella tabella di pdc
                        if  y*x < len(valori) and valori[0+y*x].lower() in lab_c:
                              break 
                        else:
                              pass         
                  const = x   
                  # Controlliamo il valore di label per capire se mancano dei dati nella tabella di pdc
                  if y*const < len(valori) and  valori[0+y*const].lower() in lab_c:
                        # I valori da salvare nella tabella 
                        label = valori[0+y*const]
                        study = valori[1+y*const]
                        spectral_count = valori[3+y*const]
                        distinct_peptides = valori[4+y*const]
                        unshared_peptides = valori[5+y*const]
                        # lo study è study name -> serve prendere prima lo study id dal db, effettuare check esistenza del progetto e fare la successiva insert nel database 
                        study = getProjectId(study, cursor)

                        if study != None  and y < len(aliquot) and checkExistAliquote(aliquot[y], cursor) and checkExistGene(gene_id, cursor):
                              if checkExistProtein_PDC(gene_id, study, aliquot[y], cursor)  :
                                    query = "UPDATE  public.protein_PDC SET label = '{}' , spectral_count = '{}' , distinct_pepetides = '{}', unshared_peptides = '{}'  WHERE gene_id = '{}' and aliquot = '{}' and project_id = '{}' ;".format(label, spectral_count, distinct_peptides, unshared_peptides,  gene_id, aliquot[y], study )
                                    #print(query)
                                    cursor.execute(query)
                                    conn.commit()
                              else:
                                    query = "INSERT INTO public.protein_PDC VALUES ('{}', '{}', '{}', '{}', '{}', '{}', '{}', '{}', '{}') ON CONFLICT (gene_id,project_id,aliquot) DO NOTHING;".format(gene_id, label, spectral_count, distinct_peptides, unshared_peptides, 0, 0, aliquot[y], study )
                                    #print(query)
                                    cursor.execute(query)
                                    conn.commit()

            try:  
                  # Scorriamo alla pagina successiva     
                  driver.find_element(By.XPATH,'//*[@id="aliquotRecordTable"]/div/p-paginator/div/a[3]/span').click()
                  # Attendiamo un tempo casuale tra 0.86s e 4.86s in modo da cercare di evitare il controllo di sicurezza anti-robot con captcha
                  sleep(0.95 + random.randint(0,4))
            except:
                  cambio = False

              

