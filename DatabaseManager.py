
import psycopg2
import params
import os 
import datetime
import pathlib


'''
Costruisce la struttura delle tabelle del database e richiama la funzione necessaria al riempimento delle tabelle di base
'''
def databaseConstruction():
    cursor,conn = databaseConnection()
    if cursor == None:
        print("FATAL ERROR")
        return cursor, conn

    errore = False
    check_table = "SELECT COUNT(*) FROM primary_site"

    try:
        cursor.execute(check_table)
        result = cursor.fetchone()
    except psycopg2.Error as db_error:
        if "non esiste" in str(db_error):
            # Procedo con la creazione delle tabelle 
                    cursor.close()
                    conn.close()
                    cursor, conn = databaseConnection()
                    # Carico le singole query da eseguire per creare le tabelle 
                    query = open(".\\load\\Schema.sql", "r").read()
                    #query = open(".\\Backup\\schema_NUOVO_DA_TESTARE.sql", "r").read()

                    query = query.split(";")
                    # Eseguo le singole query e controllo il loro risultato 
                    for x in query: 
                        try:
                            cursor.execute(x)
                            conn.commit()
                        except psycopg2.Error as db_error:
                            if "can't execute an empty query" not in str(db_error):
                                print("FATAL", db_error)
                                errore = True
                                cursor = None
                                conn = None
                                break
                    
                    if not errore:
                        print("Struttura del databse creata")    

                        
                        if fillingBasicTable(cursor, conn):
                            #print("Errore nella riempimento delle tabelle base")
                            cursor = None
                            conn = None
                        else:
                            cursor, conn = databaseConnection()

                    return cursor, conn
    return cursor,conn   
'''
Funzione per la creazione del database e la creazione delle sue tabelle 
 '''
def databaseCreation():
    try: 
        conn = psycopg2.connect(database = "postgres", user = params.user, password = params.password , host = params.host , port = params.port)
        conn.autocommit = True
        cursor = conn.cursor()
        sql = 'CREATE database "GDC"'
        cursor.execute(sql)
        print("Database creato")

        # Se esiste un backup del database lo ripristino 
        backup_file = [stringa for stringa in os.listdir(".//Backup") if stringa.startswith("Backup")]
        backup_file = []
        if len(backup_file) > 0:
            print("RECUPERO DATI DAL BACKUP: ")
            reloadData()
            print("RECUPERO DEI DATI DAL BACKUP TERMINATO")

            # Chiudo la vecchia connessione ormai non più valida
            cursor.close()
            conn.close()

            #Genero una nuova connessione al database
            connection = psycopg2.connect(**params.db_params)
            cursor = connection.cursor()
            connection.autocommit = False 

            return cursor, connection
        else:
            return databaseConstruction()


    except psycopg2.Error as db_error:
        print("Errore nella connesione al database")
        print("Dettaglio errore: ", db_error)
        return None, None

'''
Funzione che si occupa di gestire la creazione di una connessione con il database 
Se il database non esiste viene creato e vengono create tutte le tabelle necessarie 
'''
def databaseConnection():
    try: 
        # Crea una connessione al database PostgreSQL
        connection = psycopg2.connect(**params.db_params)
        # Crea un cursore per eseguire query SQL
        cursor = connection.cursor()
        # Inizia la transazione
        connection.autocommit = False 
        return cursor, connection
    except psycopg2.Error as db_error:
        # Se il database non esiste lo vado a creare 
        if "non esiste" in str(db_error) :
            print("DATABASE NON ESISTENTE" )
            return databaseCreation()
        else:
            print("Errore nella connesione al database")
            print("Dettaglio errore: ", db_error)
            return None 


''' 
Funzione che si occupa di riempire le tabelle di base necessarie al funzionamento delle altre funzioni
Sfrutta dei file contenuti nella cartella load contenenti i dati da caricare 
'''
def fillingBasicTable(cursor,conn):
    #basicTable = ["primary_site", "biospecimen_type", "data_category", "data_type", "disease", "experimental_strategy", "gene_type"]
    basicTable = [ "biospecimen_type"]

    basicQuery = "INSERT INTO public.{} VALUES (DEFAULT, '{}');"

    try: 
        for tab in basicTable:
            # dobbiamo aprire il relativo file e poi fare un for per le insert 
            path = ('.\\load\\{}.txt').format(tab)
            ins = open( path, "r").read()
            ins = ins.split('\n')
            
            for val in ins:

                quer = basicQuery.format(tab, val )
                cursor.execute(quer)
                conn.commit()
    except:
        print("ERRORE FATALE NEL RIEMPIMENTO DELLE TABELLE")
        return True

    return False

'''
Funzione che si occupa di eseguire un backup del database e gestisce le copie salvate, lasciandone 2 sempre
'''
def saveDatabase():
    timestamp = datetime.datetime.now()
    timestamp = str(timestamp.day) + "-" + str(timestamp.month) + "-" + str(timestamp.year) + "_" + str(timestamp.hour) + "-" + str(timestamp.minute) + "-" + str(timestamp.second)  
    # Questa è la riga di comando fondamentale che attraverso il tool pg_dump nativo di postgres esegue il backup del database 
    # Necessario cambiare path per pg_dump.exe e quello da cui si recupera il file (il path indicato dopo >)
    # Necessarop cambiare credenziali di accesso a postgress nel path postgresql://{username}:{password}@localhost:{port}/GDC
    comando = '" PATH\pg_dump.exe" postgresql://postgres:1234@localhost:5432/GDC > "PATH\Backup_"'+ timestamp
    
    os.system(('cmd /C {}').format(comando))
    print("BACKUP CREATO: Backup_" + timestamp)
   
    # Gestione del salvataggio di solo 2 copie del database 
    # Prendo solo i file che inziano con la stringa Backup
    backup_file = [stringa for stringa in os.listdir(".//Backup") if stringa.startswith("Backup")]
    # Se ho piu' di 2 backup allora creo un dizionario nomeFile : timestamp creazione (ottenuta da os)
    if len(backup_file) > 2:
        dict_backup_date = dict()
        for x in backup_file:
            dict_backup_date[x] = float(pathlib.Path(".//Backup//"+x).stat().st_mtime)

        # Ordino i file in base al timestamp 
        dict_backup_date = sorted(dict_backup_date, key=lambda k: dict_backup_date[k], reverse=True)

        # Eliminiamo le copie di backup più vecchie mantenendo solo le ultime 2 più recenti
        for fileName in dict_backup_date[2:]:
            os.remove(".//Backup//"+fileName)

'''
Funzione che si occupa di eseguire il SOLO ripristino del database dalla copia piu' recente disponibile 
(SE SI VUOLE RIPRISTINARE SENZA AVERE IL DB CREATO LANCIARE PRIMA LA FUNZIONE databaseConstruction() )
'''
def reloadData():

    backup_file = [stringa for stringa in os.listdir(".//Backup") if stringa.startswith("Backup")]

    # creo un dizionario nomeFile : timestamp creazione (ottenuta da os) per ordinare i file in base al timestamp ed ottenere poi il più recente 
    dict_backup_date = dict()
    for x in backup_file:
        dict_backup_date[x] = float(pathlib.Path(".//Backup//"+x).stat().st_mtime)
    # Ordino i file in base al timestamp 
    dict_backup_date = sorted(dict_backup_date, key=lambda k: dict_backup_date[k], reverse=True)

    print(dict_backup_date[0])
    # Questa è la riga di comando fondamentale che attraverso psql di postgres esegue la restore del database dal backup più recente 
    # Necessario cambiare path per psql.exe e quello da cui si recupera il file (il path indicato dopo >)
    # Necessarop cambiare credenziali di accesso a postgress nel path postgresql://{username}:{password}@localhost:{port}/GDC
    comando = '"PATH\psql.exe" postgresql://postgres:1234@localhost:5432/GDC < "PATH\Backup\\"' + dict_backup_date[0]
    os.system(('cmd /C {}').format(comando))


#reloadData()
databaseConnection()
#saveDatabase()
#databaseCreation()