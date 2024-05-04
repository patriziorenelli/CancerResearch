'''
Verifica se esiste un id nella tabella biospcimen
'''
def checkExistBiospecimen(id, cursor):
    cursor.execute("SELECT COUNT(*) FROM biospecimen WHERE id = %s", (id,))
    result = cursor.fetchone()
    if result[0] == 0: 
        return False
    return True


def insertNewBiospecimen(sample_id, case_submitter_id, type, cursor,conn):
    query = "INSERT INTO public.biospecimen VALUES ('{}', '{}', '{}') ON CONFLICT (id) DO NOTHING;".format(sample_id, case_submitter_id, type)
    cursor.execute(query)
    conn.commit()

'''
Verifica se esiste un Aliquote all'inteno del database
'''
def checkExistAliquote(id, cursor):
    cursor.execute("SELECT COUNT(*) FROM aliquote WHERE aliquote_id = %s", (id,))
    result = cursor.fetchone()
    if result[0] == 0: 
        return False
    return True


'''
Funzione che si occupa di inserire un nuovo analyte  nel database 
'''
def insertNewAnalyte(sample_id, analyte_id, concentration, cursor,conn):
    cursor.execute("INSERT INTO public.analyte VALUES (%s, %s, %s) ON CONFLICT (analyte_id) DO NOTHING;", ( sample_id, analyte_id, concentration,))
    conn.commit()


'''
Funzione che si occupa di inserire un nuovo portion  nel database 
'''
def insertNewPortion(sample_id, analyte_id, cursor,conn):
    cursor.execute("INSERT INTO public.portion VALUES (%s, %s) ON CONFLICT (portion_id) DO NOTHING;", ( sample_id, analyte_id,))
    conn.commit()

'''
Funzione che si occupa di inserire un nuovo aliquote nel database 
'''
def insertNewAliquote(sample_id, aliquote_id, concentration, cursor,conn):
    cursor.execute("INSERT INTO public.aliquote VALUES (%s, %s, %s) ON CONFLICT (aliquote_id) DO NOTHING;", ( sample_id, aliquote_id, concentration,))
    conn.commit()

'''
Verifica se esiste un progetto all'interno del database
'''
def checkExistProject(id, cursor):
    cursor.execute("SELECT COUNT(*) FROM project WHERE project_id = %s", (id,))
    result = cursor.fetchone()
    if result[0] == 0: 
        return False
    return True

'''
Funzione che ottiene il project id partendo dal nome 
'''
def getProjectId(name, cursor):
    cursor.execute("SELECT project_id FROM project WHERE name = %s", (name,))
    result = cursor.fetchall()
    if result == None or len(result) == 0:
         return None
    else:  return result[0][0]


'''
Funzione che ottiene il lo2_ratio di una riga da protein_pdc
'''
def getLog2Ratio(gene, project_id, aliq,cursor):
    cursor.execute("SELECT log2_ratio FROM protein_pdc WHERE gene_id = %s and aliquot = %s and project_id = %s", (gene, aliq, project_id,))
    result = cursor.fetchall()
    if result == None or len(result) == 0:
         return None
    else:  return result[0][0]


def getUnsharedLog2Ratio(gene, project_id, aliq,cursor):
    cursor.execute("SELECT unshared_log2_ratio FROM protein_pdc WHERE gene_id = %s and aliquot = %s and project_id = %s", (gene, aliq, project_id,))
    result = cursor.fetchall()
    if result == None or len(result) == 0:
         return None
    else:  return result[0][0]  

'''
Funzione che si occupa di inserire un nuovo progetto nel database 
'''
def insertNewProject(id, name, cursor,conn):
    query = "INSERT INTO public.project VALUES ('{}', '{}') ON CONFLICT (project_id) DO NOTHING;".format(id,name)
    cursor.execute(query)
    conn.commit()

'''
Verifica se esiste un caso all'interno del database
'''
def checkExistCase(id, cursor):
    query = "SELECT COUNT(*) FROM public.case WHERE case_id ='{}'".format(id)
    #print(query)
    cursor.execute(query)
    result = cursor.fetchone()
    if result[0] == 0: 
        return False
    return True


'''
Funzione che si occupa di inserire un nuovo case nel database 
'''
def insertNewCase(id, ethnicity, gender, race, vital_status, project, site, disease, cursor,conn):
    query = "INSERT INTO public.case VALUES ('{}', '{}', '{}', '{}', '{}', '{}', '{}', '{}') ON CONFLICT (case_id) DO NOTHING;".format(id, ethnicity, gender, race, vital_status, project, site, disease)
    #print(query)
    cursor.execute(query)
    conn.commit()

'''
Inserisce un nuovo sample all'interno del database
'''
def insertNewSample(sample_id, sample_type_id, tumor_code_id, cursor, conn):
    cursor.execute("INSERT INTO sample VALUES (%s, %s, %s) ON CONFLICT (sample_id) DO NOTHING;", (sample_id, sample_type_id, tumor_code_id))
    conn.commit()

'''
Verifica se esiste un file all'interno del database
'''
def checkExistFile(id, cursor):
    cursor.execute("SELECT COUNT(*) FROM analysis WHERE file_id = %s", (id,))
    result = cursor.fetchone()
    if result[0] == 0: 
        return False
    return True

'''
Cerca all'interno del database l'id relativo al type di un sample
'''
def searchSampleTypeId(type,cursor):
    q = "SELECT type_id FROM sample_type WHERE LOWER(type) = LOWER('{}')"
    cursor.execute(q.format(type))
    sample_type_id =  cursor.fetchone()
    return sample_type_id

'''
Inserisce un nuovo sample type all'interno del database
'''
def insertNewSampleType(type_id, type, cursor, conn):
    query = "INSERT INTO public.sample_type VALUES ('{}', '{}') ON CONFLICT (type_id) DO NOTHING;".format(type_id, type)
    cursor.execute(query)
    conn.commit()

'''
Verifica se esiste un tumor all'interno del database
'''
def checkExistTumor(id, cursor):
    cursor.execute("SELECT COUNT(*) FROM tumor WHERE tumor_code_id = %s", (id,))
    result = cursor.fetchone()
    if result[0] == 0: 
        return False
    return True

'''
Inserisce una nuova terna gene - proteina - studio all'interno del database
'''
def insertNewGeneProteinStudy(gene,protein,study,cursor,conn):
    query = "INSERT INTO public.protein_gene VALUES ('{}', '{}', '{}') ON CONFLICT (gene,study, protein) DO NOTHING;".format(gene, study, protein )
    #print(query)
    cursor.execute(query)
    conn.commit()

'''
Funzione che controlla se all'interno della table protein_pdc è già salvata una tupla uguale 
'''
def checkExistProtein_PDC(gene_id, project_id, aliquot,cursor):
    cursor.execute("SELECT COUNT(*) FROM protein_pdc WHERE gene_id = %s and aliquot = %s and project_id = %s", (gene_id,aliquot, project_id,))
    result = cursor.fetchone()
    if result[0] == 0: 
        return False
    return True


'''
Funzione che si occupa di inserire un nuovo tumore nel database 
'''
def insertNewTumor(id, tumor_code, tumor_description, cursor,conn):
    query = "INSERT INTO public.tumor VALUES ('{}', '{}', '{}') ON CONFLICT (tumor_code_id) DO NOTHING;".format(id, tumor_code, tumor_description)
    cursor.execute(query)
    conn.commit()

'''
Funzione che si occupa di verificare l'esistenza di un gene all'interno del database 
'''
def checkExistGene(id, cursor):
    cursor.execute("SELECT COUNT(*) FROM gene WHERE gene_id = %s", (id,))
    result = cursor.fetchone()
    if result[0] == 0: 
        return False
    return True

'''
funzione che si occupa di salvare un nuovo gene all'interno del database
'''
def insertNewGene(gene_id, gene_name, gene_type_id, cursor, conn):
    cursor.execute("INSERT INTO gene VALUES (%s, %s, %s) ON CONFLICT (gene_id) DO NOTHING;", (gene_id, gene_name , gene_type_id))
    conn.commit()
'''
Funzione che si occupa di controllare se un primary site è presente nel db, se non lo trova invoca la funzione che lo aggiunge, ottiene il suo nuovo id e lo restituisce 
'''
def getPrimarySite(data,cursor,conn):
    search_primary_site = "SELECT site_id FROM primary_site WHERE LOWER(site) = LOWER('{}')"
    cursor.execute(search_primary_site.format(data))
    primary_site =  cursor.fetchone()
    if primary_site == None:
            insertNewPrimarySite(data, cursor, conn)
    else:
         return primary_site[0]
    cursor.execute(search_primary_site.format(data)) 
    res = cursor.fetchone()
    return res[0]

'''
Inserisce un nuovo primary site all'interno del db 
'''
def insertNewPrimarySite(data, cursor, conn):
    quer = "INSERT INTO public.primary_site VALUES (DEFAULT, '{}')".format(data)
    cursor.execute(quer)
    conn.commit()

'''
Funzione che si occupa di controllare se un disease è presente nel db, se non lo trova invoca la funzione che lo aggiunge, ottiene il suo nuovo id e lo restituisce 
'''
def getDisease(data, cursor,conn):
    search_disease = "SELECT  disease_id FROM disease WHERE  LOWER(type) = LOWER('{}')"
    cursor.execute(search_disease.format(data))
    diase =  cursor.fetchone()
    if diase == None:
            insertNewDisease(data, cursor, conn)
    else:
         return diase[0]
    cursor.execute(search_disease.format(data)) 
    res = cursor.fetchone()
    return res[0]

'''
Inserisce un nuovo disease all'interno del db 
'''
def insertNewDisease(data, cursor, conn):
    quer = "INSERT INTO public.disease VALUES (DEFAULT, '{}')".format(data)
    cursor.execute(quer)
    conn.commit()


'''
Funzione che si occupa di controllare se una experimental strategy è presente nel db, se non la trova invoca la funzione che lo aggiunge, ottiene il suo nuovo id e lo restituisce 
'''
def getExperimentalStrategy(data, cursor,conn):
    search_exp = "SELECT strategy_id FROM experimental_strategy WHERE  LOWER(strategy) = LOWER('{}')"
    cursor.execute(search_exp.format(data))
    exp =  cursor.fetchone()
    if exp == None:
            insertNewExperimentalStrategy(data, cursor, conn)
    else:
         return exp[0]
    cursor.execute(search_exp.format(data)) 
    res = cursor.fetchone()
    return res[0]

'''
Inserisce una nuova experimental strategy all'interno del db 
'''
def insertNewExperimentalStrategy(data, cursor, conn):
    quer = "INSERT INTO public.experimental_strategy VALUES (DEFAULT, '{}')".format(data)
    cursor.execute(quer)
    conn.commit()


'''
Funzione che si occupa di controllare se una data category è presente nel db, se non lo trova invoca la funzione che la aggiunge, ottiene il suo nuovo id e lo restituisce 
'''
def getDataCategory(data, cursor,conn):
    search_cat = "SELECT category_id FROM data_category WHERE  LOWER(category) = LOWER('{}')"
    cursor.execute(search_cat.format(data))
    cat =  cursor.fetchone()
    if cat == None:
            insertNewDataCategory(data, cursor, conn)
    else:
         return cat[0]
    cursor.execute(search_cat.format(data)) 
    res = cursor.fetchone()
    return res[0]

'''
Inserisce una nuova data category all'interno del db 
'''
def insertNewDataCategory(data, cursor, conn):
    quer = "INSERT INTO public.data_category VALUES (DEFAULT, '{}')".format(data)
    cursor.execute(quer)
    conn.commit()


'''
Funzione che si occupa di controllare se una data type è presente nel db, se non lo trova invoca la funzione che la aggiunge, ottiene il suo nuovo id e lo restituisce 
'''
def getDataType(data, cursor,conn):
    search_data_type = "SELECT type_id FROM data_type WHERE  LOWER(type) = LOWER('{}')"
    cursor.execute(search_data_type.format(data))
    typ =  cursor.fetchone()
    if typ == None:
            insertNewDataType(data, cursor, conn)
    else:
         return typ[0]
    cursor.execute(search_data_type.format(data)) 
    res = cursor.fetchone()
    return res[0]

'''
Inserisce una nuova data type all'interno del db 
'''
def insertNewDataType(data, cursor, conn):
    quer = "INSERT INTO public.data_type VALUES (DEFAULT, '{}')".format(data)
    cursor.execute(quer)
    conn.commit()


'''
Funzione che si occupa di controllare se un gene type è presente nel db, se non lo trova invoca la funzione che la aggiunge, ottiene il suo nuovo id e lo restituisce 
'''
def getGeneType(data, cursor,conn):
    search_exp = "SELECT type_id FROM gene_type WHERE  LOWER(type) = LOWER('{}')"
    cursor.execute(search_exp.format(data))
    gene_type =  cursor.fetchone()
    if gene_type == None:
            insertNewGeneType(data, cursor, conn)
    else:
         return gene_type[0]
    cursor.execute(search_exp.format(data)) 
    res = cursor.fetchone()
    return res[0]

'''
Inserisce un nuovo gene type all'interno del db 
'''
def insertNewGeneType(data, cursor, conn):
    quer = "INSERT INTO public.gene_type VALUES (DEFAULT, '{}')".format(data)
    cursor.execute(quer)
    conn.commit()