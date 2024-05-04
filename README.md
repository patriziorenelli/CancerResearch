# CancerResearch

## Thesis project for a bachelor's degree in computer science at the "La Sapienza" University of Rome
The developed code allows to download, manipulate and archive data regarding various types of tumors.
All data is collected from open access repositories provided by the NIH (NATIONAL CANCER INSTITUTE) partners:
 * [GDC (GENOMIC DATA COMMONS)]( https://gdc.cancer.gov/)
 * [PDC (PROTEOMIC DATA COMMONS)](https://proteomic.datacommons.cancer.gov/pdc/)
The data can then be used to train machine learning models to promptly detect the onset of tumors

### How is the data downloaded?
To download the data, use the APIs offered by the two repositories.
The documentation of the APIs used:
 * [GDC (GENOMIC DATA COMMONS)]( https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/)
 * [PDC (PROTEOMIC DATA COMMONS)](https://pdc.cancer.gov/data-dictionary/publicapi-documentation/#!/Case/allCases)
Some data available on the PDC site cannot be downloaded via API because some APIs are still in development, 
webscraping was used to obtain all data not available via API.


### Technologies used - Requirements

> You must manually install all necessary software and libraries

<p>
    <img  src="https://github.com/patriziorenelli/CancerResearch/assets/19751923/3b487ca7-ade5-452f-aa2a-ccd92ada55b3" hspace="100" width="200" height="200">
    <img  src="https://github.com/patriziorenelli/CancerResearch/assets/19751923/13a63054-1597-45f5-8e81-5271629efb70"  width="250" height="250" >
</p>

All code was developed in Python
A Postgres database is used for local data storage.
For web scraping it is necessary to install chrome and the chromedriver compatible with the version of chrome installed

Various libraries have been used within the code including:
 *Requests: allows you to make HTTP requests to the API to obtain the necessary data,
 *Json: allows you to convert the data obtained from the API into JSON format and work with this type of data,
 *Psycopg2: allows us to establish a connection with the PostgresSQL database and subsequently carry out all the operations we need on the db,
 *Pandas: allows you to manipulate and analyze the data obtained from the API, creating specific data structures,
 *Os: allows you to interact with the operating system, for managing files and executing commands,
 *Datetime: allows you to obtain and manipulate the date and time,
 *Pathlib: allows the manipulation and construction of paths suitable for filesystems,
 *Time: allows you to obtain and manipulate date/time and provides us with the sleep() function used to insert a pause during program execution,
 *BeautifulSoup: allows you to extract data from HTML files by creating analysis trees that allow the analysis of the HTML code obtained through web scraping,
 *Selenium: allows the automation of a web browser to perform web scraping


### Code structure

The code is structured in various files based on the functionality they perform. 
We then find the files:
 *DatabaseGenericInteraction: contains the various functions that perform INSERT, UPDATE and SELECT queries in the local database
 *DatabaseManager:contains the various functions that make the connection with the database. 
  It also allows you to:
    *automatically build the local database and create all the necessary tables
    *make a local backup of the database
    *load the database from a previously performed backup
 *GDC_DownloadData:contains the various functions that allow you to download the various data available on the GDC repository
 *GetGeneProtInformation: contains the various functions that allow you to download the values of the various peptides for each portion associated with the genes through web scraping
 *GetLog2GeneAliquot: contains the various functions that allow you to download and archive the values of the various log2_ratio and unshared_log2_ratio fields
 *params: contains the parameters necessary to connect to the local Postgres database
 *PDC_DownloadData: contains functions that allow you to download all data available on the PDC repository
 
### Database structure
![image](https://github.com/patriziorenelli/CancerResearch/assets/19751923/5356e0d4-bcb0-45be-a05f-d2922dd73cd5)

The code allows you to download the data offered by the two repositories (GDC - PDC) trying to convert them into a common standard, 
the database structure is designed to contain all the various data at the same time.
There are tables that contain data from both repositories and other tables that contain specific data from a single repository.
For the specific structure, the description of all the tables and their fields, 
you must refer to the complete documentation of the project (the documentation is written in Italian and automatically translated into English )






