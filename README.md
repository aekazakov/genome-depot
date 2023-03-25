# CGCMS
Compararive genomic content management system

## Prerequisites

1. Linux-based OS with installed Python 3.8+ and conda (miniconda, anaconda etc.).

2. MySQL server

3. Apache2 web-server with mod_wsgi

4. Muscle

5. HMMER

6. NCBI-BLAST

7. Perl

8. zlib1g-dev package

In Ubuntu-based distibutions, you can install the dependencies using the APT package repository:
sudo apt install apache2 mysql-server muscle hmmer ncbi-blast+ build-essential zlib1g-dev libexpat1-dev python3-dev libmysqlclient-dev

9. Jbrowse genome browser v. 1.16.11 (https://jbrowse.org/jbrowse1.html). You can find installation instructions here: https://jbrowse.org/docs/installation.html#download-jbrowse-minified-release


## Installation

1. Create a directory where CGCMS and external tools will be installed (for example, ~/cgcms). Create "app" subdirectory. Do it only for the first CGCMS installation.

cd ~

mkdir cgcms

cd cgcms

mkdir apps

cd apps


2. Create a directory for the new CGCMS installation in cgcms/apps (for example, mygenomes) and clone the repository into it.

mkdir mygenomes 

cd mygenomes

git clone https://github.com/aekazakov/CGCMS


3. Install external tools and create virtual environment.

cd CGSMS

bash install_tools.sh


4. Create mysql user (for example, cgcmsuser) or use existing mysql account.


5. Create database.

log into mysql as root: mysql -u root -p

CREATE DATABASE cgcms CHARACTER SET utf8;

GRANT ALL PRIVILEGES ON cgcms.* TO 'cgcmsuser'@'localhost';

quit


6. Open cgcms/apps/mygenomes/CGCMS/genomebrowser/secrets.json in a text editor. Enter secret key, hostname, database name (like cgcms), database username (like cgcmsuser), database password and URL to static files directory.


7. Open cgcms/apps/mygenomes/CGCMS/genomebrowser/configs.txt in a text editor. Check:

paths to Jbrowse utilites: prepare-refseqs.pl, flatfile-to-json.pl, generate-names.pl,

path to strain metadata file.


8. Activate virtual environment cgcms-venv, change directory to cgcms/app/mygenomes/CGCMS/genomebrowser and run Django configuration commands

source cgcms/cgcms-venv/bin/activate

cd cgcms/app/mygenomes/CGCMS/genomebrowser

python manage.py collectstatic

python manage.py makemigrations

python manage.py migrate

python manage.py createsuperuser (Enter your desired username, email and password).

python manage.py configure_cgcms -i configs.txt

python manage.py createcachetable


9. Run python manage.py runserver 127.0.0.1:8000 and open in the browser http://127.0.0.1:8000/admin. You should be able to log in with the username you entered at the previous step.

If test server cannot find static files, check if www-data user can read from the cgcms/static/my_genomes directory (giving rw permissions for www-data group would work).


10. Configure web-server for CGCMS. Open apache2 site configuration file (it may be default_ssl.conf) in a text editor and add the following (with correct paths):

	WSGIDaemonProcess cgcmspy python-home=/path/to/cgcms/cgcms-venv python-path=/path/to/cgcms/app/mygenomes/CGCMS/genomebrowser

	WSGIScriptAlias /mygenomes /path/to/cgcms/app/mygenomes/CGCMS/genomebrowser/genomebrowser/wsgi.py process-group=cgcmspy application-group=%{GLOBAL}

	<Directory /path/to/cgcms/app/mygenomes/CGCMS/genomebrowser/genomebrowser/>

	    <Files wsgi.py>

            Require all granted

	    </Files>

	</Directory>

	<Directory /path/to/cgcms/static/>

		Options -Indexes +FollowSymLinks

		<IfModule mod_headers.c>

		  Header set Access-Control-Allow-Origin http://127.0.0.1:8000

		</IfModule>

		Require all granted

	</Directory>

	Alias /cgcmsstatic /path/to/cgcms/static/

You may have to add "Header always set X-Frame-Options "SAMEORIGIN"" to web server configuration if the embedded genome viewer is not properly displayed.

11. Restart apache2:

sudo systemctl restart apache2

Now you would be able to open https://your.domain.name/mygenomes in a web browser.

Note: If POEM fails to predict operons, and run_poem.sh script throws an error "AttributeError: module 'tensorflow' has no attribute 'get_default_graph'", it means the versions of keras and tensorflow are not compatible. Activate conda cgcms-poem environment and run "pip install tensorflow==1.13.1". 

## Genome import from the command line

1. Download genomes in genbank format (files may be gzipped). Make a tab-separated file (for example, genomes.txt) with six columns:

- path to Genbank file

- genome ID (no spaces)

- strain name (no spaces)

- sample ID (no spaces)

- URL (link to NCBI genome assembly etc.)

- External ID (someting like "NCBI:GCF_000006945.2")


2. Import genomes into database, if input files have been uploaded to the server:

activate virtual environment (source /path/to/virtualenv/bin/activate)

cd cgcms/CGCMS/genomebrowser

python manage.py import_genomes -i genomes.txt

Input file 
Depending on the number of genomes, this command may run from several hours to several days. After that, you should see the genomes on the web site.

In the process of genome import, CGCMS annotation pipeline runs eggnog-mapper to generate EggNOG, KEGG, GO, EC, TC, CAZy and COG mappings for all proteins annotated in the input file. The pipeline predicts operons with POEM, maps Pfam domains with hmmsearch and generates functional annotations with several annotation tools. 


## Genome import from admin panel

Alternatively, genomes can be imported from admin panel. 

Log into siteURL/admin with username and password you created during installation process.

Click "Import genomes" button and follow the instructions on the page.


# Other commands

Note: activate virtual environment (source /path/to/virtualenv/bin/activate) before running any command, then change directory to cd cgcms/CGCMS/genomebrowser


1. Command to re-create Pfam and TIGRFAM domain mappings:

python manage.py update_domain_mappings -i genomes.txt


2. Command to generate functional annotations (if annotation pipeline failed or you added a new tool):

python manage.py update_annotations -i genomes.txt

3. Command to export genomes with CGCMS annotations in genbank format

python manage.py export_genomes -g <comma-separated list of genome names> -o <output directory>

4. Command to import regulon data

python manage.py add_regulons -i <path to input file>

The input file for this command is a tab-separated text file with nine columns:

- regulon name

- genome name

- locus tag of regulator 

- locus tag of target gene 

- contig id

- site start

- site end

- site strand (1 or -1)

- site sequence




# Image Credits

All images courtesy of Unsplash (https://unsplash.com).

- Jonathan Pie
    https://unsplash.com/@r3dmax?photo=3l3RwQdHRHg
    https://unsplash.com/@r3dmax?photo=1hpE3fROU0I
    https://unsplash.com/@r3dmax?photo=3N5ccOE3wGg
    https://unsplash.com/@r3dmax?photo=-3h8OXvt4-0
    https://unsplash.com/@r3dmax?photo=7FfG8zcPcXU
    https://unsplash.com/@r3dmax?photo=EvKBHBGgaUo
    https://unsplash.com/@r3dmax?photo=8I49k45G-3A
    https://unsplash.com/@r3dmax?photo=iokiwAq05UU
