# GenomeDepot
Data management system for microbial comparative genomics.

## Prerequisites

1. Linux-based OS with installed Python 3.8+ and conda (miniconda, anaconda etc.).

2. MySQL server

3. Apache2 web-server with mod_wsgi

4. Muscle

5. HMMER

6. NCBI-BLAST (blastp and megablast required)

7. Perl

8. zlib1g-dev package

9. curl

10. pkg-config

In Ubuntu-based distibutions, you can install the dependencies using the APT package repository:
sudo apt install apache2 mysql-server muscle hmmer ncbi-blast+-legacy build-essential zlib1g-dev libexpat1-dev curl pkg-config

9. Jbrowse genome browser v. 1.16.11 (https://jbrowse.org/jbrowse1.html). You can find installation instructions here: https://jbrowse.org/docs/installation.html#download-jbrowse-minified-release


## Installation

1. Create a directory where GenomeDepot and external tools will be installed (for example, ~/genomedepot). Create "app" subdirectory. Do it only for the first GenomeDepot installation.

cd ~

mkdir genomedepot

cd genomedepot

mkdir apps

cd apps


2. Create a directory for the new GenomeDepot installation in genomedepot/apps (for example, mygenomes) and clone the repository into it.

mkdir mygenomes 

cd mygenomes

git clone https://github.com/aekazakov/genome-depot


3. Install external tools and create virtual environment.

cd genomedepot

bash install.sh


4. Create mysql user (for example, gduser) or use existing mysql account.


5. Create database.

log into mysql as root: mysql -u root -p

CREATE DATABASE genomedepot CHARACTER SET utf8;

GRANT ALL PRIVILEGES ON genomedepot.* TO 'gduser'@'localhost';

quit


6. Open genomedepot/apps/mygenomes/genome-depot/genomebrowser/secrets.json in a text editor. Enter secret key, hostname, database name (like genomedepot), database username (like gduser), database password and URL to static files directory.

7. Open genomedepot/apps/mygenomes/genome-depot/genomebrowser/static/jbrowse/jbrowse.conf in a text editor. Find documentDomain parameter, uncomment it and enter hostname.

8. Open genomedepot/apps/mygenomes/genome-depot/genomebrowser/configs.txt in a text editor. Check:

paths to Jbrowse utilites: prepare-refseqs.pl, flatfile-to-json.pl, generate-names.pl,

path to strain metadata file.


9. Activate virtual environment genomedepot-venv, change directory to genomedepot/app/mygenomes/genome-depot/genomebrowser and run Django configuration commands

source genomedepot/genomedepot-venv/bin/activate

cd genomedepot/app/mygenomes/genome-depot/genomebrowser

python manage.py collectstatic

python manage.py makemigrations

python manage.py migrate

python manage.py createsuperuser (Enter your desired username, email and password).

python manage.py import_config -i configs.txt

python manage.py createcachetable

python manage.py update_taxonomy


9. Run python manage.py runserver 127.0.0.1:8000 and open in the browser http://127.0.0.1:8000/admin. You should be able to log in with the username you entered at the previous step.

If test server cannot find static files, check if www-data user can read from the genomedepot/static/my_genomes directory (giving rw permissions for www-data group would work).


10. Configure web-server for GenomeDepot. Open apache2 site configuration file (it may be default_ssl.conf) in a text editor and add the following (with correct paths):

	WSGIDaemonProcess genomedepotpy python-home=/path/to/genomedepot/genomedepot-venv python-path=/path/to/genomedepot/app/mygenomes/genome-depot/genomebrowser

	WSGIScriptAlias /mygenomes /path/to/genomedepot/app/mygenomes/genome-depot/genomebrowser/genomebrowser/wsgi.py process-group=genomedepotpy application-group=%{GLOBAL}

	<Directory /path/to/genomedepot/app/mygenomes/genome-depot/genomebrowser/genomebrowser/>

	    <Files wsgi.py>

            Require all granted

	    </Files>

	</Directory>

	<Directory /path/to/genomedepot/static/>

		Options -Indexes +FollowSymLinks

		<IfModule mod_headers.c>

		  Header set Access-Control-Allow-Origin http://127.0.0.1:8000

		</IfModule>

		Require all granted

	</Directory>

	Alias /genomedepotstatic /path/to/genomedepot/static/

You may have to add "Header always set X-Frame-Options "SAMEORIGIN"" to web server configuration if the embedded genome viewer is not properly displayed.

11. Restart apache2:

sudo systemctl restart apache2

Now you would be able to open https://your.domain.name/mygenomes in a web browser.

Note: If POEM fails to predict operons, and run_poem.sh script throws an error "AttributeError: module 'tensorflow' has no attribute 'get_default_graph'", it means the versions of keras and tensorflow are not compatible. Activate conda genomedepot-poem environment and run "pip install tensorflow==1.13.1". 

## Genome import from the command line

1. Download genomes in genbank format (files may be gzipped). Make a tab-separated file (for example, genomes.txt) with six columns:

- path to Genbank file

- genome ID (no spaces)

- strain name (no spaces)

- sample ID (no spaces)

- URL (link to NCBI genome assembly etc.)

- External ID (someting like "NCBI:GCF_000006945.2")


2. Import genomes into database. If genome files have been uploaded to the server, they can be imported from the command line:

activate virtual environment (source /path/to/virtualenv/bin/activate)

cd genomedepot/apps/mygenomes/genome-depot/genomebrowser

python manage.py import_genomes -i genomes.txt

Input file 
Depending on the number of genomes, this command may run from several hours to several days. After that, you should see the genomes on the web site.

In the process of genome import, GenomeDepot annotation pipeline runs eggnog-mapper to generate EggNOG, KEGG, GO, EC, TC, CAZy and COG mappings for all proteins annotated in the input file. The pipeline predicts operons with POEM, maps Pfam domains with hmmsearch and generates functional annotations with several annotation tools. 


## Genome import from admin panel

Alternatively, genomes can be imported from admin panel. 

Start qcluster process:

activate virtual environment (source /path/to/virtualenv/bin/activate)

cd genomedepot/apps/mygenomes/genome-depot/genomebrowser

python manage.py qcluster

Log into siteURL/admin with superuser login and password you created during installation process.

Click "Import genomes" button and follow the instructions on the page.

Tab-separated text file is always required.

Zip-archive with genomes may be provided, if genomes haven't been uploaded to server.

For genome download from NCBI ftp, e-mail address should be provided. For security reasons, this e-mail is never stored in the database.

# Other commands

Note: activate virtual environment (source /path/to/virtualenv/bin/activate) before running any command, then change directory to  genomedepot/apps/mygenomes/genome-depot/genomebrowser


1. Command to re-create Pfam and TIGRFAM domain mappings:

python manage.py update_domain_mappings -i genomes.txt


2. Command to generate functional annotations (if annotation pipeline failed or you added a new tool):

python manage.py update_annotations -i genomes.txt

3. Command to export genomes with GenomeDepot annotations in genbank format

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

# Import test dataset of RefSeq genomes

You can test GenomDepot installation by importing a dataset of 25 RefSeq genome assemblies. 

In the GenomeDepot administration panel, go to the Genome Import page, choose option 3 (import from NCBI).
Click "Browse" and choose the testdata/demo_25genomes_import.txt file. Enter your email and click "Start import".


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
