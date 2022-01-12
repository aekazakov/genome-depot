# CGCMS
Compararive genomic content management system

## Prerequisites

1. Linux-based OS with installed Python 3.8+.

2. MySQL

3. Apache2

4. Django 

5. eggnog-mapper v. 2.0.1 (https://github.com/eggnogdb/eggnog-mapper/tree/2.0.1) with reference databases. 

6. HMMER

7. Jbrowse genome browser v. 1. (https://jbrowse.org/jbrowse1.html)

8. POEM_py3 (https://github.com/Rinoahu/POEM_py3k)


## Installation

1. Create a directory and clone repository into it;

mkdir cgcms

cd cgcms

git clone https://github.com/aekazakov/CGCMS

2. Create mysql user (for example, cgcmsuser) or use existing mysql account.

3. Create database.

log into mysql as root: mysql -u root -p

CREATE DATABASE cgcms CHARACTER SET utf8;

GRANT ALL PRIVILEGES ON cgcms.* TO 'cgcmsuser'@'localhost';

quit

4. Copy cgcms/CGCMS/genomebrowser/secrets.json.template to cgcms/CGCMS/genomebrowser/secrets.json and edit it in a text editor. Set secret key, hostname, database name (like cgcms), database username (like cgcmsuser), database password and path to directory with static files.

5. Create a directory for static files and copy static files from cgcms/CGCMS/static into it. Make it accessible for webserver (giving rw permissions for www-data group would work).

6. [Optional] Install other annotation tools:

- amrfinderplus (https://github.com/ncbi/amr)

- antismash (https://github.com/antismash/antismash)

- ecis-screen (https://github.com/ipb-jianyang/eCIS-screen)

- fama (https://github.com/novichkov-lab/fama)

- phispy (https://github.com/linsalrob/PhiSpy)

7. Copy cgcms/CGCMS/genomebrowser/configs.txt.template to cgcms/CGCMS/genomebrowser/configs.txt and edit it in a text editor. Set parameters:

cgcms.temp_dir: directory for temporary files (for example cgcms/tmp),

cgcms.eggnog-command: path to emapper.py, 

cgcms.eggnog_outdir: output directory for eggnog-mapper (for example cgcms/tmp/eggnog), 

paths to Jbrowse utilites: prepare-refseqs.pl, flatfile-to-json.pl, generate-names.pl,

paths to reference files,

paths and commands for plugins (optional). 

Delete all unnecessary plugin entries.

8. Create virtual environment for Django.

9. Activate virtual environment, change directory to cgcms/CGCMS/genomebrowser and run

python manage.py makemigrations

python manage.py migrate

python manage.py createsuperuser (Enter your desired username, email and password).

python manage.py configure_cgcms -i configs.txt

10. Run python manage.py runserver 127.0.0.1:8000 and open in the browser http://127.0.0.1:8000/admin. You should be able to log in with the username you entered at the previous step.

11. Edit apache2 site configuration file (it may be default_ssl.conf). Add the following (with correct paths):

	WSGIDaemonProcess cgcmspy python-home=/path/to/virtual/environment python-path=/path/to/cgcms/CGCMS/genomebrowser

	WSGIScriptAlias /genomes /path/to/cgcms/CGCMS/genomebrowser/genomebrowser/wsgi.py process-group=cgcmspy application-group=%{GLOBAL}

	<Directory /path/to/cgcms/CGCMS/genomebrowser/genomebrowser/>

	    <Files wsgi.py>

            Require all granted

	    </Files>

	</Directory>

	<Directory /path/to/static/files>

		Options -Indexes +FollowSymLinks

		<IfModule mod_headers.c>

		  Header set Access-Control-Allow-Origin http://127.0.0.1:8000

		</IfModule>

		Require all granted

	</Directory>

	Alias /static /path/to/static/files

	
12. Restart apache2:

sudo systemctl restart apache2

13. Test if you can open your.domain.name/genomes.


## Genome import

1. Download genomes in genbank format (files may be gzipped). Make a tab-separated file (for example, genomes.txt) with six columns:

- path to Genbank file

- genome ID (no spaces)

- strain name (no spaces)

- sample ID (no spaces)

- URL (link to NCBI genome assembly etc.)

- External ID (someting like "NCBI:GCF_000006945.2")


2. Import genomes into database:

activate virtual environment (source /path/to/virtualenv/bin/activate)

cd cgcms/CGCMS/genomebrowser

python manage.py import_genomes -i genomes.txt

Depending on the number of genomes, this command may run from several hours to several days. After that, you should see the genomes on the web site.


3. Generate Pfam and TIGRFAM domain mappings (optional):

python manage.py update_domain_mappings -i genomes.txt


4. If have any plug-ins configured, generate additional annotations (optional):

python manage.py update_annotations -i genomes.txt

