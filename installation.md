[About](introduction) | **Installation** | [User guide](user) | [Administration guide](admin) | [Developer guide](developer)

# Installation

This document uses the /opt/genomedepot directory as an example, but GenomeDepot can be installed in another directory as well. **Installation in the /opt/genomedepot directory requires sudo rights.** 

Installing GenomeDepot in the user’s home directory **is not recommended** because of possible issues with file access permissions and web app execution by the Apache web server.

In this document, https://example.com/mygenomes is used as a base URL for the GenomeDepot-based web portal. For your installation, use your own domain name instead of example.com.


## Prerequisites

* Linux-based OS (GenomeDepot was developed and tested in 64-bit Ubuntu Linux system)
* Python 3.8+ 
* conda (miniconda is recommended)
* MySQL server
* Apache2 web-server with mod_wsgi and ssl extensions
* Muscle
* HMMER
* NCBI-BLAST (blastp and megablast required)
* Perl
* zlib1g-dev package
* curl
* git
* pkg-config

In Ubuntu-based distibutions, you can install most prerequisites using the APT package manager:
```
sudo apt install apache2 mysql-server muscle hmmer ncbi-blast+-legacy build-essential zlib1g-dev libexpat1-dev python3-dev libmysqlclient-dev curl git pkg-config libapache2-mod-wsgi-py3 python3.10-venv
```

If you don't have conda, install Miniconda [as described in the Conda user guide]( https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html).


## Install dependencies

Create a directory where GenomeDepot and external tools will be installed. 

```
cd /opt
sudo mkdir genomedepot
```
If you created the genomedepot directory with sudo, change the directory ownership to a non-root user (replace USERNAME with a real user name):
```
sudo chown USERNAME:www-data genomedepot
```
Create apps subdirectory. In this directory, GenomeDepot-based portals will be installed. 
```
cd genomedepot
mkdir apps
cd apps
```
Create a directory for a new GenomeDepot portal in genomedepot/apps (for example, mygenomes) and clone the repository into it.
```
mkdir mygenomes
cd mygenomes
git clone https://github.com/aekazakov/genome-depot
```
Run the install.sh script, which installs external tools and creates a Python virtual environment for GenomeDepot. Running the script may take quite some time for the deployment of the first GenomeDepot-based portal because it will create all the conda environments and download reference data for the tools.
```
cd genome-depot
bash install.sh
```
The install.sh script installs all required Python libraries including Django framework, genome annotation tools and other dependencies. If it fails, check the error message, fix the problem and start install.sh again.
The most common problem is a name conflict with existing conda environments. In this case, the installation script does not overwrite the existing environment but stops the installation. User can rename or delete the existing environment before restarting the installation.

Another common problem is an incomplete installation of the operon prediction tool POEM. If POEM fails to predict operons, and running poem.sh script throws an error "AttributeError: module 'tensorflow' has no attribute 'get_default_graph'", it means the versions of keras and tensorflow are not compatible. Activate conda genomedepot-poem environment and run `pip install tensorflow==1.13.1`.

## Initial GenomeDepot configuration

Create MySQL user (for example, gduser), if needed, or use an existing mysql account.

Create MySQL database (for example, gdgenomes):
log into mysql as root:
```
mysql -u root -p
```

If MySQL returns error "Access denied for user 'root'@'localhost'", use sudo before the command:
```
mysql -u root -p

```

Enter commands that create the database and grant the user access to the database:
```
CREATE DATABASE gdgenomes CHARACTER SET utf8;
GRANT ALL PRIVILEGES ON gdgenomes.* TO 'gduser'@'localhost';
quit
```
Copy the genomedepot/apps/mygenomes/genome-depot/genomebrowser/.env.template file to genomedepot/apps/mygenomes/genome-depot/genomebrowser/.env. Open genomedepot/apps/mygenomes/genome-depot/genomebrowser/.env in a text editor and enter the settings:

* SECRET_KEY: Django [secret key](https://docs.djangoproject.com/en/5.1/ref/settings/#secret-key)
* ALLOWED_HOSTS: comma-separated list of host names and IP addresses for the web server (for example,  example.com,127.0.0.1,testserver)
* INTERNAL_IPS: comma-separated list of IP addresses for this site (usually, 127.0.0.1,)
* DB_USER: mysql user name (created at step 1)
* DB_PASSWORD: mysql password (created at step 1)
* DB_NAME: mysql database name (created at step 2)
* STATIC_URL: URL for static files directory (for example, https://example.com/gdstatic/mygenomes/)
* TITLE: web site title 
* BASE_URL: URL of the web site. It can be a domain name (https://example.com/) or a subdirectory (https://example.com/mygenomes/)
* ADMIN: admin name and email separated by comma (John Doe,johndoe@example.com)
* EMAIL_HOST_USER: an email address for GenomeDepot to send messages from (mygenomes@example.com)
* EMAIL_HOST_PASSWORD: a password for the EMAIL_HOST_USER email 
* EMAIL_HOST: external e-mail server address
* STATIC_ROOT: full path of a directory for static files from which the web server serves up STATIC_URL. For example, /opt/genomedepot/static/mygenomes
* STATICFILES_DIR: a directory with static files in the GenomeDepot installation. For example, /opt/genomedepot/apps/mygenomes/genome-depot/genomebrowser/static
* LOGVIEWER_LOGS: full path to the django.log file in the genomebrowser directory, followed by comma. For example, /opt/genomedepot/apps/mygenomes/genome-depot/genomebrowser/django.log,

Open /opt/genomedepot/apps/mygenomes/genome-depot/genomebrowser/configs.txt in a text editor. This file is generated by the GenomeDepot installation script and contains GenomeDepot configuration parameters to be saved in the database. Check the parameters for correctness.

Activate virtual environment genomedepot-venv, change directory to /opt/genomedepot/app/mygenomes/genome-depot/genomebrowser and run Django configuration commands
```
source /opt/genomedepot/genomedepot-venv/bin/activate
cd /opt/genomedepot/app/mygenomes/genome-depot/genomebrowser
python manage.py collectstatic
python manage.py makemigrations
python manage.py migrate
python manage.py createsuperuser 
# Enter your desired username, email and password
python manage.py import_config -i configs.txt
python manage.py createcachetable
```

Now your new GenomeDepot portal is ready for testing. Run 
```
python manage.py runserver 127.0.0.1:8000
```
and open in the browser <http://127.0.0.1:8000/admin>. You should be able to log in with the superuser username and password you just created.
 
If test server cannot find static files, check if www-data user can read from the /opt/genomedepot/static/mygenomes directory (giving rw permissions for www-data group would work).

## How to configure the Apache web server for GenomeDepot

You need sudo privileges for making changes in Apache configuration files. Open apache2 site configuration file (it may be /etc/apache2/sites_available/default_ssl.conf) in a text editor and add the following (change file paths, if needed):
```
WSGIDaemonProcess genomedepotpy python-home=/opt/genomedepot/genomedepot-venv python-path=/opt/genomedepot/app/mygenomes/genome-depot/genomebrowser
WSGIScriptAlias /mygenomes /opt/genomedepot/app/mygenomes/genome-depot/genomebrowser/genomebrowser/wsgi.py process-group=genomedepotpy application-group=%{GLOBAL}
<Directory /opt/genomedepot/app/mygenomes/genome-depot/genomebrowser/genomebrowser/>
<Files wsgi.py>

    Require all granted

</Files>
</Directory>
<Directory /opt/genomedepot/static/>
Options -Indexes +FollowSymLinks

<IfModule mod_headers.c>

  Header set Access-Control-Allow-Origin http://127.0.0.1:8000

</IfModule>

Require all granted
Alias /gdstatic /opt/genomedepot/static/
</Directory>
```
You may have to add "Header always set X-Frame-Options "SAMEORIGIN"" to web server configuration if the embedded genome viewer is not properly displayed.

Change group ownership to www-data for the genome-depot/genomebrowser/django.log file. For example:

```
sudo chown :www-data /opt/genomedepot/apps/mygenomes/genome-depot/genomebrowser/django.log
```

Restart apache2:
```
sudo systemctl restart apache2
```
Now you would be able to open https://your.domain.name/mygenomes in a web browser.

## Start Django Q cluster from the command line

A Django Q cluster must be started for the execution of GenomeDepot pipeline jobs from the task queue. To start the cluster, start a new screen session with “screen” command and enter:
```
conda deactivate
source /<GenomeDepot venv path>/bin/activate
cd /<GenomeDepot path>/genomebrowser
python manage.py qcluster
```
Optionally, the screen session can be named. Press Ctrl-A, then type “:sessionname “ and type a name.

After that, detach the screen session (press Ctrl-A, then d).

If the cluster is up and running, the number of clusters in the administration panel changes from 0 to 1.
There is no need to run more than one cluster for a GenomeDepot-based portal, since only one task can be processed at any moment. The GenomeDepot pipeline does not let parallel processing of tasks to prevent concurrent modification of the data.

To stop the Django Q cluster, attach the screen session with “screen -r” command, and press Ctrl-C.

## Start Django Q cluster from crontab

The Django Q cluster can be started on server restart as a cron job. First, create a bash script that starts the cluster:

```
#!/usr/bin/bash
conda deactivate
source /<GenomeDepot venv path>/bin/activate
cd /<GenomeDepot path>/genomebrowser
python manage.py qcluster
```
Then, open cron file (crontab -e) and add a line to the end:
```
@reboot(. ~/.profile /usr/bin/screen -dmS gdcluster <path to the shell script>)
```
Save the file and exit the editor.
To stop a cluster, run “screen -r gdcluster” command and press Ctrl-C. Warning: if the cluster is running a task, it may not stop immediately, and some data may be lost or corrupted.
If you have more than one GenomeDepot-based portal, change “gdcluster” to a unique name for each screen session. 

## How to change the background image of your web-site

There are two page background images in the repository, one for the dark mode (genomebrowser/static/images/background.jpg) and the other for the light mode (genomebrowser/static/images/background_light.jpg). 
The background image on the start page is genomebrowser/static/images/slide02.jpg.
You can replace any of these files with your own images, then run `python manage.py collectstatic` command and restart the Apache web server.

To keep git from overwriting image files during future updates, run:
```
git update-index --skip-worktree genomebrowser/static/images/slide02.jpg
git update-index --skip-worktree genomebrowser/static/images/background.jpg
git update-index --skip-worktree genomebrowser/static/images/background_light.jpg
```

## How to change links at the bottom of web pages

Edit the genomebrowser/browser/templates/footer_links.html file.
To replace icons at the links, you can choose from [Font Awesome 4.0.3 Icons set](https://btsai.github.io/font_awesome4_cheatsheet/index.html).

To keep git from  overwriting the file during future updates, run:
```
git update-index --skip-worktree genomebrowser/browser/templates/footer_links.html
```


## Tweak MySQL configuration 

As the genome database grows, MySQL performance may decrease because of excessive I/O usage. To increase MySQL performance, change innodb_buffer_pool_size parameter in MySQL configuration file. Run the following query in the mysql window to find optimal innodb_buffer_pool_size:
```
SELECT CEILING(Total_InnoDB_Bytes*1.6/POWER(1024,3)) RIBPS FROM (SELECT SUM(data_length+index_length) Total_InnoDB_Bytes FROM information_schema.tables WHERE engine='InnoDB') A;
```
It will show RIBPS value in gigabytes. If your server has enough memory, enter that number into the mysqld section of mysqld.conf (/etc/mysql/mysql.conf.d/mysqld.conf) and add “G”, for example:

```
[mysqld]
innodb_buffer_pool_size=8G
```

After making the changes, restart MySQL:
```
sudo systemctl restart mysql
```


[Continue to user guide...](user)

[About](introduction) | **Installation** | [User guide](user) | [Administration guide](admin) | [Developer guide](developer)
