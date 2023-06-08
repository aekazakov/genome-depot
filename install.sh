#!/bin/bash
WORKDIR=$(pwd)
CGCMSDIR=$(dirname "$WORKDIR")
APPNAME=$(basename "$CGCMSDIR")
CGCMSDIR=$(dirname "$CGCMSDIR")
CGCMSDIR=$(dirname "$CGCMSDIR")
# Deactivate all conda environments (including base) because python virtual environment must not use python installed in conda
conda deactivate
conda deactivate
if ! [ -d "$CGCMSDIR/cgcms-venv" ]; then
	python3 -m venv "$CGCMSDIR/cgcms-venv"
fi

CONDADIR=$(dirname "$(which conda)")
CONDADIR=$(dirname "$CONDADIR")
CONDA="$CONDADIR/etc/profile.d/conda.sh"
if ! [ -f "$CONDA" ]; then
		echo "Conda is not installed. Usage: setup.sh <conda dir>"
		exit
fi
echo "Found Conda installation: $CONDADIR"
#Installing dependencies in the virtual environment
echo "Installing python dependencies in cgcms-venv virtual environment"
source "$CGCMSDIR/cgcms-venv/bin/activate"
pip install "django==3.2.6" django_admin_shortcuts django_cors_headers django_q django_debug_toolbar openpyxl "parasail==1.2.4" biopython toytree urllib3 mysqlclient --no-cache-dir
deactivate

source $CONDA
conda config --add channels bioconda
conda config --add channels conda-forge
# Create directory structure
[ -d "$CGCMSDIR/external_tools" ] || mkdir "$CGCMSDIR/external_tools"
[ -d "$CGCMSDIR/external_refdata" ] || mkdir "$CGCMSDIR/external_refdata"
[ -d "$CGCMSDIR/static" ] || mkdir "$CGCMSDIR/static"
[ -d "$CGCMSDIR/static/$APPNAME" ] || mkdir "$CGCMSDIR/static/$APPNAME"
[ -d "$CGCMSDIR/static/$APPNAME/genomes" ] || mkdir "$CGCMSDIR/static/$APPNAME/genomes"
[ -d "$CGCMSDIR/static/$APPNAME/genomes/json" ] || mkdir "$CGCMSDIR/static/$APPNAME/genomes/json"
[ -d "$CGCMSDIR/static/$APPNAME/genomes/gbff" ] || mkdir "$CGCMSDIR/static/$APPNAME/genomes/gbff"
[ -d "$WORKDIR/appdata" ] || mkdir "$WORKDIR/appdata"
[ -d "$WORKDIR/temp" ] || mkdir "$WORKDIR/temp"
[ -d "$WORKDIR/temp/eggnog" ] || mkdir "$WORKDIR/temp/eggnog"

# Download taxonomy data
if ! [ -f "$WORKDIR/ref_data/ref_taxonomy.txt" ]; then
	echo "Downloading taxonomy data"
	cd "$WORKDIR/ref_data"
    curl -LJO -q http://iseq.lbl.gov/mydocs/cgcms_downloads/ref_taxonomy.txt.gz
	gunzip ref_taxonomy.txt.gz
fi

# Download JBrowse
cd "$CGCMSDIR/external_tools"
if ! [ -d "$CGCMSDIR/external_tools/jbrowse" ]; then
    # Install Jbrowse v.1.16.11
    curl -L -O https://github.com/GMOD/jbrowse/releases/download/1.16.11-release/JBrowse-1.16.11.zip
    unzip JBrowse-1.16.11.zip
    mv JBrowse-1.16.11 jbrowse
    rm JBrowse-1.16.11.zip
    cd jbrowse
    ./setup.sh
    cd "$CGCMSDIR/external_tools"
fi

# Install eggnog-mapper
if conda env list | grep 'cgcms-emapper' >/dev/null 2>&1; then
	conda activate cgcms-emapper
	if ! { emapper.py -v --data_dir "$CGCMSDIR/external_refdata/eggnog-mapper_v2.1.7" | grep 'emapper-2.1.7'; }>/dev/null 2>&1; then
		echo 'Conda environment cgcms-emapper exists but eggnog-mapper v2.1.7 was not properly installed. Remove the environment and restart CGCMS installation script.'
		echo 'To remove the environment, run:'
		echo '   conda remove -n cgcms-emapper --all'
		conda deactivate
		exit 1
	else
		echo 'EggNOG-mapper v.2.1.7 found'
	fi
	conda deactivate
else
	echo "Installing EggNOG mapper"
	conda create -y -n cgcms-emapper python=3.8
	conda activate cgcms-emapper
	conda install -c bioconda -y eggnog-mapper=2.1.7
	conda deactivate
	cd "$CGCMSDIR/external_tools"
fi
# Install eggnog-mapper reference data
if ! [ -d "$CGCMSDIR/external_refdata/eggnog-mapper_v2.1.7" ]; then
	conda activate cgcms-emapper
	mkdir "$CGCMSDIR/external_refdata/eggnog-mapper_v2.1.7"
	echo "Downloading reference databases for eggnog-mapper_v2.1.7. Answer \"y\" to all questions."
	download_eggnog_data.py -y --data_dir "$CGCMSDIR/external_refdata/eggnog-mapper_v2.1.7"
	conda deactivate
fi
# Install poem_py3
if conda env list | grep 'cgcms-poem' >/dev/null 2>&1; then
	echo "Found cgcms-poem environment"
	conda activate cgcms-poem
	if [ -f "$CGCMSDIR/external_tools/POEM_py3k/bin/run_poem_cgcms.sh" ]
	then
		echo "POEM_py3k found"
	else
		echo 'Conda environment cgcms-poem exists but POEM_py3k was not properly installed. Remove the environment and restart CGCMS installation script.'
		echo 'To remove the environment, run:'
		echo '   conda remove -n cgcms-poem --all'
		conda deactivate
		exit 1
	fi
	conda deactivate
else
	echo "Installing POEM_py3"
	conda create -y -n cgcms-poem python=3.7
	conda activate cgcms-poem
	git clone https://github.com/aekazakov/POEM_py3k.git
	cd ./POEM_py3k
	bash ./install.sh cgcms-poem
	conda deactivate
	cd "$CGCMSDIR/external_tools"
fi
#Install amrfinder
if conda env list | grep 'cgcms-amrfinder' >/dev/null 2>&1; then
	echo "Found cgcms-amrfinder environment"
	conda activate cgcms-amrfinder
	if amrfinder --version|grep "3.11.11" >/dev/null 2>&1; then
		echo "AMRFinder found"
	else
		echo 'Conda environment cgcms-amrfinder exists but AMRFinder was not properly installed. Remove the environment and restart CGCMS installation script.'
		echo 'To remove the environment, run:'
		echo '   conda remove -n cgcms-amrfinder --all'
		conda deactivate
		exit 1
	fi
	conda deactivate
else
	echo "Installing AMRFinder"
	conda create -y -n cgcms-amrfinder python=3.8
	conda activate cgcms-amrfinder
	conda install -y -c bioconda -c conda-forge ncbi-amrfinderplus
	amrfinder -u
	conda deactivate
	cd "$CGCMSDIR/external_tools"
fi
# Install antismash
if conda env list | grep 'cgcms-antismash' >/dev/null 2>&1; then
    echo "Found cgcms-antismash environment"
    conda activate cgcms-antismash
    if antismash -h|grep "antiSMASH 6.1.1" >/dev/null 2>&1;     then
		echo "antiSMASH 6.1.1 found"
    else
		echo 'Conda environment cgcms-antismash exists but antiSMASH was not properly installed. Remove the environment and restart CGCMS installation script.'
		echo 'To remove the environment, run:'
		echo '   conda remove -n cgcms-antismash --all'
		conda deactivate
		exit 1
	fi
	conda deactivate
else
	echo "Installing AntiSMASH"
	conda create -y -n cgcms-antismash python=3.8
	conda activate cgcms-antismash
	conda install -y antismash
	conda deactivate
	cd "$CGCMSDIR/external_tools"
fi
# Install ecis-screen
if conda env list | grep 'cgcms-ecis-screen' >/dev/null 2>&1; then
	echo "Found cgcms-ecis-screen environment"
	conda activate cgcms-ecis-screen
    if $CGCMSDIR/external_tools/eCIS-screen/HMMsearch_genomesII_fast.pl 2>&1 | grep "Pipeline for screening" >/dev/null 2>&1
	then
	echo "eCIS-screen found"
	else
	echo 'Conda environment cgcms-ecis-screen exists but eCIS-screen was not properly installed. Remove the environment and restart CGCMS installation script.'
		echo 'To remove the environment, run:'
		echo '   conda remove -n cgcms-ecis-screen --all'
	conda deactivate
	exit 1
	fi
	conda deactivate
else
	echo "Installing eCIS-screen"
	conda create -y -n cgcms-ecis-screen
	conda activate cgcms-ecis-screen
	conda install -y perl-bioperl hmmer
	git clone https://github.com/aekazakov/eCIS-screen
	cd eCIS-screen
	chmod 755 filter_hmmtab.pl
	chmod 755 gbk2IDs.pl
	chmod 755 gbk2seq.pl
	chmod 755 hmmsearch2tab.pl
	chmod 755 HMMsearch_genomesII.pl
	chmod 755 parse_hmmtab4eCIS.pl
	conda deactivate
	cd "$CGCMSDIR/external_tools"
fi
# Install Fama
if conda env list | grep 'cgcms-fama' >/dev/null 2>&1; then
    echo "Found cgcms-fama environment"
    conda activate cgcms-fama
    if python $CGCMSDIR/external_tools/fama/py/fama.py | grep "usage: fama.py" >/dev/null 2>&1
    then
	echo "Fama found"
    else
	echo 'Conda environment cgcms-fama exists but Fama was not properly installed. Remove the environment and restart CGCMS installation script.'
		echo 'To remove the environment, run:'
		echo '   conda remove -n cgcms-fama --all'
	conda deactivate
	exit 1
    fi
    conda deactivate
else
	echo "Installing Fama"
	conda create -y -n cgcms-fama python=3.8
	conda activate cgcms-fama
	conda install -y diamond krona
	git clone https://github.com/aekazakov/fama.git
	cd fama
	pip install -r requirements.txt
	/bin/sh install_reference_data_cgcms.sh "$CGCMSDIR/external_refdata/"
	conda deactivate
	cd "$CGCMSDIR/external_tools"
fi
# Install phispy
if conda env list | grep 'cgcms-phispy' >/dev/null 2>&1; then
    echo "Found cgcms-phispy environment"
    conda activate cgcms-phispy
    if phispy -v | grep "4.2.21" >/dev/null 2>&1; then
		echo "PhiSpy 4.2.21 found"
    else
		echo 'Conda environment cgcms-phispy exists but PhiSpy 4.2.21 was not properly installed. Remove the environment and restart CGCMS installation script.'
		echo 'To remove the environment, run:'
		echo '   conda remove -n cgcms-phispy --all'
		conda deactivate
		exit 1
    fi
    conda deactivate
else
	echo "Installing Phispy"
	conda create -y -n cgcms-phispy python=3.8
	conda activate cgcms-phispy
	conda install -y -c bioconda phispy 
	conda deactivate
	cd "$CGCMSDIR/external_tools"
fi
#Install GapMind
if conda env list | grep 'cgcms-gapmind' >/dev/null 2>&1; then
    echo "Found cgcms-gapmind environment"
    conda activate cgcms-gapmind
    if perl "$CGCMSDIR/external_tools/PaperBLAST/bin/gapsearch.pl" 2>&1|grep "Usage: gapsearch.pl" >/dev/null ; then
		echo "GapMind found"
    else
		echo 'Conda environment cgcms-gapmind exists but GapMind was not properly installed. Remove the environment and restart CGCMS installation script.'
		echo 'To remove the environment, run:'
		echo '   conda remove -n cgcms-gapmind --all'
		conda deactivate
		exit 1
    fi
    conda deactivate
else
	echo "Installing GapMind"
	git clone https://github.com/aekazakov/PaperBLAST.git
	cd PaperBLAST
	bash setup.sh
	cd "$CGCMSDIR/external_tools"
fi
#Install DefenseFinder
if conda env list | grep 'cgcms-defensefinder' >/dev/null 2>&1; then
    echo "Found cgcms-defensefinder environment"
    conda activate cgcms-defensefinder
    if defense-finder run --help|grep "Usage: defense-finder" >/dev/null 2>&1; then
		echo "DefenseFinder found"
    else
		echo 'Conda environment cgcms-defensefinder exists but DefenseFinder was not properly installed. Remove the environment and restart CGCMS installation script.'
		echo 'To remove the environment, run:'
		echo '   conda remove -n cgcms-defensefinder --all'
		conda deactivate
		exit 1
    fi
    conda deactivate
else
	echo "Installing DefenseFinder"
	conda create -y -n cgcms-defensefinder
	conda activate cgcms-defensefinder
	mkdir "$CGCMSDIR/external_refdata/defensefinder"
    mkdir "$CGCMSDIR/external_refdata/defensefinder/data"
    conda install -y -c bioconda hmmer
    conda install -y pip
    pip install mdmparis-defense-finder
    defense-finder update --models-dir "$CGCMSDIR/external_refdata/defensefinder/data"
	conda deactivate
	cd "$CGCMSDIR/external_tools"
fi
# Download HMM libraries
if ! [ -f "$CGCMSDIR/external_refdata/phispy/pvogs.hmm" ]; then
	echo "Downloading pVOG HMMs"
	[ -d "$CGCMSDIR/external_refdata/phispy" ] || mkdir "$CGCMSDIR/external_refdata/phispy"
	cd "$CGCMSDIR/external_refdata/phispy"
    curl -LJO -q http://iseq.lbl.gov/mydocs/cgcms_downloads/pvogs.hmm.gz
	gunzip pvogs.hmm.gz
	hmmpress pvogs.hmm
fi
if ! [ -f "$CGCMSDIR/external_refdata/pfam/Pfam-A.hmm" ]; then
	echo "Downloading PFAM HMMs"
	[ -d "$CGCMSDIR/external_refdata/pfam" ] || mkdir "$CGCMSDIR/external_refdata/pfam"
	cd "$CGCMSDIR/external_refdata/pfam"
    curl -LJO -q http://iseq.lbl.gov/mydocs/cgcms_downloads/pfam35.tar.gz
	tar xvf pfam35.tar.gz
	rm pfam35.tar.gz
	hmmpress Pfam-A.hmm
fi
if ! [ -f "$CGCMSDIR/external_refdata/tigrfam/TIGRFAM.HMM" ]; then
	echo "Downloading TIGRFAM HMMs"
	[ -d "$CGCMSDIR/external_refdata/tigrfam" ] || mkdir "$CGCMSDIR/external_refdata/tigrfam"
	cd "$CGCMSDIR/external_refdata/tigrfam"
    curl -LJO -q http://iseq.lbl.gov/mydocs/cgcms_downloads/tigrfam.tar.gz
	tar xvf tigrfam.tar.gz
	rm tigrfam.tar.gz
	hmmpress TIGRFAM.HMM
fi
echo "Activate virtual environment $CGCMSDIR/cgcms-venv and install Django before running CGCMS"
cd "$WORKDIR/genomebrowser"
if [ -f "configs.txt" ]; then
	cp configs.txt configs.txt~
	echo "Existing configs.txt copied to configs.txt~"
fi
cp configs.txt.template configs.txt
echo "cgcms.conda_path = $CONDA" >> configs.txt
echo "cgcms.temp_dir = $WORKDIR/temp" >> configs.txt
echo "cgcms.static_dir = $CGCMSDIR/static/$APPNAME/genomes" >> configs.txt
echo "cgcms.eggnog-mapper.data_dir = $CGCMSDIR/external_refdata/eggnog-mapper_v2.1.7" >> configs.txt
echo "cgcms.eggnog-mapper.dmnd_db = $CGCMSDIR/external_refdata/eggnog-mapper_v2.1.7/eggnog_proteins.dmnd" >> configs.txt
echo "cgcms.eggnog_outdir = $WORKDIR/temp/eggnog" >> configs.txt
echo "cgcms.eggnog_taxonomy = $WORKDIR/ref_data/eggnog_taxonomy_rules.txt" >> configs.txt
echo "cgcms.json_dir = $CGCMSDIR/static/$APPNAME/genomes/json" >> configs.txt
echo "cgcms.poem_command = $CGCMSDIR/external_tools/POEM_py3k/bin/run_poem_cgcms.sh" >> configs.txt
echo "cgcms.prepare_refseqs_command = $CGCMSDIR/external_tools/jbrowse/bin/prepare-refseqs.pl" >> configs.txt
echo "cgcms.flatfile_to_json_command = $CGCMSDIR/external_tools/jbrowse/bin/flatfile-to-json.pl" >> configs.txt
echo "cgcms.generate_names_command = $CGCMSDIR/external_tools/jbrowse/bin/generate-names.pl" >> configs.txt
echo "cgcms.search_db_dir = $WORKDIR/appdata" >> configs.txt
echo "cgcms.search_db_nucl = $WORKDIR/appdata/nucl.fna" >> configs.txt
echo "cgcms.search_db_prot = $WORKDIR/appdata/prot.faa" >> configs.txt
echo "plugins.antismash.antismash_ref = $WORKDIR/ref_data/ref_antismash.txt" >> configs.txt
echo "plugins.ecis_screen.ecis_hmm = $CGCMSDIR/external_tools/eCIS-screen/eCIS.hmm" >> configs.txt
echo "plugins.ecis_screen.ecis-screen_cmd = $CGCMSDIR/external_tools/eCIS-screen/HMMsearch_genomesII_fast.pl" >> configs.txt
echo "plugins.fama.fama_dir = $CGCMSDIR/external_tools/fama/py" >> configs.txt
echo "plugins.fama.fama_config = $CGCMSDIR/external_tools/fama/config.ini" >> configs.txt
echo "plugins.fama.fama_cazy_lib = $CGCMSDIR/external_refdata/fama/cazy2/cazy_v2_functions.txt" >> configs.txt
echo "plugins.fama.fama_nitrate_lib = $CGCMSDIR/external_refdata/fama/nitrogen11/fama_nitrogen-cycle_v.11.0_functions_thresholds.tsv" >> configs.txt
echo "plugins.fama.fama_universal_lib = $CGCMSDIR/external_refdata/fama/universal1.4/fama_function_thresholds_v.1.4.txt" >> configs.txt
echo "plugins.phispy.pvog_path = $CGCMSDIR/external_refdata/phispy/pvogs.hmm" >> configs.txt
echo "plugins.gapmind.gapmind_dir = $CGCMSDIR/external_tools/PaperBLAST" >> configs.txt
echo "plugins.defensefinder.defensefinder_models_dir = $CGCMSDIR/external_refdata/defensefinder/data" >> configs.txt
echo "ref.cazy_file = $WORKDIR/ref_data/ref_cazy.txt" >> configs.txt
echo "ref.cog_codes_file = $WORKDIR/ref_data/ref_cog_codes.txt" >> configs.txt
echo "ref.ec_file = $WORKDIR/ref_data/ref_ec.txt" >> configs.txt
echo "ref.go_file = $WORKDIR/ref_data/ref_go.txt" >> configs.txt
echo "ref.kegg_orthologs_file = $WORKDIR/ref_data/ref_kegg_ko.txt" >> configs.txt
echo "ref.kegg_pathways_file = $WORKDIR/ref_data/ref_kegg_pathways.txt" >> configs.txt
echo "ref.kegg_reactions_file = $WORKDIR/ref_data/ref_kegg_reactions.txt" >> configs.txt
echo "ref.tc_file = $WORKDIR/ref_data/ref_tc.txt" >> configs.txt
echo "ref.taxonomy = $WORKDIR/ref_data/ref_taxonomy.txt" >> configs.txt
echo "ref.pfam_hmm_lib = $CGCMSDIR/external_refdata/pfam/Pfam-A.hmm" >> configs.txt
echo "ref.pfam_hmm_list = $CGCMSDIR/external_refdata/pfam/ref_pfam.txt" >> configs.txt
echo "ref.tigrfam_hmm_lib = $CGCMSDIR/external_refdata/tigrfam/TIGRFAM.HMM" >> configs.txt
echo "ref.tigrfam_hmm_list = $CGCMSDIR/external_refdata/tigrfam/ref_tigrfam.txt" >> configs.txt
echo "configs.txt created."

if [ -f "secrets.json" ]; then
	cp secrets.json secrets.json~
	echo "Existing secrets.json copied to secrets.json~"
fi
cp secrets.json.template secrets.json
echo "    \"STATIC_ROOT\": \"$CGCMSDIR/static/$APPNAME\"," >> secrets.json
echo "    \"STATICFILES_DIR\": \"$WORKDIR/genomebrowser/static\"" >> secrets.json
echo "}" >> secrets.json
echo "secrets.json created."

echo "Edit secrets.json and genomebrowser/settings.py before running \"python manage.py configure_cgcsms -i configs.txt\""
