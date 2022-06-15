#!/bin/bash
WORKDIR=$(pwd)
CGCMSDIR=$(dirname "$WORKDIR")
APPNAME=$(basename "$CGCMSDIR")
CGCMSDIR=$(dirname "$CGCMSDIR")
CGCMSDIR=$(dirname "$CGCMSDIR")
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
source $CONDA
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

cd "$CGCMSDIR/external_tools"
# Install eggnog-mapper
if ! { conda env list | grep 'cgcms-emapper'; } >/dev/null 2>&1; then
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
if ! { conda env list | grep 'cgcms-poem'; } >/dev/null 2>&1; then
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
if ! { conda env list | grep 'cgcms-amrfinder'; } >/dev/null 2>&1; then
	echo "Installing AMRFinder"
	conda create -y -n cgcms-amrfinder python=3.8
	conda activate cgcms-amrfinder
	conda install -y -c bioconda -c conda-forge ncbi-amrfinderplus
	conda deactivate
	cd "$CGCMSDIR/external_tools"
fi
# Install antismash
if ! { conda env list | grep 'cgcms-antismash'; } >/dev/null 2>&1; then
	echo "Installing AntiSMASH"
	conda create -y -n cgcms-antismash python=3.8
	conda activate cgcms-antismash
	conda install -y antismash
	conda deactivate
	cd "$CGCMSDIR/external_tools"
fi
# Install ecis-screen
if ! { conda env list | grep 'cgcms-ecis-screen'; } >/dev/null 2>&1; then
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
if ! { conda env list | grep 'cgcms-fama'; } >/dev/null 2>&1; then
	echo "Installing Fama"
	conda create -y -n cgcms-fama python=3.8
	conda activate cgcms-fama
	conda install -y diamond 
	git clone https://github.com/novichkov-lab/fama.git
	cd fama
	pip install -r requirements.txt
	/bin/sh install_reference_data_cgcms.sh
	conda deactivate
	cd "$CGCMSDIR/external_tools"
fi
# Install phispy
if ! { conda env list | grep 'cgcms-phispy'; } >/dev/null 2>&1; then
	echo "Installing Phispy"
	conda create -y -n cgcms-phispy python=3.8
	conda activate cgcms-phispy
	conda install -y -c bioconda phispy 
	conda deactivate
	cd "$CGCMSDIR/external_tools"
fi
# Download HMM libraries
if ! [ -f "$CGCMSDIR/external_refdata/phispy/pVOGs.hmm" ]; then
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
    curl -LJO -q http://iseq.lbl.gov/mydocs/cgcms_downloads/pfam.tar.gz
	tar xvf pfam.tar.gz
	rm pfam.tar.gz
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

echo "Activate virtual environment $CGCMSDIR/cgcms-venv  and install Django before running CGCMS"
cd "$WORKDIR/genomebrowser"
if [ -f "configs.txt" ]; then
	cp configs.txt configs.txt~
	echo "Existing configs.txt copied to configs.txt~"
fi
cp configs.txt.template configs.txt
echo "cgcms.conda_path = $CONDA" >> configs.txt
echo "cgcms.eggnog-mapper.data_dir = $CGCMSDIR/external_refdata/eggnog-mapper_v2.1.7" >> configs.txt
echo "cgcms.eggnog-mapper.dmnd_db = $CGCMSDIR/external_refdata/eggnog-mapper_v2.1.7/eggnog_proteins.dmnd" >> configs.txt
echo "cgcms.eggnog_outdir = $WORKDIR/temp/eggnog" >> configs.txt
echo "cgcms.eggnog_taxonomy = $WORKDIR/ref_data/eggnog_taxonomy_rules.txt" >> configs.txt
echo "cgcms.json_dir = $WORKDIR/ref_data/eggnog_taxonomy_rules.txt" >> configs.txt
echo "cgcms.poem_command = $CGCMSDIR/external_tools/POEM_py3k/bin/run_poem_cgcms.sh" >> configs.txt
echo "cgcms.search_db_dir = $WORKDIR/appdata" >> configs.txt
echo "cgcms.search_db_nucl = $WORKDIR/appdata/nucl.fna" >> configs.txt
echo "cgcms.search_db_prot = $WORKDIR/appdata/prot.faa" >> configs.txt
echo "cgcms.static_dir = $CGCMSDIR/static/$APPNAME/genomes" >> configs.txt
echo "cgcms.temp_dir = $WORKDIR/temp" >> configs.txt
echo "plugins.antismash.antismash_ref = $WORKDIR/ref_data/ref_antismash.txt" >> configs.txt
echo "plugins.ecis_screen.ecis_hmm = $CGCMSDIR/external_tools/eCIS-screen/eCIS.hmm" >> configs.txt
echo "plugins.fama.fama_dir	$CGCMSDIR/external_tools/fama/py" >> configs.txt
echo "plugins.fama.fama_config = $CGCMSDIR/external_tools/fama/config.ini" >> configs.txt
echo "plugins.fama.fama_cazy_lib = $CGCMSDIR/external_refdata/fama/cazy2/collection_functions.txt" >> configs.txt
echo "plugins.fama.fama_nitrate_lib = $CGCMSDIR/external_refdata/fama/nitrogen11/fama_nitrogen-cycle_v.11.0_functions_thresholds.tsv" >> configs.txt
echo "plugins.fama.fama_universal_lib = $CGCMSDIR/external_refdata/fama/universal1.4/fama_function_thresholds_v.1.4.txt" >> configs.txt
echo "plugins.phispy.pvog_path = $CGCMSDIR/external_refdata/phispy/pvogs.hmm" >> configs.txt
echo "ref.cazy_file = $WORKDIR/ref_data/ref_cazy.txt" >> configs.txt
echo "ref.cog_codes_file = $WORKDIR/ref_data/ref_cog_codes.txt" >> configs.txt
echo "ref.ec_file = $WORKDIR/ref_data/ec_file" >> configs.txt
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

#Installing dependencies in the virtual environment
echo "Installing python dependencies in cgcms-venv virtual environment"
source "$CGCMSDIR/cgcms-venv/bin/activate"
pip install "django==3.2.6" django_admin_shortcuts django_cors_headers django_q django_debug_toolbar openpyxl parasail biopython toytree urllib3 mysqlclient --no-cache-dir
deactivate
echo "Edit secrets.json and genomebrowser/settings.py before running \"python manage.py configure_cgcsms -i configs.txt\""
