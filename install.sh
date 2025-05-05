#!/bin/bash
WORKDIR=$(pwd)
TOPDIR=$(dirname "$WORKDIR")
APPNAME=$(basename "$TOPDIR")
TOPDIR=$(dirname "$TOPDIR")
TOPDIR=$(dirname "$TOPDIR")

CONDADIR=$(dirname "$(which conda)")
CONDADIR=$(dirname "$CONDADIR")
CONDA="$CONDADIR/etc/profile.d/conda.sh"
if ! [ -f "$CONDA" ]; then
        echo "Conda is not installed. Usage: setup.sh <conda dir>"
        exit
fi
echo "Found Conda installation: $CONDADIR"

# Deactivate all conda environments (including base) because python virtual environment must use system python, not the python installed in conda
conda deactivate
conda deactivate
if ! [ -d "$TOPDIR/genomedepot-venv" ]; then
    python3 -m venv "$TOPDIR/genomedepot-venv"
fi

#Installing dependencies in the virtual environment
echo "Installing python dependencies in genomedepot-venv virtual environment"
source "$TOPDIR/genomedepot-venv/bin/activate"
pip install "django==3.2.6" django_admin_shortcuts django_cors_headers django_q django_debug_toolbar openpyxl "parasail==1.2.4" biopython "toytree==2.0.1" urllib3 PyMySQL plotly pandas python-decouple django-admin-logviewer django-cookiebanner --no-cache-dir
deactivate

source $CONDA
conda config --add channels bioconda
conda config --add channels conda-forge
# Create directory structure
[ -d "$TOPDIR/external_tools" ] || mkdir "$TOPDIR/external_tools"
[ -d "$TOPDIR/external_refdata" ] || mkdir "$TOPDIR/external_refdata"
[ -d "$TOPDIR/static" ] || mkdir "$TOPDIR/static"
[ -d "$TOPDIR/static/$APPNAME" ] || mkdir "$TOPDIR/static/$APPNAME"
[ -d "$TOPDIR/static/$APPNAME/genomes" ] || mkdir "$TOPDIR/static/$APPNAME/genomes"
[ -d "$TOPDIR/static/$APPNAME/genomes/json" ] || mkdir "$TOPDIR/static/$APPNAME/genomes/json"
[ -d "$TOPDIR/static/$APPNAME/genomes/gbff" ] || mkdir "$TOPDIR/static/$APPNAME/genomes/gbff"
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
cd "$TOPDIR/external_tools"
if ! [ -d "$TOPDIR/external_tools/jbrowse" ]; then
    # Install Jbrowse v.1.16.11
    curl -L -O https://github.com/GMOD/jbrowse/releases/download/1.16.11-release/JBrowse-1.16.11.zip
    unzip JBrowse-1.16.11.zip
    mv JBrowse-1.16.11 jbrowse
    rm JBrowse-1.16.11.zip
    cd jbrowse
    ./setup.sh
    cd "$TOPDIR/external_tools"
fi

# Install samtools and tabix
cd "$TOPDIR/external_tools"
if conda env list | grep "genomedepot-jbrowse" >/dev/null 2>&1; then
    conda activate genomedepot-jbrowse
    if ! { samtools version | grep samtools; }>/dev/null 2>&1; then
        echo "Conda environment genomedepot-jbrowse exists but samtools was not properly installed. Remove the environment and restart the GenomeDepot installation script."
        echo "To remove the environment, run:"
        echo "   conda remove -n genomedepot-jbrowse --all"
        conda deactivate
        exit 1
    else
        echo "samtools found"
    fi
    conda deactivate
else
    echo "Installing samtools and tabix"
    conda create -y -n genomedepot-jbrowse
    conda activate genomedepot-jbrowse
    conda install -y -c bioconda samtools
    conda install -y -c bioconda tabix
    conda deactivate
fi

# Install eggnog-mapper
if conda env list | grep "genomedepot-emapper" >/dev/null 2>&1; then
    conda activate genomedepot-emapper
    if ! { emapper.py -v --data_dir "$TOPDIR/external_refdata/eggnog-mapper_v2.1.7" | grep "emapper-2.1.7"; }>/dev/null 2>&1; then
        echo "Conda environment genomedepot-emapper exists but eggnog-mapper v2.1.7 was not properly installed. Remove the environment and restart the GenomeDepot installation script."
        echo "To remove the environment, run:"
        echo "   conda remove -n genomedepot-emapper --all"
        conda deactivate
        exit 1
    else
        echo "EggNOG-mapper v.2.1.7 found"
    fi
    conda deactivate
else
    echo "Installing EggNOG mapper"
    conda create -y -n genomedepot-emapper python=3.8
    conda activate genomedepot-emapper
    conda install -c bioconda -y eggnog-mapper=2.1.7
    conda deactivate
    cd "$TOPDIR/external_tools"
fi
# Install eggnog-mapper reference data
if ! [ -d "$TOPDIR/external_refdata/eggnog-mapper_v2.1.7" ]; then
    conda activate genomedepot-emapper
    mkdir "$TOPDIR/external_refdata/eggnog-mapper_v2.1.7"
    echo "Downloading reference databases for eggnog-mapper_v2.1.7. Answer \"y\" to all questions."
    download_eggnog_data.py -y --data_dir "$TOPDIR/external_refdata/eggnog-mapper_v2.1.7"
    conda deactivate
fi
# Install poem_py3
if conda env list | grep "genomedepot-poem" >/dev/null 2>&1; then
    echo "Found genomedepot-poem environment"
    conda activate genomedepot-poem
    if [ -f "$TOPDIR/external_tools/POEM_py3k/bin/run_poem_cgcms.sh" ]
    then
        echo "POEM_py3k found"
    else
        echo "Conda environment genomedepot-poem exists but POEM_py3k was not properly installed. Remove the environment and restart the GenomeDepot installation script."
        echo "To remove the environment, run:"
        echo "   conda remove -n genomedepot-poem --all"
        conda deactivate
        exit 1
    fi
    conda deactivate
else
    echo "Installing POEM_py3"
    conda create -y -n genomedepot-poem python=3.7
    conda activate genomedepot-poem
    git clone https://github.com/aekazakov/POEM_py3k.git
    cd ./POEM_py3k
    bash ./install.sh genomedepot-poem
    conda deactivate
    cd "$TOPDIR/external_tools"
fi
#Install amrfinder
if conda env list | grep "genomedepot-amrfinder" >/dev/null 2>&1; then
        echo "Found genomedepot-amrfinder environment"
        conda activate genomedepot-amrfinder
        if amrfinder --version|grep "3.11.11" >/dev/null 2>&1; then
                echo "AMRFinder found"
        else
                echo "Conda environment genomedepot-amrfinder exists but AMRFinder was not properly installed. Remove the environment and restart the GenomeDepot installation script."
                echo "To remove the environment, run:"
                echo "   conda remove -n genomedepot-amrfinder --all"
                conda deactivate
                exit 1
        fi
        conda deactivate
else
        echo "Installing AMRFinder"
        conda create -y -n genomedepot-amrfinder python=3.8
        conda activate genomedepot-amrfinder
        conda install -y -c bioconda -c conda-forge ncbi-amrfinderplus=3.11.11
        amrfinder -u
        conda deactivate
        cd "$TOPDIR/external_tools"
fi
# Install antismash
if conda env list | grep "genomedepot-antismash" >/dev/null 2>&1; then
    echo "Found genomedepot-antismash environment"
    conda activate genomedepot-antismash
    if antismash -h|grep "antiSMASH 6.1.1" >/dev/null 2>&1;     then
                echo "antiSMASH 6.1.1 found"
    else
                echo "Conda environment genomedepot-antismash exists but antiSMASH was not properly installed. Remove the environment and restart the GenomeDepot installation script."
                echo "To remove the environment, run:"
                echo "   conda remove -n genomedepot-antismash --all"
                conda deactivate
                exit 1
        fi
        conda deactivate
else
        echo "Installing AntiSMASH"
        conda create -y -n genomedepot-antismash python=3.8
        conda activate genomedepot-antismash
        conda install -y antismash
        conda deactivate
        cd "$TOPDIR/external_tools"
fi
# Install ecis-screen
if conda env list | grep "genomedepot-ecis-screen" >/dev/null 2>&1; then
        echo "Found genomedepot-ecis-screen environment"
        conda activate genomedepot-ecis-screen
    if $TOPDIR/external_tools/eCIS-screen/HMMsearch_genomesII_fast.pl 2>&1 | grep "Pipeline for screening" >/dev/null 2>&1
        then
        echo "eCIS-screen found"
        else
        echo "Conda environment genomedepot-ecis-screen exists but eCIS-screen was not properly installed. Remove the environment and restart GenomeDepot installation script."
                echo "To remove the environment, run:"
                echo "   conda remove -n genomedepot-ecis-screen --all"
        conda deactivate
        exit 1
        fi
        conda deactivate
else
        echo "Installing eCIS-screen"
        conda create -y -n genomedepot-ecis-screen
        conda activate genomedepot-ecis-screen
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
        cd "$TOPDIR/external_tools"
fi
# Install Fama
if conda env list | grep "genomedepot-fama" >/dev/null 2>&1; then
    echo "Found genomedepot-fama environment"
    conda activate genomedepot-fama
    if python $TOPDIR/external_tools/fama/py/fama.py | grep "usage: fama.py" >/dev/null 2>&1
    then
        echo "Fama found"
    else
        echo "Conda environment genomedepot-fama exists but Fama was not properly installed. Remove the environment and restart GenomeDepot installation script."
                echo "To remove the environment, run:"
                echo "   conda remove -n genomedepot-fama --all"
        conda deactivate
        exit 1
    fi
    conda deactivate
else
        echo "Installing Fama"
        conda create -y -n genomedepot-fama python=3.8
        conda activate genomedepot-fama
        conda install -y diamond krona
        git clone https://github.com/aekazakov/fama.git
        cd fama
        pip install -r requirements.txt
        /bin/sh install_reference_data_cgcms.sh "$TOPDIR/external_refdata/"
        conda deactivate
        cd "$TOPDIR/external_tools"
fi
# Install phispy
if conda env list | grep "genomedepot-phispy" >/dev/null 2>&1; then
    echo "Found genomedepot-phispy environment"
    conda activate genomedepot-phispy
    if phispy -v | grep "4.2.21" >/dev/null 2>&1; then
                echo "PhiSpy 4.2.21 found"
    else
                echo "Conda environment genomedepot-phispy exists but PhiSpy 4.2.21 was not properly installed. Remove the environment and restart GenomeDepot installation script."
                echo "To remove the environment, run:"
                echo "   conda remove -n genomedepot-phispy --all"
                conda deactivate
                exit 1
    fi
    conda deactivate
else
    echo "Installing Phispy"
    conda create -y -n genomedepot-phispy python=3.8
    conda activate genomedepot-phispy
    conda install -y -c bioconda biopython==1.80 
    conda install -y -c bioconda phispy 
    conda deactivate
    cd "$TOPDIR/external_tools"
fi
#Install GapMind
if conda env list | grep "genomedepot-gapmind" >/dev/null 2>&1; then
    echo "Found genomedepot-gapmind environment"
    conda activate genomedepot-gapmind
    if perl "$TOPDIR/external_tools/PaperBLAST/bin/gapsearch.pl" 2>&1|grep "Usage: gapsearch.pl" >/dev/null ; then
        echo "GapMind found"
    else
        echo "Conda environment genomedepot-gapmind exists but GapMind was not properly installed. Remove the environment and restart GenomeDepot installation script."
        echo "To remove the environment, run:"
        echo "   conda remove -n genomedepot-gapmind --all"
        conda deactivate
        exit 1
    fi
    conda deactivate
else
    echo "Installing GapMind"
    git clone https://github.com/aekazakov/PaperBLAST.git
    cd PaperBLAST
    bash setup.sh genomedepot-gapmind
    cd "$TOPDIR/external_tools"
fi
#Install DefenseFinder
if conda env list | grep "genomedepot-defensefinder" >/dev/null 2>&1; then
    echo "Found genomedepot-defensefinder environment"
    conda activate genomedepot-defensefinder
    if defense-finder run --help|grep "Usage: defense-finder" >/dev/null 2>&1; then
        echo "DefenseFinder found"
    else
        echo "Conda environment genomedepot-defensefinder exists but DefenseFinder was not properly installed. Remove the environment and restart GenomeDepot installation script."
        echo "To remove the environment, run:"
        echo "   conda remove -n genomedepot-defensefinder --all"
        conda deactivate
        exit 1
    fi
    conda deactivate
else
    echo "Installing DefenseFinder"
    conda create -y -n genomedepot-defensefinder
    conda activate genomedepot-defensefinder
    mkdir "$TOPDIR/external_refdata/defensefinder"
    mkdir "$TOPDIR/external_refdata/defensefinder/data"
    conda install -y -c bioconda hmmer
    conda install -y pip
    pip install mdmparis-defense-finder
    defense-finder update --models-dir "$TOPDIR/external_refdata/defensefinder/data"
    conda deactivate
    cd "$TOPDIR/external_tools"
fi
#Install MacsyFinder
if conda env list | grep "genomedepot-macsyfinder" >/dev/null 2>&1; then
    echo "Found genomedepot-macsyfinder environment"
    conda activate genomedepot-macsyfinder
    if macsyfinder --version|grep "Macsyfinder" >/dev/null 2>&1; then
        echo "MacsyFinder found"
    else
        echo "Conda environment genomedepot-macsyfinder exists but MacsyFinder was not properly installed. Remove the environment and restart GenomeDepot installation script."
        echo "To remove the environment, run:"
        echo "   conda remove -n genomedepot-macsyfinder --all"
        conda deactivate
        exit 1
    fi
    conda deactivate
else
    echo "Installing MacsyFinder"
    conda create -y -n genomedepot-macsyfinder
    conda activate genomedepot-macsyfinder
    conda install -y -c bioconda macsyfinder
    mkdir "$TOPDIR/external_refdata/macsyfinder"
    mkdir "$TOPDIR/external_refdata/macsyfinder/data"
    macsydata install --target "$TOPDIR/external_refdata/macsyfinder/data" TXSScan
    conda deactivate
    cd "$TOPDIR/external_tools"
fi
#Install geNomad
if conda env list | grep "genomedepot-genomad" >/dev/null 2>&1; then
    echo "Found genomedepot-genomad environment"
    conda activate genomedepot-genomad
    if genomad | grep geNomad >/dev/null 2>&1; then
                echo "geNomad found"
    else
                echo "Conda environment genomedepot-genomad exists but geNomad was not properly installed. Remove the environment and restart GenomeDepot installation script."
                echo "To remove the environment, run:"
                echo "   conda remove -n genomedepot-genomad --all"
                conda deactivate
                exit 1
    fi
    conda deactivate
else
    echo "Installing geNomad"
    conda create -y -n genomedepot-genomad -c conda-forge -c bioconda genomad
    conda activate genomedepot-genomad
    mkdir "$TOPDIR/external_refdata/geNomad"
    cd "$TOPDIR/external_refdata/geNomad"
    genomad download-database .
    conda deactivate
    cd "$TOPDIR/external_tools"
fi

# Install HMMER for Pfam and TIGRFAM plugins
if conda env list | grep "genomedepot-hmmsearch" >/dev/null 2>&1; then
    echo "Found genomedepot-hmmsearch environment"
    conda activate genomedepot-hmmsearch
    if hmmsearch -h|grep "HMMER" >/dev/null 2>&1; then
        echo "hmmsearch found"
    else
        echo "Conda environment genomedepot-hmmsearch exists but HMMER was not properly installed. Remove the environment and restart GenomeDepot installation script."
        echo "To remove the environment, run:"
        echo "   conda remove -n genomedepot-hmmsearch --all"
        conda deactivate
        exit 1
    fi
    conda deactivate
else
    echo "Installing HMMER"
    conda create -y -n genomedepot-hmmsearch
    conda activate genomedepot-hmmsearch
    conda install -y -c bioconda hmmer
    conda deactivate
    cd "$TOPDIR/external_tools"
fi

# Download HMM libraries
if ! [ -f "$TOPDIR/external_refdata/phispy/pvogs.hmm" ]; then
    echo "Downloading pVOG HMMs"
    [ -d "$TOPDIR/external_refdata/phispy" ] || mkdir "$TOPDIR/external_refdata/phispy"
    cd "$TOPDIR/external_refdata/phispy"
    curl -LJO -q http://iseq.lbl.gov/mydocs/cgcms_downloads/pvogs.hmm.gz
    gunzip pvogs.hmm.gz
    conda activate genomedepot-phispy
    hmmpress pvogs.hmm
    conda deactivate
fi
if ! [ -f "$TOPDIR/external_refdata/pfam/Pfam-A.hmm" ]; then
    echo "Downloading PFAM HMMs"
    [ -d "$TOPDIR/external_refdata/pfam" ] || mkdir "$TOPDIR/external_refdata/pfam"
    cd "$TOPDIR/external_refdata/pfam"
    curl -LJO -q http://iseq.lbl.gov/mydocs/cgcms_downloads/pfam35.tar.gz
    tar xvf pfam35.tar.gz
    rm pfam35.tar.gz
    conda activate genomedepot-hmmsearch
    hmmpress Pfam-A.hmm
    conda deactivate
fi
if ! [ -f "$TOPDIR/external_refdata/tigrfam/TIGRFAM.HMM" ]; then
    echo "Downloading TIGRFAM HMMs"
    [ -d "$TOPDIR/external_refdata/tigrfam" ] || mkdir "$TOPDIR/external_refdata/tigrfam"
    cd "$TOPDIR/external_refdata/tigrfam"
    curl -LJO -q http://iseq.lbl.gov/mydocs/cgcms_downloads/tigrfam.tar.gz
    tar xvf tigrfam.tar.gz
    rm tigrfam.tar.gz
    conda activate genomedepot-hmmsearch
    hmmpress TIGRFAM.HMM
    conda deactivate
fi
echo "Activate virtual environment $TOPDIR/genomedepot-venv and install Django before running GenomeDepot"
cd "$WORKDIR/genomebrowser"
if [ -f "configs.txt" ]; then
    cp configs.txt configs.txt~
    echo "Existing configs.txt copied to configs.txt~"
fi
cp configs.txt.template configs.txt
echo "core.conda_path = $CONDA" >> configs.txt
echo "core.temp_dir = $WORKDIR/temp" >> configs.txt
echo "core.static_dir = $TOPDIR/static/$APPNAME/genomes" >> configs.txt
echo "core.eggnog-mapper.data_dir = $TOPDIR/external_refdata/eggnog-mapper_v2.1.7" >> configs.txt
echo "core.eggnog-mapper.dmnd_db = $TOPDIR/external_refdata/eggnog-mapper_v2.1.7/eggnog_proteins.dmnd" >> configs.txt
echo "core.eggnog_outdir = $WORKDIR/temp/eggnog" >> configs.txt
echo "core.eggnog_taxonomy = $WORKDIR/ref_data/eggnog_taxonomy_rules.txt" >> configs.txt
echo "core.json_dir = $TOPDIR/static/$APPNAME/genomes/json" >> configs.txt
echo "core.poem_dir = $TOPDIR/external_tools/POEM_py3k" >> configs.txt
echo "core.prepare_refseqs_command = $TOPDIR/external_tools/jbrowse/bin/prepare-refseqs.pl" >> configs.txt
echo "core.flatfile_to_json_command = $TOPDIR/external_tools/jbrowse/bin/flatfile-to-json.pl" >> configs.txt
echo "core.generate_names_command = $TOPDIR/external_tools/jbrowse/bin/generate-names.pl" >> configs.txt
echo "core.search_db_dir = $WORKDIR/appdata" >> configs.txt
echo "core.search_db_nucl = $WORKDIR/appdata/nucl.fna" >> configs.txt
echo "core.search_db_prot = $WORKDIR/appdata/prot.faa" >> configs.txt
echo "plugins.antismash.antismash_ref = $WORKDIR/ref_data/ref_antismash.txt" >> configs.txt
echo "plugins.ecis_screen.ecis_hmm = $TOPDIR/external_tools/eCIS-screen/eCIS.hmm" >> configs.txt
echo "plugins.ecis_screen.ecis-screen_cmd = $TOPDIR/external_tools/eCIS-screen/HMMsearch_genomesII_fast.pl" >> configs.txt
echo "plugins.fama.fama_dir = $TOPDIR/external_tools/fama/py" >> configs.txt
echo "plugins.fama.fama_config = $TOPDIR/external_tools/fama/config.ini" >> configs.txt
echo "plugins.fama.fama_cazy_lib = $TOPDIR/external_refdata/fama/cazy2/cazy_v2_functions.txt" >> configs.txt
echo "plugins.fama.fama_nitrate_lib = $TOPDIR/external_refdata/fama/nitrogen11/fama_nitrogen-cycle_v.11.0_functions_thresholds.tsv" >> configs.txt
echo "plugins.fama.fama_universal_lib = $TOPDIR/external_refdata/fama/universal1.4/fama_function_thresholds_v.1.4.txt" >> configs.txt
echo "plugins.phispy.pvog_path = $TOPDIR/external_refdata/phispy/pvogs.hmm" >> configs.txt
echo "plugins.gapmind.gapmind_dir = $TOPDIR/external_tools/PaperBLAST" >> configs.txt
echo "plugins.defensefinder.defensefinder_models_dir = $TOPDIR/external_refdata/defensefinder/data" >> configs.txt
echo "plugins.macsyfinder.models_dir = $TOPDIR/external_refdata/macsyfinder/data" >> configs.txt
echo "plugins.genomad.ref_db = $TOPDIR/external_refdata/geNomad/genomad_db" >> configs.txt
echo "ref.cazy_file = $WORKDIR/ref_data/ref_cazy.txt" >> configs.txt
echo "ref.cog_codes_file = $WORKDIR/ref_data/ref_cog_codes.txt" >> configs.txt
echo "ref.ec_file = $WORKDIR/ref_data/ref_ec.txt" >> configs.txt
echo "ref.go_file = $WORKDIR/ref_data/ref_go.txt" >> configs.txt
echo "ref.kegg_orthologs_file = $WORKDIR/ref_data/ref_kegg_ko.txt" >> configs.txt
echo "ref.kegg_pathways_file = $WORKDIR/ref_data/ref_kegg_pathways.txt" >> configs.txt
echo "ref.kegg_reactions_file = $WORKDIR/ref_data/ref_kegg_reactions.txt" >> configs.txt
echo "ref.tc_file = $WORKDIR/ref_data/ref_tc.txt" >> configs.txt
echo "ref.taxonomy = $WORKDIR/ref_data/ref_taxonomy.txt" >> configs.txt
echo "plugins.hmmsearch_pfam.hmm_lib = $TOPDIR/external_refdata/pfam/Pfam-A.hmm" >> configs.txt
echo "plugins.hmmsearch_pfam.ref_data = $TOPDIR/external_refdata/pfam/ref_pfam.txt" >> configs.txt
echo "plugins.hmmsearch_tigrfam.hmm_lib = $TOPDIR/external_refdata/tigrfam/TIGRFAM.HMM" >> configs.txt
echo "plugins.hmmsearch_tigrfam.ref_data = $TOPDIR/external_refdata/tigrfam/ref_tigrfam.txt" >> configs.txt
echo "configs.txt created."

if [ -f ".env" ]; then
        cp .env .env~
        echo "Existing .env file copied to .env~"
fi
cp .env.template .env
echo "STATIC_ROOT=$TOPDIR/static/$APPNAME" >> .env
echo "STATICFILES_DIR=$WORKDIR/genomebrowser/static" >> .env
echo "LOGVIEWER_LOGS=$WORKDIR/genomebrowser/django.log," >> .env

touch django.log
chmod 664 django.log

echo "Edit .env file before running \"python manage.py configure_cgcsms -i configs.txt\""
