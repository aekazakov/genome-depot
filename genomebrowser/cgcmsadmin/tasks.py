# Code implementing admin tasks
# Functions should be called from async_tasks.py
from browser.dataimport.importer import Importer

def test_task_impl(request, genome_names):
    print(request)
    print(genome_names)
    return 'Genomes:' + genome_names

def delete_genome(modeladmin, request, queryset):
    importer = Importer()
    genome_names = [genome.name for genome in queryset]
    queryset.delete()
    Protein.objects.filter(gene=None).delete()
    Strain.objects.filter(genome=None).delete()
    Sample.objects.filter(genome=None).delete()
    importer.export_proteins()
    importer.export_contigs()
    importer.delete_search_databases()
    shutil.copyfile(os.path.join(importer.config['cgcms.temp_dir'], os.path.basename(importer.config['cgcms.search_db_nucl'])), importer.config['cgcms.search_db_nucl'])
    shutil.copyfile(os.path.join(importer.config['cgcms.temp_dir'], os.path.basename(importer.config['cgcms.search_db_prot'])), importer.config['cgcms.search_db_prot'])
    importer.create_search_databases()
    for genome in genome_names:
        shutil.rmtree(os.path.join(importer.config['cgcms.json_dir'], genome))
    os.remove(importer.config['cgcms.search_db_nucl'])
    os.remove(importer.config['cgcms.search_db_prot'])
delete_genome.short_description = 'Delete genome and static content'

def import_genomes_impl(lines):
    print ('Asynchronous task received. Starting import.')
    print(lines)
    importer = Importer()
    result = importer.import_genomes(lines)
    return result
