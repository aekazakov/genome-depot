from django.db import models

# Create your models here.
class ChangeLog(models.Model):
    data = models.TextField(null=False, blank=True)
    action = models.CharField(max_length=16, null=False, blank=True)  # saved or deleted
    timestamp = models.DateTimeField(null=False, blank=True)

    class Meta:
        app_label = "browser"
        db_table = "model_change_logs"


class Config(models.Model):
    param = models.CharField(max_length=255, unique=True)
    value = models.CharField(max_length=255)

    def __str__(self):
        return self.param + ':' + self.value


class Taxon(models.Model):
    taxonomy_id = models.CharField(max_length=10, unique=True, db_index=True)
    eggnog_taxid = models.CharField(max_length=10, blank=True, null=True)
    rank = models.CharField(max_length=20)
    parent_id = models.CharField(max_length=10)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name + ' (' + self.taxonomy_id + ')'
        
    @property
    def admin_name(self):
        return self.name + ' (' + self.taxonomy_id + '): ' + self.rank


class Strain(models.Model):
    strain_id = models.CharField(max_length=100, unique=True, db_index=True)
    full_name = models.CharField(max_length=200)
    order = models.CharField(max_length=100)
    taxon = models.ForeignKey(Taxon, on_delete=models.CASCADE, blank = False, null = False)

    def __str__(self):
        return self.strain_id + ' (' + self.order + ')'


class Sample(models.Model):
    sample_id = models.CharField(max_length=100, unique=True)
    full_name = models.CharField(max_length=200)
    description = models.CharField(max_length=250)

    def __str__(self):
        return self.sample_id + ' (' + self.full_name + ')'


class Strain_metadata(models.Model):
    strain = models.ForeignKey(Strain, on_delete=models.CASCADE)
    source = models.CharField(max_length=30)
    url = models.CharField(max_length=250)
    key = models.CharField(max_length=250)
    value = models.TextField()

    def __str__(self):
        return self.source + ': ' + self.value


class Sample_metadata(models.Model):
    sample = models.ForeignKey(Sample, on_delete=models.CASCADE)
    source = models.CharField(max_length=30)
    url = models.CharField(max_length=250)
    key = models.CharField(max_length=250)
    value = models.TextField()

    def __str__(self):
        return self.source + ': ' + self.value


class Genome(models.Model):
    name = models.CharField(max_length=200, unique=True, db_index=True)
    description = models.TextField()
    strain = models.ForeignKey(Strain, on_delete=models.CASCADE, blank = True, null = True)
    sample = models.ForeignKey(Sample, on_delete=models.CASCADE, blank = True, null = True)
    contigs = models.PositiveIntegerField()
    size = models.PositiveIntegerField()
    genes = models.PositiveIntegerField()
    json_url = models.CharField(max_length=200)
    pub_date = models.DateTimeField('date published')
    external_url = models.CharField(max_length=200)
    external_id = models.CharField(max_length=40)
    gbk_filepath = models.CharField(max_length=200)
    taxon = models.ForeignKey(Taxon, on_delete=models.CASCADE, blank = False, null = False)

    def __str__(self):
        return self.name


class Contig(models.Model):
    contig_id = models.CharField(max_length=100)
    name = models.CharField(max_length=250)
    size = models.IntegerField()
    genome = models.ForeignKey(Genome, on_delete=models.CASCADE, blank = True, null = True)
    
    def __str__(self):
        return self.contig_id


class Operon(models.Model):
    name = models.CharField(max_length=250)
    start = models.IntegerField()
    end = models.IntegerField()
    strand = models.IntegerField()
    genome = models.ForeignKey(Genome, on_delete=models.CASCADE, blank = True, null = True)
    contig = models.ForeignKey(Contig, on_delete=models.SET_NULL, blank = True, null = True)
    
    def __str__(self):
        return self.name

    class Meta:
        unique_together = ('genome', 'name',)


class Cog_class(models.Model):
    cog_id = models.CharField(max_length=1, unique=True)
    description = models.CharField(max_length=80)

    def __str__(self):
        return self.cog_id


class Kegg_reaction(models.Model):
    kegg_id = models.CharField(max_length=10, unique=True)
    description = models.TextField()

    def __str__(self):
        return self.kegg_id


class Kegg_pathway(models.Model):
    kegg_id = models.CharField(max_length=10, unique=True)
    description = models.CharField(max_length=200)

    def __str__(self):
        return self.kegg_id


class Kegg_ortholog(models.Model):
    kegg_id = models.CharField(max_length=10, unique=True)
    description = models.TextField()

    def __str__(self):
        return self.kegg_id


class Go_term(models.Model):
    go_id = models.CharField(max_length=12, unique=True)
    go_namespace = models.CharField(max_length=50)
    description = models.TextField()

    def __str__(self):
        return self.go_id


class Cazy_family(models.Model):
    cazy_id = models.CharField(max_length=12, unique=True)
    description = models.TextField()

    def __str__(self):
        return self.cazy_id


class Ec_number(models.Model):
    ec_number = models.CharField(max_length=12, unique=True)
    description = models.TextField()

    def __str__(self):
        return self.ec_number


class Tc_family(models.Model):
    tc_id = models.CharField(max_length=15, unique=True)
    description = models.TextField()

    def __str__(self):
        return self.tc_id


class Ortholog_group(models.Model):
    eggnog_id = models.CharField(max_length=15, db_index=True)
    taxon = models.ForeignKey(Taxon, on_delete=models.SET_NULL, blank = True, null = True)

    def __str__(self):
        return self.eggnog_id

    class Meta:
        unique_together = ('eggnog_id', 'taxon',)


class Eggnog_description(models.Model):
    fingerprint = models.CharField(max_length=32, unique=True)
    description = models.TextField()

    def __str__(self):
        return self.description


class Protein(models.Model):
    name = models.CharField(max_length=100)
    length = models.IntegerField()
    protein_hash = models.CharField(max_length=32, unique=True, db_index=True)
    sequence = models.TextField()
    cog_classes = models.ManyToManyField(Cog_class)
    kegg_reactions = models.ManyToManyField(Kegg_reaction)
    kegg_pathways = models.ManyToManyField(Kegg_pathway)
    kegg_orthologs = models.ManyToManyField(Kegg_ortholog)
    go_terms = models.ManyToManyField(Go_term)
    cazy_families = models.ManyToManyField(Cazy_family)
    ec_numbers = models.ManyToManyField(Ec_number)
    tc_families = models.ManyToManyField(Tc_family)
    ortholog_groups = models.ManyToManyField(Ortholog_group)
    taxonomy_id = models.ForeignKey(Taxon, on_delete=models.SET_NULL, blank = True, null = True)
    eggnog_description = models.ForeignKey(Eggnog_description, on_delete=models.SET_NULL, blank = True, null = True)

    def __str__(self):
        return self.protein_hash


class Gene(models.Model):
    name = models.CharField(max_length=50, db_index=True)
    locus_tag = models.CharField(max_length=50, db_index=True)
    contig = models.ForeignKey(Contig, on_delete=models.SET_NULL, blank = True, null = True)
    type = models.CharField(max_length=20)
    start = models.IntegerField()
    end = models.IntegerField()
    strand = models.IntegerField()
    genome = models.ForeignKey(Genome, on_delete=models.CASCADE)
    function = models.CharField(max_length=250, db_index=True) # prokka function
    protein = models.ForeignKey(Protein, on_delete=models.SET_NULL, blank = True, null = True)
    operon = models.ForeignKey(Operon, on_delete=models.SET_NULL, blank = True, null = True, related_name='genes')
    
    def __str__(self):
        return self.locus_tag


class Regulon(models.Model):
    name = models.CharField(max_length=50, db_index=True)
    regulators = models.ManyToManyField(Gene)
    genome = models.ForeignKey(Genome, on_delete=models.CASCADE)
    description = models.TextField()
    
    def __str__(self):
        return self.name + '(' + self.genome.name + ')'


class Site(models.Model):
    name = models.CharField(max_length=50)
    type = models.CharField(max_length=20)
    start = models.IntegerField()
    end = models.IntegerField()
    strand = models.IntegerField()
    sequence = models.TextField()
    regulon = models.ForeignKey(Regulon, on_delete=models.CASCADE)
    genome = models.ForeignKey(Genome, on_delete=models.CASCADE)
    contig = models.ForeignKey(Contig, on_delete=models.SET_NULL, blank = True, null = True)
    operons = models.ManyToManyField(Operon, blank = True)
    genes = models.ManyToManyField(Gene, blank = True)
    
    def __str__(self):
        return self.name
    
    class Meta:
        unique_together = ('regulon', 'name',)


class Annotation(models.Model):
    gene_id = models.ForeignKey(Gene, on_delete=models.CASCADE)
    source = models.CharField(max_length=30, db_index=True)
    url = models.CharField(max_length=300)
    key = models.CharField(max_length=30, db_index=True)
    value = models.CharField(max_length=50, db_index=True)
    note = models.TextField()

    def __str__(self):
        return self.source + ': ' + self.value


