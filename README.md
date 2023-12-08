# Introduction

The comparative genomic content management system (CGCMS) is an open-source web-based platform for annotation, management and comparative analysis of microbial genomic sequences and associated data including ortholog families, protein domains, operons, regulatory interactions, metagenomic samples, strains taxonomy and metadata. 

CGCMS is a tool developed to create web portals for microbial genome collections each containing hundreds and thousands of genomes. The web portals are built on the Django framework and backed by a MySQL database that aggregates gene annotations generated by various bioinformatic tools. The genome annotation tools are installed in separate Conda environments and run by the CGCMS annotation pipeline. CGCMS employs Django Q, a multiprocessing task queue, for scheduling and executing pipeline jobs. In addition to the pipeline-generated data, administrators of a CGCMS-based portal can import gene annotations from text files or enter them manually in the site administration panel.

[Demo CGCMS-based genome collection portal](https://iseq.lbl.gov/genomes)


# Contents

[CGCMS installation and configuration](installation)

[User manual](user)

[Administrator manual](admin)

[For developer](developer)


# Image Credits

*Streptococcus pneumoniae* bacterial colonies that were grown on primary isolation medium. Photo by [Dr. Richard Facklam, USCDCP](https://pixnio.com/science/microscopy-images/streptococcus-pneumoniae/streptococcus-pneumoniae-bacterial-colonies-that-were-grown-on-primary-isolation-medium) on [Pixnio](https://pixnio.com/)

Earth spinning rotating animation. [Amirabbaszakavi](https://commons.wikimedia.org/wiki/File:Earth-spinning-rotating-animation-40.gif), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0), via Wikimedia Commons
