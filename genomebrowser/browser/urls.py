from django.urls import path

from . import views

urlpatterns = [
    path('', views.startpage, name='index'),
    path('textsearch/', views.textsearch, name='textsearch'),
    path('genome/<str:name>', views.genome_detail, name='genomedetails'),
    path('sample/<str:sample_id>', views.sample_detail, name='sampledetails'),
    path('strain/<str:strain_id>/', views.strain_detail, name='straindetails'),
    path('taxonomy/<str:taxonomy_id>', views.taxon_detail, name='taxondetails'),
    path('gene/<str:genome>/<str:locus_tag>/', views.gene_detail, name='genedetails'),
    path('operon/<str:genome>/<str:name>/', views.operon_detail, name='operondetails'),
    path('operon/<str:name>/', views.operon_detail, name='operondetails'),
    path('operons/<str:genome>/', views.OperonListView.as_view(), name='operonlist'),
    path('site/<str:genome>/<str:name>/', views.site_detail, name='sitedetails'),
    path('sites/<str:genome>/', views.SiteListView.as_view(), name='sitelist'),
    path('regulon/<str:genome>/<str:name>/', views.regulon_detail, name='regulondetails'),
    path('regulons/<str:genome>/', views.RegulonListView.as_view(), name='regulonlist'),
    path('getgene/', views.gene_byname, name='genebyname'),
    path('strains/', views.StrainListView.as_view(), name='strain_list'),
    path('samples/', views.SampleListView.as_view(), name='sample_list'),
    path('genomes/', views.GenomeListView.as_view(), name='genome_list'),
    path('genes/', views.GeneListView.as_view(), name='gene_list'),
    path('searchgene/', views.GeneSearchResultsAjaxView.as_view(), name='searchgene'),
    path('loadinggenesearch/', views.GeneSearchResultsAjaxView.ajax_view, name='loadinggenesearch'),
    path('searchannotation/', views.AnnotationSearchResultsAjaxView.as_view(), name='searchannotation'),
    path('loadingtextsearch/',views.AnnotationSearchResultsAjaxView.ajax_view,name="loadingtextsearch"),
    path('export/', views.export_csv, name='export'),
    path('exportfasta/', views.export_fasta, name='exportfasta'),
    path('searchgenome/', views.GenomeSearchResultsView.as_view(), name='searchgenome'),
    path('searchstrain/', views.StrainSearchResultsView.as_view(), name='searchstrain'),
    path('searchsample/', views.StrainSearchResultsView.as_view(), name='searchsample'),
    path('seqsearch', views.protein_search_external, name='proteinsearchexternal'),
    path('gos/', views.GoSearchResultsView.as_view(), name='gos'),
    path('kos/', views.KoSearchResultsView.as_view(), name='kos'),
    path('pathways/', views.KpSearchResultsView.as_view(), name='kps'),
    path('reactions/', views.KrSearchResultsView.as_view(), name='krs'),
    path('enzymes/', views.EcSearchResultsView.as_view(), name='enzymes'),
    path('transporters/', views.TcSearchResultsView.as_view(), name='transporters'),
    path('cazy/', views.CazySearchResultsView.as_view(), name='cazy'),
    path('cogs/', views.CogSearchResultsView.as_view(), name='cogs'),
    path('tag/<str:name>', views.tag_detail, name='tagdetails'),
    path('comparative/', views.ComparativeView.as_view(), name='comparative'),
    path('loadingscribl/',views.ComparativeView.ajax_view,name="loadingscribl"),
    path('cregulon/', views.cregulon_view, name='cregulon'),
    path('help/', views.show_help, name='help'),
    path('nuclsearchform/',views.nucleotidesearchform,name="nucleotidesearchform"),
    path('nuclsearch/',views.NsearchResultView.as_view(),name="nucleotidesearch"),
    path('protsearchform/',views.proteinsearchform,name="proteinsearchform"),
    path('protsearch/',views.PsearchResultView.as_view(),name="proteinsearch"),
    path('loadingnuclsearch/',views.NsearchResultView.ajax_view,name="loadingnuclsearch"),
    path('loadingprotsearch/',views.PsearchResultView.ajax_view,name="loadingprotsearch"),

]
handler404 = views.handler404
