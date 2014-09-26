#!perl
use strict;
use warnings FATAL => qw( all );

use SOAP::Lite;
use HTTP::Cookies;
use Data::Dumper;
use Getopt::Long;
use Carp;

my ($deseq_file, $direction, $foldchange, $min_padj);
GetOptions
(
	"deseq-file=s" => \$deseq_file,
	"direction=s" => \$direction,  # "up", "down", "up_or_down" 
	"foldchange=f" => \$foldchange, # e.g. absolute fold change (without sign), e.g. 1.0
	"min-padj=f" => \$min_padj # minimum adjusted p-value
);

# get gene IDs
my (@diffgenes);
my %ensembl2hgnc;
open(G, "$deseq_file") or die "ERROR: Could not read file $deseq_file\n";
<G>;
while(<G>)
{
	my ($ensembl, $hgnc_symbol, $description, $baseMean, $log2FoldChange, $lfcSE, $stat, $pvalue, $padj) = split("\t");
	$ensembl2hgnc{$ensembl} = $hgnc_symbol;
	
	next if ($padj ne "NA" and $padj > $min_padj);
	next if ($log2FoldChange eq "NA");
	
	push(@diffgenes, $ensembl) if ($direction eq "up" and $log2FoldChange >= $foldchange);
	push(@diffgenes, $ensembl) if ($direction eq "down" and $log2FoldChange <= -$foldchange);
	push(@diffgenes, $ensembl) if ($direction eq "up_or_down" and abs($log2FoldChange) >= $foldchange);
}
close(G);

# connect to DAVID Web Service
my $soap = SOAP::Lite                             
	-> uri('http://service.session.sample')                
	-> proxy('http://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService',
			cookie_jar => HTTP::Cookies->new(ignore_discard=>1));

# user authentication by email address
# For new user registration, go to http://david.abcc.ncifcrf.gov/webservice/register.htm
my $check = $soap->authenticate('christian.frech@ccri.at')->result;
print STDERR "User authentication: $check\n";

die "ERROR: Could not authenticate at http://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService\n"
	if (lc($check) ne "true"); 


# list conversion types
# AFFYMETRIX_3PRIME_IVT_ID,AFFYMETRIX_EXON_GENE_ID,AFFYMETRIX_SNP_ID,AGILENT_CHIP_ID,AGILENT_ID,AGILENT_OLIGO_ID,ENSEMBL_GENE_ID,ENSEMBL_TRANSCRIPT_ID,
# ENTREZ_GENE_ID,FLYBASE_GENE_ID,FLYBASE_TRANSCRIPT_ID,GENBANK_ACCESSION,GENOMIC_GI_ACCESSION,GENPEPT_ACCESSION,ILLUMINA_ID,IPI_ID,MGI_ID,PFAM_ID,PIR_ID,
# PROTEIN_GI_ACCESSION,REFSEQ_GENOMIC,REFSEQ_MRNA,REFSEQ_PROTEIN,REFSEQ_RNA,RGD_ID,SGD_ID,TAIR_ID,UCSC_GENE_ID,UNIGENE,UNIPROT_ACCESSION,UNIPROT_ID,
# UNIREF100_ID,WORMBASE_GENE_ID,WORMPEP_ID,ZFIN_ID
# my $conversionTypes = $soap ->getConversionTypes()->result;
# print STDERR "\nConversion Types: \n$conversionTypes\n"; 
	 
# list all annotation category names
# BBID,BIND,BIOCARTA,BLOCKS,CGAP_EST_QUARTILE,CGAP_SAGE_QUARTILE,CHROMOSOME,COG_NAME,COG_ONTOLOGY,CYTOBAND,DIP,EC_NUMBER,ENSEMBL_GENE_ID,ENTREZ_GENE_ID,
# ENTREZ_GENE_SUMMARY,GENETIC_ASSOCIATION_DB_DISEASE,GENERIF_SUMMARY,GNF_U133A_QUARTILE,GENETIC_ASSOCIATION_DB_DISEASE_CLASS,GOTERM_BP_2,GOTERM_BP_1,
# GOTERM_BP_4,GOTERM_BP_3,GOTERM_BP_FAT,GOTERM_BP_5,GOTERM_CC_1,GOTERM_BP_ALL,GOTERM_CC_3,GOTERM_CC_2,GOTERM_CC_5,GOTERM_CC_4,GOTERM_MF_1,GOTERM_MF_2,
# GOTERM_CC_FAT,GOTERM_CC_ALL,GOTERM_MF_5,GOTERM_MF_FAT,GOTERM_MF_3,GOTERM_MF_4,HIV_INTERACTION_CATEGORY,HIV_INTERACTION_PUBMED_ID,GOTERM_MF_ALL,
# HIV_INTERACTION,KEGG_PATHWAY,HOMOLOGOUS_GENE,INTERPRO,OFFICIAL_GENE_SYMBOL,NCICB_CAPATHWAY_INTERACTION,MINT,PANTHER_MF_ALL,PANTHER_FAMILY,
# PANTHER_BP_ALL,OMIM_DISEASE,PFAM,PANTHER_SUBFAMILY,PANTHER_PATHWAY,PIR_SUPERFAMILY,PIR_SUMMARY,PIR_SEQ_FEATURE,PROSITE,PUBMED_ID,REACTOME_INTERACTION,
# REACTOME_PATHWAY,PIR_TISSUE_SPECIFICITY,PRINTS,PRODOM,PROFILE,SMART,SP_COMMENT,SP_COMMENT_TYPE,SP_PIR_KEYWORDS,SCOP_CLASS,SCOP_FAMILY,
# SCOP_FOLD,SCOP_SUPERFAMILY,UP_SEQ_FEATURE,UNIGENE_EST_QUARTILE,ZFIN_ANATOMY,UP_TISSUE,TIGRFAMS,SSF,UCSC_TFBS
# my $allCategoryNames= $soap ->getAllAnnotationCategoryNames()->result;	 	  	
# print STDERR  "\nAll available annotation category names: \n$allCategoryNames\n";
 
my $default_categories = "KEGG_PATHWAY,OMIM_DISEASE,COG_ONTOLOGY,SP_PIR_KEYWORDS,GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,BBID,BIOCARTA,INTERPRO,PIR_SUPERFAMILY,SMART,NCICB_CAPATHWAY_INTERACTION";
my $additional_categories = "";
$soap->setCategories("$default_categories,$additional_categories")
 	or croak "ERROR: could not set categories\n";
# print "$validated_categories[0]\n";
 
 #addList
print STDERR "Testing following genes for enrichment: ", join(",", @diffgenes), "\n";

 my $inputIds = join(",", @diffgenes);
 my $idType = 'ENSEMBL_GENE_ID';
 my $listName = 'up';
 my $listType=0;
 #to add background list, set listType=1
 my $list = $soap->addList($inputIds, $idType, $listName, $listType)->result;
 print STDERR "Percentage of mapped genes: $list\n"; 
 
#list all species  names
my $allSpecies = $soap->getSpecies()->result;	 	  	
print STDERR  "\nAll species: \n$allSpecies\n"; 
#list current species  names
my $currentSpecies= $soap->getCurrentSpecies()->result;	 	  	
print STDERR  "Current species: $currentSpecies\n"; 

#set user defined species 
#my $species = $soap ->setCurrentSpecies("1")->result;

#print STDERR "\nCurrent species: \n$species\n"; 
 
#set user defined categories 
#my $categories = $soap ->setCategories("BBID,BIOCARTA,COG_ONTOLOGY,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART,UP_SEQ_FEATURE")->result;
#to use DAVID default categories, send empty string to setCategories():
my $categories = $soap->setCategories("")->result;
print STDERR "Valid categories: $categories\n";  

my $thd = 0.1;
my $ct = 2;
my $chartReport = $soap->getChartReport($thd, $ct);

my @chartRecords = $chartReport->paramsout;

#shift(@chartRecords,($chartReport->result));
#print $chartReport->result."\n";
print STDERR "Total chart records: ".(@chartRecords+1)."\n";
#my $retval = %{$chartReport->result};

print "Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\tFDR\n";

my @chartRecordKeys = keys %{$chartReport->result};
	
#print "@chartRecordKeys\n";
	
my @chartRecordValues = values %{$chartReport->result};
		
for my $j (0 .. (@chartRecords-1))
{			
	my %chartRecord = %{$chartRecords[$j]};
	my $categoryName = $chartRecord{"categoryName"};
	my $termName = $chartRecord{"termName"};
	my $listHits = $chartRecord{"listHits"};
	my $percent = $chartRecord{"percent"};
	my $ease = $chartRecord{"ease"};
	my $genes = $chartRecord{"geneIds"};
	my @hgncs;
	map { push(@hgncs, $ensembl2hgnc{$_} ? $ensembl2hgnc{$_} : $_) } split(", ", $genes);
	my $listTotals = $chartRecord{"listTotals"};
	my $popHits = $chartRecord{"popHits"};
	my $popTotals = $chartRecord{"popTotals"};
	my $foldEnrichment = $chartRecord{"foldEnrichment"};
	my $bonferroni = $chartRecord{"bonferroni"};
	my $benjamini = $chartRecord{"benjamini"};
	my $FDR = $chartRecord{"afdr"};
					
	print "$categoryName\t$termName\t$listHits\t$percent\t$ease\t".join(",", @hgncs)."\t$listTotals\t$popHits\t$popTotals\t$foldEnrichment\t$bonferroni\t$benjamini\t$FDR\n";				 
}		  	

print STDERR "Done.\n";
