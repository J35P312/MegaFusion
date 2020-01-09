# MegaFusion
Convert RNA fusion files to SV VCF. MegaFusion should function on most recent versions of python.
MegaFusion accepts a fusion transcript file produced from any tool (such as Ariba), as well as a JSON file, specifying which columns to put in the output vcf file.

# Install
MegaFusion requires some recent version of python.
MegaFusion do not require any particular modules.

# Run
MegaFusion is run using the following command:

	python MegaFusion.py --json input.json --fusion fusion.tab > output.vcf
	
Specify a sample name (defualt = Bob)

	python MegaFusion.py --json input.json --fusion fusion.tab --sample Sven > output.notbob.vcf

The json parameter specify the path of a JSON config file, and the fusion parameter specify the path to a fusion file (produced using Ariba, star-fusion or similar tool).
A JSON files for converting ariba and star-fusion files are provided in the json folder.

# The JSON config file

A config file needs to be created for each fusion transcript caller. Examples are provided in the json folder.

header:
	The header key indicate a letter or word for recognizing the header of the input fusion transcript file.
	
source:
	This key indicate the source of the input fusion transcript file; The source may be a pipeline or fusion-caller transcript caller.
 
delimiter:
	The delimiter of the fusion transcript file, such as "\t" if tav delimited or "," if comma separated.

required:
	A dictionary specifying some required information, these entries must be set properly, or the vcf may be invalid.
	These entries specify which column (counting from 0) to retrieve such required information, and if the information in that column needs to be split into multiple entries. 

		split_reads:
			The number of reads supporting the fusion event
			example:
				"split_reads": 12,

			indicating that the split_read information is found in column 12 of the input fusion transcript file.

		discordant_reads:
			The number of fragments/pairs supporting the fusion event
		geneA:
			one of the genes that are part of the fusion	
		geneB:
			The other gene that is part of the fusion

		ChromosomeA:
			The chromosome of geneA. 
			example:
				ChromosomeA:{"column":4,
                                	"delimiter": ":",
                                	"element": 0
                                },
			The information of a single column may be split to allow the extraction of the reuqired entry.
			This is done by specifying a column of the input fusion file, a delimiter (i.e which charachter to split at), as well as which element to extract.

		ChromosomeB:
			The chromsome of geneB

		posA:
			The position of the breakpoint within gene A

		posB:
			The position of the breakpoint within gene B

		strand1:
			The orientation/strand of the breakpoint at gene A

		strand2:
			The orientation/strand of the breakpoint at gene B

Filter:
	Set the FILTER column of the vcf file. 

	"filter":{
		"filter":1,
		"column":19,
		"delimiter": ",",
		"pass":"."
	},

	filter - enable filtering by setting this entry to any number except 0.
	column - the column describing which filter to apply (counting from 0 in the input fusion transcript file).
	delimiter - some callers may specify multiple filters, delimited by some charachter (commonly ",")
	pass - set the filter to PASS if the following string is found (in the input fusion transcript file)

	Disable filtering by setting the filter entry to 0;

	"filter":{
		"filter":0,
		"column":19,
		"delimiter": ",",
		"pass":"."
	},

	When disabled, all fusions are set to PASS.

	All fusion events are written to the vcf, regardless of filter.

Custom information:
	Custom information may be added to the INFO or FORMAT column of the vcf file.
	The custom dictionary may also be left empty if not custom information is desired:

		custom:{}

	A custom entry is specified through the following syntax:
	"custom": {
		"PEPTIDE":{
		"column":-2,
		"description":"peptide sequence",
		"type":"String",
		"entry":"INFO",
		"none":"."
		}
	}

	Here PEPTIDE specify the name of the entry (as written in the output vcf).
	column:
		 specify which column to extract (from the input fusion transcript file). The column may be counted from 0, or backwards (-1 is the last column, -2 second last etc)
	description:
		The description written in the vcf header
	type:
		a type (String,Integer, or Float), should be set in accordance to the vcf format.
	entry:
		The column of the custom information (INFO or FORMAT).
	none:
		(optional) Do not add the entry if equal to the following string

	Any number of custom entries may be selected.
	
