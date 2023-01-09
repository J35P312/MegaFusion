import argparse
import json

def retrieve_required_entry(data,tag,content):
	if type(data["required"][tag]) is dict:
		info=content[ data["required"][tag]["column"] ].split(data["required"][tag]["delimiter"])[data["required"][tag]["element"] ]
	else:
		info=content[ data["required"][tag] ]
	return(info)


version = "0.0.1"
parser = argparse.ArgumentParser("""Fusion2VCF-{}: convert gene fusion tab files to SV vcf""".format(version))
parser.add_argument('--fusion'        ,required=True , type=str,  help="path to fusion tab file")
parser.add_argument('--json'        ,required=True , type=str,  help="path to a config json file")
parser.add_argument('--sample'        , type=str,default="Bob",  help="Sample name (default=Bob")
args = parser.parse_args()
args.version=version

with open(args.json) as json_file:
    data = json.load(json_file)

print("##fileformat=VCFv4.1")
print("##source={}".format(data["source"]))
print("##ALT=<ID=BND,Description=\"Break end\">")
print("##INFO=<ID=CHRA,Number=1,Type=String,Description=\"Chromosome A\">")
print("##INFO=<ID=CHRB,Number=1,Type=String,Description=\"Chromosome B\">")
print("##INFO=<ID=GENEA,Number=.,Type=String,Description=\"Gene A\">")
print("##INFO=<ID=GENEB,Number=.,Type=String,Description=\"Gene B\">")
print("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">")
print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
print("##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"Number of paired-ends that support the event\">")
print("##FORMAT=<ID=RV,Number=1,Type=Integer,Description=\"Number of split reads that support the event\">")

required_format="GT:DV:RV"
for entry in sorted(data["custom"].keys()):
	try:
		print ("##{}=<ID={},Number=.,Type={},Description=\"{}\">".format(data["custom"][entry]["entry"],entry,data["custom"][entry]["type"],data["custom"][entry]["description"]) )
		if "FORMAT" == data["custom"][entry]["entry"]:
			required_format+= ":{}".format(entry)
	except:
		print ("invalid config")
		quit()

filter_string="##FILTER=<ID={},Description=\"Not up to snuff\">"

filter_set=set([])

if data["filter"]["filter"]:
	for line in open(args.fusion):
		if line.startswith(data["header"]):
			continue
		content=line.strip().split(data["delimiter"])
		if data["filter"]["pass"] == content[data["filter"]["column"]]:
			continue

		if not "delimiter" in data["filter"]:
			if not content[data["filter"]["column"]] in filter_set:
				print (filter_string.format(content[data["filter"]["column"]]) )
			filter_set.add( content[data["filter"]["column"]])
		else:
			for entry in content[data["filter"]["column"]].split(data["filter"]["delimiter"]):
				if not entry in filter_set:
					print (filter_string.format(entry) )
			filter_set.add(entry)

print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format(args.sample))

required_INFO="SVTYPE=BND;CHRA={};CHRB={};GENEA={};GENEB={};ORIENTATION={},{}"
required_columns="{}\t{}\t{}\tN\t{}\t{}\t{}\t{}\t{}\t{}"
required_sample="./1:{}:{}"
i=1
for line in open(args.fusion):
	if line.startswith(data["header"]):
		continue

	content=line.strip().split(data["delimiter"])

	chrA=retrieve_required_entry(data,"chromosomeA",content)
	chrB=retrieve_required_entry(data,"chromosomeB",content)
	posA=retrieve_required_entry(data,"posA",content)
	posB=retrieve_required_entry(data,"posB",content)
	geneA=retrieve_required_entry(data,"geneA",content)
	geneB=retrieve_required_entry(data,"geneB",content)

	strand1=retrieve_required_entry(data,"strand1",content)
	strand2=retrieve_required_entry(data,"strand2",content)


	if not strand1 in ["+","-"] or not strand2 in ["+","-"]:
		altA="N[{}:{}[".format(chrB,posB)
		#altB="]{}:{}]N".format(chrB,posB)
	elif strand1 == "-" and strand2 == "-":
		altA="[{}:{}[N".format(chrB,posB)
		#altB="[{}:{}[N".format(chrB,posB)
	elif strand1 == "+" and strand2 == "-":
		altA="N]{}:{}]".format(chrB,posB)
		#altB="N]{}:{}]".format(chrB,posB)
	elif strand1 == "-" and strand2 == "+":
		altA="N]{}:{}]".format(chrB,posB)
		#altB="N]{}:{}]".format(chrB,posB)
	else:
		altA="N[{}:{}[".format(chrB,posB)
		#altB="]{}:{}]N".format(chrB,posB)



	INFO=required_INFO.format(chrA,chrB,geneA,geneB,strand1,strand2)
	for entry in data["custom"]:
		if data["custom"][entry]["entry"] == "INFO":
			if "none" in data["custom"][entry]:
				if data["custom"][entry]["none"] ==  content[ data["custom"][entry]["column"] ]:
					continue
			if "remove" in data["custom"][entry]:
				for r in data["custom"][entry]["remove"]:
					content[ data["custom"][entry]["column"] ]=content[ data["custom"][entry]["column"] ].replace(r,"")

			INFO+= ";{}={}".format(entry,content[ data["custom"][entry]["column"] ])


	pairs=retrieve_required_entry(data,"discordant_pairs",content)
	reads=retrieve_required_entry(data,"split_reads",content)
	FORMAT=required_sample.format(pairs,reads)

	for entry in sorted(data["custom"].keys()):
		if data["custom"][entry]["entry"] == "FORMAT":
			if "none" in data["custom"][entry]:
				if data["custom"][entry]["none"] ==  content[ data["custom"][entry]["column"] ]:
					content[ data["custom"][entry]["column"] ] = "."
			if "remove" in data["custom"][entry]:
				for r in data["custom"][entry]["remove"]:
					entry=entry.replace(r,"")
			FORMAT+= ":{}".format(content[ data["custom"][entry]["column"] ])

	ID="{}_Fusion_{}".format(data["source"],i)
	qual="."
	filt="PASS"
	if data["filter"]["filter"]:
		if content[ data["filter"]["column"] ] == data["filter"]["pass"]:
			pass
		else:
			if not "delimiter" in data["filter"]:
				qual=content[ data["filter"]["column"] ]
			else:
				qual=content[ data["filter"]["column"] ].replace(data["filter"]["delimiter"],",")

	print( required_columns.format(chrA,posA,ID+"_1",altA,qual,filt,INFO,required_format,FORMAT) )

	i+=1
