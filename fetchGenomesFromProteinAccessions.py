from Bio import Entrez
import time
import xml.etree.ElementTree as ET

# setting email for Entrez
Entrez.email = "idamico@andrew.cmu.edu"

# get_genome_assembly takes as input an NCBI protein accession number.
# It returns the corresponding NCBI nucleotide accession number for the first genome in NCBI that encodes this protein. If there are no genomes that encode this protein on NCBI, it will return an error.
def get_genome_assembly(protein_accession):
    try:
        # search for corresponding nucleotide entries for this protein accession with Entrez
        handle = Entrez.elink(dbfrom="protein",db="nucleotide",id=protein_accession)

        # read the Entrez data manually (Entrez.read was bugging out)
        xml_data = handle.read()
        root = ET.fromstring(xml_data)

        # close our Entrez handle
        handle.close()

        # grab the first linked nucleotide 
        for linksetdb in root.findall(".//LinkSetDb"):
            if not linksetdb.findall("Link"):
                return (protein_accession, "No Link to Database, Error Level 0.25")

            for link in linksetdb.findall("Link"):
                nucleotide_id = link.findtext("Id")

                # if we have a nucleotide ID, we're going to find its NCBI accession 
                if nucleotide_id:
                    # look for a RefSeq accession 
                    handle = Entrez.efetch(db="nucleotide", id=nucleotide_id, rettype="refseq", retmode="xml")
                    records = Entrez.read(handle)

                    # close our handle
                    handle.close()

                    # loop through all our returned RefSeq accessions (can be a lot)
                    for entries in records:
                        # we need to find ones with "GBSeq_xrefs" dictionaries, as these have the actual GCF accession # in them
                        if "GBSeq_xrefs" in entries:
                            # loop over items in the "GBSeq_xrefs" dictionary
                            for types in entries["GBSeq_xrefs"]:
                                if "GBXref_dbname" in types:
                                    # we specifically want the Assembly accession ID
                                    if types["GBXref_dbname"]=='Assembly':
                                        # the first one we find, we will return
                                        return (protein_accession,types["GBXref_id"])

                                    # no Assembly accessions in our "GBSeq_xrefs" dictionary
                                    else:
                                        return (protein_accession, "No RefSeq Genome Accession, Error Level 3")

                                # no "GBXref_dbname" in our "GBSeq_xrefs" dictionary
                                else:
                                    return (protein_accession, "No RefSeq Genome Accession, Error Level 2")

                        # no "GBSeq_xrefs" dictionaries in this nucleotide accession
                        else:
                            return (protein_accession, "No RefSeq Genome Accession, Error Level 1")
                
                else:
                    return (protein_accession, "No Nucleotide ID, Error Level 0.5")

    except:
        # if we hit a wall anywhere in that process, we'll return an error
        return (protein_accession, "No Nucleotide Information, Error Level 0")

# list of our protein IDs
protein_ids = [
    "MFO7165040.1",
    "MFN8226402.1",
    "MFN8070568.1",
    "MFN8068980.1",
    "MFN8043739.1",
    "MFM9032779.1",
    "MEX0579565.1",
    "MET0756685.1",
    "MET0699879.1",
    "MET0455862.1",
    "MEO9220447.1",
    "MEN3319087.1",
    "MEI7546272.1",
    "MEI7518123.1",
    "MEE9243538.1",
    "MDT7760835.1",
    "MDT5391486.1",
    "MDT5387419.1",
    "MDT5355714.1",
    "MDT5348027.1",
    "MDT5342622.1",
    "MDT5336884.1",
    "MDT5318321.1",
    "MDT5257240.1",
    "MDT5236085.1",
    "MDT5226017.1",
    "MDT5222019.1",
    "MDT5221117.1",
    "MDT5186294.1",
    "MDT5159697.1",
    "MDT5155044.1",
    "MDT5142647.1",
    "MDT5138530.1",
    "MDT5114213.1",
    "MDT5093397.1",
    "MDT5088035.1",
    "MDT5083807.1",
    "MDT5082889.1",
    "MDT5075045.1",
    "MDT5064946.1",
    "MDT5044703.1",
    "MDT5020888.1",
    "MDR3656772.1",
    "MDD4868433.1",
    "MCX6488706.1",
    "MCX6484175.1",
    "MCX6479031.1",
    "MCW2731434.1",
    "MCW2653947.1",
    "MCW2592405.1",
    "MCW2557663.1",
    "MCW2553412.1",
    "MCW2521450.1",
    "MCU1694367.1",
    "MCE9515400.1",
    "MCB1265451.1",
    "MCB0948228.1",
    "MCB0941027.1",
    "MCB0934384.1",
    "MCB0926215.1",
    "MBY0286386.1",
    "MBX9640984.1",
    "MBX7450248.1",
    "MBV9722919.1",
    "MBV9350010.1",
    "MBV9318471.1",
    "MBV9090688.1",
    "MBV8964618.1",
    "MBV8928018.1",
    "MBV8860621.1",
    "MBV8179382.1",
    "MBI3689161.1",
    "MBI3223522.1",
    "MBI3216428.1",
    "HZQ32820.1",
    "HZN79198.1",
    "HZC91290.1",
    "HZC51810.1",
    "HZC09767.1",
    "HYR14823.1",
    "HYQ35909.1",
    "HYO03946.1",
    "HYB37475.1",
    "HXO51188.1",
    "HXO49300.1",
    "HXO11067.1",
    "HXA89225.1",
    "HXA11572.1",
    "HWX99089.1",
    "HWT49419.1",
    "HWF30324.1",
    "HVQ83177.1",
    "HUL98416.1",
    "HUE33280.1",
    "HTZ16453.1",
    "HTX96141.1",
    "HTM84663.1",
    "HTI74448.1",
    "HPY26259.1",
    "HOB49349.1",
    "HME13920.1",
    "HKV18549.1",
    "HKH52309.1",
    "HEY6819389.1",
    "HEY6645620.1",
    "HEY5844673.1",
    "HEY5843588.1",
    "HEY5841728.1",
    "HEY5150782.1",
    "HEY2501814.1",
    "HEY2448391.1",
    "HEY2196988.1",
    "HEY1441860.1",
    "HEY0224811.1",
    "HEX9176916.1",
    "HEX7426095.1",
    "HEX5252895.1",
    "HEX5144392.1",
    "HEX4587596.1",
    "HEX4558481.1",
    "HEX4391829.1",
    "HEX3546307.1",
    "HEX3288378.1",
    "HEX2399371.1",
    "HEX2286788.1",
    "HEX2214097.1",
    "HEV7851787.1",
    "HEV7583870.1",
    "HEV7421491.1",
    "HEV7418343.1",
    "HEV7361846.1",
    "HET9889137.1",
    "HET9563769.1",
    "HET7073938.1",
    "HCA51568.1",
    "GJF10822.1",
    "CAM4325576.1",
    "AQT80856.1"

]

# initialize our results
results = []

# loop over all our protein IDs
for pid in protein_ids:
    # kick off the genome accession finding process for each protein ID
    result = get_genome_assembly(pid)

    # add the results to our result list
    results.append(result)

    # so we don't piss NCBI off
    time.sleep(0.34)  
    
# print out our results
for protein, genome in results:
    print(f"{protein}\t{genome}")
