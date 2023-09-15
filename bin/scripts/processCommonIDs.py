"""Process for the commond IDS

By using the extracted TaxIds from the species, this code returns
the LCA of the species by using the ete3 tools library.

"""


from ete3 import NCBITaxa
import argparse


def get_species_id(id_taxa, ncbi):
    lineage = ncbi.get_rank(ncbi.get_lineage(id_taxa)).items()
    key_genus = [id[0] for id in lineage if id[1] == "genus"]
    return key_genus[0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("text_ids", nargs="+")
    args = parser.parse_args()
    ncbi = NCBITaxa()

    if len(args.text_ids) > 0:
        list_species_ids = []
        for gen_id in args.text_ids:
            species_id = get_species_id(gen_id, ncbi)
            list_species_ids.append(species_id)
        if len(set(list_species_ids)) == 1:
            with open("common_species_id.txt", "w") as f:
                f.write(str(list_species_ids[0]))
        else:
            raise ValueError("Not common genus found for species")
    else:
        with open("common_species_id.txt", "w") as f:
                f.write(str(get_species_id(args.text_ids)))
