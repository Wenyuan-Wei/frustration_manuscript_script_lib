import json
from Bio import PDB
from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio.PDB import PDBParser, PPBuilder
import sys 
import glob
import os 
import shutil

def is_standard_residue(residue):
    try:
        PDB.Polypeptide.three_to_one(residue.resname)
        return True
    except KeyError:
        return False

def align_sequences(pdb_file, fasta_file, chain_id):
    # Parse the PDB file
    parser = PDBParser()
    structure = parser.get_structure('PDB', pdb_file)
    model = structure[0]  # Assuming single model
    chain = model[chain_id]

    # Extract the sequence from the PDB chain, considering the actual residue IDs
    ppb = PPBuilder()
    pdb_seq = ''
    pdb_residue_ids = []  # This will store the actual residue IDs

    for pp in ppb.build_peptides(chain):
        for residue in pp:
            pdb_seq += PDB.Polypeptide.three_to_one(residue.get_resname().strip())  # Get the residue name
            pdb_residue_ids.append(residue.get_id()[1])  # Get the residue ID (resSeq)

    pdb_id=os.path.basename(pdb_file).split(".")[0]
    
    # Write the PDB sequence to a temporary FASTA file
    pdb_fasta_file = f"temp_pdb.fasta"
    with open(pdb_fasta_file, "w") as f:
        f.write(f">{pdb_id}\n{pdb_seq}\n")
    
    # Read the reference FASTA file
    reference_sequence = SeqIO.parse(fasta_file, "fasta")
    for record in reference_sequence:
        if record.id != pdb_id:
            fasta_seq = record.seq
            fasta_id = record.id
    
    # Append the reference sequence to the temporary FASTA file
    with open(pdb_fasta_file, "a") as f:
        f.write(f">{fasta_id}\n{fasta_seq}\n")

    # Align the sequences using ClustalW
    clustalw_cline = ClustalwCommandline("clustalw2", infile=pdb_fasta_file)
    clustalw_cline()

    # Read the alignment result
    alignment_file = pdb_fasta_file.replace(".fasta", ".aln")
    if os.path.exists(alignment_file) and os.path.getsize(alignment_file) > 0:
        alignment = AlignIO.read(alignment_file, "clustal")
    else:
        raise ValueError(f"Alignment file {alignment_file} not found or is empty")

    # Create the mapping dictionary
    mapping = {}
    pdb_aligned = None
    fasta_aligned = None

    for rec in alignment:
        if rec.id == pdb_id:
            pdb_aligned = rec.seq
        else:
            fasta_aligned = rec.seq

    if pdb_aligned and fasta_aligned:
        pdb_index = 0
        fasta_index = 0
        pdb_len = len(pdb_residue_ids)

        for pdb_res, fasta_res in zip(pdb_aligned, fasta_aligned):
            if pdb_res != '-' and fasta_res != '-':
                if pdb_index < pdb_len:
                    mapping[pdb_residue_ids[pdb_index]] = fasta_index + 1
                else:
                    print(f"Warning: pdb_index {pdb_index} out of bounds for pdb_residue_ids with length {pdb_len}")
                pdb_index += 1
                fasta_index += 1
            elif pdb_res != '-':
                pdb_index += 1
            elif fasta_res != '-':
                fasta_index += 1

    # Clean up temporary files
    os.remove(alignment_file.replace(".aln", ".dnd"))

    shutil.move(alignment_file, f"./clustalw_files/{pdb_id}.aln")

    return mapping


def identical(mapping):
    # Calculate the percentage of identical key-value pairs
    identical_pairs = sum(1 for k, v in mapping.items() if int(k) == int(v))
    total_pairs = len(mapping)
    
    if total_pairs == 0:
        return False
    
    percentage_identical = (identical_pairs / total_pairs) * 100
    return percentage_identical > 0

def main():

    pdb_file_path = sys.argv[1]
    fasta_file_path = sys.argv[2]
    ref_json_file = sys.argv[3]
    output_file = sys.argv[4]

    with open(ref_json_file, 'r') as f:
        ref_dict = json.load(f)

    if os.path.exists(output_file):
        print(f"Loading existing mapping from {output_file}")
        with open(output_file, 'r') as out:
            pdb_to_fasta_mapping = json.load(out)
    else:
        pdb_to_fasta_mapping = {}
    
    for key in ref_dict.keys():
        for k in ref_dict[key].keys():
            if k in pdb_to_fasta_mapping.keys():
                continue
            print(k)
            pdb_file = os.path.join(pdb_file_path, f"{k}.pdb")
            chain_id = ref_dict[key][k]["G_alpha_chain"].split(",")[0]
            mapping = align_sequences(pdb_file, os.path.join(fasta_file_path, f"{k}_cat.fasta"), chain_id)

            pdb_to_fasta_mapping[k] = mapping

            with open(output_file, 'w') as out:
                json.dump(pdb_to_fasta_mapping, out, indent=4)


if __name__ == '__main__':
    main()
