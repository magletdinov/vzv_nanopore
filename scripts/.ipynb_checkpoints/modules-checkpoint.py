from pathlib import Path
import numpy as np
from scipy.stats import mode
import editdistance
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pybedtools

ROOT = Path().cwd()
PREPROCESSING = Path("../preprocessing")
ASSEMBLING = Path("../reference_genomes")
DATA = Path("../data")
CLADE_FASTA = DATA / "multifasta_for_clades"
ALIGNMENTS = DATA / "alignments"
CONSENSUSES = DATA / "consensuses"
CONSENSUSES_ALIGNMENTS = DATA / "consensuses_alignments"

def alignment_to_np(to_align):
    align_np = [np.frombuffer(seq_record.seq.encode('utf-8'), dtype=np.int8) for seq_record in SeqIO.parse(to_align, "fasta")]
    align_np = np.stack(align_np)
    return align_np

def delete_insertions(align_np):
    mask = np.isin(align_np[0], np.frombuffer("-".encode('utf-8'), dtype=np.int8))
    return align_np[:, ~mask]

def align_to_consensus(align_np):
    consensus_np = mode(align_np[1:], axis=0, keepdims=False)
    return consensus_np

def consensus_to_seq(consensus_np):
    #consensus = Seq(consensus_np.mode.tobytes().decode("utf-8"))
    consensus = Seq(consensus_np.tobytes().decode("utf-8"))
    return consensus

def are_all_rows_unique(align_np):
    unique_rows = np.unique(align_np, axis=0)
    return unique_rows.shape[0] == align_np.shape[0]

def create_distance_matrix(array):
    # Создание массива для хранения расстояний между строками
    num_rows, num_columns = array.shape
    distances = np.zeros((num_rows, num_rows), dtype="int8")

    # Вычисление расстояний между строками
    for i in range(num_rows):
        for j in range(i + 1, num_rows):
            distances[i, j] = editdistance.eval(array[i], array[j])
            distances[j, i] = distances[i, j]  # расстояние Левенштейна симметрично
    return distances

def check_region(array): 
    number_of_unique_rows = np.unique(array, axis=0, return_counts=True)[1].shape[0]
    return number_of_unique_rows

def insertion_coordinats_finder(array):
    #поиск реальных коордиант инсерций после их вырезания из выравнивания
    insetion_coordinats = []
    for i in np.where(array == 45)[0]:
        #sub_array = a[0][:i]
        insetion_coordinats.append(i-len(insetion_coordinats))
    insetion_coordinats = np.unique(np.asarray(insetion_coordinats))
    return insetion_coordinats

def job_1(to_align):
    align_np = alignment_to_np(to_align)
    align_np_modif = delete_insertions(align_np)
    consensus_np = align_to_consensus(align_np_modif)
    #consensus = consensus_to_seq(consensus_np)
    to_save = DATA / "consensuses" / to_align.parent.name / f"{to_align.name.split('-')[0]}_consensus.npy"
    if to_save.exists():
        print(f"{to_save.name} exists")
        return
    to_save.parent.mkdir(exist_ok=True, parents=True)
    return np.save(to_save, consensus_np.mode)

def job_2(start, window, primary_array):
    if are_all_rows_unique(primary_array[:, start+20:start+window-20]):
        primers_area = pybedtools.BedTool(f'''
                                NC_001348.1 {start} {start+20}
                                NC_001348.1 {start+window-20} {start+window}
                               ''', from_string=True)
        number_of_insertion_problems = primers_area.intersect(insertion_regions).count()
        if number_of_insertion_problems > 0:
            return None

        left_area = primary_array[:, start:start+window][:, :20]
        right_area = primary_array[:, start:start+window][:, -20:]
        assert left_area.shape[1] == right_area.shape[1] == 20
        left_score, right_score = check_region(left_area), check_region(right_area)

        total_distance = create_distance_matrix(primary_array[:, start:start+window]).sum()
        record = SeqRecord(seq=ref_NC.seq[start:start+window], id=f"region {start}:{start+window} in NC_001348.1")
        to_save = PREPROCESSING / "selection_of_unique_region" / f"window_{window}" / f"{number_of_insertion_problems}_insertion_problems" / f"l_score_{left_score}-r_score_{right_score}" / f"total_dist-{total_distance}" / f"region-{start}:{start+window}.fasta"
        #to_save = PREPROCESSING / "selection_of_unique_region" / f"window_{window}" / f"{number_of_insertion_problems}_insertion_problems" / f"region-{start}:{start+window}.fasta"    
        #to_save = PREPROCESSING / "selection_of_unique_region" / f"window_{window}" / f"{number_of_insertion_problems}_insertion_problems" / f"l_score_{left_score}-r_score_{right_score}" / f"region-{start}:{start+window}.fasta"
        if to_save.exists():
            print(f"{to_save.name} exists")
            return None

        to_save.parent.mkdir(exist_ok=True, parents=True)
        SeqIO.write(record, to_save, "fasta")