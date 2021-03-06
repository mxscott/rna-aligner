import utils
import evaluation
import shared
from project import *
from timeit import default_timer as timer


toy_genome = 'AGAAGCAACCTGAGAAACGTCTTTGGGATGTGTTCATTCACTTCACAATGATGAACGTTTCTTTTGATTGAGAAGTTTGTAAGGATAACTTTTGTAGAATCTGCAAAGGGATATATGTGAGCCCCTTGATTCCTATGGCAAAATAGGAATTATCTTGAGATAAAAGCCAGACAGAAGATTTCTGAGAAACTTTTTTGTGATGTGTACTTTCATCTCACAGAGTTGAAAAATTCTTTTGATTGAGCAGTTTGGAAACAGTCTTTTCGTATCATCTGCAAATGGATGTTTGGGGCGCTTTGTGGCCTAAGGTGAAAATGGAAACACCTTCACATAAAAACTAGACAGAAGAATTCTGAGGAACCTCTTTATGATGTGTGCATTCATCTCAGATGGGTGAAATTTTCTTTTGATGGAGCAGTTTGGAAACAGTCTTTTTCTAGTATCTGCAGAAGGATATTTGTGAGCGGTGTAAGGCCTATGGTGAAAAAGGAAATATCTTCACATAAAAAACAGACAGAAGCTTTCTGAGGAACTTTTTGTGAGGTGTGCATTCATCTCACCGTGTTGAAACTTTATTTTATTTGAGCAGTTTAGAGACAGTCTTTCTCTGCAATCTGCAAAGGTCTAATTCTGAGCCCTTTGAGGTCTATGGTGAAAAAGAAATATCTTCCCATTTAAACTAGACAGAAGCATTCTGAGGAACTTCGTTGTGATGCCTCTCCATTCATCGGACAGAGTTGAAGGTTTCTTTTAATTCAGCACTTTGGAAAGCATATTTTTGTAGAATCTGCAAAGGGATATTTTTGAGACATTTGAAGCCTAGAGTGAAATAGTAAATATCTTCCCATGAAAACTAGACAGGAGAATTCTGAGAAACTTCATTCTGATGTGTGCATTAACCTCACAGAATTTAACCTTTCTTTTGATTGAGAAGTATGGAAATGGTGGTCTTTTAGAACCTGGAAAGGGATATTTCTTAGCCCTTTGAGGCCTATGGTGAGACTGGAAATATCATCACATGAAAACTAGTCCGAAGCTTTCGGAGAAACTTCTTTGAGATGTGTGCTTTCACCTCACAGAGTTAATCACTTTCTTTTGATTGAGCAGTTTGGAAACACTCTTTCTGTGACATCTGTAAATGGATATTAGGAGTGCTTTGAGGCCAATG'
iso_1 = shared.Isoform('1', [shared.Exon('1', 100, 200), shared.Exon('2', 250, 400)])
iso_2 = shared.Isoform('2', [shared.Exon('1', 100, 200), shared.Exon('3', 700, 850)])
iso_3 = shared.Isoform('3', [shared.Exon('2', 250, 400), shared.Exon('3', 700, 850)])
iso_4 = shared.Isoform('4', [shared.Exon('4', 900, 1000)])

gene1 = shared.Gene('1', [iso_1, iso_2, iso_3])
gene2 = shared.Gene('2', [iso_1, iso_3])
gene3 = shared.Gene('3', [iso_4])

toy_genes = {gene1, gene2, gene3}

read = 'GATTTCTGAGAAACTTTTTTGTGA'

'''
with open('genome.fa') as f:
    f.readline()
    genome = f.readline()

genes = utils.parse_tab_file('genes.tab')

'''
start_construct = timer()
aligner = Aligner(toy_genome, toy_genes)
end_construct = timer()


print('Initialization took: ' + str(end_construct - start_construct) + ' seconds')

#print(aligner.index_dict[100])
#print(aligner.index_dict[350])
#print(aligner.index_dict[1000])

#print(aligner.transcriptome)


start_align = timer()
alignment = aligner.align(read)
end_align = timer()

print('Alignment took: ' + str(end_align - start_align) + ' seconds')

print(alignment)

