# Author: Tam Ho tamho@uoregon.edu

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

# Last updated: 08/10/2023

__version__ = "0.6"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = set("atgcnATGCN")
RNA_bases = set("augcnAUGCN")
base_subs = {"A":"T","T":"A","G":"C","C":"G","N":"N"}

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter) - 33


def qual_score(phred_score: str) -> float:
    '''Calculate the average quality score of the whole phred score string'''
    sum: int = 0
    for i in range (len(phred_score)):
        score = phred_score[i]
        converted_score = convert_phred (score)
        sum += converted_score
    return sum / len(phred_score)

def validate_base_seq(seq,RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNA_bases if RNAflag else DNA_bases)

def gc_content(seq):
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq (seq) == True , "non-DNA"
    seq = seq.upper()
    Gs = seq.count("G")
    Cs = seq.count("C")
    return (Gs+Cs)/len(seq)

def oneline_fasta(file):
    '''Takes a FASTA file with wrapped sequence line and output same FASTA file but with no sequence wrapping'''
    with open (file,"r") as f_in, open (f'{file}_oneline.fa',"w") as f_out:
        aaseq = ""
        line = f_in.readline()
        f_out.write(line)
        for line in f_in:
            if line[0] == ">":
                f_out.write(aaseq+"\n")
                f_out.write(line)
                aaseq = ""
            else: 
                aaseq += line.strip()
        f_out.write(aaseq+"\n")

def calc_median (sorted_list: list) -> float:
    '''Takes in a sorted list and returns the median of all values in the list as decimal'''
    n = len(sorted_list) // 2
    if len(sorted_list)%2 == 0:
        median = ((sorted_list[n] + sorted_list[n-1])/2.0)
    else: 
        median = (sorted_list[n])
    return median

def rev_comp_DNA (DNA:str) -> str:
    '''Takes in a string of DNA sequence, flip this entire sequence, loop through it and replace each nucleotide with its complement, then return the reverse complemented DNA sequence'''
    rc_DNA = ""
    r_DNA = DNA [::-1]
    for i in range (len(DNA)):
        rc_DNA += base_subs[r_DNA[i]]
    return rc_DNA

if __name__ == "__main__":
    # Test validate_base_seq()
    assert validate_base_seq("AAUAGAU", True), "RNA test failed"
    print("RNA test passed!")
    assert validate_base_seq("AATAGAT"), "DNA test failed"
    print("DNA test passed!")
    assert validate_base_seq("R is the best!")==False, "R sux"
    print("non-nucleic test passsed!")
    # Test gc_content()
    assert gc_content("GCGCGC") == 1
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATGCAT") == 0.5
    print("correctly calculated GC content")
    # Test convert_phred()
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")
    # Test qual_score()
    assert qual_score("EEE") == 36
    assert qual_score("#I") == 21
    assert qual_score("EJ") == 38.5
    print("You calculated the correct average phred score")
    # Test calc_median()
    assert calc_median([1,2,3]) == 2, "calc_median function does not work for odd length list"
    assert calc_median([1,2]) == 1.5, "calc_median function does not work for even length list"
