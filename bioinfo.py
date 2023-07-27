# Author: Tam Ho tamho@uoregon.edu

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

# Last updated: 07/21/2023

__version__ = "0.4"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = set("atgcnATGCN")
RNA_bases = set("augcnAUGCN")

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
    seq = seq.upper()
    Gs = seq.count("G")
    Cs = seq.count("C")
    return (Gs+Cs)/len(seq)

def oneline_fasta():
    '''Takes a FASTA file with wrapped sequence line and output same FASTA file but with no sequence wrapping'''
    import argparse
    def get_args():
        parser = argparse.ArgumentParser(description="Reads FASTA file and concatenate all sequence lines into one single line")
        parser.add_argument("-f", "--filename", help="Path to file", required=True)
        parser.add_argument("-o", "--outputfile", help="Path to output file", type=str)
        return parser.parse_args()

    args=get_args()	

    f = args.filename  # input file name
    o = args.outputfile # output file name


    with open (f,"r") as f_in, open (o,"w") as f_out:
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
