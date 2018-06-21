import bootcamp_utils

def is_valid_sequence(seq):
    for aa in seq:
        if aa not in bootcamp_utils.aa.keys():
            raise RuntimeError(aa + ' is not a valid amino acid.')

def number_negatives(seq, pH = 7):
    """Number of negative residues a protein sequence"""
    # Convert sequence to upper case
    seq = seq.upper()

    # Check for a valid sequence
    is_valid_sequence(seq)
    
    # Count E's and D's, since these are the negative residues, includes      pH dependence
    
    if pH < 4:
        return 0
    else:
        return seq.count('E') + seq.count('D')

def number_positives(seq, pH = 7):
    """Number of positive residues a protein sequence"""
    # Convert sequence to upper case
    seq = seq.upper()

    # Check for a valid sequence
    is_valid_sequence(seq)
    
    # Count R's and K's and H's, since these are the positive residues,     includes pH dependence
    
    if pH < 8:
        return seq.count('R') + seq.count('K') + seq.count('H')
    elif pH >= 8 and pH < 10:
        return seq.count('R') + seq.count('K')
    else:
        return 0

def net_charge(seq, pH = 7):
    """Calculates the net charge of a protein
    Parameters
    ----------
    seq : string
    Protein sequence
    
    pH : int or float
    pH of solution
    """
    
    # Convert sequence to upper case
    seq = seq.upper()

    # Check for a valid sequence
    is_valid_sequence(seq)
    
    # Calculates net charge
    net_charge = number_positives(seq, pH) - number_negatives(seq, pH)
    
    return net_charge
    
    