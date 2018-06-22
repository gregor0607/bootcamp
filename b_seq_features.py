import bootcamp_utils

pKas = {'D' : 3.9,
       'E' : 4.3,
       'R' : 12.0,
       'K' : 10.5,
       'H' : 6.08}

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

def fraction_negatives(seq, pH = 7):
    """Number of negative residues a protein sequence"""
    # Convert sequence to upper case
    seq = seq.upper()

    # Check for a valid sequence
    is_valid_sequence(seq)
    
    # Counts fraction of 'E'
    pKaE = pKas.get('E')

    f_E = (10**-pKaE / 10**-pH) / (1 + (10**-pKaE / 10**-pH))
    
    # Counts fraction of 'D'
    pKaD = pKas.get('D')

    f_D = (10**-pKaD / 10**-pH) / (1 + (10**-pKaD / 10**-pH))
    
    # return fractional total negative charge
    
    return seq.count('E')*f_E + seq.count('D')*f_D
    
    
def fraction_positives(seq, pH = 7):
    """Number of positive residues a protein sequence"""
    # Convert sequence to upper case
    seq = seq.upper()

    # Check for a valid sequence
    is_valid_sequence(seq)
    
      # Counts fraction of 'R'
    pKaR = pKas.get('R')

    f_R = (10**-pKaR / 10**-pH) / (1 + (10**-pKaR / 10**-pH))
    
    # Counts fraction of 'K'
    pKaK = pKas.get('K')

    f_K = (10**-pKaK / 10**-pH) / (1 + (10**-pKaK / 10**-pH))
    
    # Counts fraction of 'H'
    pKaH = pKas.get('H')

    f_H = (10**-pKaH / 10**-pH) / (1 + (10**-pKaH / 10**-pH))
    
    # return fractional total positive charge
    
    if seq.count('R') + seq.count('K') + seq.count('H') == 0:
        return 0
    else:
        return (1 - seq.count('R') * f_R) + (1 - seq.count('K') * f_K) + (1 - seq.count('H') * f_H)
    

def fractional_net_charge(seq, pH = 7):
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
    net_charge = fraction_positives(seq, pH) - fraction_negatives(seq, pH)
    
    return net_charge
    
    