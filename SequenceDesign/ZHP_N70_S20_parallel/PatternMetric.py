# -*- coding: utf-8 -*-
"""
Created on Wed May 15 15:29:31 2024

@author: chatt
"""
from numpy import array, log, diff, ma, mean, arange   

import numpy as np


class StickerPattern:
    def __init__(self, sequence, amino_acid):
        # amino_acid type is the sticker within sequence 
        self.seq = sequence 
        self.aa = amino_acid
    
    @staticmethod
    def calc_H(counts):
        p_arr = array(counts)/sum(counts)
        H = sum(-1 * p_arr*ma.log(p_arr).filled(0))
        return H 
    
    
    def process_patch(self):
        # seq = 'XXPPXXXPX'; P = polar, X = other
        # {'S', 'T', 'Y'} : polar residues 
        N = self.seq.count(self.aa) 
        E = 0 
        for i in range(len(self.seq) - 1): 
            if (self.seq[i] == self.aa and self.seq[i+1] == self.aa):
                E += 1 
        
        return E/(N-1) 
    
    def get_H_frameshifted(self,  bL=3):
        sequence = self.seq.upper()
        amino_acid = self.aa.upper()
        N_seq = len(sequence)
        
        boxLength = bL  # aa per box 
        
        i = 0
        H_seq = []
        nBox = 0

        while i<boxLength:
            j = 0
            myseq = []
            counts = []
            while j < N_seq:
                seq = sequence[j+i: j+i+boxLength]
                if len(seq) != boxLength:
                    idx_periodic = boxLength - len(seq)
                    seq_extra = sequence[:idx_periodic]
                    seq += seq_extra
                counts.append(seq.count(amino_acid))
                myseq.append(seq)
                    
                j += boxLength
            
            nBox = len(counts)
            
            H = self.calc_H(counts)
            
            H_seq.append(H/log(nBox))
            
            i += 1 
        
        
        return mean(array(H_seq))
    
    def calc_H_ave(self):
        H_stat = []
        for bL in [4,5]:
            H = self.get_H_frameshifted(bL=bL)
            H_stat.append(H)
        return mean(array(H_stat))
    

    
    








    
   









